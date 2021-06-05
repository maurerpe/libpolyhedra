/*
  Copyright (C) 2020 Paul Maurer

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  
  1. Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.
  
  2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
  
  3. Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <limits.h>
#include <math.h>
#include <string.h>

#include "libpolyhedra.h"
#include "util.h"

struct lp_transform {
  float wxyz[4];
  float trans[3];
  int mat_valid;
  float mat[9];
};

struct lp_transform *LP_Transform_New(void) {
  struct lp_transform *trans;

  if ((trans = malloc(sizeof(*trans))) == NULL)
    goto err;
  memset(trans, 0, sizeof(*trans));
  trans->wxyz[0] = 1;
  
  return trans;
  
 err:
  return NULL;
}

void LP_Transform_Free(struct lp_transform *trans) {
  free(trans);
}

void LP_Transform_Copy(struct lp_transform *dest, const struct lp_transform *src) {
  memcpy(dest, src, sizeof(*dest));
}

void LP_Transform_SetToIdentity(struct lp_transform *trans) {
  memset(trans, 0, sizeof(*trans));
  trans->wxyz[0] = 1;
}

static inline float maxf(float a, float b) {
  if (a >= b)
    return a;
  
  return b;
}

void LP_Transform_SetToMatrix4x4(struct lp_transform *trans, const float *m) {
  float i33, m00, m11, m22, x, y, z;
  
  i33 = 1 / m[4*3+3];
  m00 = m[4*0+0] * i33;
  m11 = m[4*1+1] * i33;
  m22 = m[4*2+2] * i33;

  /* https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/ */
  trans->wxyz[0] = 0.5f * sqrtf(maxf(0, 1 + m00 + m11 + m22));
  x = 0.5f * sqrtf(maxf(0, 1 + m00 - m11 - m22));
  y = 0.5f * sqrtf(maxf(0, 1 - m00 + m11 - m22));
  z = 0.5f * sqrtf(maxf(0, 1 - m00 - m11 + m22));
  trans->wxyz[1] = copysignf(x, (m[4*2+1] - m[4*1+2]) * i33);
  trans->wxyz[2] = copysignf(y, (m[4*0+2] - m[4*2+0]) * i33);
  trans->wxyz[3] = copysignf(z, (m[4*1+0] - m[4*0+1]) * i33);
  Normalize4d(trans->wxyz);
  
  trans->trans[0] = m[4*0+3] * i33;
  trans->trans[1] = m[4*1+3] * i33;
  trans->trans[2] = m[4*2+3] * i33;
}

void LP_Transform_Translate(struct lp_transform *trans, float dx, float dy, float dz) {
  trans->trans[0] += dx;
  trans->trans[1] += dy;
  trans->trans[2] += dz;
}

void LP_Transform_Rotate(struct lp_transform *trans, float angle_rad, float axis_x, float axis_y, float axis_z) {
  struct lp_transform tt;
  float x, y, z, inv_mag, sina;
  
  if (angle_rad == 0)
    return;
  
  LP_Transform_SetToIdentity(&tt);
  
  inv_mag = 1.0 / sqrtf(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
  x = axis_x * inv_mag;
  y = axis_y * inv_mag;
  z = axis_z * inv_mag;
  
  sina = sinf(0.5 * angle_rad);
  tt.wxyz[0] = cosf(0.5 * angle_rad);
  tt.wxyz[1] = sina * x;
  tt.wxyz[2] = sina * y;
  tt.wxyz[3] = sina * z;
  
  LP_Transform_Combine(trans, &tt, trans);
}

void LP_Transform_ApplyQauternion(struct lp_transform *trans, const float *wxyz) {
  struct lp_transform tt;

  LP_Transform_SetToIdentity(&tt);
  tt.wxyz[0] = wxyz[0];
  tt.wxyz[1] = wxyz[1];
  tt.wxyz[2] = wxyz[2];
  tt.wxyz[3] = wxyz[3];
  
  LP_Transform_Combine(trans, &tt, trans);
}

static void MultQuat(float *dest_w, float *dest_xyz, float a_w, const float *a_xyz, float b_w, const float *b_xyz) {
  float cross[3];
  
  Cross(cross, a_xyz, b_xyz);
  
  if (dest_w)
    *dest_w = a_w * b_w - Dot(a_xyz, b_xyz);

  if (dest_xyz) {
    dest_xyz[0] = cross[0] + a_w * b_xyz[0] + b_w * a_xyz[0];
    dest_xyz[1] = cross[1] + a_w * b_xyz[1] + b_w * a_xyz[1];
    dest_xyz[2] = cross[2] + a_w * b_xyz[2] + b_w * a_xyz[2];
  }
}

static void RotPoint(float *dest, const float *wxyz, const float *src) {
  float ww, xyz[3], conj[3];
  
  MultQuat(&ww, xyz, wxyz[0], wxyz + 1, 0, src);
  
  conj[0] = -wxyz[1];
  conj[1] = -wxyz[2];
  conj[2] = -wxyz[3];
  MultQuat(NULL, dest, ww, xyz, wxyz[0], conj);
}

void LP_Transform_Combine(struct lp_transform *dest, const struct lp_transform *a, const struct lp_transform *b) {
  float trans[3];
  
  /* Translation */
  RotPoint(trans, a->wxyz, b->trans);
  dest->trans[0] = trans[0] + a->trans[0];
  dest->trans[1] = trans[1] + a->trans[1];
  dest->trans[2] = trans[2] + a->trans[2];
  
  /* Rotation */
  MultQuat(dest->wxyz, dest->wxyz + 1, a->wxyz[0], a->wxyz + 1, b->wxyz[0], b->wxyz + 1);
  Normalize4d(dest->wxyz);
  dest->mat_valid = 0;
}

void LP_Transform_Invert(struct lp_transform *trans) {
  trans->wxyz[1] = -trans->wxyz[1];
  trans->wxyz[2] = -trans->wxyz[2];
  trans->wxyz[3] = -trans->wxyz[3];
  
  trans->trans[0] = -trans->trans[0];
  trans->trans[1] = -trans->trans[1];
  trans->trans[2] = -trans->trans[2];
  
  RotPoint(trans->trans, trans->wxyz, trans->trans);
}

static void BuildMat(struct lp_transform *trans) {
  float q0, q1, q2, q3, q00, q01, q02, q03, q11, q12, q13, q22, q23, q33;
  
  /* https://automaticaddison.com/how-to-convert-a-quaternion-to-a-rotation-matrix/ */
  q0 = trans->wxyz[0];
  q1 = trans->wxyz[1];
  q2 = trans->wxyz[2];
  q3 = trans->wxyz[3];
  
  q00 = q0 * q0;
  q01 = q0 * q1;
  q02 = q0 * q2;
  q03 = q0 * q3;
  q11 = q1 * q1;
  q12 = q1 * q2;
  q13 = q1 * q3;
  q22 = q2 * q2;
  q23 = q2 * q3;
  q33 = q3 * q3;
  
  trans->mat[3*0 + 0] = 2 * (q00 + q11) - 1;
  trans->mat[3*0 + 1] = 2 * (q12 - q03);
  trans->mat[3*0 + 2] = 2 * (q13 + q02);
  
  trans->mat[3*1 + 0] = 2 * (q12 + q03);
  trans->mat[3*1 + 1] = 2 * (q00 + q22) - 1;
  trans->mat[3*1 + 2] = 2 * (q23 - q01);
  
  trans->mat[3*2 + 0] = 2 * (q13 - q02);
  trans->mat[3*2 + 1] = 2 * (q23 + q01);
  trans->mat[3*2 + 2] = 2 * (q00 + q33) - 1;
  trans->mat_valid = 1;
}

void LP_Transform_Point(const struct lp_transform *trans, float *dest, const float *src, int options) {
  float off[3];
  int count;

  if (!trans->mat_valid)
    BuildMat((struct lp_transform *) trans);
  
  if (options & LP_TRANSFORM_INVERT) {
    if (!(options & LP_TRANSFORM_NO_OFFSET)) {
      off[0] = src[0] - trans->trans[0];
      off[1] = src[1] - trans->trans[1];
      off[2] = src[2] - trans->trans[2];
      src = off;
    }
  
    for (count = 0; count < 3; count++) {
      dest[count] =
	trans->mat[3*0 + count] * src[0] +
	trans->mat[3*1 + count] * src[1] +
	trans->mat[3*2 + count] * src[2];
    }
    
    return;
  }
  
  for (count = 0; count < 3; count++) {
    dest[count] =
      trans->mat[3*count + 0] * src[0] +
      trans->mat[3*count + 1] * src[1] +
      trans->mat[3*count + 2] * src[2];
  }

  if (!(options & LP_TRANSFORM_NO_OFFSET)) {
    dest[0] += trans->trans[0];
    dest[1] += trans->trans[1];
    dest[2] += trans->trans[2];
  }
}

struct lp_vertex_list *LP_Transform_VertexList(const struct lp_transform *trans, const struct lp_vertex_list *src, int options) {
  struct lp_vertex_list *vl;
  size_t fpv, num_vert, num_ind, count;
  float *ff, *vert;
  unsigned int *ind;

  if ((fpv = LP_VertexList_FloatsPerVert(src)) < 3) {
    fprintf(stderr, "Too few floats per vertext to transform\n");
    goto err;
  }
  num_vert = LP_VertexList_NumVert(src);
  num_ind  = LP_VertexList_NumInd(src);
  
  if ((vl = LP_VertexList_New(3, LP_VertexList_PrimativeType(src))) == NULL)
    goto err;
  
  if ((ff = calloc(num_vert, 3 * sizeof(float))) == NULL)
    goto err2;
  
  vert = LP_VertexList_GetVert(src);
  for (count = 0; count < num_vert; count++)
    LP_Transform_Point(trans, &ff[3 * count], &vert[fpv * count], options);

  ind = LP_VertexList_GetInd(src);
  for (count = 0; count < num_ind; count++)
    if (LP_VertexList_Add(vl, &ff[3 * ind[count]]) == UINT_MAX)
      goto err3;
  
  free(ff);
  return vl;

 err3:
  free(ff);
 err2:
  LP_VertexList_Free(vl);
 err:
  return NULL;
}
