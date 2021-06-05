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

#include <math.h>
#include <string.h>

#include "libpolyhedra.h"

/* Reference:
 * Fast and Accurate Computation of Polyhedral Mass Properties
 * Brian Mirtich
 */

/* Based on volInt.c public domain code by Brian Mirtich */

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

struct face {
  double norm[3];
  double w;
  size_t vert[3];
};

struct axis {
  int A;
  int B;
  int C;
};

struct proj_int {
  double P1;
  double Pa;
  double Pb;
  double Paa;
  double Pab;
  double Pbb;
  double Paaa;
  double Paab;
  double Pabb;
  double Pbbb;
};

struct face_int {
  double Fa;
  double Fb;
  double Fc;
  double Faa;
  double Fbb;
  double Fcc;
  double Faaa;
  double Fbbb;
  double Fccc;
  double Faab;
  double Fbbc;
  double Fcca;
};

static void InitFace(struct face *face, const unsigned int *verts, const double *offset, const float *data, size_t fpv) {
  double dx1, dy1, dz1, dx2, dy2, dz2, nx, ny, nz, len, inv;
  
  memset(face, 0, sizeof(*face));
  
  face->vert[0] = fpv * verts[0];
  face->vert[1] = fpv * verts[1];
  face->vert[2] = fpv * verts[2];
  
  dx1 = (double) data[face->vert[1] + 0] - (double) data[face->vert[0] + 0];
  dy1 = (double) data[face->vert[1] + 1] - (double) data[face->vert[0] + 1];
  dz1 = (double) data[face->vert[1] + 2] - (double) data[face->vert[0] + 2];
  dx2 = (double) data[face->vert[2] + 0] - (double) data[face->vert[1] + 0];
  dy2 = (double) data[face->vert[2] + 1] - (double) data[face->vert[1] + 1];
  dz2 = (double) data[face->vert[2] + 2] - (double) data[face->vert[1] + 2];
  nx = dy1 * dz2 - dy2 * dz1;
  ny = dz1 * dx2 - dz2 * dx1;
  nz = dx1 * dy2 - dx2 * dy1;
  len = sqrt(nx * nx + ny * ny + nz * nz);
  if (len == 0)
    inv = 0;
  else
    inv = 1.0 / len;
  face->norm[0] = nx * inv;
  face->norm[1] = ny * inv;
  face->norm[2] = nz * inv;
  face->w =
    -face->norm[0] * (data[face->vert[0] + 0] - offset[0]) +
    -face->norm[1] * (data[face->vert[0] + 1] - offset[1]) +
    -face->norm[2] * (data[face->vert[0] + 2] - offset[2]);
}

static void ProjInt(struct proj_int *pi, const struct axis *axis, const struct face *face, const double *offset, const float *data) {
  double a0, a1, da;
  double b0, b1, db;
  double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
  double a1_2, a1_3, b1_2, b1_3;
  double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
  double Cab, Kab, Caab, Kaab, Cabb, Kabb;
  int count;
  
  memset(pi, 0, sizeof(*pi));
  
  for (count = 0; count < 3; count++) {
    a0 = data[face->vert[count] + axis->A] - offset[axis->A];
    b0 = data[face->vert[count] + axis->B] - offset[axis->B];
    a1 = data[face->vert[(count + 1) % 3] + axis->A] - offset[axis->A];
    b1 = data[face->vert[(count + 1) % 3] + axis->B] - offset[axis->B];
    /* printf("a0 = %f, b0 = %f, a1 = %f, b1 = %f\n", a0, b0, a1, b1); */
    da = a1 - a0;
    db = b1 - b0;
    a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
    b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
    a1_2 = a1 * a1; a1_3 = a1_2 * a1;
    b1_2 = b1 * b1; b1_3 = b1_2 * b1;
    
    C1 = a1 + a0;
    Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
    Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
    Cab = 3*a1_2 + 2*a1*a0 + a0_2; Kab = a1_2 + 2*a1*a0 + 3*a0_2;
    Caab = a0*Cab + 4*a1_3; Kaab = a1*Kab + 4*a0_3;
    Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
    Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;
    
    pi->P1 += db*C1;
    pi->Pa += db*Ca;
    pi->Paa += db*Caa;
    pi->Paaa += db*Caaa;
    pi->Pb += da*Cb;
    pi->Pbb += da*Cbb;
    pi->Pbbb += da*Cbbb;
    pi->Pab += db*(b1*Cab + b0*Kab);
    pi->Paab += db*(b1*Caab + b0*Kaab);
    pi->Pabb += da*(a1*Cabb + a0*Kabb);
  }

  pi->P1 /= 2.0;
  pi->Pa /= 6.0;
  pi->Paa /= 12.0;
  pi->Paaa /= 20.0;
  pi->Pb /= -6.0;
  pi->Pbb /= -12.0;
  pi->Pbbb /= -20.0;
  pi->Pab /= 24.0;
  pi->Paab /= 60.0;
  pi->Pabb /= -60.0;
}

static void FaceInt(struct face_int *fi, const struct axis *axis, const struct face *face, const double *offset, const float *data) {
  double w, na, nb, nc;
  double k1, k2, k3, k4;
  struct proj_int pi;

  ProjInt(&pi, axis, face, offset, data);

  w = face->w;
  na = face->norm[axis->A];
  nb = face->norm[axis->B];
  nc = face->norm[axis->C];
  k1 = 1 / nc; k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

  fi->Fa = k1 * pi.Pa;
  fi->Fb = k1 * pi.Pb;
  fi->Fc = -k2 * (na*pi.Pa + nb*pi.Pb + w*pi.P1);

  fi->Faa = k1 * pi.Paa;
  fi->Fbb = k1 * pi.Pbb;
  fi->Fcc = k3 * (SQR(na)*pi.Paa + 2*na*nb*pi.Pab + SQR(nb)*pi.Pbb
		  + w*(2*(na*pi.Pa + nb*pi.Pb) + w*pi.P1));
  
  fi->Faaa = k1 * pi.Paaa;
  fi->Fbbb = k1 * pi.Pbbb;
  fi->Fccc = -k4 * (CUBE(na)*pi.Paaa + 3*SQR(na)*nb*pi.Paab 
		    + 3*na*SQR(nb)*pi.Pabb + CUBE(nb)*pi.Pbbb
		    + 3*w*(SQR(na)*pi.Paa + 2*na*nb*pi.Pab + SQR(nb)*pi.Pbb)
		    + w*w*(3*(na*pi.Pa + nb*pi.Pb) + w*pi.P1));

  fi->Faab = k1 * pi.Paab;
  fi->Fbbc = -k2 * (na*pi.Pabb + nb*pi.Pbbb + w*pi.Pbb);
  fi->Fcca = k3 * (SQR(na)*pi.Paaa + 2*na*nb*pi.Paab + SQR(nb)*pi.Pabb
		   + w*(2*(na*pi.Paa + nb*pi.Pab) + w*pi.Pa));
}

void LP_MassProperties(const struct lp_vertex_list *in, struct lp_mass_properties *properties) {
  float *data;
  unsigned int *idx;
  size_t fpv, num_vert, num_idx, count;
  double nx, ny, nz, T0, T1[3], T2[3], TP[3], offset[3], r[3];
  struct face face;
  struct axis axis;
  struct face_int fi;
  
  memset(properties, 0, sizeof(*properties));
  
  if ((fpv = LP_VertexList_FloatsPerVert(in)) < 3) {
    printf("Cannot determine mass properties: too few floats per vertex\n");
    return;
  }
  
  data = LP_VertexList_GetVert(in);
  num_vert = LP_VertexList_NumVert(in);
  memset(offset, 0, sizeof(offset));
  for (count = 0; count < num_vert; count++) {
    offset[0] += data[0];
    offset[1] += data[1];
    offset[2] += data[2];
    data += fpv;
  }
  offset[0] /= num_vert;
  offset[1] /= num_vert;
  offset[2] /= num_vert;
  
  T0 = T1[0] = T1[1] = T1[2] 
     = T2[0] = T2[1] = T2[2] 
     = TP[0] = TP[1] = TP[2] = 0;
  
  num_idx = LP_VertexList_NumInd(in) - 2;
  idx = LP_VertexList_GetInd(in);
  data = LP_VertexList_GetVert(in);
  for (count = 0; count < num_idx; count += 3) {
    InitFace(&face, idx + count, offset, data, fpv);
    
    nx = fabs(face.norm[0]);
    ny = fabs(face.norm[1]);
    nz = fabs(face.norm[2]);
    if (nx > ny && nx > nz) axis.C = 0;
    else axis.C = (ny > nz) ? 1 : 2;
    axis.A = (axis.C + 1) % 3;
    axis.B = (axis.A + 1) % 3;
    
    FaceInt(&fi, &axis, &face, offset, data);
    
    T0 += face.norm[0] * ((axis.A == 0 ? fi.Fa : ((axis.B == 0) ? fi.Fb : fi.Fc)));
    T1[axis.A] += face.norm[axis.A] * fi.Faa;
    T1[axis.B] += face.norm[axis.B] * fi.Fbb;
    T1[axis.C] += face.norm[axis.C] * fi.Fcc;
    T2[axis.A] += face.norm[axis.A] * fi.Faaa;
    T2[axis.B] += face.norm[axis.B] * fi.Fbbb;
    T2[axis.C] += face.norm[axis.C] * fi.Fccc;
    TP[axis.A] += face.norm[axis.A] * fi.Faab;
    TP[axis.B] += face.norm[axis.B] * fi.Fbbc;
    TP[axis.C] += face.norm[axis.C] * fi.Fcca;
  }
  
  T1[0] /= 2; T1[1] /= 2; T1[2] /= 2;
  T2[0] /= 3; T2[1] /= 3; T2[2] /= 3;
  TP[0] /= 2; TP[1] /= 2; TP[2] /= 2;
  
  /* Volume */
  properties->volume = T0;
  
  /* Center of mass */
  r[0] = T1[0] / T0;
  r[1] = T1[1] / T0;
  r[2] = T1[2] / T0;
  
  properties->center_of_mass[0] = r[0] + offset[0];
  properties->center_of_mass[1] = r[1] + offset[1];
  properties->center_of_mass[2] = r[2] + offset[2];
  
  /* Inertia tensor */
  properties->inertia_tensor[0] = T2[1] + T2[2];
  properties->inertia_tensor[4] = T2[2] + T2[0];
  properties->inertia_tensor[8] = T2[0] + T2[1];
  
  properties->inertia_tensor[1] = -TP[0];
  properties->inertia_tensor[5] = -TP[1];
  properties->inertia_tensor[2] = -TP[2];
  
  /* Offset inertia tensor to center of mass */
  properties->inertia_tensor[0] -= T0 * (r[1]*r[1] + r[2]*r[2]);
  properties->inertia_tensor[4] -= T0 * (r[2]*r[2] + r[0]*r[0]);
  properties->inertia_tensor[8] -= T0 * (r[0]*r[0] + r[1]*r[1]);
  
  properties->inertia_tensor[1] += T0 * (r[0] * r[1]);
  properties->inertia_tensor[5] += T0 * (r[1] * r[2]);
  properties->inertia_tensor[2] += T0 * (r[2] * r[0]);
  
  properties->inertia_tensor[3] = properties->inertia_tensor[1];
  properties->inertia_tensor[7] = properties->inertia_tensor[5];
  properties->inertia_tensor[6] = properties->inertia_tensor[2];
}
