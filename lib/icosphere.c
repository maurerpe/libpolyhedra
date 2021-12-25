/*
  Copyright (C) 2021 Paul Maurer

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

/* Alorithm:
 * http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <limits.h>
#include <math.h>
#include <string.h>

#include "libpolyhedra.h"
#include "util.h"

static void FindMid(float *mid, float *aa, float *bb, float radius) {
  mid[0] = aa[0] + bb[0];
  mid[1] = aa[1] + bb[1];
  mid[2] = aa[2] + bb[2];
  
  Normalize(mid);
  mid[0] *= radius;
  mid[1] *= radius;
  mid[2] *= radius;
}

static int AddTri(struct lp_vertex_list *vl, float *v1, float *v2, float *v3) {
  if (LP_VertexList_Add(vl, v1) == UINT_MAX)
    return -1;
  if (LP_VertexList_Add(vl, v2) == UINT_MAX)
    return -1;
  if (LP_VertexList_Add(vl, v3) == UINT_MAX)
    return -1;
  
  return 0;
}

static struct lp_vertex_list *SubDivide(const struct lp_vertex_list *in, float radius) {
  struct lp_vertex_list *out;
  size_t num_vert, count;
  float *v1, *v2, *v3, aa[3], bb[3], cc[3];
  
  if ((out = LP_VertexList_New(3, lp_pt_triangle)) == NULL)
    goto err;
  
  num_vert = LP_VertexList_NumInd(in);
  for (count = 0; count < num_vert; count += 3) {
    v1 = LP_VertexList_LookupVert(in, count);
    v2 = LP_VertexList_LookupVert(in, count + 1);
    v3 = LP_VertexList_LookupVert(in, count + 2);
    
    FindMid(aa, v1, v2, radius);
    FindMid(bb, v1, v3, radius);
    FindMid(cc, v2, v3, radius);
    
    if (AddTri(out, v1, aa, bb) < 0)
      goto err2;
    if (AddTri(out, v2, cc, aa) < 0)
      goto err2;
    if (AddTri(out, v3, bb, cc) < 0)
      goto err2;
    if (AddTri(out, aa, cc, bb) < 0)
      goto err2;
  }
  
  return out;

 err2:
  LP_VertexList_Free(out);
 err:
  return NULL;
}

#define ADD_PT(xx, yy, zz)	\
  do {				\
    *cur++ = (xx);		\
    *cur++ = (yy);		\
    *cur++ = (zz);		\
  } while (0)

#define ADD_TRI(aa, bb, cc)						\
  do {									\
    if (AddTri(out, vv + 3 * (aa), vv + 3 * (bb), vv + 3 * (cc)) < 0)	\
      goto err2;							\
  } while (0)

static struct lp_vertex_list *MakeIcohedron(float radius) {
  struct lp_vertex_list *out;
  float s, t, scale, vv[36], *cur;
  
  t = (1.0f + sqrtf(5.0f)) / 2.0f;
  scale = radius / sqrtf(1 + t * t);
  t *= scale;
  s  = scale;

  /* Verticies */
  cur = vv;
  ADD_PT(-s,  t,  0);
  ADD_PT( s,  t,  0);
  ADD_PT(-s, -t,  0);
  ADD_PT( s, -t,  0);

  ADD_PT( 0, -s,  t);
  ADD_PT( 0,  s,  t);
  ADD_PT( 0, -s, -t);
  ADD_PT( 0,  s, -t);

  ADD_PT( t,  0, -s);
  ADD_PT( t,  0,  s);
  ADD_PT(-t,  0, -s);
  ADD_PT(-t,  0,  s);
  
  if ((out = LP_VertexList_New(3, lp_pt_triangle)) == NULL)
    goto err;
  
  /* Faces */
  ADD_TRI(0, 11, 5);
  ADD_TRI(0, 5, 1);
  ADD_TRI(0, 1, 7);
  ADD_TRI(0, 7, 10);
  ADD_TRI(0, 10, 11);
  
  ADD_TRI(1, 5, 9);
  ADD_TRI(5, 11, 4);
  ADD_TRI(11, 10, 2);
  ADD_TRI(10, 7, 6);
  ADD_TRI(7, 1, 8);

  ADD_TRI(3, 9, 4);
  ADD_TRI(3, 4, 2);
  ADD_TRI(3, 2, 6);
  ADD_TRI(3, 6, 8);
  ADD_TRI(3, 8, 9);

  ADD_TRI(4, 9, 5);
  ADD_TRI(2, 4, 11);
  ADD_TRI(6, 2, 10);
  ADD_TRI(8, 6, 7);
  ADD_TRI(9, 8, 1);

  return out;

 err2:
  LP_VertexList_Free(out);
 err:
  return NULL;
}

struct lp_vertex_list *LP_IcoSphere(float radius, int num_subdiv) {
  struct lp_vertex_list *cur, *next;
  int count;
  
  if ((cur = MakeIcohedron(radius)) == NULL)
    goto err;
  
  for (count = 0; count < num_subdiv; count++) {
    if ((next = SubDivide(cur, radius)) == NULL)
      goto err2;
    LP_VertexList_Free(cur);
    cur = next;
  }
  
  return cur;
  
 err2:
  LP_VertexList_Free(cur);
 err:
  return NULL;
}
