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

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <limits.h>
#include <math.h>
#include <string.h>

#include "libpolyhedra.h"

static int AddVert(struct lp_vertex_list *vl, float x, float y, float z) {
  float vert[3];

  vert[0] = x;
  vert[1] = y;
  vert[2] = z;

  if (LP_VertexList_Add(vl, vert) == UINT_MAX)
    return -1;
  
  return 0;
}

struct lp_vertex_list *LP_Cube(float x, float y, float z) {
  struct lp_vertex_list *pts, *hull;
  
  if ((pts = LP_VertexList_New(3, lp_pt_point)) == NULL)
    goto err;
  
  if (AddVert(pts,  x,  y,  z) < 0)
    goto err2;
  if (AddVert(pts,  x,  y, -z) < 0)
    goto err2;
  if (AddVert(pts,  x, -y,  z) < 0)
    goto err2;
  if (AddVert(pts,  x, -y, -z) < 0)
    goto err2;
  if (AddVert(pts, -x,  y,  z) < 0)
    goto err2;
  if (AddVert(pts, -x,  y, -z) < 0)
    goto err2;
  if (AddVert(pts, -x, -y,  z) < 0)
    goto err2;
  if (AddVert(pts, -x, -y, -z) < 0)
    goto err2;
  
  if ((hull = LP_ConvexHull(pts)) == NULL)
    goto err2;
  
  return hull;

 err2:
  LP_VertexList_Free(pts);
 err:
  return NULL;
}
