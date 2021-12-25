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

struct lp_vertex_list *LP_UVSphere(float radius, int segs, int rings) {
  struct lp_vertex_list *pts, *hull;
  int ang_count, azi_count;
  float ang, azi, ang_incr, azi_incr, rr, xx, yy, zz;

  if (segs < 3)
    segs = 3;

  if (rings < 2)
    rings = 2;
  
  if ((pts = LP_VertexList_New(3, lp_pt_point)) == NULL)
    goto err;
  
  AddVert(pts, 0, 0,  radius);
  AddVert(pts, 0, 0, -radius);
  ang_incr = 2 * M_PI / segs;
  azi_incr = M_PI / rings;
  ang = 0;
  azi = azi_incr - M_PI_2;
  for (azi_count = 1; azi_count < rings; azi_count++, azi += azi_incr) {
    rr = radius * cos(azi);
    zz = radius * sin(azi);
    for (ang_count = 0; ang_count < segs; ang_count++, ang += ang_incr) {
      xx = rr * cos(ang);
      yy = rr * sin(ang);
      AddVert(pts, xx, yy, zz);
    }
  }
  
  if ((hull = LP_ConvexHull(pts)) == NULL)
    goto err2;
  
  return hull;

 err2:
  LP_VertexList_Free(pts);
 err:
  return NULL;
}
