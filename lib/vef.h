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

#ifndef VEH_H
#define VEH_H

#include "libpolyhedra.h"

#include "hash.h"

struct face;

struct vert {
  float point[3];
  struct hash *edges;
};

struct edge {
  struct vert *vert[2];
  struct face *face[2];

  int info_vld;
  float  x_vec[3];
  float  z_vec[3];
  float  ang;
};

struct face {
  struct vert *vert[3];
  struct edge *edge[3];
  float norm[3];
  float dist;
  
  int basis_vld;
  float basis_x[3];
  float basis_y[3];

  int coord_2d_vld;
  float v1_x_len;
  float v2_pos[2];
};

struct vef {
  struct hash *verts;
  struct hash *edges;
  struct hash *faces;
  float min[3];
  float max[3];
};

struct vef *Vef_New(const struct lp_vertex_list *vl);
void Vef_Free(struct vef *vef);

void Vef_CalcInfo(struct edge *edge);
void Vef_CalcBasis(struct face *face);
void Vef_CalcCoord2D(struct face *face);

/* vef must be convex, pt must be interior */
/* returns -INFINITY on error */
/* start can be used to speed up queries for nearby points */
float Vef_ConvexInteriorDist(const struct vef *vef, const float *pt, struct face **start);

/* vef must be convex */
/* returns -INFINITY on error */
/* start can be used to speed up queries for nearby points */
float Vef_ConvexRayDist(const struct vef *vef, const float *pt, const float *dir, struct face **start);

#endif
