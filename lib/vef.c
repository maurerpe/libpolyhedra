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

#include "hash.h"
#include "queue.h"
#include "SipHash/siphash.h"
#include "util.h"
#include "vef.h"

#define PRESENT ((void *) 1)

static struct vert *Vert_New(const float *pt, struct vef *vef) {
  struct vert *vert;
  int count;
  
  if ((vert = Hash_Lookup(vef->verts, pt, NULL)))
    return vert;
  
  if ((vert = malloc(sizeof(*vert))) == NULL) {
    fprintf(stderr, "Error: Could not allocate memory for vertex\n");
    goto err;
  }
  memset(vert, 0, sizeof(*vert));
  vert->point[0] = pt[0];
  vert->point[1] = pt[1];
  vert->point[2] = pt[2];
  
  for (count = 0; count < 3; count++) {
    if (pt[count] < vef->min[count])
      vef->min[count] = pt[count];
    if (pt[count] > vef->max[count])
      vef->max[count] = pt[count];
  }
  
  if ((vert->edges = Hash_NewPtr(NULL, NULL, NULL, NULL, NULL)) == NULL)
    goto err2;
  
  if (Hash_Insert(vef->verts, vert->point, vert, NULL) < 0)
    goto err3;
  
  return vert;

 err3:
  Hash_Free(vert->edges);
 err2:
  free(vert);
 err:
  return NULL;
}

static void Vert_Free_Func(void *user, void *data) {
  struct vert *vert = (struct vert *) data;

  if (vert == NULL)
    return;

  Hash_Free(vert->edges);
  free(vert);
}

static struct edge *Edge_New(struct vert *v1, struct vert *v2, struct vef *vef) {
  struct edge *edge;
  
  if ((edge = Hash_Lookup(v1->edges, v2, NULL)))
    return edge;
  
  if ((edge = malloc(sizeof(*edge))) == NULL) {
    fprintf(stderr, "Error: Could not allocate memory for edge\n");
    goto err;
  }
  memset(edge, 0, sizeof(*edge));
  
  edge->vert[0] = v1;
  edge->vert[1] = v2;
  
  if (Hash_Insert(v1->edges, v2, edge, NULL) < 0)
    goto err2;

  if (Hash_Insert(v2->edges, v1, edge, NULL) < 0)
    goto err2;
  
  if (Hash_Insert(vef->edges, edge, PRESENT, NULL) < 0)
    goto err2;
  
  return edge;
  
 err2:
  free(edge);
 err:
  return NULL;
}

static void Edge_Free_Func(void *user, void *data) {
  free(data);
}

struct face *Face_New(const float *p1, const float *p2, const float *p3, struct vef *vef) {
  struct face *face;
  struct edge *edge;
  int count;

  if ((face = malloc(sizeof(*face))) == NULL) {
    fprintf(stderr, "Error: Could not allocate memory for face\n");
    goto err;
  }
  memset(face, 0, sizeof(*face));

  PlaneNorm(face->norm, p1, p2, p3);
  face->dist = Dot(face->norm, p1);

  if ((face->vert[0] = Vert_New(p1, vef)) == NULL)
    goto err2;
  if ((face->vert[1] = Vert_New(p2, vef)) == NULL)
    goto err2;
  if ((face->vert[2] = Vert_New(p3, vef)) == NULL)
    goto err2;
  
  for (count = 0; count < 3; count++) {
    if ((edge = face->edge[count] = Edge_New(face->vert[count], face->vert[(count + 1) % 3], vef)) == NULL)
      goto err2;
    
    edge->face[edge->face[0] == NULL ? 0 : 1] = face;
  }
  
  if (Hash_Insert(vef->faces, face, PRESENT, NULL) < 0)
    goto err2;
  
#ifdef DEBUG
  printf("Face w/ norm (%g,%g,%g)\n  (%g,%g,%g)\n  (%g, %g, %g)\n  (%g,%g,%g)\n\n",
	 face->norm[0],
	 face->norm[1],
	 face->norm[2],
	 face->vert[0]->point[0],
	 face->vert[0]->point[1],
	 face->vert[0]->point[2],
	 face->vert[1]->point[0],
	 face->vert[1]->point[1],
	 face->vert[1]->point[2],
	 face->vert[2]->point[0],
	 face->vert[2]->point[1],
	 face->vert[2]->point[2]);
#endif
  
  return face;
  
 err2:
  free(face);
 err:
  return NULL;
}

static void Face_Free_Func(void *user, void *data) {
  free(data);
}

struct face *Face_Adj(struct face *face, int count) {
  struct edge *edge = face->edge[count];
  
  return edge->face[edge->face[0] == face ? 1 : 0];
}

struct vef *Vef_New(const struct lp_vertex_list *vl) {
  struct vef *vef;
  size_t count, num;
  int mcount;
  
  if ((vef = malloc(sizeof(*vef))) == NULL) {
    fprintf(stderr, "Error: Could not allcoate memory for vef\n");
    goto err;
  }
  memset(vef, 0, sizeof(*vef));

  for (mcount = 0; mcount < 3; mcount++) {
    vef->min[mcount] =  INFINITY;
    vef->max[mcount] = -INFINITY;
  }
  
  if ((vef->verts = Hash_NewFixed(3 * sizeof(float), NULL, NULL, Vert_Free_Func, NULL)) == NULL)
    goto err2;
  
  if ((vef->edges = Hash_NewPtr(NULL, Edge_Free_Func, NULL, NULL, NULL)) == NULL)
    goto err3;
  
  if ((vef->faces = Hash_NewPtr(NULL, Face_Free_Func, NULL, NULL, NULL)) == NULL)
    goto err4;
  
  num = LP_VertexList_NumInd(vl);
  for (count = 0; count < num - 2; count += 3) {
    if (Face_New(LP_VertexList_LookupVert(vl, count),
		 LP_VertexList_LookupVert(vl, count + 1),
		 LP_VertexList_LookupVert(vl, count + 2),
		 vef) == NULL)
      goto err5;
  }
  
  return vef;

 err5:
  Hash_Free(vef->faces);
 err4:
  Hash_Free(vef->edges);
 err3:
  Hash_Free(vef->verts);
 err2:
  free(vef);
 err:
  return NULL;
}

void Vef_Free(struct vef *vef) {
  if (vef == NULL)
    return;
  
  Hash_Free(vef->faces);
  Hash_Free(vef->edges);
  Hash_Free(vef->verts);
  free(vef);
}

void Vef_CalcInfo(struct edge *edge) {
  float dx, dy, *y_vec, *norm;
  
  if (edge->info_vld)
    return;
  
  edge->z_vec[0] = edge->vert[1]->point[0] - edge->vert[0]->point[0];
  edge->z_vec[1] = edge->vert[1]->point[1] - edge->vert[0]->point[1];
  edge->z_vec[2] = edge->vert[1]->point[2] - edge->vert[0]->point[2];
  Normalize(edge->z_vec);
  
  y_vec = edge->face[0]->norm;
  Cross(edge->x_vec, y_vec, edge->z_vec);
  Normalize(edge->x_vec);
  
  norm = edge->face[1]->norm;
  dx = -Dot(norm, y_vec);
  dy =  Dot(norm, edge->x_vec);

  edge->ang = atan2(dy, dx);

  if (edge->ang < 0)
    edge->ang += 2 * M_PI;
  
  edge->info_vld = 1;
}

void Vef_CalcBasis(struct face *face) {
  struct vert *v0, *v1;
  
  if (face->basis_vld)
    return;
  
  v0 = face->vert[0];
  v1 = face->vert[1];
  
  face->basis_x[0] = v1->point[0] - v0->point[0];
  face->basis_x[1] = v1->point[1] - v0->point[1];
  face->basis_x[2] = v1->point[2] - v0->point[2];
  Normalize(face->basis_x);
  
  Cross(face->basis_y, face->norm, face->basis_x);
  Normalize(face->basis_y);
  
  face->basis_vld = 1;
}

void Vef_CalcCoord2D(struct face *face) {
  struct vert *v0, *v1, *v2;
  float d1[3], d2[3];

  if (face->coord_2d_vld)
    return;
  
  Vef_CalcBasis(face);
  
  v0 = face->vert[0];
  v1 = face->vert[1];
  v2 = face->vert[2];

  d1[0] = v1->point[0] - v0->point[0];
  d1[1] = v1->point[1] - v0->point[1];
  d1[2] = v1->point[2] - v0->point[2];
  
  d2[0] = v2->point[0] - v0->point[0];
  d2[1] = v2->point[1] - v0->point[1];
  d2[2] = v2->point[2] - v0->point[2];
  
  face->v1_x_len  = Dot(d1, face->basis_x);
  face->v2_pos[0] = Dot(d2, face->basis_x);
  face->v2_pos[1] = Dot(d2, face->basis_y);

#ifdef DEBUG
  if (fabsf(Dot(d1, face->basis_y)) > face->v1_x_len * 1e-4)
    printf("Warning: v1_x_len has y component\n");
  if (fabsf(Dot(d1, face->norm)) > face->v1_x_len * 1e-4)
    printf("Warning: v1_x_len has z component\n");
  if (fabsf(Dot(d2, face->norm)) > face->v1_x_len * 1e-4)
    printf("Warning: v2_pos has z component\n");
#endif
  
  face->coord_2d_vld = 1;
}

float Vef_ConvexInteriorDist(const struct vef *vef, const float *pt, struct face **start) {
  struct face *face, *adj_face, *min_face = NULL;
  struct hash *hash;
  struct queue *queue;
  float min = INFINITY, dist, tol;
  int count;
  
  tol = 1e-6 * Dist(vef->max, vef->min);
  
  if (Hash_Lookup(vef->verts, pt, NULL))
    return 0;
  
  if ((hash = Hash_NewPtr(NULL, NULL, NULL, NULL, NULL)) == NULL)
    goto err;
  if ((queue = Queue_New()) == NULL)
    goto err2;
  
  if (start && *start)
    face = *start;
  else if ((face = (struct face *) Hash_GetFirstKey(vef->faces)) == NULL)
    goto err3;
  
  if (Hash_Insert(hash, face, PRESENT, NULL) < 0)
    goto err3;
  if (Queue_PushBack(queue, face) < 0)
    goto err3;
  
  while (Queue_Length(queue) > 0) {
    if ((face = (struct face *) Queue_Pop(queue)) == NULL)
      goto err3;
    
    if ((dist = face->dist - Dot(face->norm, pt)) < -tol) {
      min = dist;
      min_face = face;
      break;
    }
    if (dist > min + tol)
      continue;
    if (dist < min) {
      min = dist;
      min_face = face;
    }
    
    for (count = 0; count < 3; count++) {
      adj_face = Face_Adj(face, count);
      if (Hash_Lookup(hash, adj_face, NULL))
	continue;
      if (Hash_Insert(hash, adj_face, PRESENT, NULL) < 0)
	goto err3;
      if (Queue_PushBack(queue, adj_face) < 0)
	goto err3;
    }
  }
  
  if (start)
    *start = min_face;
  
  Queue_Free(queue);
  Hash_Free(hash);
  return min;

 err3:
  Queue_Free(queue);
 err2:
  Hash_Free(hash);
 err:
  return -INFINITY;
}

static int Edge2D(const float *pt, float v1_x_len, const float *v2_pos, float tol) {
  float dist, max, norm[2], delta[2];
  int edge;
#ifdef DEBUG
  float d0, d1;
#endif
  
  if (v1_x_len <= 0)
    fprintf(stderr, "Warning: v1_x_len <= 0: %g\n", v1_x_len);
  if (v2_pos[1] <= 0)
    fprintf(stderr, "Warning: v2_pos[1] >= 0: %g\n", v2_pos[1]);
  
  max = -pt[1];
  edge = 0;
#ifdef DEBUG
  d0 = max;
#endif

  norm[0] = v2_pos[1];
  norm[1] = -(v2_pos[0] - v1_x_len);
  Normalize2d(norm);
  delta[0] = pt[0] - v1_x_len;
  delta[1] = pt[1];
  dist = Dot2d(delta, norm);
  if (dist > max) {
    max = dist;
    edge = 1;
  }
#ifdef DEBUG
  d1 = dist;
#endif
  
  norm[0] = -v2_pos[1];
  norm[1] = v2_pos[0];
  Normalize2d(norm);
  dist = Dot2d(pt, norm);
  if (dist > max) {
    max = dist;
    edge = 2;
  }

#ifdef DEBUG
  printf("Edge2D: Dist = (%g,%g,%g)\n", d0, d1, dist);
#endif
  
  if (max < tol)
    return 3;
  
  return edge;
}

static struct face *OtherFace(struct edge *edge, struct face *face) {
  if (edge->face[0] == face)
    return edge->face[1];
  
  return edge->face[0];
}

float Vef_ConvexRayDist(const struct vef *vef, const float *pt, const float *dir, struct face **start) {
  struct face *face;
  struct hash *hash;
  float tol, div, dist, *ref, pt3d[3], pt2d[2], com[2], scale;
  int edge;
  
#ifdef DEBUG
  printf("\n******************************************************************\nFinding dist of ray (%g,%g,%g) -> (%g,%g,%g)\n",
	 pt[0], pt[1], pt[2],
	 dir[0], dir[1], dir[2]);
#endif
  
  tol = 2e-6 * Dist(vef->max, vef->min);
  
  if ((hash = Hash_NewPtr(NULL, NULL, NULL, NULL, NULL)) == NULL)
    goto err;
  
  if (start && *start)
    face = *start;
  else if ((face = (struct face *) Hash_GetFirstKey(vef->faces)) == NULL)
    goto err2;
  
  while (1) {
#ifdef DEBUG
    float norm[3];
    PlaneNorm(norm, face->vert[0]->point, face->vert[1]->point, face->vert[2]->point);
    printf("\nTrying face w/ norm (%g,%g,%g) -> (%g,%g,%g):\n  (%g,%g,%g)\n  (%g,%g,%g)\n  (%g,%g,%g)\n",
	   face->norm[0],
	   face->norm[1],
	   face->norm[2],
	   norm[0],
	   norm[1],
	   norm[2],
	   face->vert[0]->point[0],
	   face->vert[0]->point[1],
	   face->vert[0]->point[2],
	   face->vert[1]->point[0],
	   face->vert[1]->point[1],
	   face->vert[1]->point[2],
	   face->vert[2]->point[0],
	   face->vert[2]->point[1],
	   face->vert[2]->point[2]);
#endif
    
    if (Hash_Lookup(hash, face, NULL)) {
      fprintf(stderr, "Error: Going around in circles finding ray distance\n");
      goto err2;
    }
    if (Hash_Insert(hash, face, PRESENT, NULL) < 0)
      goto err2;
    
    Vef_CalcCoord2D(face);
    com[0] = (face->v2_pos[0] + face->v1_x_len) / 3;
    com[1] =  face->v2_pos[1] / 3;
    scale = 2 * (Norm2d(face->v2_pos) + fabsf(face->v1_x_len));
    div = Dot(dir, face->norm);
    edge = -1;
    if (div < -0.5 || div >= 1e-6) {
      dist = (face->dist - Dot(pt, face->norm)) / div;
      ref = face->vert[0]->point;
      pt3d[0] = pt[0] + dir[0] * dist - ref[0];
      pt3d[1] = pt[1] + dir[1] * dist - ref[1];
      pt3d[2] = pt[2] + dir[2] * dist - ref[2];
      pt2d[0] = Dot(pt3d, face->basis_x);
      pt2d[1] = Dot(pt3d, face->basis_y);
      if (div < 0) {
#ifdef DEBUG
	printf("Raw pt2d: (%g,%g), com = (%g,%g), scale = %g\n",
	       pt2d[0], pt2d[1],
	       com[0], com[1],
	       scale);
#endif
	pt2d[0] -= com[0];
	pt2d[1] -= com[1];
	Normalize2d(pt2d);
	pt2d[0] *= -scale;
	pt2d[1] *= -scale;
	pt2d[0] += com[0];
	pt2d[1] += com[1];
      }
      edge = Edge2D(pt2d, face->v1_x_len, face->v2_pos, tol);
#ifdef DEBUG
      printf("Edge2D A: dist = %g, div = %g, pt2d = %g,%g, v1x = %g, v2pos = %g,%g, tol = %g => %d\n",
	     dist,
	     div,
	     pt2d[0],
	     pt2d[1],
	     face->v1_x_len,
	     face->v2_pos[0],
	     face->v2_pos[1],
	     tol,
	     edge);
#endif
      if (div > 0 && edge >= 3)
	break;
    }
    
    if (edge < 0) {
      pt2d[0] = Dot(dir, face->basis_x);
      pt2d[1] = Dot(dir, face->basis_y);
      Normalize2d(pt2d);
      pt2d[0] *= scale;
      pt2d[1] *= scale;
      pt2d[0] += com[0];
      pt2d[1] += com[1];
      edge = Edge2D(pt2d, face->v1_x_len, face->v2_pos, tol);
#ifdef DEBUG
      printf("Edge2D B: div = %g, pt2d = %g,%g, v1x = %g, v2pos = %g,%g, tol = %g => %d\n",
	     div,
	     pt2d[0],
	     pt2d[1],
	     face->v1_x_len,
	     face->v2_pos[0],
	     face->v2_pos[1],
	     tol,
	     edge);
#endif
    }
    
    if (edge < 0 || edge >= 3)
      goto err2;
    
    face = OtherFace(face->edge[edge], face);
  }
  
  if (start)
    *start = face;
  
#ifdef DEBUG
  printf("Found face: dist = %g\n", dist);
#endif
  Hash_Free(hash);
  return dist;

 err2:
  Hash_Free(hash);
 err:
  return -INFINITY;
}
