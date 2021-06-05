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

#include "hash.h"
#include "queue.h"
#include "SipHash/siphash.h"
#include "util.h"

#define PRESENT ((void *) 1)

struct plane {
  float norm[3];
  float x_axis[3];
  float y_axis[3];
  float dist;
};

struct face;

struct vert {
  float point[3];
  struct hash *edges;
  float dist;
};

struct edge {
  struct vert *vert[2];
  struct face *face[2];
  float inter[3];
  int inter_vld;
};

struct face {
  struct vert *vert[3];
  struct edge *edge[3];
  int visited;
};

struct shape {
  struct hash *verts;
  struct hash *edges;
  struct hash *faces;
  struct hash *pt2d;
  struct hash *edge2d;
  struct lp_vertex_list *poly2d;
};

static struct vert *Vert_New(const float *point, const struct plane *plane, struct shape *shape) {
  struct vert *vert;
  float tol;
  
  if ((vert = Hash_Lookup(shape->verts, point, NULL)))
    return vert;
  
  if ((vert = malloc(sizeof(*vert))) == NULL) {
    fprintf(stderr, "Could not allocate vertex for plane cut");
    goto err;
  }
  memset(vert, 0, sizeof(*vert));
  
  vert->point[0] = point[0];
  vert->point[1] = point[1];
  vert->point[2] = point[2];

  if (plane) {
    tol = Norm(point);
    if (fabsf(plane->dist) > tol)
      tol = fabsf(plane->dist);
    tol *= 1e-5;
    
    vert->dist = Dot(point, plane->norm) - plane->dist;
    if (fabsf(vert->dist) < tol)
      vert->dist = 0;
  }
  
  if ((vert->edges = Hash_NewPtr(NULL, NULL, NULL, NULL, NULL)) == NULL)
    goto err2;
  
  if (Hash_Insert(shape->verts, point, vert, NULL) < 0)
    goto err3;

  return vert;

 err3:
  Hash_Free(vert->edges);
 err2:
  free(vert);
 err:
  fprintf(stderr, "Error: Could not create vertex\n");
  return NULL;
}

static void Vert_Free(struct vert *vert) {
  if (vert == NULL)
    return;

  Hash_Free(vert->edges);
  free(vert);
}

static void Vert_Free_Func(void *user, void *data) {
  Vert_Free((struct vert *) data);
}

static struct edge *Edge_New(struct vert *v1, struct vert *v2, const struct plane *plane, struct shape *shape) {
  struct edge *edge;
  float x, y;
  
  if ((edge = Hash_Lookup(v1->edges, v2, NULL)))
    return edge;
  
  if ((edge = malloc(sizeof(*edge))) == NULL) {
    perror("Could not allocate edge for plane cut");
    goto err;
  }
  memset(edge, 0, sizeof(*edge));
  
  edge->vert[0] = v1;
  edge->vert[1] = v2;
  
  if (Hash_Insert(v1->edges, v2, edge, NULL) < 0)
    goto err2;
  if (Hash_Insert(v2->edges, v1, edge, NULL) < 0)
    goto err2;
  
  if (plane &&
      ((v1->dist > 0 && v2->dist < 0) ||
       (v1->dist < 0 && v2->dist > 0))) {
    /* Intersects plane */
    x = -v1->dist / (v2->dist - v1->dist);
    y = 1 - x;
    edge->inter[0] = y * v1->point[0] + x * v2->point[0];
    edge->inter[1] = y * v1->point[1] + x * v2->point[1];
    edge->inter[2] = y * v1->point[2] + x * v2->point[2];
    edge->inter_vld = 1;
  }
  
  if (Hash_Insert(shape->edges, edge, PRESENT, NULL) < 0)
    goto err2;
  
  return edge;
  
 err2:
  free(edge);
 err:
  fprintf(stderr, "Error: Could not create edge\n");
  return NULL;
}

static void Edge_Free_Func(void *user, void *data) {
  free(data);
}

static struct face *Face_New(float *p1, float *p2, float *p3, struct shape *shape) {
  struct face *face;
  struct edge *edge;
  float *pt[3];
  int count;
  
  pt[0] = p1;
  pt[1] = p2;
  pt[2] = p3;
  
  if ((face = malloc(sizeof(*face))) == NULL) {
    perror("Error: Could not allocate face for plane cut");
    goto err;
  }
  memset(face, 0, sizeof(*face));

  for (count = 0; count < 3; count++)
    if ((face->vert[count] = Vert_New(pt[count], NULL, shape)) == NULL)
      goto err2;

  for (count = 0; count < 3; count++) {
    if ((edge = face->edge[count] = Edge_New(face->vert[count],
					     face->vert[(count + 1) % 3],
					     NULL,
					     shape)) == NULL)
      goto err2;
    
    edge->face[edge->face[0] == NULL ? 0 : 1] = face;
  }
  
  if (Hash_Insert(shape->faces, face, PRESENT, NULL) < 0)
    goto err2;
  
  return face;
  
 err2:
  free(face);
 err:
  fprintf(stderr, "Error: Could not create face\n");
  return NULL;
}

static void Face_Free_Func(void *user, void *data) {
  free(data);
}

static int Make_Quad(float *p1, float *p2, float *p3, float *p4, struct shape *shape) {
  /* Split into triangles along shortest diagonal */
  if (Dist2(p1, p3) > Dist2(p2, p4)) {
    if (Face_New(p2, p3, p4, shape) == NULL)
      return -1;
    if (Face_New(p1, p2, p4, shape) == NULL)
      return -1;
  } else {
    if (Face_New(p1, p3, p4, shape) == NULL)
      return -1;
    if (Face_New(p1, p2, p3, shape) == NULL)
      return -1;
  }
  
  return 0;
}

static int Add2dPoint(float *pt, const struct plane *plane, struct shape *shape) {
  struct vert *v;
  float ff[2];
  
  if ((v = Vert_New(pt, NULL, shape)) == NULL)
    return -1;
  
  ff[0] = Dot(pt, plane->x_axis);
  ff[1] = Dot(pt, plane->y_axis);  
  if (LP_VertexList_Add(shape->poly2d, ff) == UINT_MAX)
    return -1;
  if (Hash_Insert(shape->pt2d, ff, v, NULL) < 0)
    return -1;
  
  return 0;
}

static int Make_Faces(float *p1, float *p2, float *p3, const struct plane *plane, struct shape **shape) {
  struct vert *v[3];
  struct edge *e[3], *edge;
  struct face *face;
  struct shape *ss;
  float *pt[3];
  int count, non1, non2, i1, i2;
  
  pt[0] = p1;
  pt[1] = p2;
  pt[2] = p3;
  
  for (count = 0; count < 3; count++)
    if ((v[count] = Vert_New(pt[count], plane, shape[2])) == NULL)
      return -1;
  
  for (count = 0; count < 3; count++)
    if ((e[count] = Edge_New(v[count], v[(count + 1) % 3], plane, shape[2])) == NULL)
      return -1;
  
  switch (e[0]->inter_vld +
	  e[1]->inter_vld +
	  e[2]->inter_vld) {
  case 0:
    /* Doesn't intersect plane */
    switch ((v[0]->dist == 0) +
	    (v[1]->dist == 0) +
	    (v[2]->dist == 0)) {
    case 0:
    case 1:
      /* At most 1 point on plane */
      if (v[0]->dist != 0)
	non1 = 0;
      else
	non1 = 1;
      
      if (Face_New(p1, p2, p3, shape[v[non1]->dist > 0]) == NULL)
	return -1;
      break;
      
    case 2:
      /* Two points lie on plane: exactly one edge on plane */
      if (v[0]->dist != 0)
	non1 = 0;
      else if (v[1]->dist != 0)
	non1 = 1;
      else
	non1 = 2;
      
      i1 = (non1 + 1) % 3;
      i2 = (non1 + 2) % 3;
      ss = shape[v[non1]->dist > 0];
      
      if ((face = Face_New(p1, p2, p3, ss)) == NULL)
	return -1;
      
      edge = face->edge[i1];
      if (Hash_Lookup(ss->edge2d, edge, NULL)) {
	Hash_Remove(ss->edge2d, edge);
      } else {
	if (Hash_Insert(ss->edge2d, edge, PRESENT, NULL) < 0)
	  return -1;
      }
      break;
      
    case 3:
      /* All points on plane, ignore face */
      break;
    }
    break;
    
  case 1:
    /* One edge intersects plane */
    /* One point exactly on plane, other two straddle plane */
    if (e[0]->inter_vld)
      i1 = 0;
    else if (e[1]->inter_vld)
      i1 = 1;
    else
      i1 = 2;
    
    non1 = (i1 + 1) % 3;
    non2 = (i1 + 2) % 3;
    
    if (v[non2]->dist != 0) {
      fprintf(stderr, "Internal Error: plane_cut.c: Expected point to be on plane\n");
      return -1;
    }
    
    if (Add2dPoint(e[i1]->inter, plane, shape[0]) < 0)
      return -1;
    if (Add2dPoint(e[i1]->inter, plane, shape[1]) < 0)
      return -1;
    if (Add2dPoint(v[non2]->point, plane, shape[0]) < 0)
      return -1;
    if (Add2dPoint(v[non2]->point, plane, shape[1]) < 0)
      return -1;
    
    if (Face_New(v[non1]->point, v[non2]->point, e[i1]->inter, shape[v[non1]->dist > 0]) == NULL)
      return -1;
    if (Face_New(v[non2]->point, v[i1]->point,   e[i1]->inter, shape[v[i1]->dist > 0]) == NULL)
      return -1;
    break;

  case 2:
    /* Two edges intersect plane */
    /* Two points on one side, one on the other */
    if (!e[0]->inter_vld)
      non1 = 0;
    else if (!e[1]->inter_vld)
      non1 = 1;
    else
      non1 = 2;
  
    i1 = (non1 + 1) % 3;
    i2 = (non1 + 2) % 3;
    
    if (Add2dPoint(e[i1]->inter, plane, shape[0]) < 0)
      return -1;
    if (Add2dPoint(e[i1]->inter, plane, shape[1]) < 0)
      return -1;
    if (Add2dPoint(e[i2]->inter, plane, shape[0]) < 0)
      return -1;
    if (Add2dPoint(e[i2]->inter, plane, shape[1]) < 0)
      return -1;
    
    if (Face_New(v[i2]->point, e[i2]->inter, e[i1]->inter, shape[v[i2]->dist > 0]) == NULL)
      return -1;
    if (Make_Quad(v[non1]->point, v[i1]->point, e[i1]->inter, e[i2]->inter, shape[v[i1]->dist > 0]) < 0)
      return -1;
    break;

  default:
    fprintf(stderr, "Internal Error: plane_cut.c: Invalid number of edges intersects plane\n");
    return -1;
  }  
  
  return 0;
}

static int AddEdge2d(struct shape *shape, const struct plane *plane) {
  struct hash_iterator *hi;
  struct edge *edge;

  if ((hi = Hash_IteratorNew(shape->edge2d)) == NULL)
    goto err;
  while (Hash_IteratorNext(hi)) {
    edge = (struct edge *) Hash_IteratorGetKey(hi);
    if (Add2dPoint(edge->vert[0]->point, plane, shape) < 0)
      goto err2;
    if (Add2dPoint(edge->vert[1]->point, plane, shape) < 0)
      goto err2;
  }
  Hash_IteratorFree(hi);
  
  return 0;
  
 err2:
  Hash_IteratorFree(hi);
 err:
  return -1;
}

static int Make_Face_From_2d(const float *p1, const float *p2, const float *p3, struct shape *shape, int sense) {
  struct vert *v1, *v2, *v3;
  
  if ((v1 = (struct vert *) Hash_Lookup(shape->pt2d, p1, NULL)) == NULL ||
      (v2 = (struct vert *) Hash_Lookup(shape->pt2d, p2, NULL)) == NULL ||
      (v3 = (struct vert *) Hash_Lookup(shape->pt2d, p3, NULL)) == NULL) {
    fprintf(stderr, "Error: Unexpected 2d point when slicing polyhedron\n");
    return -1;
  }
  
  if (sense) {
    if (Face_New(v1->point, v2->point, v3->point, shape) == NULL)
      return -1;
  } else {
    if (Face_New(v1->point, v3->point, v2->point, shape) == NULL)
      return -1;
  }
  
  return 0;
}

static struct face *Get_Next(struct face *face, struct edge *edge) {
  return edge->face[edge->face[0] == face ? 1 : 0];
}

static int Build_Poly3d(struct lp_vertex_list *poly3d, struct face *face) {
  struct queue *queue;
  struct face *next;
  int count;
  
  if ((queue = Queue_New()) == NULL)
    goto err;
  
  face->visited = 1;
  if (Queue_Push(queue, face) < 0)
    goto err2;
  
  while (Queue_Length(queue) > 0) {
    face = (struct face *) Queue_Pop(queue);
    
    for (count = 0; count < 3; count++) {
      if (LP_VertexList_Add(poly3d, face->vert[count]->point) < 0)
	goto err2;
      
      if ((next = Get_Next(face, face->edge[count])) == NULL) {
	fprintf(stderr, "Warning: Could not find adjacent face\n");
	continue;
      }
      if (next->visited)
	continue;
      next->visited = 1;
      if (Queue_Push(queue, next) < 0)
	goto err2;
    }
  }
  
  Queue_Free(queue);
  return 0;

 err2:
  Queue_Free(queue);
 err:
  fprintf(stderr, "Error: Could not build 3d polyhedron\n");
  return -1;
}

static struct shape *Shape_New(void) {
  struct shape *shape;

  if ((shape = malloc(sizeof(*shape))) == NULL) {
    fprintf(stderr, "Could not allocate shape for plane cut");
    goto err;
  }
  
  if ((shape->verts = Hash_NewFixed(3 * sizeof(float), NULL, NULL, Vert_Free_Func, NULL)) == NULL)
    goto err2;
  
  if ((shape->edges = Hash_NewPtr(NULL, Edge_Free_Func, NULL, NULL, NULL)) == NULL)
    goto err3;

  if ((shape->faces = Hash_NewPtr(NULL, Face_Free_Func, NULL, NULL, NULL)) == NULL)
    goto err4;
  
  if ((shape->pt2d = Hash_NewFixed(2 * sizeof(float), NULL, NULL, NULL, NULL)) == NULL)
    goto err5;

  if ((shape->edge2d = Hash_NewPtr(NULL, NULL, NULL, NULL, NULL)) == NULL)
    goto err6;
  
  if ((shape->poly2d = LP_VertexList_New(2, lp_pt_line)) == NULL)
    goto err7;
  
  return shape;

 err7:
  Hash_Free(shape->edge2d);
 err6:
  Hash_Free(shape->pt2d);
 err5:
  Hash_Free(shape->faces);
 err4:
  Hash_Free(shape->edges);
 err3:
  Hash_Free(shape->verts);
 err2:
  free(shape);
 err:
  return NULL;
}

static void Shape_Free(struct shape *shape) {
  if (shape == NULL)
    return;
  
  LP_VertexList_Free(shape->poly2d);
  Hash_Free(shape->edge2d);
  Hash_Free(shape->pt2d);
  Hash_Free(shape->faces);
  Hash_Free(shape->edges);
  Hash_Free(shape->verts);
  free(shape);
}

struct lp_vl_list *LP_PlaneCut(const struct lp_vertex_list *in, const float *norm, float dist) {
  struct plane plane;
  struct shape *shape[3];
  struct lp_vertex_list *tri, *poly3d;
  struct lp_vl_list *out = NULL, *temp;
  struct hash_iterator *hi;
  struct face *face;
  size_t count, num;
  int s_count;
  int valid = 1;
  
  plane.norm[0] = norm[0];
  plane.norm[1] = norm[1];
  plane.norm[2] = norm[2];
  plane.dist    = dist;
  Normalize(plane.norm);
  
  BasisVectors(plane.x_axis, plane.y_axis, plane.norm);
  
  if (LP_VertexList_FloatsPerVert(in) < 3) {
    fprintf(stderr, "Error: Insufficent floats per vert for plane cut: %zu\n", LP_VertexList_FloatsPerVert(in));
    goto err;
  }
  
  if (LP_VertexList_PrimativeType(in) != lp_pt_triangle) {
    fprintf(stderr, "Error: Can only plane cut triangular shapes\n");
    goto err;
  }
  
  /* Shape 0 is side 0 (dist < 0) */
  if ((shape[0] = Shape_New()) == NULL)
    goto err;
  /* Shape 1 is side 1 (dist > 0) */
  if ((shape[1] = Shape_New()) == NULL)
    goto err2;
  /* Shape 2 holds orginal vertices and edges for determining
   * sides and intersections
   */
  if ((shape[2] = Shape_New()) == NULL)
    goto err3;

  num = LP_VertexList_NumInd(in);
#ifdef DEBUG
  printf("Making %zu faces\n", num / 3);
#endif
  for (count = 0; count < num; count += 3) {
    if (Make_Faces(LP_VertexList_LookupVert(in, count),
		   LP_VertexList_LookupVert(in, count + 1),
		   LP_VertexList_LookupVert(in, count + 2),
		   &plane, shape) < 0)
      goto err4;
  }
  
  for (s_count = 0; s_count < 2; s_count++)
    if (AddEdge2d(shape[s_count], &plane) < 0)
      goto err4;

#ifdef DEBUG
  struct lp_vl_list list;
  list.vl = shape[0]->poly2d;
  list.next = NULL;
  LP_VertexList_Write("test0.svg", &list, 1000);
  list.vl = shape[1]->poly2d;
  list.next = NULL;
  LP_VertexList_Write("test1.svg", &list, 1000);
#endif
  
  for (s_count = 0; s_count < 2; s_count++) {
#ifdef DEBUG
    printf("Triangulating\n");
#endif
    if ((tri = LP_Triangulate2D(shape[s_count]->poly2d)) == NULL) {
      valid = 0;
    } else {
      num = LP_VertexList_NumInd(tri);
#ifdef DEBUG
      printf("Adding %zu triangles\n", num / 3);
      list.vl = tri;
      list.next = NULL;
      LP_VertexList_Write(s_count == 0 ? "tri0.svg" : "tri1.svg", &list, 1000);
#endif
      for (count = 0; count < num; count += 3) {
	if (Make_Face_From_2d(LP_VertexList_LookupVert(tri, count),
			      LP_VertexList_LookupVert(tri, count + 1),
			      LP_VertexList_LookupVert(tri, count + 2),
			      shape[s_count],
			      s_count) < 0)
	  goto err5;
      }
    }
    
    if ((hi = Hash_IteratorNew(shape[s_count]->faces)) == NULL)
      goto err5;
    while (Hash_IteratorNext(hi)) {
      face = Hash_IteratorGetKey(hi);
      if (face->visited)
	continue;
      if ((poly3d = LP_VertexList_New(3, lp_pt_triangle)) == NULL)
	goto err6;
      if (Build_Poly3d(poly3d, face) < 0)
	goto err7;
      if ((temp = LP_VertexList_ListAppend(out, poly3d)) == NULL)
	goto err7;
      out = temp;
    }
    Hash_IteratorFree(hi);
    LP_VertexList_Free(tri);
  }
  
  if (!valid) {
#ifdef DEBUG
    LP_VertexList_Write("plane_cut.obj", out, 1);
#endif
    goto err4;
  }
  
  Shape_Free(shape[2]);
  Shape_Free(shape[1]);
  Shape_Free(shape[0]);
#ifdef DEBUG
  printf("Returning with %zu polyhedrons\n", LP_VertexList_ListLength(out));
#endif
  return out;
  
 err7:
  LP_VertexList_Free(poly3d);
 err6:
  Hash_IteratorFree(hi);
 err5:
  LP_VertexList_Free(tri);
 err4:
  Shape_Free(shape[2]);
 err3:
  Shape_Free(shape[1]);
 err2:
  Shape_Free(shape[0]);
 err:
  LP_VertexList_ListFree(out);
  fprintf(stderr, "Error: Could not cut polyhedron with a plane\n");
  return NULL;
}
