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

#include "array.h"
#include "ftree.h"
#include "hash.h"
#include "queue.h"
#include "util.h"

#define PRESENT ((void *) 1)
#define X       0
#define Y       1
#define LEFT    0
#define RIGHT   1

#ifdef DEBUG
size_t mono_poly_count = 0;
#endif

struct mono_poly;

struct vert {
  float point[2];
  struct hash *edges;
  struct ftree_node *node;
};

struct edge {
  struct vert *verts[2];
  struct mono_poly *mp;
};

struct mono_poly {
  struct queue *stack[2];
  struct vert *top;
  struct edge *active_edge[2];
  struct ftree_node *node;
  int top_side;
#ifdef DEBUG
  size_t count;
#endif
};

struct poly {
  struct ftree *mtree;
  struct array *edges;
  struct ftree *vtree;
  struct vert *verts;
  size_t num_verts;
};

#define HAS_CUSP(mp) ((mp)->stack[1] != NULL)

static size_t Vert_GetIdx(struct vert *vert, struct vert *verts) {
  return vert - verts;
}

static int Vert_Init(struct vert *vert, const float *point, struct ftree *vtree) {
  memset(vert, 0, sizeof(*vert));

  vert->point[X] = point[X];
  vert->point[Y] = point[Y];
  
  if ((vert->edges = Hash_NewPtr(NULL, NULL, NULL, NULL, NULL)) == NULL)
    goto err;
  
  if ((vert->node = FTree_Insert(vtree, point[Y], vert, NULL)) == NULL)
    goto err2;
  
  return 0;

 err2:
  Hash_Free(vert->edges);
 err:
  return -1;
}

static void Vert_Destroy(struct vert *vert) {
  if (vert == NULL)
    return;
  
  Hash_Free(vert->edges);
}

static int Edge_New(struct vert *v1, struct vert *v2, struct array *edges) {
  struct edge *edge;
  
  if (v1 == v2)
    /* Degenerate edge */
    return 0;
  
  if (Hash_Lookup(v1->edges, v2, NULL)) {
    /* Duplicate edge, remove */
    Hash_Remove(v1->edges, v2);
    Hash_Remove(v2->edges, v1);
    return 0;
  }
  
  if ((edge = malloc(sizeof(*edge))) == NULL)
    goto err;
  memset(edge, 0, sizeof(*edge));
  
  edge->verts[0] = v1;
  edge->verts[1] = v2;
  
  if (Hash_Insert(v1->edges, v2, edge, NULL) < 0)
    goto err2;
  if (Hash_Insert(v2->edges, v1, edge, NULL) < 0)
    goto err2;
  
  if (Array_Add(edges, edge) < 0)
    goto err2;
  
  return 0;

 err2:
  free(edge);
 err:
  return -1;
}

static void Edge_Free_Func(void *data) {
  free(data);
}

static float Edge_Ang(struct edge *edge, struct vert *ref) {
  float *p1, *p2, ang;
  
  p1 = ref->point;
  if (edge->verts[0] == ref) {
    p2 = edge->verts[1]->point;
  } else if (edge->verts[1] == ref) {
    p2 = edge->verts[0]->point;
  } else {
    fprintf(stderr, "Internal Error: triangulate2d.c: Could not find reference on edge\n");
    return 0;
  }
  
  ang = atan2f(p2[Y] - p1[Y], p2[X] - p1[X]);
  
  if (fabsf(M_PI - ang) < 1e-5 && edge->mp == NULL)
    return -ang;
  return ang;
}

static float Edge_Pos(struct edge *edge, float yy) {
  float *aa = edge->verts[0]->point, *bb = edge->verts[1]->point;
  
  if (aa[Y] == bb[Y])
    return 0.5 * (aa[X] + bb[X]);
  
  return (bb[X] - aa[X]) * (yy - aa[Y]) / (bb[Y] - aa[Y]) + aa[X];
}

static float MonoPoly_Key(struct mono_poly *mp, float *yy) {
  return Edge_Pos(mp->active_edge[LEFT], *yy);
}

static float MonoPoly_Key_Func(void *data, void *user) {
  return MonoPoly_Key((struct mono_poly *) data, (float *) user);
}

static int Edge_Orient(struct edge *edge, struct vert *top) {
  if (edge->verts[0] != top) {
    if (edge->verts[1] != top) {
      fprintf(stderr, "Internal Error: triangulate2d.c: Edge does not contain top vertex\n");
      return -1;
    }
    edge->verts[1] = edge->verts[0];
    edge->verts[0] = top;
  }
  
  return 0;
}

static struct mono_poly *MonoPoly_New(struct edge *left, struct edge *right, struct vert *start, struct ftree *mtree) {
  struct mono_poly *mp;
  int count;
  
  if ((mp = malloc(sizeof(*mp))) == NULL)
    goto err;
  memset(mp, 0, sizeof(*mp));
  
  mp->active_edge[0] = left;
  mp->active_edge[1] = right;
#ifdef DEBUG
  mp->count = mono_poly_count++;
  printf("Creating MP #%zu\n", mp->count);
#endif
  
  for (count = 0; count < 2; count++)
    if (Edge_Orient(mp->active_edge[count], start) < 0)
      goto err2;
  
  mp->top = start;
  
  if ((mp->node = FTree_Insert(mtree, 0, mp, &start->point[Y])) == NULL)
    goto err2;
  
  left->mp = mp;
  right->mp = mp;
  
  return mp;

 err2:
  free(mp);
 err:
  return NULL;
}

static void MonoPoly_Free(struct mono_poly *mp) {
  if (mp == NULL)
    return;

  Queue_Free(mp->stack[0]);
  Queue_Free(mp->stack[1]);
  free(mp);
}

static void MonoPoly_Free_Func(void *data) {
  MonoPoly_Free((struct mono_poly *) data);
}

static int MonoPoly_AddTriangle(struct lp_vertex_list *out, struct vert *p1, struct vert *p2, struct vert *p3, int is_opp) {
  float v1[2], v2[2], det, d1, d2, d3, temp, tol;
  
#ifdef DEBUG
  printf("Trying to add triangle: %g,%g %g,%g %g,%g... ",
	 p1->point[X],
	 p1->point[Y],
	 p2->point[X],
	 p2->point[Y],
	 p3->point[X],
	 p3->point[Y]);
  det = INFINITY;
  tol = 0;
#endif
  
  if (!is_opp) {
    v1[X] = p2->point[X] - p1->point[X];
    v1[Y] = p2->point[Y] - p1->point[Y];
    
    v2[X] = p3->point[X] - p2->point[X];
    v2[Y] = p3->point[Y] - p2->point[Y];

    det = v2[X] * v1[Y] - v2[Y] * v1[X];
    d1 = Dist2d2(p1->point, p2->point);
    d2 = Dist2d2(p1->point, p3->point);
    d3 = Dist2d2(p2->point, p3->point);
    if (d2 > d1) {
      temp = d1;
      d1 = d2;
      d2 = temp;
    }
    if (d3 > d2) {
      d2 = d3;
    }
    tol = 1e-6 * sqrtf(d1) * sqrtf(d2);
  
    if (det <= tol) {
#ifdef DEBUG
      printf("failure %g <= %g\n", det, tol);
#endif
      return 0;
    }
  }
  
  if (LP_VertexList_Add(out, p1->point) == UINT_MAX)
    return -1;
  if (LP_VertexList_Add(out, p2->point) == UINT_MAX)
    return -1;
  if (LP_VertexList_Add(out, p3->point) == UINT_MAX)
    return -1;
  
#ifdef DEBUG
  printf("success %g > %g\n", det, tol);
#endif
  return 1;
}

static int MonoPoly_AddVertSimple(struct lp_vertex_list *out, struct mono_poly *mp, struct vert *vert, int side) {
  struct vert *prev, *prev2, *hold;
  int ret;
  
  if (mp->stack[0] == NULL) {
    if ((mp->stack[0] = Queue_New()) == NULL)
      return -1;
    if (Queue_Push(mp->stack[0], mp->top) < 0)
      return -1;
    
    mp->top = vert;
    mp->top_side = side;
    return 0;
  }
  
  prev = mp->top;
  
  if (side == mp->top_side) {
    /* Same Side */
#ifdef DEBUG
    printf("Same side: %d\n", side);
#endif
    while (Queue_Length(mp->stack[0]) > 0) {
      prev2 = (struct vert *) Queue_Pop(mp->stack[0]);
      if (side == LEFT)
	ret = MonoPoly_AddTriangle(out, vert, prev, prev2, 0);
      else
	ret = MonoPoly_AddTriangle(out, vert, prev2, prev, 0);
      if (ret < 0)
	return -1;
      if (ret == 0) {
	if (Queue_Push(mp->stack[0], prev2) < 0)
	  return -1;
	break;
      }
      prev = prev2;
    }
  } else {
    /* Opposite Side */
#ifdef DEBUG
    printf("Opposite side: %d\n", side);
#endif
    hold = prev;
    while (Queue_Length(mp->stack[0]) > 0) {
      prev2 = (struct vert *) Queue_Pop(mp->stack[0]);
      if (side == LEFT)
	ret = MonoPoly_AddTriangle(out, vert, prev2, hold, 1);
      else
	ret = MonoPoly_AddTriangle(out, vert, hold, prev2, 1);
      if (ret < 0)
	return -1;
      if (ret == 0) {
	if (Queue_Push(mp->stack[0], prev2) < 0)
	  return -1;
	break;
      }
      hold = prev2;
    }
  }
  if (Queue_Push(mp->stack[0], prev) < 0)
    return -1;
  mp->top = vert;
  mp->top_side = side;
  
  return 0;
}

static int MonoPoly_AddVert(struct lp_vertex_list *out, struct mono_poly *mp, struct vert *vert, int side) {
  struct mono_poly temp_mp;
  
  if (HAS_CUSP(mp)) {
#ifdef DEBUG
    printf("MP #%zu: Resolving cusp\n", mp->count);
#endif
    memset(&temp_mp, 0, sizeof(temp_mp));
    temp_mp.stack[0] = mp->stack[side];
    temp_mp.top = mp->top;
    temp_mp.top_side = !side;
    if (MonoPoly_AddVertSimple(out, &temp_mp, vert, side) < 0)
      return -1;
    Queue_Free(temp_mp.stack[0]);
    if (side == LEFT)
      mp->stack[0] = mp->stack[1];
    mp->stack[1] = NULL;
    mp->top_side = side;
  }
  
  return MonoPoly_AddVertSimple(out, mp, vert, side);
}

static int MonoPoly_AdvEdge(struct lp_vertex_list *out, struct mono_poly *mp, struct edge *edge, struct vert *vert) {
  int side;
  
  if (Edge_Orient(edge, vert) < 0)
    return -1;
  
  if (mp->active_edge[0]->verts[1] == vert)
    side = 0;
  else if (mp->active_edge[1]->verts[1] == vert)
    side = 1;
  else {
    fprintf(stderr, "Internal Error: triangulate2d.c: Vertex not found when advancing edge\n");
    return -1;
  }
  
  mp->active_edge[side] = edge;
  edge->mp = mp;
  
#ifdef DEBUG
  printf("MP #%zu: AdvEdge: %g,%g\n", mp->count, vert->point[0], vert->point[1]);
#endif
  
  MonoPoly_AddVert(out, mp, vert, side);
  return 0;
}

static int MonoPoly_Merge(struct lp_vertex_list *out, struct ftree *mtree, struct mono_poly *left, struct mono_poly *right, struct vert *vert) {
  struct edge *edge;
  
#ifdef DEBUG
  printf("Merging: MP #%zu & #%zu @ %g,%g vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n", left->count, right->count, vert->point[0], vert->point[1]);
#endif
  
  if (left->active_edge[RIGHT]->verts[1] != vert ||
      right->active_edge[LEFT]->verts[1] != vert) {

    /* Check for swapped left and right */
    if (left->active_edge[LEFT]->verts[1] == vert &&
	right->active_edge[RIGHT]->verts[1] == vert) {
      fprintf(stderr, "Warning: swapped left and right in merge\n");
      return MonoPoly_Merge(out, mtree, right, left, vert);
    }
    
    /* Check for backwards left */
    if (left->active_edge[LEFT]->verts[1] == vert) {
      fprintf(stderr, "Warning: polynominal crossing detected\n");
      edge = left->active_edge[LEFT];
      left->active_edge[LEFT] = left->active_edge[RIGHT];
      left->active_edge[RIGHT] = edge;
    }
    
    /* Check for backwards right */
    if (right->active_edge[RIGHT]->verts[1] == vert) {
      fprintf(stderr, "Warning: polynominal crossing detected\n");
      edge = right->active_edge[LEFT];
      right->active_edge[LEFT] = right->active_edge[RIGHT];
      right->active_edge[RIGHT] = edge;
    }
    
    /* Check orginal test again */
    if (left->active_edge[RIGHT]->verts[1] != vert ||
	right->active_edge[LEFT]->verts[1] != vert) {
      fprintf(stderr, "Internal Error: triangulate2d.c: Incorrect vertex when merging\n");
      return -1;
    }
  }
  
  if (MonoPoly_AddVert(out, left, vert, RIGHT) < 0)
    return -1;
  if (MonoPoly_AddVert(out, right, vert, LEFT) < 0)
    return -1;
  left->stack[RIGHT] = right->stack[0];
  right->stack[0] = NULL;
  left->active_edge[RIGHT] = right->active_edge[RIGHT];
  left->active_edge[RIGHT]->mp = left;
  FTree_Delete(mtree, right->node);
  
  return 0;
}

static int MonoPoly_Split(struct lp_vertex_list *out, struct ftree *mtree, struct mono_poly *mp, struct mono_poly *mp_new) {
  struct edge *left, *right;
  struct vert *vert;
  
#ifdef DEBUG
  printf("Splitting: MP #%zu -> #%zu: %g,%g %d\n", mp->count, mp_new->count, mp_new->top->point[0], mp_new->top->point[1], mp->top_side);
#endif
  
  left  = mp_new->active_edge[LEFT];
  right = mp_new->active_edge[RIGHT];
  vert  = mp_new->top;

  mp_new->top = mp->top;
  mp_new->top_side = mp->top_side;
  mp_new->active_edge[RIGHT] = mp->active_edge[RIGHT];
  mp_new->active_edge[RIGHT]->mp = mp_new;
  mp_new->active_edge[LEFT] = right;
  mp->active_edge[RIGHT] = left;
  mp->active_edge[RIGHT]->mp = mp;
  
  if (HAS_CUSP(mp)) {
    mp_new->stack[0] = mp->stack[RIGHT];
    mp_new->top_side = LEFT;
    mp->stack[1] = NULL;
    mp->top_side = RIGHT;
  } else if (mp->top_side == LEFT) {
    mp_new->stack[0] = mp->stack[0];
    mp->stack[0] = NULL;
  }
  
  if (MonoPoly_AddVertSimple(out, mp_new, vert, LEFT) < 0)
    return -1;
  if (MonoPoly_AddVertSimple(out, mp, vert, RIGHT) < 0)
    return -1;
  
  return 0;
}

static int MonoPoly_NewSmart(struct lp_vertex_list *out, struct edge *left, struct edge *right, struct vert *start, struct ftree *mtree) {
  struct mono_poly *mp, *mp_new;
  struct ftree_node *prev;
  
  if ((mp_new = MonoPoly_New(left, right, start, mtree)) == NULL)
    goto err;
#ifdef DEBUG
  printf("NewSmart: %g,%g L=%g,%g R=%g,%g ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n",
	 start->point[X],
	 start->point[Y],
	 left->verts[1]->point[X],
	 left->verts[1]->point[Y],
	 right->verts[1]->point[X],
	 right->verts[1]->point[Y]);
#endif
  
  prev = FTree_Prev(mtree, mp_new->node);
  if (prev == NULL)
    return 0;
  mp = FTree_GetData(prev);
  if (Edge_Pos(mp->active_edge[RIGHT], start->point[Y]) > start->point[X])
    return MonoPoly_Split(out, mtree, mp, mp_new);
  
  return 0;
  
 err:
  return -1;
}

static int MonoPoly_Finish(struct lp_vertex_list *out, struct mono_poly *mp, struct vert *vert) {
#ifdef DEBUG
  printf("MP #%zu: Finish: %g,%g\n", mp->count, vert->point[0], vert->point[1]);
#endif
    
  if (HAS_CUSP(mp)) {
    if (MonoPoly_AddVert(out, mp, vert, RIGHT) < 0)
      return -1;
    mp->top = Queue_Pop(mp->stack[0]);
    mp->top_side = 1;
  }
  
  return MonoPoly_AddVertSimple(out, mp, vert, !mp->top_side);
}

static struct poly *Poly_New(size_t num_verts) {
  struct poly *poly;

  if ((poly = malloc(sizeof(*poly))) == NULL)
    goto err;

  poly->num_verts = num_verts;
  
  if ((poly->mtree = FTree_New(NULL, MonoPoly_Free_Func, MonoPoly_Key_Func)) == NULL)
    goto err2;
  
  if ((poly->edges = Array_New(num_verts + (num_verts >> 2), Edge_Free_Func)) == NULL)
    goto err3;
  
  if ((poly->vtree = FTree_New(NULL, NULL, NULL)) == NULL)
    goto err4;
  
  if ((poly->verts = calloc(num_verts, sizeof(*poly->verts))) == NULL)
    goto err5;
  
  return poly;

 err5:
  FTree_Free(poly->vtree);
 err4:
  Array_Free(poly->edges);
 err3:
  FTree_Free(poly->mtree);
 err2:
  free(poly);
 err:
  return NULL;
}

static void Poly_Free(struct poly *poly) {
  size_t count;
  
  if (poly == NULL)
    return;
  
  for (count = 0; count < poly->num_verts; count++)
    Vert_Destroy(&poly->verts[count]);
  free(poly->verts);
  FTree_Free(poly->vtree);
  Array_Free(poly->edges);
  FTree_Free(poly->mtree);
  free(poly);
}

static int Poly_Setup(struct poly *poly, const struct lp_vertex_list *in) {
  size_t count, num;
  float *data;
  unsigned int *ind;

  num = LP_VertexList_NumVert(in);
  data = LP_VertexList_GetVert(in);
  for (count = 0; count < num; count++)
    if (Vert_Init(&poly->verts[count], data + 2 * count, poly->vtree) < 0)
      return -1;
  
  num = LP_VertexList_NumInd(in) & -2;
  ind = LP_VertexList_GetInd(in);
  for (count = 0; count < num; count += 2)
    if (Edge_New(&poly->verts[ind[count]], &poly->verts[ind[count + 1]], poly->edges) < 0)
      return -1;
  
  return 0;
}

static int Poly_Triangulate(struct lp_vertex_list *out, struct poly *poly) {
  struct ftree_node *node;
  struct vert *vert;
  struct edge *edge;
  size_t num_edges;
  struct hash_iterator *hi;
  struct ftree *top, *bot;
  struct ftree_node *top_node, *next_node, *bot_node;
  float ang;
  
  if ((top = FTree_New(NULL, NULL, NULL)) == NULL)
    goto err;
  if ((bot = FTree_New(NULL, NULL, NULL)) == NULL)
    goto err2;
  
  for (node = FTree_Highest(poly->vtree); node; node = FTree_Prev(poly->vtree, node)) {
    vert = FTree_GetData(node);
    num_edges = Hash_NumEntries(vert->edges);
#ifdef DEBUG
    printf("\nProcessing vert %g,%g with %zu edges\n", vert->point[X], vert->point[Y], num_edges);
#endif
    
    if (num_edges == 0)
      continue;
    
    if (num_edges & 1) {
      fprintf(stderr, "Error: Vertex %zu has odd number of edges: %zu\n", Vert_GetIdx(vert, poly->verts), num_edges);
      goto err3;
    }
    
    if ((hi = Hash_IteratorNew(vert->edges)) == NULL)
      goto err3;
    while (Hash_IteratorNext(hi)) {
      edge = Hash_IteratorGetData(hi);
      ang = Edge_Ang(edge, vert);
      if (FTree_Insert(edge->mp ? top : bot, ang, edge, NULL) == NULL)
	goto err4;
    }
    Hash_IteratorFree(hi);
    
    bot_node = FTree_Lowest(bot);
    for (top_node = FTree_Highest(top); top_node; top_node = next_node) {
      next_node = FTree_Prev(top, top_node);
      edge = (struct edge *) FTree_GetData(top_node);
      if (next_node && ((struct edge *) FTree_GetData(next_node))->mp == edge->mp) {
	if (MonoPoly_Finish(out, edge->mp, vert) < 0)
	  goto err3;
	FTree_Delete(poly->mtree, edge->mp->node);
	next_node = FTree_Prev(top, next_node);
	continue;
      }
      
      if (bot_node) {
	if (MonoPoly_AdvEdge(out, edge->mp, (struct edge *) FTree_GetData(bot_node), vert) < 0)
	  goto err3;
	bot_node = FTree_Next(bot, bot_node);
	continue;
      }
      
      if (MonoPoly_Merge(out, poly->mtree, edge->mp, ((struct edge *) FTree_GetData(next_node))->mp, vert) < 0)
	goto err3;
      next_node = FTree_Prev(top, next_node);
    }
    
    for (; bot_node; bot_node = FTree_Next(bot, next_node)) {
      next_node = FTree_Next(bot, bot_node);
      edge = (struct edge *) FTree_GetData(bot_node);
      if (MonoPoly_NewSmart(out, edge, (struct edge *) FTree_GetData(next_node), vert, poly->mtree) < 0)
	goto err3;
    }

    FTree_Clear(top);
    FTree_Clear(bot);
  }
  
  FTree_Free(bot);
  FTree_Free(top);
  return 0;

 err4:
  Hash_IteratorFree(hi);
 err3:
  FTree_Free(bot);
 err2:
  FTree_Free(top);
 err:
  return -1;
}

struct lp_vertex_list *LP_Triangulate2D(const struct lp_vertex_list *in) {
  struct lp_vertex_list *out;
  struct poly *poly;
  size_t num_verts;
  
  if (LP_VertexList_FloatsPerVert(in) != 2) {
    fprintf(stderr, "Error: Incorrect number of floats per vert for triangulate 2D\n");
    goto err;
  }
  
  if (LP_VertexList_PrimativeType(in) != lp_pt_line) {
    fprintf(stderr, "Error: Incorrect primative type for triangulate 2D\n");
    goto err;
  }

#ifdef DEBUG
  mono_poly_count = 0;
  printf("Triangulating %zu edges\n", LP_VertexList_NumInd(in) / 2);
#endif
  
  num_verts = LP_VertexList_NumVert(in);
  if ((poly = Poly_New(num_verts)) == NULL)
    goto err;
  
  if (Poly_Setup(poly, in) < 0)
    goto err2;
  
  if ((out = LP_VertexList_New(2, lp_pt_triangle)) == NULL)
    goto err2;
  
  if (Poly_Triangulate(out, poly) < 0)
    goto err3;

#ifdef DEBUG
  printf("Returning with %zu triangles\n", LP_VertexList_NumInd(out) / 3);
#endif
  Poly_Free(poly);
  return out;

 err3:
  LP_VertexList_Free(out);
 err2:
  Poly_Free(poly);
 err:
  return NULL;
}
