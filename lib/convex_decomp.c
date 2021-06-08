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

#include "ftree.h"
#include "queue.h"
#include "util.h"
#include "vef.h"

/* Based on:
 * Thul, Daniel et ala
 * Approximate Convex Decomposition and Transfer for Animated Meshes
 * 2018
 */

#define PRESENT ((void *) 1)

#define NUM_EDGES  16
#define NUM_ANGLES 9

struct vlh_list {
  struct lp_vertex_list *vl;
  struct lp_vertex_list *hull;
  struct vlh_list *next;
  float err;
};

static float ConvexError(const struct lp_vertex_list *vl, const struct lp_vertex_list *hull) {
  struct lp_mass_properties mp, mpc;
  
  LP_MassProperties(vl,   &mp);
  LP_MassProperties(hull, &mpc);
  
  return mpc.volume - mp.volume;
}

static struct vlh_list *VhlList_New(struct lp_vertex_list *vl) {
  struct vlh_list *vlh;

  if ((vlh = malloc(sizeof(*vlh))) == NULL) {
    fprintf(stderr, "Could not allocate memeory for vlh list\n");
    goto err;
  }
  memset(vlh, 0, sizeof(*vlh));

  vlh->vl = vl;

  if ((vlh->hull = LP_ConvexHull(vl)) == NULL)
    goto err2;
  
  vlh->err = ConvexError(vl, vlh->hull);
  
  return vlh;
  
 err2:
  free(vlh);
 err:
  return NULL;
}

static void VlhList_Free(struct vlh_list *vlh) {
  struct vlh_list *next;
  
  while (vlh != NULL) {
    LP_VertexList_Free(vlh->vl);
    LP_VertexList_Free(vlh->hull);
    next = vlh->next;
    free(vlh);
    vlh = next;
  }
}

static float VlhList_TotalError(struct vlh_list *vlh) {
  float err = 0;

  while (vlh != NULL) {
    err += vlh->err;
    vlh = vlh->next;
  }
  
  return err;
}

static float VlhList_TotalSqrError(struct vlh_list *vlh) {
  float err = 0;

  while (vlh != NULL) {
    err += vlh->err * vlh->err;
    vlh = vlh->next;
  }
  
  return err;
}

static size_t VlhList_Len(struct vlh_list *vlh) {
  size_t len = 0;

  while (vlh != NULL) {
    len++;
    vlh = vlh->next;
  }
  
  return len;
}

static struct vlh_list *VlhList_Convert(struct lp_vl_list *list, float *err) {
  struct vlh_list *head = NULL, **tail = &head;
  struct lp_vl_list *cur = list;
  
  if (err)
    *err = 0;
  
  while (cur) {
    if (LP_VertexList_NumVert(cur->vl) > 4) {
      if ((*tail = VhlList_New(cur->vl)) == NULL)
	goto err;
      cur->vl = NULL;
      if (err)
	*err += (*tail)->err;
      tail = &(*tail)->next;
    } else {
      fprintf(stderr, "Warning: only %u points in polyhedron, skipping\n",
	      LP_VertexList_NumVert(cur->vl));
    }
    
    cur = cur->next;
  }
  
  LP_VertexList_ListFree(list);
  return head;
  
 err:
  VlhList_Free(head);
  LP_VertexList_ListFree(list);
  return NULL;
}

static struct lp_vl_list *Vl_Convert(struct vlh_list *vlh, int hull, int free) {
  struct lp_vl_list *head = NULL, **tail = &head;
  struct vlh_list *cur = vlh;
  
  while (cur) {
    if ((*tail = LP_VertexList_ListAppend(NULL, hull ? cur->hull : cur->vl)) == NULL)
      goto err;
    if (free) {
      if (hull)
	cur->hull = NULL;
      else
	cur->vl = NULL;
    }
    cur = cur->next;
    tail = &(*tail)->next;
  }
  
  if (free)
    VlhList_Free(vlh);
  return head;
  
 err:
  LP_VertexList_ListFree(head);
  if (free)
    VlhList_Free(vlh);
  return NULL;
}

static struct vlh_list **WorstPart(struct vlh_list **list) {
  float worst_err = -INFINITY, cur_err;
  struct vlh_list **worst = NULL;
  
  while (*list != NULL) {
    cur_err = (*list)->err;
    if (cur_err > worst_err) {
      worst_err = cur_err;
      worst = list;
    }

    list = &(*list)->next;
  }
  
  return worst;
}

static struct ftree *FurthestEdges(struct vef *full, struct vef *hull) {
  struct ftree *ftree;
  struct hash *hash;
  struct queue *queue;
  struct edge *edge, *ee;
  struct face *start = NULL;
  struct hash_iterator *hi;
  struct lp_transform *trans;
  float mid[3], dir[3], dist;
  int count;
  
  if ((ftree = FTree_New(NULL, NULL, NULL)) == NULL)
    goto err;
  if ((hash = Hash_NewPtr(NULL, NULL, NULL, NULL, NULL)) == NULL)
    goto err2;
  if ((queue = Queue_New()) == NULL)
    goto err3;
  if ((trans = LP_Transform_New()) == NULL)
    goto err4;
  
  if ((edge = Hash_GetFirstKey(full->edges)) == NULL)
    goto err5;
  if (Queue_PushBack(queue, edge) < 0)
    goto err5;
  if (Hash_Insert(hash, edge, PRESENT, NULL) < 0)
    goto err5;
  
  while (Queue_Length(queue) > 0) {
    if ((edge = (struct edge *) Queue_Pop(queue)) == NULL)
      goto err5;
    
    Vef_CalcInfo(edge);
    
#ifdef DEBUG_EDGE
    printf("Finding distance of edge: (%g,%g,%g) - (%g,%g,%g): %g deg\n",
	   edge->vert[0]->point[0],
	   edge->vert[0]->point[1],
	   edge->vert[0]->point[2],
	   edge->vert[1]->point[0],
	   edge->vert[1]->point[1],
	   edge->vert[1]->point[2],
	   edge->ang * 180 / M_PI);
#endif
    
    mid[0] = 0.5 * (edge->vert[0]->point[0] + edge->vert[1]->point[0]);
    mid[1] = 0.5 * (edge->vert[0]->point[1] + edge->vert[1]->point[1]);
    mid[2] = 0.5 * (edge->vert[0]->point[2] + edge->vert[1]->point[2]);
    
    LP_Transform_SetToIdentity(trans);
    LP_Transform_Rotate(trans,
			edge->ang / 2,
			edge->z_vec[0],
			edge->z_vec[1],
			edge->z_vec[2]);
    LP_Transform_Point(trans, dir, edge->x_vec, LP_TRANSFORM_NO_OFFSET);

#ifdef DEBUG_EDGE
    float avg[3];
    avg[0] = edge->face[0]->norm[0] + edge->face[1]->norm[0];
    avg[1] = edge->face[0]->norm[1] + edge->face[1]->norm[1];
    avg[2] = edge->face[0]->norm[2] + edge->face[1]->norm[2];
    Normalize(avg);
    printf("dir = (%g,%g,%g), norm avg = (%g,%g,%g)\n",
	   dir[0], dir[1], dir[2],
	   avg[0], avg[1], avg[2]);
#endif
    
    if (isinf(dist = Vef_ConvexRayDist(hull, mid, dir, &start)))
      goto err5;

    if (FTree_Insert(ftree, dist, edge, NULL) < 0)
      goto err5;
    
    for (count = 0; count < 2; count++) {
      if ((hi = Hash_IteratorNew(edge->vert[count]->edges)) == NULL)
	goto err5;
      while (Hash_IteratorNext(hi)) {
	ee = (struct edge *) Hash_IteratorGetData(hi);
	if (Hash_Lookup(hash, ee, NULL))
	  continue;
	if (Hash_Insert(hash, ee, PRESENT, NULL) < 0)
	  goto err6;
	if (Queue_PushBack(queue, ee) < 0)
	  goto err6;
      }
      Hash_IteratorFree(hi);
    }
  }

  LP_Transform_Free(trans);
  Queue_Free(queue);
  Hash_Free(hash);
  return ftree;

 err6:
  Hash_IteratorFree(hi);
 err5:
  LP_Transform_Free(trans);
 err4:
  Queue_Free(queue);
 err3:
  Hash_Free(hash);
 err2:
  FTree_Free(ftree);
 err:
  return NULL;
}

static int CutPart(struct vlh_list **vlh) {
  struct vef *full, *hull;
  struct ftree *ftree;
  struct ftree_node *node;
  struct lp_transform *trans;
  struct edge *edge;
  struct vlh_list *cut, *min = NULL, *last;
  int count, ang_count;
  float norm[3], *nn, err, min_err = INFINITY, dist;

#ifdef DEBUG
  struct lp_vl_list list_a, list_b;
  list_a.vl = (*vlh)->vl;
  list_a.next = &list_b;
  list_b.vl = (*vlh)->hull;
  list_b.next = NULL;
  LP_VertexList_Write("furthest.obj", &list_a, 1);
#endif
  
  if ((full = Vef_New((*vlh)->vl)) == NULL)
    goto err;
  if ((hull = Vef_New((*vlh)->hull)) == NULL)
    goto err2;
  if ((ftree = FurthestEdges(full, hull)) == NULL)
    goto err3;
  if ((trans = LP_Transform_New()) == NULL)
    goto err4;

#ifdef DEBUG
  printf("Cutting part with %zu vertices, %zu edges, and %zu faces\n",
	 Hash_NumEntries(full->verts),
	 Hash_NumEntries(full->edges),
	 Hash_NumEntries(full->faces));
#endif
  
  node = FTree_Highest(ftree);
  for (count = 0; count < NUM_EDGES && node; count++, node = FTree_Prev(ftree, node)) {
    edge = (struct edge *) FTree_GetData(node);
#ifdef DEBUG
    printf("  along edge from (%g, %g, %g) to (%g, %g, %g)\n",
	   edge->vert[0]->point[0],
	   edge->vert[0]->point[1],
	   edge->vert[0]->point[2],
	   edge->vert[1]->point[0],
	   edge->vert[1]->point[1],
	   edge->vert[1]->point[2]);
#endif
    
    norm[0] = edge->face[0]->norm[0];
    norm[1] = edge->face[0]->norm[1];
    norm[2] = edge->face[0]->norm[2];
    LP_Transform_SetToIdentity(trans);
    LP_Transform_Rotate(trans,
			edge->ang / NUM_ANGLES,
			edge->z_vec[0],
			edge->z_vec[1],
			edge->z_vec[2]);
    for (ang_count = NUM_ANGLES - 1; ang_count > 0; ang_count--) {
      nn = ang_count == 0 ? edge->face[1]->norm : norm;
      dist = Dot(nn, edge->vert[0]->point);
      if ((cut = VlhList_Convert(LP_PlaneCut((*vlh)->vl, nn, dist), &err)) == NULL)
	goto err4;
      err = VlhList_TotalSqrError(cut);
      err *= (1 + 1e-3 * fabsf(count - (NUM_EDGES - 1) / 2));
      printf("Error after cut %g\n", err);
      if (err < min_err) {
#ifdef DEBUG
	printf("************ New min found *************\n");
	LP_VertexList_Write("cut.obj", Vl_Convert(cut, 0, 0), 1);
#endif
	min_err = err;
	VlhList_Free(min);
	min = cut;
      } else {
	VlhList_Free(cut);
      }
      
      LP_Transform_Point(trans, norm, norm, LP_TRANSFORM_NO_OFFSET);
      Normalize(norm);
    }
  }
  
  if (min == NULL) {
    LP_Transform_Free(trans);
    FTree_Free(ftree);
    Vef_Free(hull);
    Vef_Free(full);
    return 1;
  }

#ifdef DEBUG
  printf("*****************************************************************\n");
  printf("* Writing min cut to 'min.obj'                                  *\n");
  printf("*****************************************************************\n");
  LP_VertexList_Write("min.obj", Vl_Convert(min, 0, 0), 1);
#endif
  
  last = min;
  while (last->next)
    last = last->next;
  last->next = (*vlh)->next;
  (*vlh)->next = NULL;
  VlhList_Free(*vlh);
  *vlh = min;
  
  LP_Transform_Free(trans);
  FTree_Free(ftree);
  Vef_Free(hull);
  Vef_Free(full);
  return 0;

 err4:
  FTree_Free(ftree);
 err3:
  Vef_Free(hull);
 err2:
  Vef_Free(full);
 err:
  VlhList_Free(min);
  return -1;
}

static const float x_axis[3] = {1, 0, 0};

struct lp_vl_list *LP_ConvexDecomp(const struct lp_vertex_list *in, float threshold) {
  struct vlh_list *vlh;
  struct lp_mass_properties mp;
  float err, thresh;
  int ret;
  
  LP_MassProperties(in, &mp);
  thresh = threshold * mp.volume;
  
  if ((vlh = VlhList_Convert(LP_PlaneCut(in, x_axis, INFINITY), &err)) == NULL)
    goto err;

#ifdef DEBUG
  printf("Init err = %g, thresh = %g, %zu parts\n", err, thresh, VlhList_Len(vlh));
#endif
  while (err > thresh) {
    if ((ret = CutPart(WorstPart(&vlh))) < 0)
      goto err2;
    if (ret == 1)
      break;
    err = VlhList_TotalError(vlh);
#ifdef DEBUG
    printf("err = %g, thresh = %g, %zu parts\n", err, thresh, VlhList_Len(vlh));

    LP_VertexList_Write("decomp.obj", Vl_Convert(vlh, 0, 0), 1.0);
#endif
  }
  
  return Vl_Convert(vlh, 1, 1);

 err2:
  VlhList_Free(vlh);
 err:
  return NULL;
}
