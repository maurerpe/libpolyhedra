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

#include "array.h"
#include "bvh_vl.h"
#include "ftree.h"
#include "util.h"

enum bvh_axis {
  x_axis,
  y_axis,
  z_axis
};

struct bvh_node {
  enum bvh_axis axis;
  float min[3];
  float max[3];
  struct bvh_node *a;
  struct bvh_node *b;
  struct array *points;
};

struct bvh_vl {
  struct lp_vertex_list *vl;
  struct bvh_node *root;
};

static struct bvh_node *BNode_New(size_t init_alloc) {
  struct bvh_node *node;

  if ((node = malloc(sizeof(*node))) == NULL)
    goto err;
  memset(node, 0, sizeof(*node));

  node->min[0] = node->min[1] = node->min[2] =  INFINITY;
  node->max[0] = node->max[1] = node->max[2] = -INFINITY;
  
  if ((node->points = Array_New(init_alloc, NULL)) == NULL)
    goto err2;
  
  return node;

 err2:
  free(node);
 err:
  return NULL;
}

static void BNode_Free(struct bvh_node *node) {
  if (node == NULL)
    return;
  
  BNode_Free(node->a);
  BNode_Free(node->b);
  Array_Free(node->points);
}

static int Split_BNode(struct bvh_node *node, float dist) {
  size_t len;
  enum bvh_axis axis;
  float range[3], median;
  struct ftree *ftree;
  float **data, **stop, *vert;
  struct bvh_node *nn;
  int count;
  
  if ((len = Array_Length(node->points)) < 4)
    return 0;
  
  range[0] = node->max[0] - node->min[0];
  range[1] = node->max[1] - node->min[1];
  range[2] = node->max[2] - node->min[2];
  
  if (range[0] >= range[1] && range[0] >= range[2])
    axis = x_axis;
  else
    axis = range[1] >= range[2] ? y_axis : z_axis;
  node->axis = axis;
  
  if (range[axis] < dist)
    return 0;
  
  if ((ftree = FTree_New(NULL, NULL, NULL)) == NULL)
    goto err;
  
  data = (float **) Array_Data(node->points);
  stop = data + len;
  for (; data < stop; data++)
    if (FTree_Insert(ftree, (*data)[axis], *data, NULL) == NULL)
      goto err2;
  
  median = FTree_GetKey(FTree_Median(ftree));
  
  if (median == node->max[axis] || median == node->min[axis])
    median = 0.5 * (node->max[axis] + node->min[axis]);
  
  if ((node->a = BNode_New(len)) == NULL)
    goto err2;
  if ((node->b = BNode_New(len)) == NULL)
    goto err3;
  
  data = (float **) Array_Data(node->points);
  for (; data < stop; data++) {
    vert = *data;
    nn = vert[axis] <= median ? node->a : node->b;
    Array_Add(nn->points, vert);

    for (count = 0; count < 3; count++) {
      if (vert[count] < nn->min[count])
	nn->min[count] = vert[count];
      
      if (vert[axis] > nn->max[count])
	nn->max[count] = vert[count];
    }
  }
  
  FTree_Free(ftree);
  Array_Free(node->points);
  node->points = NULL;
  
  Split_BNode(node->a, dist);
  Split_BNode(node->b, dist);
  
  return 0;
  
 err3:
  BNode_Free(node->a);
 err2:
  FTree_Free(ftree);
 err:
  return -1;
}

struct bvh_vl *BvhVl_New(struct lp_vertex_list *vl, float dist) {
  struct bvh_vl *bvh;
  struct bvh_node *node;
  size_t fpv;
  float *vert, *stop;
  int count;

  if ((bvh = malloc(sizeof(*bvh))) == NULL)
    goto err;
  memset(bvh, 0, sizeof(*bvh));

  bvh->vl = vl;
  
  if ((bvh->root = node = BNode_New(LP_VertexList_NumVert(vl))) == NULL)
    goto err2;
  
  fpv = LP_VertexList_FloatsPerVert(vl);
  vert = LP_VertexList_GetVert(vl);
  stop = vert + fpv * LP_VertexList_NumVert(vl);
  
  for (; vert < stop; vert += fpv) {
    for (count = 0; count < 3; count++) {
      if (vert[count] < node->min[count])
	node->min[count] = vert[count];
      
      if (vert[count] > node->max[count])
	node->max[count] = vert[count];
    }
    
    Array_Add(node->points, vert);
  }
  
  if (Split_BNode(node, dist) < 0)
    goto err2;
  
  return bvh;

 err2:
  BvhVl_Free(bvh);
 err:
  return NULL;
}

void BvhVl_Free(struct bvh_vl *bvh) {
  if (bvh == NULL)
    return;
  
  BNode_Free(bvh->root);
  
  free(bvh);
}

static float BDist2(struct bvh_node *a, struct bvh_node *b) {
  int count;
  float range[3];
  
  for (count = 0; count < 3; count++) {
    if (a->min[count] < b->min[count]) {
      if (a->max[count] >= b->min[count])
	range[count] = 0;
      else
	range[count] = b->min[count] - a->max[count];
    } else {
      if (b->max[count] >= a->min[count])
	range[count] = 0;
      else
	range[count] = a->min[count] - b->max[count];
    }
  }
  
  return Norm2(range);
}

struct pair_data {
  float dist2;
  struct bvh_vl *bvh;
  bvh_vl_pair_func_t func;
  void *user;
};

static void BNode_Pair_Search(struct bvh_node *node, struct bvh_node *base, const struct pair_data *pd) {
  if (node == NULL || node == base || BDist2(node, base) > pd->dist2)
    return;
  
  BNode_Pair_Search(node->a, base, pd);
  BNode_Pair_Search(node->b, base, pd);
  
  if (node->points) {
    float **data1, **start2, **data2, **stop1, **stop2;
    
    data1 = (float **) Array_Data(node->points);
    stop1 = data1 + Array_Length(node->points);
    start2 = (float **) Array_Data(base->points);
    stop2  = start2 + Array_Length(base->points);
    for (; data1 < stop1; data1++) {
      for (data2 = start2; data2 < stop2; data2++) {
	if (Dist2(*data1, *data2) < pd->dist2)
	  pd->func(pd->user, pd->bvh->vl, *data1, *data2);
      }
    }
  }
}

static void BNode_Pair(struct bvh_node *node, const struct pair_data *pd) {
  if (node == NULL)
    return;
  
  BNode_Pair(node->a, pd);
  BNode_Pair(node->b, pd);

  if (node->points) {
    float **data1, **data2, **stop;
    
    data1 = (float **) Array_Data(node->points);
    stop  = data1 + Array_Length(node->points);
    for (; data1 < stop; data1++) {
      for (data2 = data1 + 1; data2 < stop; data2++) {
	if (Dist2(*data1, *data2) < pd->dist2)
	  pd->func(pd->user, pd->bvh->vl, *data1, *data2);
      }
    }
    
    BNode_Pair_Search(pd->bvh->root, node, pd);
  }
}

void BvhVl_Pairs(struct bvh_vl *bvh, float dist, bvh_vl_pair_func_t func, void *user) {
  struct pair_data pd;

  pd.dist2 = dist * dist;
  pd.bvh   = bvh;
  pd.func  = func;
  pd.user  = user;

  BNode_Pair(bvh->root, &pd);
}
