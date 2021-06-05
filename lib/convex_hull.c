/*
  Copyright (C) 2020-2021 Paul Maurer

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

#include "ftree.h"
#include "hash.h"
#include "libpolyhedra.h"
#include "unique_queue.h"
#include "util.h"

/* Reference: The Quickhull algorithm for convex hulls
 * C. Bradford Barber, David P. Dobkin, Hannu Huhdanpaa
 */

#define PRESENT ((void *) 1)
#define EXTEND  ((void *) 2)
#define DELETE  ((void *) 3)

struct point_list_elem {
  size_t idx;
  struct point_list_elem *next;
};

struct point_list {
  struct point_list_elem *head;
  struct point_list_elem *tail;
  float max_dist;
};

/* Neighbor is joined by idx and next->idx */
struct face_vert {
  size_t idx;
  struct face *neighbor;
  struct face_vert *prev;
  struct face_vert *next;
};

struct face {
  struct face_vert *verts;
  float norm[3];
  float xx[3];
  float yy[3];
  struct point_list *pts;
  struct ftree_node *node;
};

struct ridge_list_elem {
  size_t idx;
  struct face *neighbor;
  int extend;
  struct ridge_list_elem *next;
};

struct ridge_list {
  struct ridge_list_elem *head;
  struct ridge_list_elem *tail;
};

#ifdef DEBUG
static void PrintPoint(FILE *out, const char *msg, size_t idx, const float *data) {
  fprintf(out, "%s %zu = (%g, %g, %g)\n",
	  msg,
	  idx,
	  (double) data[3 * idx + 0],
	  (double) data[3 * idx + 1],
	  (double) data[3 * idx + 2]);
}

static void PrintFace(FILE *out, const struct face_vert *vert, const float *data) {
  const struct face_vert *cur;

  cur = vert;
  do {
    PrintPoint(out, " ", cur->idx, data);
    cur = cur->next;
  } while (cur != vert);
}
#endif

static struct point_list_elem *PointListElem_New(size_t idx) {
  struct point_list_elem *ple;

  if ((ple = malloc(sizeof(*ple))) == NULL)
    goto err;
  memset(ple, 0, sizeof(*ple));

  ple->idx = idx;

  return ple;
  
 err:
  return NULL;
}

static void PointListElem_Free(struct point_list_elem *ple) {
  free(ple);
}

static struct point_list *PointList_New(void) {
  struct point_list *pl;

  if ((pl = malloc(sizeof(*pl))) == NULL)
    goto err;
  memset(pl, 0, sizeof(*pl));
  
  return pl;
  
 err:
  return NULL;
}

static void PointList_Clear(struct point_list *pl) {
  struct point_list_elem *cur, *next;
  
  cur = pl->head;
  while (cur) {
    next = cur->next;
    PointListElem_Free(cur);
    cur = next;
  }
  
  memset(pl, 0, sizeof(*pl));
}

static void PointList_Free(struct point_list *pl) {
  if (pl == NULL)
    return;
  
  PointList_Clear(pl);
  
  free(pl);
}

static void PointList_Add(struct point_list *pl, struct point_list_elem *elem, float dist) {
  if (pl->head == NULL) {
    pl->max_dist = dist;
    pl->head = elem;
    pl->tail = elem;
    elem->next = NULL;
    return;
  }
  
  if (dist > pl->max_dist) {
    pl->max_dist = dist;
    elem->next = pl->head;
    pl->head = elem;
    return;
  }
  
  pl->tail->next = elem;
  pl->tail = elem;
  elem->next = NULL;
}

static void PointList_Join(struct point_list *dest, struct point_list *src) {
  struct point_list_elem *elem;
  
  if (src->head == NULL)
    return;

  elem = src->head->next;
  PointList_Add(dest, src->head, src->max_dist);

  if (src->head == src->tail) {
    memset(src, 0, sizeof(*src));
    return;
  }
  
  dest->tail->next = elem;
  dest->tail = src->tail;
  memset(src, 0, sizeof(*src));
}

static struct face_vert *FaceVert_New(size_t idx, struct face_vert *prev) {
  struct face_vert *fv;

  if ((fv = malloc(sizeof(*fv))) == NULL)
    goto err;
  memset(fv, 0, sizeof(*fv));
  
  fv->idx = idx;
  if (prev) {
    if (prev->next != prev) {
      fv->next = prev->next;
      fv->prev = prev;
      prev->next->prev = fv;
      prev->next = fv;
    } else {
      prev->next = fv;
      prev->prev = fv;
      fv->next = prev;
      fv->prev = prev;
    }
  } else {
    fv->next = fv;
    fv->prev = fv;
  }
  
  return fv;
  
 err:
  return NULL;
}

static void FaceVert_Free(struct face_vert *fv) {
  struct face_vert *cur, *next;
  
  if (fv == NULL)
    return;

  cur = fv;
  do {
    next = cur->next;
    free(cur);
    cur = next;
  } while (cur != fv);
}

static struct face_vert *FaceVert_FindVert(struct face_vert *fv, size_t pt) {
  struct face_vert *cur;
  
  cur = fv;
  while (1) {
    if (cur->idx == pt)
      break;
    
    cur = cur->next;
    if (cur == fv) {
      fprintf(stderr, "Internal Error: convex_hull.c: Face does not contain requested vert\n");
      return NULL;
    }
  }
  
  return cur;
}

static struct face_vert *FaceVert_FindEdge(struct face_vert *fv, size_t pt1, size_t pt2) {
  struct face_vert *cur;
  
  if ((cur = FaceVert_FindVert(fv, pt1)) == NULL)
    return NULL;
  
  if (cur->next->idx != pt2) {
    fprintf(stderr, "Internal Error: convex_hull.c: Face does not contain requested edge\n");
    return NULL;
  }
  
  return cur;
}

static void FaceVert_PrepForRetention(struct face_vert **fv, struct hash *visited) {
  struct face_vert *cur = *fv;
  
  while (Hash_Lookup(visited, cur->neighbor, NULL) != DELETE)
    cur = cur->next;
  
  *fv = cur;
}

static void FaceVert_PrepForExtend(struct face_vert **fv, struct hash *visited) {
  struct face_vert *del, *cur = *fv;
  void *cat;

  while ((cat = Hash_Lookup(visited, cur->neighbor, NULL)) == DELETE || cat == EXTEND)
    cur = cur->prev;
  
  while ((cat = Hash_Lookup(visited, cur->neighbor, NULL)) != DELETE && cat != EXTEND)
    cur = cur->next;
  
  while ((cat = Hash_Lookup(visited, cur->next->neighbor, NULL)) == DELETE || cat == EXTEND) {
    del = cur->next;
#ifdef DEBUG
    printf("Deleting vert %zu\n", del->idx);
#endif
    del->prev->next = del->next;
    del->next->prev = del->prev;
    del->next = del;
    del->prev = del;
    FaceVert_Free(del);
  }

  *fv = cur;
}

static struct face_vert *FaceVert_Extend(struct face_vert **fv, size_t pt) {
  
  return *fv = FaceVert_New(pt, *fv);
}

static struct face *Face_New(size_t idx0, size_t idx1, size_t idx2, struct hash *faces, const float *data) {
  struct face *face;
  
  if ((face = malloc(sizeof(*face))) == NULL)
    goto err;
  memset(face, 0, sizeof(*face));

  if ((face->verts = FaceVert_New(idx0, NULL)) == NULL)
    goto err2;
  if (FaceVert_New(idx1, face->verts) == NULL)
    goto err3;
  if (FaceVert_New(idx2, face->verts->next) == NULL)
    goto err3;

  PlaneNorm(face->norm, data + 3 * idx0, data + 3 * idx1, data + 3 * idx2);
  BasisVectors(face->xx, face->yy, face->norm);
  
  if ((face->pts = PointList_New()) == NULL)
    goto err3;
  
  if (Hash_Insert(faces, face, PRESENT, NULL) < 0)
    goto err4;
  
#ifdef DEBUG
  printf("New face:\n");
  PrintFace(stdout, face->verts, data);
#endif
  return face;

 err4:
  PointList_Free(face->pts);
 err3:
  FaceVert_Free(face->verts);
 err2:
  free(face);
 err:
  return NULL;
}

static void Face_Free(struct face *face) {
  if (face == NULL)
    return;
  
  PointList_Free(face->pts);
  FaceVert_Free(face->verts);
  free(face);
}

static void Face_Free_Func(void *user, void *data) {
  Face_Free((struct face *) data);
}

static int Face_Update(struct face *face, struct ftree *ftree) {
  if (face->pts->head == NULL) {
    if (face->node) {
      FTree_Delete(ftree, face->node);
      face->node = NULL;
    }
    return 0;
  }
  
  if (face->node == NULL) {
    face->node = FTree_Insert(ftree, face->pts->max_dist, face, NULL);
    return face->node ? 0 : -1;
  }
  
  if (face->pts->max_dist != FTree_GetKey(face->node))
    FTree_Rekey(ftree, face->node, face->pts->max_dist, NULL);
  
  return 0;
}

static void *Categorize(const struct face *face, size_t idx, const float *data, float *dist_out) {
  const float *pt, *vert;
  float delta[3], dist, x1, x2, y1, y2, dx, dy, dd, max, area, tol, dpt;
  struct face_vert *fv;
  
  pt   = data + 3 * idx;
  fv = face->verts;
  vert = data + 3 * fv->prev->idx;
  delta[0] = vert[0] - pt[0];
  delta[1] = vert[1] - pt[1];
  delta[2] = vert[2] - pt[2];
  dist = Dot(delta, face->norm);
  if (dist_out)
    *dist_out = dist;

  area = 0;
  max = -INFINITY;
  x2 = Dot(delta, face->xx);
  y2 = Dot(delta, face->yy);
  do {
    x1 = x2;
    y1 = y2;
    vert = data + 3 * fv->idx;
    delta[0] = vert[0] - pt[0];
    delta[1] = vert[1] - pt[1];
    delta[2] = vert[2] - pt[2];
    x2 = Dot(delta, face->xx);
    y2 = Dot(delta, face->yy);

    area += x1 * y2 - y1 * x2;
    
    dx = x2 - x1;
    dy = y2 - y1;
    
    dd = (dy * x1 - dx * y1) / sqrtf(dx * dx + dy * dy);
    if (dd > max)
      max = dd;
    
    fv = fv->next;
  } while (fv != face->verts);
  
  tol = 1e-5 * sqrtf(fabsf(area));
  
#ifdef DEBUG
  printf("Max = %g, dist = %g, tol = %g\n", max, dist, tol);
#endif
  
  if (max > 0) {
    if (fabsf(dist) < tol || fabsf(dist) < 1e-5 * max)
      return EXTEND;
    if (dist > 0)
      return DELETE;
    return PRESENT;
  }
  
  if (dist > tol)
    return DELETE;
  
  dpt = dist + tol;
  if (dpt * dpt + max * max < 4 * tol * tol)
    return EXTEND;
  
  return PRESENT;
}

static void Face_AssignPoints(struct face *face, struct point_list *pool, const float *data) {
  struct point_list_elem *prev, *cur;
  float dist;
  
  prev = pool->head;
  while ((cur = prev->next)) {
    if (Categorize(face, cur->idx, data, &dist) == DELETE) {
      prev->next = cur->next;
      PointList_Add(face->pts, cur, dist);
    } else {
      prev = cur;
    }
  }
  pool->tail = prev;
}

static struct ridge_list_elem *RidgeListElem_New(size_t idx, int extend, struct face *neighbor) {
  struct ridge_list_elem *rle;

  if ((rle = malloc(sizeof(*rle))) == NULL)
    goto err;
  memset(rle, 0, sizeof(*rle));

  rle->idx = idx;
  rle->extend = extend;
  rle->neighbor = neighbor;
  
  return rle;
  
 err:
  return NULL;
}

static struct ridge_list_elem *RidgeListElem_NewV(struct face *neighbor, struct hash *visited) {
  int extend;

  if ((extend = Hash_Lookup(visited, neighbor, NULL) == EXTEND))
    FaceVert_PrepForExtend(&neighbor->verts, visited);
  else
    FaceVert_PrepForRetention(&neighbor->verts, visited);

  return RidgeListElem_New(neighbor->verts->next->idx, extend, neighbor);
}

static void RidgeListElem_Free(struct ridge_list_elem *rle) {
  free(rle);
}

static struct ridge_list *RidgeList_New(void) {
  struct ridge_list *rl;

  if ((rl = malloc(sizeof(*rl))) == NULL)
    goto err;
  memset(rl, 0, sizeof(*rl));

  return rl;
  
 err:
  return NULL;
}

static void RidgeList_Clear(struct ridge_list *rl) {
  struct ridge_list_elem *cur, *next;
  
  cur = rl->head;
  while (cur) {
    next = cur->next;
    RidgeListElem_Free(cur);
    cur = next;
  }
  
  memset(rl, 0, sizeof(*rl));
}

static void RidgeList_Free(struct ridge_list *rl) {
  if (rl == NULL)
    return;
  
  RidgeList_Clear(rl);
  
  free(rl);
}

static void RidgeList_AppendRle(struct ridge_list *rl, struct ridge_list_elem *rle) {
  rle->next = NULL;
  if (rl->head == NULL) {
    rl->head = rle;
  } else {
    rl->tail->next = rle;
  }
  rl->tail = rle;
}

static int RidgeList_Append(struct ridge_list *rl, size_t idx, int extend, struct face *neighbor) {
  struct ridge_list_elem *rle;
  
  if ((rle = RidgeListElem_New(idx, extend, neighbor)) == NULL)
    return -1;
  
  RidgeList_AppendRle(rl, rle);
  
  return 0;
}

static int RidgeList_AppendV(struct ridge_list *rl, struct face *neighbor, struct hash *visited) {
  struct ridge_list_elem *rle;
  
  if ((rle = RidgeListElem_NewV(neighbor, visited)) == NULL)
    return -1;
  
  RidgeList_AppendRle(rl, rle);
  
  return 0;
}

static int BuildVl_AddFace(struct lp_vertex_list *out, struct face *face, const float *data) {
  struct face_vert *cur;
  
  cur  = face->verts->next->next;
  while (cur != face->verts) {
    if (LP_VertexList_Add(out, data + 3 * face->verts->idx) == UINT_MAX)
      return -1;
    if (LP_VertexList_Add(out, data + 3 * cur->idx) == UINT_MAX)
      return -1;
    if (LP_VertexList_Add(out, data + 3 * cur->prev->idx) == UINT_MAX)
      return -1;
    cur = cur->next;
  }
  
  return 0;
}

static struct lp_vertex_list *BuildVl(struct hash *faces, const float *data) {
  struct lp_vertex_list *out;
  struct hash_iterator *hi;
  
  if ((out = LP_VertexList_New(3, lp_pt_triangle)) == NULL)
    goto err;
  
  if ((hi = Hash_IteratorNew(faces)) == NULL)
    goto err2;
  while (Hash_IteratorNext(hi)) {
    if (BuildVl_AddFace(out, (struct face *) Hash_IteratorGetKey(hi), data) < 0)
      goto err3;
  }
  Hash_IteratorFree(hi);

  return out;

 err3:
  Hash_IteratorFree(hi);
 err2:
  LP_VertexList_Free(out);
 err:
  return NULL;
}

static int BuildNewFaces(struct ridge_list *rl, struct point_list *pool, struct hash *faces, struct ftree *faces_with_pts, const float *data) {
  struct ridge_list_elem *rle;
  struct face *face, *face_prev, **neighbor_prev, *first_face, **first_neighbor;
  struct face_vert *vcur;
  size_t prev_idx, idx;

  first_neighbor = NULL;
  neighbor_prev = NULL;
  first_face = NULL;
  idx = pool->head->idx;
  prev_idx = rl->tail->idx;
  face_prev = NULL;
  for (rle = rl->head; rle; prev_idx = rle->idx, rle = rle->next, face_prev = face) {
    if (rle->extend) {
      face = rle->neighbor;
      if (FaceVert_Extend(&face->verts, idx) == NULL)
	goto err;
#ifdef DEBUG
      printf("After extension:\n");
      PrintFace(stdout, face->verts, data);
      if (face->verts->idx == face->verts->prev->idx)
	exit(1);
#endif
      face->verts->prev->neighbor = face_prev;
      if (first_neighbor == NULL)
	first_neighbor = &face->verts->prev->neighbor;
    } else {
      if ((face = Face_New(idx, rle->idx, prev_idx, faces, data)) == NULL)
	goto err;
      face->verts->prev->neighbor = face_prev;
      face->verts->next->neighbor = rle->neighbor;
      if (first_neighbor == NULL)
	first_neighbor = &face->verts->prev->neighbor;

      if ((vcur = FaceVert_FindEdge(rle->neighbor->verts, prev_idx, rle->idx)) == NULL)
	goto err;
      vcur->neighbor = face;
    }
    
    if (neighbor_prev)
      *neighbor_prev = face;
    if (rle->extend)
      neighbor_prev = &rle->neighbor->verts->neighbor;
    else
      neighbor_prev = &face->verts->neighbor;
    if (first_face == NULL)
      first_face = face;

    Face_AssignPoints(face, pool, data);

    if (Face_Update(face, faces_with_pts) < 0)
      goto err;
  }
  
  *first_neighbor = face;
  *neighbor_prev = first_face;
  
#ifdef DEBUG
  struct point_list_elem *cur;
  for (cur = pool->head->next; cur; cur = cur->next)
    PrintPoint(stdout, "Dropping interior point", cur->idx, data);

  struct hash_iterator *hi;
  struct face_vert *fv;
  int found_err = 0;
  if ((hi = Hash_IteratorNew(faces)) == NULL)
    goto err;
  while (Hash_IteratorNext(hi)) {
    face = (struct face *) Hash_IteratorGetKey(hi);
    vcur = face->verts;
    do {
      if ((fv = FaceVert_FindEdge(vcur->neighbor->verts, vcur->next->idx, vcur->idx)) == NULL || fv->neighbor != face) {
	printf("Internal Error: convex_hull.c: Incorrect neighbor after building faces\n");
	found_err = 1;
      }
      vcur = vcur->next;
    } while (vcur != face->verts);
  }
  Hash_IteratorFree(hi);
  if (found_err)
    goto err;
#endif
  
  return 0;
  
 err:
  return -1;
}

static int FindHull(struct hash *faces, struct ftree *faces_with_pts, const float *data) {
  struct point_list *pool;
  struct ridge_list *rl;
  struct hash *visited;
  struct unique_queue *queued;
  struct hash_iterator *hi;
  struct face *face, *neighbor, *no_view;
  struct face_vert *cur;
  struct ftree_node *node;
  size_t idx, first_idx;
  void *cat;
  int found;
#ifdef DEBUG
  int count = 0;
  char buf[256];
#endif
  
  if ((pool = PointList_New()) == NULL)
    goto err;
  
  if ((rl = RidgeList_New()) == NULL)
    goto err2;
  
  if ((visited = Hash_NewPtr(NULL, NULL, NULL, NULL, NULL)) == NULL)
    goto err3;
  
  if ((queued = UniqueQueue_New()) == NULL)
    goto err4;
  
  while ((node = FTree_Highest(faces_with_pts))) {
    /* Found face with points above */
    face = (struct face *) FTree_GetData(node);
    if (face->pts->head == NULL) {
      Face_Update(face, faces_with_pts);
      continue;
    }
    PointList_Join(pool, face->pts);
    idx = pool->head->idx;

#ifdef DEBUG
    struct lp_vl_list list;
    memset(&list, 0, sizeof(list));
    if ((list.vl = BuildVl(faces, data))) {
      snprintf(buf, sizeof(buf), "hull%03d.obj", count++);
      LP_VertexList_Write(buf, &list, 1);
      LP_VertexList_Free(list.vl);
    }
    PrintPoint(stdout, "\nAdding point", idx, data);
#endif
    
    /* Find deletion face */
    if (Categorize(face, idx, data, NULL) != DELETE) {
      found = 0;
      cur = face->verts;
      do {
	if (Categorize(cur->neighbor, idx, data, NULL) == DELETE) {
	  found = 1;
	  face = cur->neighbor;
	  break;
	}

	cur = cur->next;
      } while (cur != face->verts);
      
      if (!found) {
	/* No deletion face, reassign points in the pool and try next point */
	Face_AssignPoints(face, pool, data);
	Face_Update(face, faces_with_pts);
	do {
	  Face_AssignPoints(cur->neighbor, pool, data);
	  Face_Update(face, faces_with_pts);
	  
	  cur = cur->next;
	} while (cur != face->verts);
	
#ifdef DEBUG
	printf("Could not find deletion face\n");
	
	struct point_list_elem *cur;
	for (cur = pool->head; cur != NULL; cur = cur->next)
	  PrintPoint(stdout, "Dropping point", cur->idx, data);
#endif
	PointList_Clear(pool);
	continue;
      }
    }
    
    /* Identify all faces with view of point */
    no_view = NULL;
    do {
      cat = Categorize(face, idx, data, NULL);
#ifdef DEBUG
      printf("Face marked for %s\n", cat == DELETE ? "deletion" : cat == EXTEND ? "extension" : "retention");
      PrintFace(stdout, face->verts, data);
#endif
      if (Hash_Insert(visited, face, cat, NULL) < 0)
	goto err5;
      if (cat != DELETE) {
	no_view = face;
	continue;
      }
      
      face->pts->max_dist = 0;
      PointList_Join(pool, face->pts);
      
      cur = face->verts;
      do {
	if (!Hash_Lookup(visited, cur->neighbor, NULL))
	  if (UniqueQueue_PushBack(queued, cur->neighbor) < 0)
	    goto err5;
	cur = cur->next;
      } while (cur != face->verts);
    } while ((face = UniqueQueue_Pop(queued)));
    if (no_view == NULL) {
      fprintf(stderr, "Internal error: convex_hull.c: All faces can view point\n");
      goto err5;
    }
    
    /* Found face that cannot see point: must be on ridge */
    face = no_view;
#ifdef DEBUG
    printf("Before first append\n");
    PrintFace(stdout, face->verts, data);
#endif
    if (RidgeList_AppendV(rl, face, visited) < 0)
      goto err5;
#ifdef DEBUG
    printf("After first append\n");
    PrintFace(stdout, face->verts, data);
#endif
    first_idx = face->verts->idx;
    face = rl->tail->neighbor;
    
    /* Trace around ridge */
#ifdef DEBUG
    PrintPoint(stdout, "Ridge: First point", first_idx, data);
#endif
    while (rl->tail->idx != first_idx) {
#ifdef DEBUG
      PrintPoint(stdout, "Ridge: Current point", rl->tail->idx, data);
      PrintFace(stdout, face->verts, data);
#endif
      if ((cur = FaceVert_FindVert(face->verts, rl->tail->idx)) == NULL) {
#ifdef DEBUG
	printf("Could not find vert %zu on face\n", rl->tail->idx);
#endif
	goto err5;
      }
      
      neighbor = cur->neighbor;
      if ((cat = Hash_Lookup(visited, neighbor, NULL)) == DELETE) {
	RidgeList_Append(rl, cur->next->idx, 0, face);
      } else if (cat == EXTEND) {
	RidgeList_AppendV(rl, neighbor, visited);
	face = neighbor;
      } else {
	face = neighbor;
      }
    }
    
    /* Delete old faces */
    if ((hi = Hash_IteratorNew(visited)) == NULL)
      goto err5;
    while (Hash_IteratorNext(hi)) {
      if ((cat = Hash_IteratorGetData(hi)) == DELETE || cat == EXTEND) {
	face = (struct face *) Hash_IteratorGetKey(hi);
	if (Face_Update(face, faces_with_pts) < 0)
	  goto err6;
	if (cat == DELETE)
	  Hash_Remove(faces, face);
      }
    }
    Hash_IteratorFree(hi);

    /* Build new faces */
    if (BuildNewFaces(rl, pool, faces, faces_with_pts, data) < 0)
      goto err5;
    
    if (pool->head->idx != idx)
      fprintf(stderr, "Internal error: convex_hull.c: pool corruption\n");
    
    Hash_Clear(visited);
    PointList_Clear(pool);
    RidgeList_Clear(rl);
  }

  UniqueQueue_Free(queued);
  Hash_Free(visited);
  RidgeList_Free(rl);
  PointList_Free(pool);
  return 0;

 err6:
  Hash_IteratorFree(hi);
 err5:
  UniqueQueue_Free(queued);
 err4:
  Hash_Free(visited);
 err3:
  RidgeList_Free(rl);
 err2:
  PointList_Free(pool);
 err:
  return -1;
}

static int InitSimplex(size_t len, const float *data, struct hash *faces, struct ftree *faces_with_pts) {
  float ff, min_f, max_f, dd_f, dist;
  const float *max_p, *min_p;
  size_t idx, min_idx, max_idx, dd_idx, temp_idx;
  struct face *face;
  struct point_list *pool, *below, *temp_pl;
  struct point_list_elem *elem;
  struct ridge_list *rl;
  struct face_vert *cur;
  void *cat;
  
  if (len < 4) {
    fprintf(stderr, "Cannot build convex hull from less than 4 points: %zu unique points found\n", len);
    goto err;
  }
  
  /* Find points with largest and smallest x values */
  min_f = INFINITY;
  max_f = -INFINITY;
  min_idx = 0;
  max_idx = 0;
  for (idx = 0; idx < len; idx++) {
    ff = data[3 * idx];
    if (ff > max_f) {
      max_f = ff;
      max_idx = idx;
    }
    if (ff < min_f) {
      min_f = ff;
      min_idx = idx;
    }
  }

  /* Find largest combined distance from min/max x, build initial face */
  dd_f = 0;
  dd_idx = 0;
  min_p = data + 3 * min_idx;
  max_p = data + 3 * max_idx;
  for (idx = 0; idx < len; idx++) {
    dist = Dist(data + 3 * idx, min_p) + Dist(data + 3 * idx, max_p);
    if (dist > dd_f) {
      dd_f = dist;
      dd_idx = idx;
    }
  }
  
  if ((face = Face_New(min_idx, max_idx, dd_idx, faces, data)) == NULL)
    goto err;
  
  if (Norm2(face->norm) == 0) {
    fprintf(stderr, "Cannot create convex hull: All points are colinear\n");
    goto err;
  }
  
  /* Sort all points as above, below, or coplaner */
  if ((pool = PointList_New()) == NULL)
    goto err;
  
  if ((below = PointList_New()) == NULL)
    goto err2;
  
  for (idx = 0; idx < len; idx++) {
    if (idx == min_idx || idx == max_idx || idx == dd_idx)
      continue;
    
    if ((elem = PointListElem_New(idx)) == NULL)
      goto err2;

    if ((cat = Categorize(face, idx, data, &dist)) == DELETE) {
      PointList_Add(face->pts, elem, dist);
    } else if (cat == EXTEND) {
      PointList_Add(pool, elem, fabsf(dist));
    } else {
      PointList_Add(below, elem, -dist);
    }
  }

  /* If furthest point was above, flip face */
  if (face->pts->max_dist > below->max_dist) {
    temp_pl = below;
    below = face->pts;
    face->pts = temp_pl;
    temp_idx = face->verts->idx;
    face->verts->idx = face->verts->next->idx;
    face->verts->next->idx = temp_idx;
    face->norm[0] = -face->norm[0];
    face->norm[1] = -face->norm[1];
    face->norm[2] = -face->norm[2];
  }
  
  if (below->head == NULL) {
    fprintf(stderr, "Cannot create convex hull: All points coplaner\n");
    goto err3;
  }

  if (Face_Update(face, faces_with_pts) < 0)
    goto err3;
  
  /* Build remaining faces */
  PointList_Join(pool, below);
  if ((rl = RidgeList_New()) == NULL)
    goto err3;
  cur = face->verts;
  for (idx = 0; idx < 3; idx++, cur = cur->next) {
    if (RidgeList_Append(rl, cur->idx, 0, face) < 0)
      goto err4;
  }
  if (BuildNewFaces(rl, pool, faces, faces_with_pts, data) < 0)
    goto err4;
  
  /* Sucess: cleanup and return */
  RidgeList_Free(rl);
  PointList_Free(below);
  PointList_Free(pool);
  return 0;

 err4:
  RidgeList_Free(rl);
 err3:
  PointList_Free(below);
 err2:
  PointList_Free(pool);
 err:
  return -1;
}

struct lp_vertex_list *LP_ConvexHull(const struct lp_vertex_list *in) {
  struct lp_vertex_list *in3, *out;
  struct hash *faces;
  struct ftree *faces_with_pts;
  const float *data;
  size_t fpv, idx, len;
  
  if ((in3 = LP_VertexList_New(3, lp_pt_triangle)) == NULL)
    goto err;
  
  if ((fpv = LP_VertexList_FloatsPerVert(in)) != 3) {
    if (fpv < 3) {
      fprintf(stderr, "Error: Need at least 3 floats per vertex for convex hull\n");
      goto err2;
    }

    data = LP_VertexList_GetVert(in);
    len = LP_VertexList_NumVert(in);
    for (idx = 0; idx < len; idx++, data += fpv) {
      if (LP_VertexList_Add(in3, data) == UINT_MAX)
	goto err2;
    }
    in = in3;
  }
  
  data = LP_VertexList_GetVert(in);
  len  = LP_VertexList_NumVert(in);

#ifdef DEBUG
  printf("Finding convex hull of %zu points\n", len);
#endif
  
  if ((faces = Hash_NewPtr(NULL, Face_Free_Func, NULL, NULL, NULL)) == NULL)
    goto err2;
  
  if ((faces_with_pts = FTree_New(NULL, NULL, NULL)) == NULL)
    goto err3;
  
  if (InitSimplex(len, data, faces, faces_with_pts) < 0)
    goto err4;
  
  if (FindHull(faces, faces_with_pts, data) < 0)
    goto err4;
  
  if ((out = BuildVl(faces, data)) == NULL)
    goto err4;
  
  FTree_Free(faces_with_pts);
  Hash_Free(faces);
  LP_VertexList_Free(in3);
#ifdef DEBUG
  printf("Returning convex hull with %zu faces\n", LP_VertexList_NumInd(out) / 3);
#endif
  return out;

 err4:
  FTree_Free(faces_with_pts);
 err3:
  Hash_Free(faces);
 err2:
  LP_VertexList_Free(in3);
 err:
  fprintf(stderr, "Error: Could not build convex hull\n");
  return NULL;
}
