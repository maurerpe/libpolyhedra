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

#include "libpolyhedra.h"

#include "array.h"
#include "bvh_vl.h"
#include "ftree.h"
#include "hash.h"
#include "SipHash/siphash.h"
#include "util.h"

/* Reference: Surface Simplification Using Quadratic Error Metrics
 * Michael Garland, Paul S. Heckbert
 */

#define PRESENT ((void *) 1)

struct vert {
  float v[3];
  float Q[10];
  struct hash *pair_hash;
  struct array *face_arr;
};

struct pair {
  struct vert *vert[2];
  float vbar[3];
  struct ftree_node *node;
};

struct face {
  struct vert *vert[3];
};

static void CalcKp(float *kp, struct vert **vert) {
  float norm[3], a, b, c, d;
  
  PlaneNorm(norm, vert[0]->v, vert[1]->v, vert[2]->v);
  
  a = norm[0];
  b = norm[1];
  c = norm[2];
  d = -Dot(norm, vert[0]->v);
  
  kp[0] = a * a;
  kp[1] = a * b;
  kp[2] = a * c;
  kp[3] = a * d;
  kp[4] = b * b;
  kp[5] = b * c;
  kp[6] = b * d;
  kp[7] = c * c;
  kp[8] = c * d;
  kp[9] = d * d;
}

static float CalcCost(const float *vbar, const float *Qbar) {
  float prod[4];
  
  prod[0] =
    Qbar[0] * vbar[0] +
    Qbar[1] * vbar[1] +
    Qbar[2] * vbar[2] +
    Qbar[3];

  prod[1] =
    Qbar[1] * vbar[0] +
    Qbar[4] * vbar[1] +
    Qbar[5] * vbar[2] +
    Qbar[6];
  
  prod[2] =
    Qbar[2] * vbar[0] +
    Qbar[5] * vbar[1] +
    Qbar[7] * vbar[2] +
    Qbar[8];
  
  prod[3] =
    Qbar[3] * vbar[0] +
    Qbar[6] * vbar[1] +
    Qbar[8] * vbar[2] +
    Qbar[9];

  return Dot(vbar, prod) + prod[3];
}

static float CalcLowestCost(struct pair *pair) {
  float Qbar[10], mat[9], bb[3], mid[3], a, b, c;
  int count;
  
  for (count = 0; count < 10; count++)
    Qbar[count] = pair->vert[0]->Q[count] + pair->vert[1]->Q[count];
  
  mat[0] = Qbar[0];
  mat[1] = Qbar[1];
  mat[2] = Qbar[2];
  
  mat[3] = Qbar[1];
  mat[4] = Qbar[4];
  mat[5] = Qbar[5];
  
  mat[6] = Qbar[2];
  mat[7] = Qbar[5];
  mat[8] = Qbar[7];
  
  bb[0] = -Qbar[3];
  bb[1] = -Qbar[6];
  bb[2] = -Qbar[8];
  
  if (Solve3x3(pair->vbar, mat, bb) == 0)
    return CalcCost(pair->vbar, Qbar);
  
  mid[0] = 0.5 * (pair->vert[0]->v[0] + pair->vert[1]->v[0]);
  mid[1] = 0.5 * (pair->vert[0]->v[1] + pair->vert[1]->v[1]);
  mid[2] = 0.5 * (pair->vert[0]->v[2] + pair->vert[1]->v[2]);
  
  a = CalcCost(pair->vert[0]->v, Qbar);
  b = CalcCost(pair->vert[1]->v, Qbar);
  c = CalcCost(mid, Qbar);
  
  if (a <= b) {
    if (c <= a) {
      pair->vbar[0] = mid[0];
      pair->vbar[1] = mid[1];
      pair->vbar[2] = mid[2];
      return c;
    } else {
      pair->vbar[0] = pair->vert[0]->v[0];
      pair->vbar[1] = pair->vert[0]->v[1];
      pair->vbar[2] = pair->vert[0]->v[2];
      return a;
    }
  } else {
    if (c <= b) {
      pair->vbar[0] = mid[0];
      pair->vbar[1] = mid[1];
      pair->vbar[2] = mid[2];
      return c;
    } else {
      pair->vbar[0] = pair->vert[1]->v[0];
      pair->vbar[1] = pair->vert[1]->v[1];
      pair->vbar[2] = pair->vert[1]->v[2];
      return b;
    }
  }
}

static struct pair *Pair_New(struct ftree *pairs, struct vert *a, struct vert *b) {
  struct pair *pair;
  float cost;
  
  if ((pair = malloc(sizeof(*pair))) == NULL)
    goto err;
  memset(pair, 0, sizeof(*pair));
  
  pair->vert[0] = a;
  pair->vert[1] = b;
  
  if (Hash_Insert(a->pair_hash, b, pair, NULL) < 0)
    goto err2;
  
  if (Hash_Insert(b->pair_hash, a, pair, NULL) < 0)
    goto err3;
  
  cost = CalcLowestCost(pair);
  if ((pair->node = FTree_Insert(pairs, cost, pair, NULL)) == NULL)
    goto err4;
  
  return pair;

 err4:
  Hash_Remove(b->pair_hash, a);
 err3:
  Hash_Remove(a->pair_hash, b);
 err2:
  free(pair);
 err:
  return NULL;
}

static void Pair_Free(struct pair *pair) {
  if (pair == NULL)
    return;
  
  free(pair);
}

static void Pair_Free_Func(void *data) {
  Pair_Free((struct pair *) data);
}

static void Face_Cannonize(struct face *face) {
  struct vert *temp;
  
  if (face->vert[0] < face->vert[1] && face->vert[0] < face->vert[2])
    return;

  temp = face->vert[0];
  if (face->vert[1] < face->vert[2]) {
    face->vert[0] = face->vert[1];
    face->vert[1] = face->vert[2];
    face->vert[2] = temp;
    return;
  }

  face->vert[0] = face->vert[2];
  face->vert[2] = face->vert[1];
  face->vert[1] = temp;
}

static struct face *Face_New(struct hash *faces, struct vert *a, struct vert *b, struct vert *c) {
  struct face *face;
  float Kp[10];
  int count, idx;
  
  if ((face = malloc(sizeof(*face))) == NULL)
    goto err;
  memset(face, 0, sizeof(*face));
  face->vert[0] = a;
  face->vert[1] = b;
  face->vert[2] = c;
  Face_Cannonize(face);
  
  CalcKp(Kp, face->vert);
  for (count = 0; count < 3; count++) {
    for (idx = 0; idx < 10; idx++)
      face->vert[count]->Q[idx] += Kp[idx];
    if (Array_Add(face->vert[count]->face_arr, face) < 0)
      goto err2;
  }
  
  if (Hash_Insert(faces, face, PRESENT, NULL) < 0)
    goto err2;
  
  return face;
  
 err2:
  while (count-- > 0)
    Array_Remove(face->vert[count]->face_arr, -1);
  free(face);
 err:
  return NULL;
}

static void Face_Free(struct face *face) {
  if (face == NULL)
    return;
  
  free(face);
}

static void Face_Free_Func(void *user, void *data) {
  Face_Free((struct face *) data);
}

static struct vert *Vert_New(struct hash *verts, const float *v) {
  struct vert *vv;
  
  if ((vv = malloc(sizeof(*vv))) == NULL)
    goto err;
  memset(vv, 0, sizeof(*vv));
  vv->v[0] = v[0];
  vv->v[1] = v[1];
  vv->v[2] = v[2];
  
  if ((vv->face_arr = Array_New(8, NULL)) == NULL)
    goto err2;
  
  if ((vv->pair_hash = Hash_NewPtr(NULL, NULL, NULL, NULL, NULL)) == NULL)
    goto err3;
  
  if (Hash_Insert(verts, vv, PRESENT, NULL) < 0)
    goto err4;
  
  return vv;
  
 err4:
  Hash_Free(vv->pair_hash);
 err3:
  Array_Free(vv->face_arr);
 err2:
  free(vv);
 err:
  return NULL;
}

static void Vert_Free(struct vert *vert) {
  Hash_Free(vert->pair_hash);
  Array_Free(vert->face_arr);
  free(vert);
}

static void Vert_Free_Func(void *user, void *data) {
  Vert_Free((struct vert *) data);
}

static int Add_Pairs(struct ftree *pairs, struct hash *faces) {
  struct hash_iterator *hi;
  struct face *face;
  int count, cp1;
  
  if ((hi = Hash_IteratorNew(faces)) == NULL)
    goto err;
  
  while (Hash_IteratorNext(hi)) {
    face = (struct face *) Hash_IteratorGetKey(hi);
    
    for (count = 0; count < 3; count++) {
      cp1 = (count + 1) % 3;
      if (Hash_Lookup(face->vert[count]->pair_hash, face->vert[cp1], NULL) == NULL)
	if (Pair_New(pairs, face->vert[count], face->vert[cp1]) == NULL)
	  goto err2;
    }
  }
  
  Hash_IteratorFree(hi);
  return 0;
  
 err2:
  Hash_IteratorFree(hi);
 err:
  return -1;
}

struct agg_bvh_pair {
  struct ftree *pairs;
  struct vert **vert_arr;
  int err;
};

static void Agg_Bvh_Pair(void *user, struct lp_vertex_list *vl, float *fa, float *fb) {
  struct agg_bvh_pair *abp = (struct agg_bvh_pair *) user;
  struct vert *a, *b;
  float *arr;
  size_t fpv;
  
  fpv = LP_VertexList_FloatsPerVert(vl);
  arr = LP_VertexList_GetVert(vl);
  a = abp->vert_arr[(fa - arr) / fpv];
  b = abp->vert_arr[(fb - arr) / fpv];

  if (Hash_Lookup(a->pair_hash, b, NULL))
    return;
  
  if (Pair_New(abp->pairs, a, b) == NULL) {
    fprintf(stderr, "Could not create agg pair\n");
    abp->err = 1;
  }
}

static int Add_Agg_Pairs(struct ftree *pairs, struct vert **vert_arr, struct lp_vertex_list *vl, float aggregation_thresh) {
  struct bvh_vl *bvh;
  struct agg_bvh_pair abp;

  if ((bvh = BvhVl_New(vl, aggregation_thresh)) == NULL) {
    fprintf(stderr, "Could not create aggregation boundary volume hierarchy\n");
    goto err;
  }

  memset(&abp, 0, sizeof(abp));
  abp.pairs    = pairs;
  abp.vert_arr = vert_arr;
  
  BvhVl_Pairs(bvh, aggregation_thresh, Agg_Bvh_Pair, &abp);
  if (abp.err)
    goto err2;
  
  BvhVl_Free(bvh);
  return 0;
  
 err2:
  BvhVl_Free(bvh);
 err:
  return -1;
}

static int AllowedContraction(const struct pair *pair) {
  struct vert *a, *b;
  struct face **face_arr, *face;
  size_t flen, fcount;
  float orig[3], nnew[3];
  int vert;
  
  for (vert = 0; vert < 2; vert++) {
    a = pair->vert[vert];
    b = pair->vert[1 - vert];
    
    face_arr = (struct face **) Array_Data(a->face_arr);
    flen = Array_Length(a->face_arr);
    for (fcount = 0; fcount < flen; fcount++) {
      face = face_arr[fcount];
      if (face->vert[0] == b ||
	  face->vert[1] == b ||
	  face->vert[2] == b)
	/* Face will be removed.  Therefore, it cannot be inverted */
	continue;
      
      /* Check for normal inversion */
      PlaneNorm(orig, face->vert[0]->v, face->vert[1]->v, face->vert[2]->v);
      if (face->vert[0] == a) {
	PlaneNorm(nnew, pair->vbar, face->vert[1]->v, face->vert[2]->v);
      } else if (face->vert[1] == a) {
	PlaneNorm(nnew, face->vert[0]->v, pair->vbar, face->vert[2]->v);
      } else if (face->vert[2] == a) {
	PlaneNorm(nnew, face->vert[0]->v, face->vert[1]->v, pair->vbar);
      } else {
	printf("Warning: vertex not in face\n");
	continue;
      }
      if (Dot(nnew, orig) < 0)
	return 0;
    }
  }
  
  return 1;
}

static int Contract_Pair(struct ftree *pairs, struct hash *verts, struct hash *faces) {
  struct pair *pair, *pp;
  struct ftree_node *node;
  struct vert *a, *b, *c, *vv;
  struct hash_iterator *hi;
  struct face *face, **face_arr, **arr2;
  int count;
  float cost;
  size_t fcount, flen, fcount2, flen2;

  while (1) {
    if ((node = FTree_Lowest(pairs)) == NULL)
      return -1;
    
    pair = FTree_GetData(node);
    if (isinf(FTree_GetKey(node))) {
      fprintf(stderr, "Failure: All remianing pairs are disallowed\n");
      return -1;
    }
    
    if (AllowedContraction(pair))
      break;
    
    FTree_Rekey(pairs, node, INFINITY, NULL);
  }
  
  //printf("Contracting (%f, %f, %f) and (%f, %f, %f) to (%f, %f, %f)\n",
  //	   pair->vert[0]->v[0], pair->vert[0]->v[1], pair->vert[0]->v[2],
  //	   pair->vert[1]->v[0], pair->vert[1]->v[1], pair->vert[1]->v[2],
  //	   pair->vbar[0], pair->vbar[1], pair->vbar[2]);
  
  a = pair->vert[0];
  b = pair->vert[1];
  for (count = 0; count < 10; count++)
    a->Q[count] += b->Q[count];
  a->v[0] = pair->vbar[0];
  a->v[1] = pair->vbar[1];
  a->v[2] = pair->vbar[2];
  Hash_Remove(a->pair_hash, b);
  Hash_Remove(b->pair_hash, a);
  
  if ((hi = Hash_IteratorNew(a->pair_hash)) == NULL)
    return -1;
  while (Hash_IteratorNext(hi)) {
    pp = Hash_IteratorGetData(hi);
    cost = CalcLowestCost(pp);
    FTree_Rekey(pairs, pp->node, cost, NULL);
  }
  Hash_IteratorFree(hi);
  
  if ((hi = Hash_IteratorNew(b->pair_hash)) == NULL)
    return -1;
  while (Hash_IteratorNext(hi)) {
    pp = Hash_IteratorGetData(hi);
    vv = Hash_IteratorGetKey(hi);
    Hash_Remove(vv->pair_hash, b);
    if (Hash_Lookup(a->pair_hash, vv, NULL)) {
      FTree_Delete(pairs, pp->node);
      continue;
    }
    pp->vert[pp->vert[0] == b ? 0 : 1] = a;
    Hash_Insert(a->pair_hash, vv, pp, NULL);
    Hash_Insert(vv->pair_hash, a, pp, NULL);
    cost = CalcLowestCost(pp);
    FTree_Rekey(pairs, pp->node, cost, NULL);
  }
  Hash_IteratorFree(hi);
  
  face_arr = (struct face **) Array_Data(a->face_arr);
  flen = Array_Length(a->face_arr);
  fcount = 0;
  while (fcount < flen) {
    face = face_arr[fcount];
    if (face->vert[0] == b ||
	face->vert[1] == b ||
	face->vert[2] == b) {
      Array_Remove(a->face_arr, fcount);
      flen--;
    } else {
      fcount++;
    }
  }
  
  face_arr = (struct face **) Array_Data(b->face_arr);
  flen = Array_Length(b->face_arr);
  for (fcount = 0; fcount < flen; fcount++) {
    face = face_arr[fcount];
    if (face->vert[0] == a ||
	face->vert[1] == a ||
	face->vert[2] == a) {
      if (face->vert[0] != a && face->vert[0] != b)
	c = face->vert[0];
      else if (face->vert[1] != a && face->vert[1] != b)
	c = face->vert[1];
      else
	c = face->vert[2];
      
      arr2 = (struct face **) Array_Data(c->face_arr);
      flen2 = Array_Length(c->face_arr);
      for (fcount2 = 0; fcount2 < flen2; fcount2++) {
	if (arr2[fcount2] == face) {
	  Array_Remove(c->face_arr, fcount2);
	  break;
	}
      }
      
      Hash_Remove(faces, face);      
    } else {
      if (face->vert[0] == b) face->vert[0] = a;
      if (face->vert[1] == b) face->vert[1] = a;
      if (face->vert[2] == b) face->vert[2] = a;
      Face_Cannonize(face);
      Array_Add(a->face_arr, face);
    }
  }
  
  FTree_Delete(pairs, node);
  Hash_Remove(verts, b);
  
  return 0;
}

struct lp_vertex_list *LP_Simplify(const struct lp_vertex_list *in, size_t num_faces_out, float aggregation_thresh) {
  struct hash *faces, *verts;
  struct ftree *pairs;
  struct lp_vertex_list *vl, *out;
  struct hash_iterator *hi;
  struct face *face;
  struct vert **vert_arr, *vert[3];
  float *vv, *ii;
  int count;
  size_t cc, num, fpv;
  unsigned int idx, *arr;
  
  if (LP_VertexList_FloatsPerVert(in) < 3) {
    fprintf(stderr, "Error: Too few floats per vert to simplify\n");
    goto err;
  }
  
  if (LP_VertexList_PrimativeType(in) != lp_pt_triangle) {
    fprintf(stderr, "Error: Can only simplify triangular polyhedra\n");
    goto err;
  }
  
  if ((faces = Hash_NewPtr(NULL, Face_Free_Func, NULL, NULL, NULL)) == NULL)
    goto err;
  
  if ((verts = Hash_NewPtr(NULL, Vert_Free_Func, NULL, NULL, NULL)) == NULL)
    goto err2;
  
  if ((pairs = FTree_New(NULL, Pair_Free_Func, NULL)) == NULL)
    goto err3;
  
  if ((vert_arr = calloc(sizeof(*vert_arr), LP_VertexList_NumVert(in))) == NULL)
    goto err4;
  
  if ((vl = LP_VertexList_New(3, lp_pt_triangle)) == NULL)
    goto err5;
  
  fpv = LP_VertexList_FloatsPerVert(in);
  vv  = LP_VertexList_GetVert(in);
  arr = LP_VertexList_GetInd(in);
  num = LP_VertexList_NumInd(in) / 3;
  for (cc = 0; cc < num; cc++) {
    for (count = 0; count < 3; count++) {
      ii = &vv[fpv * arr[3 * cc + count]];
      if ((idx = LP_VertexList_Add(vl, ii)) == UINT_MAX)
	goto err6;
      if ((vert[count] = vert_arr[idx]) == NULL) {
	if ((vert[count] = Vert_New(verts, ii)) == NULL)
	  goto err6;
	vert_arr[idx] = vert[count];
      }
    }
    
    if (Face_New(faces, vert[0], vert[1], vert[2]) == NULL)
      goto err6;
  }
  
  if (Add_Pairs(pairs, faces) < 0)
    goto err6;
  
  if (aggregation_thresh > 0 && Add_Agg_Pairs(pairs, vert_arr, vl, aggregation_thresh) < 0) {
    fprintf(stderr, "Aggregation failed\n");
    goto err6;
  }
  
  printf("Simplifing polyhedron with %zu faces\n", Hash_NumEntries(faces));
  while (Hash_NumEntries(faces) > num_faces_out) {
    if (Contract_Pair(pairs, verts, faces) < 0) {
      fprintf(stderr, "Error: Unable to contract pair with %zu faces remaining\n", Hash_NumEntries(faces));
      break;
    }
  }
  
  if ((out = LP_VertexList_New(3, lp_pt_triangle)) == NULL)
    goto err6;
  
  if ((hi = Hash_IteratorNew(faces)) == NULL)
    goto err7;
  
  while (Hash_IteratorNext(hi)) {
    face = (struct face *) Hash_IteratorGetKey(hi);
    for (count = 0; count < 3; count++) {
      LP_VertexList_Add(out, face->vert[count]->v);
    }
  }
  Hash_IteratorFree(hi);
  
  LP_VertexList_Free(vl);
  free(vert_arr);
  FTree_Free(pairs);
  Hash_Free(verts);
  Hash_Free(faces);
  return out;
  
 err7:
  LP_VertexList_Free(out);
 err6:
  LP_VertexList_Free(vl);
 err5:
  free(vert_arr);
 err4:
  FTree_Free(pairs);
 err3:
  Hash_Free(verts);
 err2:
  Hash_Free(faces);
 err:
  return NULL;
}
