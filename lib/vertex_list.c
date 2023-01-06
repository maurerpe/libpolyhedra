/*
  Copyright (C) 2020-21 Paul Maurer

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
#include <stddef.h>
#include <string.h>

#include "libpolyhedra.h"

#include "file_obj.h"
#include "file_stl.h"
#include "file_svg.h"
#include "hash.h"
#include "SipHash/siphash.h"

#define PRESENT ((void *) 1)

struct lp_vertex_list {
  size_t floats_per_vert;
  size_t vert_size;
  enum primative_type primative_type;

  unsigned int vert_alloc;
  unsigned int vert_used;
  float *vert;

  size_t ind_alloc;
  size_t ind_used;
  unsigned int *ind;

  struct hash *vert_hash;
};

/* Key stored in hash is (void *) index + 1 */
/* Note that index cannot be used because it can be zero, which is NULL */
static uint64_t VlHash(const void *user, const unsigned char secret[16], const void *key) {
  struct lp_vertex_list *vl = (struct lp_vertex_list *) user;
  
  return siphash(secret, (const unsigned char *) key, vl->vert_size);
}

static int VlCmp(const void *user, const void *key_a, const void *key_b) {
  struct lp_vertex_list *vl = (struct lp_vertex_list *) user;
  float *kk = vl->vert + (((unsigned int) (ptrdiff_t) key_a) - 1) * vl->floats_per_vert;

  return memcmp(kk, key_b, vl->vert_size);
}

static void *VlCopy(void *user, const void *key) {
  struct lp_vertex_list *vl = (struct lp_vertex_list *) user;
  size_t new_alloc;
  float *new_mem, *vv;
  
  if (vl->vert_used >= vl->vert_alloc) {
    if (vl->vert_alloc == UINT_MAX) {
      fprintf(stderr, "Error: Too many vertices in a single vertex list\n");
      goto err;
    }

    if (vl->vert_alloc > (UINT_MAX >> 1))
      new_alloc = UINT_MAX;
    else
      new_alloc = vl->vert_alloc << 1;
    
    if (SIZE_MAX / vl->vert_size <= new_alloc) {
      fprintf(stderr, "Error: Out of memory adding vertex\n");
      goto err;
    }

    if ((new_mem = realloc(vl->vert, new_alloc * vl->vert_size)) == NULL) {
      fprintf(stderr, "Error: Out of memory adding vertex\n");
      goto err;
    }
    
    vl->vert = new_mem;
    vl->vert_alloc = new_alloc;
  }

  vv = vl->vert + vl->vert_used * vl->floats_per_vert;
  memcpy(vv, key, vl->vert_size);
  return (void *) (ptrdiff_t) ++vl->vert_used;

 err:
  return NULL;
}

struct lp_vertex_list *LP_VertexList_New(size_t floats_per_vert, enum primative_type pt) {
  struct lp_vertex_list *vl;

  if (SIZE_MAX / sizeof(float) < floats_per_vert) {
    fprintf(stderr, "Error: Too many floats in a vertex\n");
    goto err;
  }
  
  if ((vl = malloc(sizeof(*vl))) == NULL) {
    perror("Error: Could not allocate vertex list");
    goto err;
  }
  memset(vl, 0, sizeof(*vl));

  vl->floats_per_vert = floats_per_vert;
  vl->primative_type  = pt;
  vl->vert_size = floats_per_vert * sizeof(float);
  vl->vert_alloc = 64;
  vl->ind_alloc  = 64;

  if ((vl->vert = calloc(vl->vert_alloc, vl->vert_size)) == NULL) {
    perror("Error: Could not allocate vertexes");
    goto err2;
  }

  if ((vl->ind = calloc(vl->ind_alloc, sizeof(unsigned int))) == NULL) {
    perror("Error: Could not allocate vertex indices");
    goto err3;
  }

  if ((vl->vert_hash = Hash_New(vl, VlHash, VlCmp, VlCopy, NULL, NULL, NULL, NULL)) == NULL)
    goto err4;
  
  return vl;

 err4:
  free(vl->ind);
 err3:
  free(vl->vert);
 err2:
  free(vl);
 err:
  fprintf(stderr, "Error: Could not allocate memory for vertex list\n");
  return NULL;
}

void LP_VertexList_Free(struct lp_vertex_list *vl) {
  if (vl == NULL)
    return;

  Hash_Free(vl->vert_hash);
  free(vl->ind);
  free(vl->vert);
  free(vl);
}

void LP_VertexList_Clear(struct lp_vertex_list *vl) {
  vl->vert_used = 0;
  vl->ind_used = 0;
  
  Hash_Clear(vl->vert_hash);
}

struct lp_vertex_list *LP_VertexList_Copy(const struct lp_vertex_list *vl, size_t new_floats_per_vert) {
  struct lp_vertex_list *out;
  size_t count, num, fpv;
  
  if ((fpv = LP_VertexList_FloatsPerVert(vl)) < new_floats_per_vert) {
    fprintf(stderr, "Error too few vertices to copy\n");
    goto err;
  }
  
  if (new_floats_per_vert == SIZE_MAX)
    new_floats_per_vert = fpv;
  
  if ((out = LP_VertexList_New(new_floats_per_vert, LP_VertexList_PrimativeType(vl))) == NULL)
    goto err;
  
  num = LP_VertexList_NumInd(vl);
  for (count = 0; count < num; count++)
    if (LP_VertexList_Add(out, LP_VertexList_LookupVert(vl, count)) == UINT_MAX)
      goto err2;
  
  return out;
  
 err2:
  LP_VertexList_Free(out);
 err:
  return NULL;
}

static unsigned int AddVert(struct lp_vertex_list *vl, const float *vert) {
  void *key_out;
  
  if (Hash_Insert(vl->vert_hash, vert, PRESENT, &key_out) < 0) {
    fprintf(stderr, "Error: Could not add vertex to hash\n");
    return UINT_MAX;
  }
  
  return ((ptrdiff_t) key_out) - 1;
}

static unsigned int AddInd(struct lp_vertex_list *vl, unsigned int ind) {
  size_t new_alloc;
  unsigned int *new_mem;
  
  if (vl->ind_used >= vl->ind_alloc) {
    if (((SIZE_MAX / sizeof(unsigned int)) >> 1) < vl->ind_alloc) {
      fprintf(stderr, "Error: Too many indices in a single vertex list\n");
      goto err;
    }
    
    new_alloc = vl->ind_alloc << 1;
    if ((new_mem = realloc(vl->ind, new_alloc * sizeof(unsigned int))) == NULL) {
      fprintf(stderr, "Error: Out of memory adding vertex index\n");
      goto err;
    }
    
    vl->ind = new_mem;
    vl->ind_alloc = new_alloc;
  }
  
  vl->ind[vl->ind_used] = ind;
  vl->ind_used++;
  return ind;
  
 err:
  return UINT_MAX;
}

unsigned int LP_VertexList_Add(struct lp_vertex_list *vl, const float *vert) {
  unsigned int ind;
  
  if ((ind = AddVert(vl, vert)) == UINT_MAX)
    return UINT_MAX;
  
  return LP_VertexList_AddIndex(vl, ind);
}

unsigned int LP_VertexList_AddIndex(struct lp_vertex_list *vl, unsigned int index) {
  if (index > vl->vert_used) {
    printf("Error: Vertex index is out of range: %u, %u\n", index, vl->vert_used);
    return UINT_MAX;
  }
  
  return AddInd(vl, index);
}

void LP_VertexList_Finalize(struct lp_vertex_list *vl) {
  Hash_Free(vl->vert_hash);
  vl->vert_hash = NULL;
}

size_t LP_VertexList_FloatsPerVert(const struct lp_vertex_list *vl) {
  return vl->floats_per_vert;
}

enum primative_type LP_VertexList_PrimativeType(const struct lp_vertex_list *vl) {
  return vl->primative_type;
}

unsigned int LP_VertexList_NumVert(const struct lp_vertex_list *vl) {
  return vl->vert_used;
}

size_t LP_VertexList_NumInd(const struct lp_vertex_list *vl) {
  return vl->ind_used;
}

float *LP_VertexList_GetVert(const struct lp_vertex_list *vl) {
  return vl->vert;
}

unsigned int *LP_VertexList_GetInd(const struct lp_vertex_list *vl) {
  return vl->ind;
}

float *LP_VertexList_LookupVert(const struct lp_vertex_list *vl, unsigned int index) {
  return vl->vert + vl->floats_per_vert * vl->ind[index];
}

/********************** VertexList lists ***************************/
struct lp_vl_list *LP_VertexList_ListAppend(struct lp_vl_list *list, struct lp_vertex_list *vl) {
  struct lp_vl_list **end, *new;
  
  if (vl == NULL)
    return list;
  
  end = &list;
  while (*end)
    end = &(*end)->next;
  
  if ((new = malloc(sizeof(*list))) == NULL)
    return NULL;
  memset(new, 0, sizeof(*list));
  
  new->vl = vl;
  *end = new;
  
  return list;
}

struct lp_vl_list *LP_VertexList_ListJoin(struct lp_vl_list *list1, struct lp_vl_list *list2) {
  struct lp_vl_list **tail;
  
  for (tail = &list1; *tail != NULL; tail = &(*tail)->next)
    ;
  
  *tail = list2;
  
  return list1;
}

size_t LP_VertexList_ListLength(struct lp_vl_list *list) {
  size_t count = 0;
  
  while (list != NULL) {
    list = list->next;
    count++;
  }
  
  return count;
}

void LP_VertexList_ListFree(struct lp_vl_list *list) {
  struct lp_vl_list *next;

  while (list) {
    next = list->next;
    LP_VertexList_Free(list->vl);
    free(list);
    list = next;
  }
}

/********************** Read & Write *******************************/

enum file_type {
  ft_obj,
  ft_stl,
  ft_svg
};

static enum file_type FileType(const char *filename) {
  size_t len;

  len = strlen(filename);
  
  if (len > 4 && strcasecmp(filename + len - 4, ".obj") == 0)
    return ft_obj;

  if (len > 4 && strcasecmp(filename + len - 4, ".stl") == 0)
    return ft_stl;

  if (len > 4 && strcasecmp(filename + len - 4, ".svg") == 0)
    return ft_svg;

  return UINT_MAX;
}

struct lp_vl_list *LP_VertexList_Read(const char *filename, float scale) {
  FILE *in;
  struct lp_vl_list *list, *cur;
  enum file_type ft;

  if ((ft = FileType(filename)) == UINT_MAX) {
    fprintf(stderr, "Error: Unkown mesh format '%s', must be .obj, .stl, or .svg\n", filename);
    goto err;
  }
  
  printf("Reading mesh(es) from '%s'\n", filename);
  if ((in = fopen(filename, "r")) == NULL) {
    perror("Error: Could not open file for reading");
    goto err;
  }

  switch (ft) {
  case ft_obj: list = FileObj_Read(in, scale); break;
  case ft_stl: list = FileStl_Read(in, scale); break;
  case ft_svg: list = FileSvg_Read(in, scale); break;
  }
  
  if (list == NULL) {
    fprintf(stderr, "Error: No polyhedra returned from file read\n");
    goto err2;
  }
  
  fclose(in);
  for (cur = list; cur != NULL; cur = cur->next)
    LP_VertexList_Finalize(cur->vl);
  return list;

 err2:
  fclose(in);
 err:
  return NULL;
}

int LP_VertexList_Write(const char *filename, struct lp_vl_list *list, float scale) {
  FILE *out;
  enum file_type ft;
  int ret;

  if ((ft = FileType(filename)) == UINT_MAX) {
    fprintf(stderr, "Error: Unkown mesh format '%s', must be .obj or .stl\n", filename);
    goto err;
  }
  
  printf("Writing %zu mesh(es) to '%s'\n", LP_VertexList_ListLength(list), filename);
  if ((out = fopen(filename, "w")) == NULL) {
    perror("Error: Could not open file for writing");
    goto err;
  }
  
  switch (ft) {
  case ft_obj: ret = FileObj_Write(out, list, scale); break;
  case ft_stl: ret = FileStl_Write(out, list, scale); break;
  case ft_svg: ret = FileSvg_Write(out, list, scale); break;
  }
  
  if (ret < 0) {
    fprintf(stderr, "Error: Could not write polyhedra to file\n");
    goto err2;
  }
  
  fclose(out);
  return 0;

 err2:
  fclose(out);
 err:
  return -1;
}
