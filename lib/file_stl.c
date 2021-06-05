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
#include <string.h>

#include "file_stl.h"
#include "util.h"

static inline int IsLittleEndian(void) {
  const union {uint16_t i; unsigned char c[2];} one = {1};
  
  return one.c[0];
}

static void MakeLittleInt16(uint16_t *val) {
  if (IsLittleEndian())
    return;
  
  *val = (*val >> 16) | (*val << 16);
}

static void MakeLittleInt32(uint32_t *val) {
  if (IsLittleEndian())
    return;
  
  *val =
    ((*val >> 24)) |
    ((*val >>  8) & 0xFF00) |
    ((*val <<  8) & 0xFF0000) |
    ((*val << 24));
}

static void MakeLittleFloat(float *val) {
  MakeLittleInt32((uint32_t *) val);
}

struct face {
  float norm[3];
  float v[9];
};

static void MakeLittleFace(struct face *face) {
  MakeLittleFloat(&face->norm[0]);
  MakeLittleFloat(&face->norm[1]);
  MakeLittleFloat(&face->norm[2]);
  MakeLittleFloat(&face->v[0]);
  MakeLittleFloat(&face->v[1]);
  MakeLittleFloat(&face->v[2]);
  MakeLittleFloat(&face->v[3]);
  MakeLittleFloat(&face->v[4]);
  MakeLittleFloat(&face->v[5]);
  MakeLittleFloat(&face->v[6]);
  MakeLittleFloat(&face->v[7]);
  MakeLittleFloat(&face->v[8]);
}

static void FixWindingOrder(struct face *face) {
  float norm[3], temp[3];
  
  PlaneNorm(norm, &face->v[0], &face->v[3], &face->v[6]);
  if (Dot(norm, face->norm) >= 0)
    return;
  
  temp[0] = face->v[3];
  temp[1] = face->v[4];
  temp[2] = face->v[5];
  face->v[3] = face->v[6];
  face->v[4] = face->v[7];
  face->v[5] = face->v[8];
  face->v[6] = temp[0];
  face->v[7] = temp[1];
  face->v[8] = temp[2];
}

static int ReadBinaryStl(FILE *in, struct lp_vertex_list *vl, float scale) {
  char head[74];
  uint32_t num_faces, count;
  uint16_t attr_bytes, attr_count;
  struct face face;
  float ff[6];
  int vert;
  
  if (fread(head, sizeof(head), 1, in) != 1) {
    fprintf(stderr, "Error: Unable to read stl header(2)\n");
    return -1;
  }
  
  if (fread(&num_faces, sizeof(num_faces), 1, in) != 1) {
    fprintf(stderr, "Error: Unable to read number of faces\n");
    return -1;
  }
  
  MakeLittleInt32(&num_faces);
  
  for (count = 0; count < num_faces; count++) {
    if (fread(&face, sizeof(face), 1, in) != 1) {
      fprintf(stderr, "Error: Unable to read face %lu\n", (unsigned long) count);
      return -1;
    }
    
    MakeLittleFace(&face);
    FixWindingOrder(&face);
    
    for (vert = 0; vert < 3; vert++) {
      ff[0] = face.v[3 * vert    ] * scale;
      ff[1] = face.v[3 * vert + 1] * scale;
      ff[2] = face.v[3 * vert + 2] * scale;
      ff[3] = face.norm[0];
      ff[4] = face.norm[1];
      ff[5] = face.norm[2];
      
      if (LP_VertexList_Add(vl, ff) == UINT_MAX)
	return -1;
    }
    
    if (fread(&attr_bytes, sizeof(attr_bytes), 1, in) != 1) {
      fprintf(stderr, "Error: Unable to read face %lu attribute size\n", (unsigned long) count);
      return -1;
    }
    
    if (attr_bytes) {
      MakeLittleInt16(&attr_bytes);
      for (attr_count = 0; attr_count < attr_bytes; attr_count++) {
	if (fread(head, 1, 1, in) != 1) {
	  fprintf(stderr, "Error: Unable to read face %lu attribute byte %u\n", (unsigned long) count, (int) attr_count);
	}
      }
    }
  }
  
  return 0;
}

static int ReadAsciiStl(FILE *in, struct lp_vertex_list *vl, float scale) {
  fprintf(stderr, "Error: ASCII .stl not yet supported\n");
  return -1;
}

static struct lp_vertex_list *FileStl_ReadSingle(FILE *in, float scale) {
  struct lp_vertex_list *vl;
  char head[6];
  
  if ((vl = LP_VertexList_New(6, lp_pt_triangle)) == NULL) {
    fprintf(stderr, "Error: Could not allocate memory for vertex list");
    goto err;
  }
  
  if (fread(head, sizeof(head), 1, in) != 1) {
    fprintf(stderr, "Error: Unable to read stl header\n");
    goto err2;
  }
  
  if (memcmp(head, "solid ", 6) == 0) {
    if (ReadAsciiStl(in, vl, scale) < 0)
      goto err2;
  } else {
    if (ReadBinaryStl(in, vl, scale) < 0)
      goto err2;
  }
  
  return vl;
  
 err2:
  LP_VertexList_Free(vl);
 err:
  return NULL;
}

struct lp_vl_list *FileStl_Read(FILE *in, float scale) {
  return LP_VertexList_ListAppend(NULL, FileStl_ReadSingle(in, scale));
}

const uint16_t zero_attr = 0;

static int FileStl_WriteSingle(FILE *out, const struct lp_vertex_list *vl, float scale) {
  struct face face;
  size_t count, num;
  uint32_t num_tri;
  char head[80];
  float *ff;

  if (LP_VertexList_FloatsPerVert(vl) < 3) {
    fprintf(stderr, "Error: Too few floats per vert for .stl file\n");
    return -1;
  }
  if (LP_VertexList_PrimativeType(vl) != lp_pt_triangle) {
    fprintf(stderr, "Error: wrong primative type for .stl file\n");
    return -1;
  }
  num = LP_VertexList_NumInd(vl);
  
  memset(head, 0, sizeof(head));
  strncpy(head, "binary stl libpolyhedra\n", sizeof(head));
  if (fwrite(head, sizeof(head), 1, out) != 1)
    return -1;
  
  num_tri = num / 3;
  MakeLittleInt32(&num_tri);
  if (fwrite(&num_tri, sizeof(num_tri), 1, out) != 1)
    return -1;
  
  for (count = 0; count < num / 3; count++) {
    ff = LP_VertexList_LookupVert(vl, 3 * count);
    face.v[0] = ff[0] * scale;
    face.v[1] = ff[1] * scale;
    face.v[2] = ff[2] * scale;
    ff = LP_VertexList_LookupVert(vl, 3 * count + 1);
    face.v[3] = ff[0] * scale;
    face.v[4] = ff[1] * scale;
    face.v[5] = ff[2] * scale;
    ff = LP_VertexList_LookupVert(vl, 3 * count + 2);
    face.v[6] = ff[0] * scale;
    face.v[7] = ff[1] * scale;
    face.v[8] = ff[2] * scale;
    PlaneNorm(face.norm, &face.v[0], &face.v[3], &face.v[6]);
    MakeLittleFace(&face);
    if (fwrite(&face, sizeof(face), 1, out) != 1)
      return -1;
    if (fwrite(&zero_attr, sizeof(zero_attr), 1, out) != 1)
      return -1;
  }
  
  return 0;
}

int FileStl_Write(FILE *out, const struct lp_vl_list *list, float scale) {
  if (list == NULL || list->next != NULL) {
    fprintf(stderr, "Error: STL supports exactly one mesh per file\n");
    return -1;
  }
  
  return FileStl_WriteSingle(out, list->vl, scale);
}
