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

#include "file_obj.h"

enum obj_state {
  o_firstword,
  o_comment,
  o_floatspace,
  o_float,
  o_intspace,
  o_int,
  o_name
};

enum obj_type {
  t_v,
  t_vt,
  t_vn,
  t_f
};

static int ObjPushFloat(int type, size_t count, struct lp_vertex_list *v, struct lp_vertex_list *vn, struct lp_vertex_list *vt, float *ff, size_t line, size_t col) {
  if (type == t_v) {
    if (count < 3) {
      fprintf(stderr, "Error: Line %zu, column %zu: too few floating point numbers, expected 3\n", line, col);
      return -1;
    }
    
    if (LP_VertexList_Add(v, ff) == UINT_MAX)
      return -1;
    
    return 0;
  }
  
  if (type == t_vn) {
    if (count < 3) {
      fprintf(stderr, "Error: Line %zu, column %zu: too few floating point numbers, expected 3\n", line, col);
      return -1;
    }

    if (LP_VertexList_Add(vn, ff) == UINT_MAX)
      return -1;
    
    return 0;
  }
  
  if (type == t_vt) {
    if (count < 2) {
      fprintf(stderr, "Error: Line %zu, column %zu: too few floating point numbers, expected 2\n", line, col);
      return -1;
    }
    
    if (LP_VertexList_Add(vt, ff) == UINT_MAX)
      return -1;
    
    return 0;
  }

  fprintf(stderr, "Internal Error: file_obj.c: Invalid float category\n");
  return -1;
}

static int ObjPushVert(struct lp_vertex_list *vl, struct lp_vertex_list *v, struct lp_vertex_list *vn, struct lp_vertex_list *vt, unsigned long long *ii, size_t subcount, int has_n, int has_t, size_t line, size_t col) {
  float ff[8], *f, *cur = ff;
  int count, exp_sub = 0;

  if (has_n)
    exp_sub++;

  if (has_t)
    exp_sub++;
  
  if (subcount != exp_sub) {
    if (has_t) {
      if (has_n)
	fprintf(stderr, "Error: Line %zu, column %zu: each face vertex needs a vertex, a normal, and a uv\n", line, col);
      else
	fprintf(stderr, "Error: Line %zu, column %zu: each face vertex needs a vertex and a uv\n", line, col);
    } else {
      if (has_n)
	fprintf(stderr, "Error: Line %zu, column %zu: each face vertex needs a vertex and a normal\n", line, col);
      else
	fprintf(stderr, "Error: Line %zu, column %zu: each face vertex needs a vertex and no other values\n", line, col);
    }
    
    return -1;
  }
  
  if (ii[0] == 0 || ii[0] > LP_VertexList_NumInd(v)) {
    fprintf(stderr, "Error: Line %zu, column %zu: Vertex index out of range (1 - %zu): %llu\n", line, col, LP_VertexList_NumInd(v), ii[0]);
    return -1;
  }
  
  f = LP_VertexList_LookupVert(v, ii[0] - 1);
  *cur++ = *f++;
  *cur++ = *f++;
  *cur++ = *f;

  if (has_n) {
    count = has_t ? 2 : 1;
    if (ii[count] == 0 || ii[count] > LP_VertexList_NumInd(vn)) {
      fprintf(stderr, "Error: Line %zu, column %zu: Normal index out of range (1 - %zu): %llu\n", line, col, LP_VertexList_NumInd(vn), ii[count]);
      return -1;
    }
  
    f = LP_VertexList_LookupVert(vn, ii[count] - 1);
    *cur++ = *f++;
    *cur++ = *f++;
    *cur++ = *f;
  }

  if (has_t) {
    if (ii[1] == 0 || ii[1] > LP_VertexList_NumInd(vt)) {
      fprintf(stderr, "Error: Line %zu, column %zu: UV index out of range (1 - %zu): %llu\n", line, col, LP_VertexList_NumInd(vt), ii[1]);
      return -1;
    }
    
    f = LP_VertexList_LookupVert(vt, ii[1] - 1);
    *cur++ = *f++;
    *cur++ = *f;
  }
  
  if (LP_VertexList_Add(vl, ff) == UINT_MAX)
    return -1;
  
  return 0;
}

struct file_data {
  char buf[1024];
  char *cur;
  char *end;
  size_t line;
  size_t col;
  int prev_was_cr;
  int err;
};

static struct lp_vertex_list *FileObj_ReadSingle(FILE *in, float scale, struct lp_vertex_list *v, struct lp_vertex_list *vn, struct lp_vertex_list *vt, struct file_data *fd) {
  char str[80], *curst, ch;
  size_t len, count = 0, subcount = 0;
  int state, type = t_v, has_n = 0, has_t = 0, done = 0;
  struct lp_vertex_list *vl = NULL;
  float ff[3];
  unsigned long long ii[3];
  
  state = o_firstword;
  curst = str;
  while (!done) {
    if (++fd->cur >= fd->end) {
      if ((len = fread(fd->buf, 1, sizeof(fd->buf), in)) == 0) {
	if (ferror(in)) {
	  fprintf(stderr, "Error: Cannot read from file\n");
	  goto err;
	}
	break;
      }
      fd->cur = fd->buf;
      fd->end = fd->buf + len;
    }

    ch = *fd->cur;
    fd->col++;
    
    switch (state) {
    case o_firstword:
      if (ch == '\n' || ch == '\r') {
	curst = str;
	break;
      }
      if (ch == '#') {
	state = o_comment;
	break;
      }
      if (ch == ' ') {
	if (curst == str)
	  break;
	*curst = '\0';
	count = 0;
	subcount = 0;
	if (strcmp(str, "v") == 0) {
	  if (vl) {
	    fprintf(stderr, "Error: v entries must be before f entries\n");
	    goto err;
	  }
	  type = t_v;
	  state = o_floatspace;
	} else if (strcmp(str, "vt") == 0) {
	  if (vl) {
	    fprintf(stderr, "Error: vt entries must be before f entries\n");
	    goto err;
	  }
	  type = t_vt;
	  state = o_floatspace;
	  has_t = 1;
	} else if (strcmp(str, "vn") == 0) {
	  if (vl) {
	    fprintf(stderr, "Error: vn entries must be before f entries\n");
	    goto err;
	  }
	  type = t_vn;
	  state = o_floatspace;
	  has_n = 1;
	} else if (strcmp(str, "f") == 0) {
	  if (vl == NULL && (vl = LP_VertexList_New(3 + (has_n ? 3 : 0) + (has_t ? 2 : 0), lp_pt_triangle)) == NULL)
	    goto err;
	  type = t_f;
	  state = o_intspace;
	} else if (strcmp(str, "o") == 0) {
	  state = o_name;
	} else {/* unsupported type */
	  state = o_comment;
	}
	break;
      }
      if (curst >= str + sizeof(str) - 1) {
	state = o_comment;
	break;
      }
      
      *curst++ = ch;
      break;

    case o_comment:
      if (ch == '\n' || ch == '\r') {
	state = o_firstword;
	curst = str;
      }
      break;

    case o_floatspace:
      if (ch == ' ')
	break;

      if (ch == '\n' || ch == '\r') {
	if (ObjPushFloat(type, count, v, vn, vt, ff, fd->line, fd->col) < 0) {
	  fprintf(stderr, "Error: Line %zu, column %zu: Could not push float at end of line\n", fd->line, fd->col);
	  goto err;
	}
	state = o_firstword;
	curst = str;
	break;
      }
      
      if ((type == t_v || type == t_vn) && count > 3) {
	fprintf(stderr, "Error: Line %zu, column %zu: too many floating point numbers, expected 3\n", fd->line, fd->col);
	goto err;
      }
      if (type == t_vt && count > 2) {
	fprintf(stderr, "Error: Line %zu, column %zu: too many floating point numbers, expected 2\n", fd->line, fd->col);
	goto err;
      }
      state = o_float;
      curst = str;
      *curst++ = ch;
      break;
      
    case o_float:
      if (ch == '\n' || ch == '\r' || ch == ' ') {
	*curst = '\0';
	ff[count] = strtod(str, &curst);
	if (type == t_v)
	  ff[count] *= scale;
	if (type == t_vt)
	  ff[count] = 1.0 - ff[count];
	count++;
	if (*curst != '\0') {
	  fprintf(stderr, "Error: Line %zu, column %zu: invalid floating point number: %s\n", fd->line, fd->col, str);
	  goto err;
	}
	if (ch == ' ') {
	  state = o_floatspace;
	  break;
	}

	if (ObjPushFloat(type, count, v, vn, vt, ff, fd->line, fd->col) < 0) {
	  fprintf(stderr, "Error: Line %zu, column %zu: Could not push float\n", fd->line, fd->col);
	  goto err;
	}
	state = o_firstword;
	curst = str;
	break;
      }
      
      if (curst >= str + sizeof(str) - 1) {
	fprintf(stderr, "Error: Line %zu, column %zu: floating point number too long\n", fd->line, fd->col);
	goto err;
      }
      
      *curst++ = ch;
      break;
      
    case o_intspace:
      if (ch == ' ')
	break;

      if (ch == '\n' || ch == '\r') {
	if (count != 3) {
	  fprintf(stderr, "Error: Line %zu, column %zu: incorrect number of vertices for face, expected 3\n", fd->line, fd->col);
	  goto err;
	}

	state = o_firstword;
	curst = str;
	break;
      }

      subcount = 0;
      count++;
      if (count > 3) {
	fprintf(stderr, "Error: Line %zu, column %zu: incorrect number of vertices for face, expected 3 (only triangular faces are supported)\n", fd->line, fd->col);
	goto err;
      }
      state = o_int;
      curst = str;
      *curst++ = ch;
      break;
      
    case o_int:
      if (ch == '\n' || ch == '\r' || ch == ' ' || ch == '/') {
	*curst = '\0';
	ii[subcount] = strtoull(str, &curst, 0);
	if (*curst != '\0') {
	  fprintf(stderr, "Error: Line %zu, column %zu: invalid integer\n", fd->line, fd->col);
	  goto err;
	}

	if (ch == '/') {
	  subcount++;
	  curst = str;
	  break;
	}
	
	if (ObjPushVert(vl, v, vn, vt, ii, subcount, has_n, has_t, fd->line, fd->col) < 0) {
	  fprintf(stderr, "Error: Line %zu, column %zu: Could not push vertex\n", fd->line, fd->col);
	  goto err;
	}
	
	if (ch == ' ') {
	  state = o_intspace;
	  break;
	}
	
	if (count != 3) {
	  fprintf(stderr, "Error: Line %zu, column %zu: incorrect number of vertices for face, expected 3\n", fd->line, fd->col);
	  goto err;
	}
	state = o_firstword;
	curst = str;
	break;
      }
      
      if (curst >= str + sizeof(str) - 1) {
	fprintf(stderr, "Error: Line %zu, column %zu: integer too long\n", fd->line, fd->col);
	goto err;
      }
      
      *curst++ = ch;
      break;
          
    case o_name:
      if (ch == '\n' || ch == '\r')
	done = 1;
      break;
    }
    
    if (ch == '\r') {
      fd->line++;
      fd->col = 0;
      fd->prev_was_cr = 1;
    } else {
      if (ch == '\n' && !fd->prev_was_cr) {
	fd->line++;
	fd->col = 0;
      }
      fd->prev_was_cr = 0;
    }
  }
  
  /* End of file */
  switch (state) {
  case o_firstword:
  case o_comment:
  case o_floatspace:
  case o_float:
  case o_name:
    break;
    
  case o_intspace:
    if (count != 3) {
      fprintf(stderr, "Error: Line %zu, column %zu: incorrect number of vertices for face, expected 3\n", fd->line, fd->col);
      goto err;
    }
    break;
    
  case o_int:
    *curst = '\0';
    ii[subcount] = strtoull(str, &curst, 0);
    if (*curst != '\0') {
      fprintf(stderr, "Error: Line %zu, column %zu: invalid integer\n", fd->line, fd->col);
      goto err;
    }
    
    if (ObjPushVert(vl, v, vn, vt, ii, subcount, has_n, has_t, fd->line, fd->col) < 0) {
      fprintf(stderr, "Error: Line %zu, column %zu: Could not push final vertex\n", fd->line, fd->col);
      goto err;
    }
    break;
  }

  return vl;
  
 err:
  LP_VertexList_Free(vl);
  fd->err = 1;
  fprintf(stderr, "Error: Line %zu, column %zu: Could not parse .obj file\n", fd->line, fd->col);
  return NULL;
}

struct lp_vl_list *FileObj_Read(FILE *in, float scale) {
  struct lp_vl_list *list = NULL;
  struct lp_vertex_list *v, *vn, *vt;
  struct file_data fd;
  
  if ((v = LP_VertexList_New(3, lp_pt_point)) == NULL)
    goto err;

  if ((vn = LP_VertexList_New(3, lp_pt_unspecified)) == NULL)
    goto err2;

  if ((vt = LP_VertexList_New(2, lp_pt_unspecified)) == NULL)
    goto err3;
  
  memset(&fd, 0, sizeof(fd));
  fd.cur = fd.end = fd.buf;
  fd.line = 1;
  
  while (!feof(in) || fd.cur < fd.end) {
    list = LP_VertexList_ListAppend(list, FileObj_ReadSingle(in, scale, v, vn, vt, &fd));
    if (fd.err)
      goto err4;
  }
  
  LP_VertexList_Free(vt);
  LP_VertexList_Free(vn);
  LP_VertexList_Free(v);
  return list;

 err4:
  LP_VertexList_Free(vt);
 err3:
  LP_VertexList_Free(vn);
 err2:
  LP_VertexList_Free(v);
 err:
  LP_VertexList_ListFree(list);
  return NULL;
}

struct wface {
  size_t v;
  size_t vn;
  size_t vt;
};

static int FileObj_WriteSingle(FILE *out, size_t poly_count, const struct lp_vertex_list *vl, float scale, size_t *v_off, size_t *vn_off, size_t *vt_off) {
  struct lp_vertex_list *v, *vn, *vt;
  struct wface *wf, *wfmem;
  size_t fpv, count, num, num_verts;
  float *ff;
  int has_vn, has_vt, face;

  fpv = LP_VertexList_FloatsPerVert(vl);
  num = LP_VertexList_NumInd(vl);
  
  if (fpv < 3) {
    fprintf(stderr, "Error: Two few floats per vert to write .obj file\n");
    goto err;
  }
  
  if (LP_VertexList_PrimativeType(vl) != lp_pt_triangle) {
    fprintf(stderr, "Error: Incorrect primative type for .obj\n");
    goto err;
  }
  
  has_vn = fpv == 5 || fpv == 8;
  has_vt = fpv == 6 || fpv == 8;
  
  if ((wf = calloc(num, sizeof(*wf))) == NULL)
    goto err;
  wfmem = wf;
  
  if ((v = LP_VertexList_New(3, lp_pt_point)) == NULL)
    goto err2;
  
  if ((vn = LP_VertexList_New(3, lp_pt_unspecified)) == NULL)
    goto err3;

  if ((vt = LP_VertexList_New(2, lp_pt_unspecified)) == NULL)
    goto err4;
  
  for (count = 0; count < num; count++) {
    ff = LP_VertexList_LookupVert(vl, count);
    if ((wf[count].v = LP_VertexList_Add(v, &ff[0])) == UINT_MAX)
      goto err5;
    if (has_vn)
      if ((wf[count].vn = LP_VertexList_Add(vn, &ff[3])) == UINT_MAX)
	goto err5;
    if (has_vt)
      if ((wf[count].vt = LP_VertexList_Add(vt, &ff[has_vn ? 6 : 3])) == UINT_MAX)
	goto err5;
  }
  
  fprintf(out, "o polyhedra.%03zu\n", poly_count);
  
  num_verts = LP_VertexList_NumVert(v);
  ff = LP_VertexList_GetVert(v);
  for (count = 0; count < num_verts; count++, ff += 3)
    fprintf(out, "v %f %f %f\n", ff[0] * scale, ff[1] * scale, ff[2] * scale);
  
  num_verts = LP_VertexList_NumVert(vt);
  ff = LP_VertexList_GetVert(vt);
  for (count = 0; count < num_verts; count++, ff += 2)
    fprintf(out, "vt %f %f\n", ff[0], ff[1]);
  
  num_verts = LP_VertexList_NumVert(vn);
  ff = LP_VertexList_GetVert(vn);
  for (count = 0; count < num_verts; count++, ff += 3)
    fprintf(out, "vn %f %f %f\n", ff[0], ff[1], ff[2]);
  
  for (count = 0; count < num / 3; count++) {
    fprintf(out, "f");
    for (face = 0; face < 3; face++, wf++) {
      fprintf(out, " %zu", wf->v + *v_off);
      if (has_vt)
	fprintf(out, "/%zu", wf->vt + *vt_off);
      if (has_vn)
	fprintf(out, "%s/%zu", has_vt ? "" : "/", wf->vn + *vn_off);
    }
    fprintf(out, "\n");
  }
  
  *v_off  += LP_VertexList_NumVert(v);
  *vn_off += LP_VertexList_NumVert(vn);
  *vt_off += LP_VertexList_NumVert(vt);
  
  LP_VertexList_Free(vt);
  LP_VertexList_Free(vn);
  LP_VertexList_Free(v);
  free(wfmem);
  return 0;

 err5:
  LP_VertexList_Free(vt);
 err4:
  LP_VertexList_Free(vn);
 err3:
  LP_VertexList_Free(v);
 err2:
  free(wfmem);
 err:
  return -1;
}

int FileObj_Write(FILE *out, const struct lp_vl_list *list, float scale) {
  size_t v_off = 1, vn_off = 1, vt_off = 1, count = 0;

  fprintf(out, "# libpolyhedra\n\n");
  
  while (list) {
    if (FileObj_WriteSingle(out, count++, list->vl, scale, &v_off, &vn_off, &vt_off) < 0)
      return -1;
    
    list = list->next;
  }
  
  return 0;
}
