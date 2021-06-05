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

#include "file_svg.h"
#include "util.h"

struct lp_vl_list *FileSvg_Read(FILE *in, float scale) {
  fprintf(stderr, "Reading .svg files not yet supported\n");
  return NULL;
}

static int FileSvg_WriteSingle(FILE *out, const struct lp_vertex_list *vl, float scale) {
  size_t count, num, num_lines;
  float *ff1, *ff2, *ff3;
  unsigned int *ind;

  if (LP_VertexList_FloatsPerVert(vl) < 2) {
    fprintf(stderr, "Error: Too few floats per vert for .svg\n");
    return -1;
  }
  num = LP_VertexList_NumInd(vl);
  ind = LP_VertexList_GetInd(vl);

  switch (LP_VertexList_PrimativeType(vl)) {
  case lp_pt_line:
    num_lines = num / 2;
    for (count = 0; count < num_lines; count++) {
      ff1 = LP_VertexList_LookupVert(vl, 2 * count);
      ff2 = LP_VertexList_LookupVert(vl, 2 * count + 1);
      fprintf(out, "    <!-- %04u,%04u --><line x1=\"%g\" y1=\"%g\" x2=\"%g\" y2=\"%g\"/>\n",
	      ind[2 * count],
	      ind[2 * count + 1],
	      ff1[0] * scale,
	      ff1[1] * scale,
	      ff2[0] * scale,
	      ff2[1] * scale);
    }
    break;
    
  case lp_pt_triangle:
    num_lines = num / 3;
    for (count = 0; count < num_lines; count++) {
      ff1 = LP_VertexList_LookupVert(vl, 3 * count);
      ff2 = LP_VertexList_LookupVert(vl, 3 * count + 1);
      ff3 = LP_VertexList_LookupVert(vl, 3 * count + 2);
      fprintf(out, "    <!-- %04u,%04u,%04u --><polygon points=\"%g,%g %g,%g %g,%g\"/>\n",
	      ind[3 * count],
	      ind[3 * count + 1],
	      ind[3 * count + 2],
	      ff1[0] * scale,
	      ff1[1] * scale,
	      ff2[0] * scale,
	      ff2[1] * scale,
	      ff3[0] * scale,
	      ff3[1] * scale);
    }
    break;
    
  default:
    fprintf(stderr, "Error: Incorrect primative type for .svg\n");
    return -1;
  }
  
  return 0;
}

int FileSvg_Write(FILE *out, const struct lp_vl_list *list, float scale) {
  float min[2], max[2];
  const struct lp_vl_list *cur;
  const float *ff;
  size_t fpv, count, num, ii;

  min[0] =  INFINITY;
  min[1] =  INFINITY;
  max[0] = -INFINITY;
  max[1] = -INFINITY;
  for (cur = list; cur != NULL; cur = cur->next) {
    fpv = LP_VertexList_FloatsPerVert(cur->vl);
    if (fpv < 2) {
      fprintf(stderr, "Error: Too few floats per vert for .svg\n");
      return -1;
    }
    
    num = LP_VertexList_NumVert(cur->vl);
    ff  = LP_VertexList_GetVert(cur->vl);
    for (count = 0; count < num; count++, ff += fpv) {
      for (ii = 0; ii < 2; ii++) {
	if (ff[ii] < min[ii])
	  min[ii] = ff[ii];
	if (ff[ii] > max[ii])
	  max[ii] = ff[ii];
      }
    }
  }
  
  fprintf(out, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
  fprintf(out, "<svg viewBox=\"%g %g %g %g\" xmlns=\"http://www.w3.org/2000/svg\">\n\n",
	  min[0] * scale,
	  min[1] * scale,
 	  (max[0] - min[0]) * scale,
	  (max[1] - min[1]) * scale);
  
  for (cur = list, count = 0; cur != NULL; cur = cur->next, count++) {
    switch (LP_VertexList_PrimativeType(cur->vl)) {
    case lp_pt_line:
      fprintf(out, "  <g id=\"polyline_%03zu\" stroke=\"black\" stroke-width=\"1\" fill=\"none\">\n", count);
      break;
    case lp_pt_triangle:
      fprintf(out, "  <g id=\"polyline_%03zu\" fill=\"blue\" stroke=\"none\">\n", count);
      break;
    default:
      fprintf(stderr, "Error: Incorrect primative type for .svg\n");
      return -1;
    }
    if (FileSvg_WriteSingle(out, cur->vl, scale) < 0)
      return -1;
    fprintf(out, "  </g>\n\n");
  }
  
  fprintf(out, "</svg>\n");
  
  return 0;
}
