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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_LOCALE_H
#include <locale.h>
#endif

#include <limits.h>
#include <math.h>
#include <string.h>

#include "libpolyhedra.h"

void help(FILE *out) {
  fprintf(out, "%s: convert and operate on polyhedra with triangular faces\n", PACKAGE_STRING);
  fprintf(out, "  polyhedra [-c] [-d t] [-h] [-m] [-o <outfile>] [-p <x,y,z,d>]\n");
  fprintf(out, "    [-q] [-s <faces>] [-x <scale>] <infile>...\n\n");
  fprintf(out, "  Reads in the polyhedra contained in input files and optionally performs\n");
  fprintf(out, "  operations on them.  The operations, when selected, are always performed\n");
  fprintf(out, "  this order, regardless of the order of the options in the command:\n");
  fprintf(out, "    1. Scale (enabled with -x)\n");
  fprintf(out, "    2. Simplify (enabled with -s)\n");
  fprintf(out, "    3. Convex Hull (enabled with -c)\n");
  fprintf(out, "    4. Approximate Surface decompisition (enabled with -d)\n");
  fprintf(out, "    5. Mass properities (enabled with -m)\n\n");
  fprintf(out, "  -c\n");
  fprintf(out, "    Calculate the convex hull\n\n");
  fprintf(out, "  -d threshold\n");
  fprintf(out, "    Perform approximate surface decomposition into convex polyhedra\n\n");
  fprintf(out, "  -h\n");
  fprintf(out, "    Print this help screen and exit\n\n");
  fprintf(out, "  -m\n");
  fprintf(out, "    Calculate mass properties of each polyhedra individually:\n");
  fprintf(out, "      * volume,\n");
  fprintf(out, "      * center of mass, and\n");
  fprintf(out, "      * interia tensor\n\n");
  fprintf(out, "  -o <outfile>\n");
  fprintf(out, "    Save resulting polyhedra to <outfile>.  Default: out.obj\n");
  fprintf(out, "    To omit saving output pass an empty string as <outfile>\n\n");
  fprintf(out, "  -p <x,y,z,d>\n");
  fprintf(out, "    Cut the polyhedra along a plane define by the normal: (x, y, z) and is\n");
  fprintf(out, "    d units from the origin.\n\n");
  fprintf(out, "  -q\n");
  fprintf(out, "    Quiet.  Supress status outputs\n\n");
  fprintf(out, "  -s <faces>\n");
  fprintf(out, "    Simplify each polyhedra to no more than <faces> faces.\n\n");
  fprintf(out, "  -x <scale>\n");
  fprintf(out, "    Scale each each polyhedra by a factor of <scale>\n\n");
}

void Parse_Floats(float *out, size_t num, const char *str) {
  const char *ss = str;
  char *end;
  size_t count;
  char exp;
  
  for (count = 0; count < num; count++, ss = end + 1) {
    out[count] = strtof(ss, &end);
    exp = count == num - 1 ? '\0' : ',';
    if (*end != exp) {
      fprintf(stderr, "Expecting comma seperated list of floats: '%s'\n", str);
      help(stderr);
      exit(1);
    }
  }
}

int main(int argc, char *argv[]) {
  int mass_prop = 0;
  unsigned long long simplify = 0;
  int convex = 0;
  int decomp = 0;
  int verbose = 1;
  int plane = 0;
  float scale = 1.0, dval[1], pval[4];
  const char *outfile = "out.obj";
  char *end;
  struct lp_vl_list *data = NULL, *list, *list2, *out;
  struct lp_mass_properties mp;
  size_t count;
  int opt;

#ifdef HAVE_SETLOCALE
  setlocale(LC_NUMERIC, "C");
#endif
  
  while ((opt = getopt(argc, argv, "cd:hmo:p:qs:x:")) >= 0) {
    switch (opt) {
    case 'c':
      convex = 1;
      break;

    case 'd':
      decomp = 1;
      Parse_Floats(dval, 1, optarg);
      break;
      
    case 'h':
      help(stdout);
      exit(0);
      
    case 'm':
      mass_prop = 1;
      break;
      
    case 'o':
      outfile = strdup(optarg);
      break;
      
    case 'p':
      plane = 1;
      Parse_Floats(pval, 4, optarg);
      break;
      
    case 'q':
      verbose = 0;
      break;
      
    case 's':
      simplify = strtoull(optarg, &end, 0);
      if (*end != '\0') {
	fprintf(stderr, "Error: expected non-negative integer for -s argument: %s\n", optarg);
	help(stderr);
	exit(1);
      }
      break;

    case 'x':
      scale = strtof(optarg, &end);
      if (*end != '\0') {
	fprintf(stderr, "Error: expected floating point number for -x argument: %s\n", optarg);
	help(stderr);
	exit(1);
      }
      break;
      
    default:
      fprintf(stderr, "Error unknown option '%c'\n", (char) opt);
      help(stderr);
      exit(1);
    }
  }
  
  if (optind >= argc) {
    fprintf(stderr, "Error: At least one input file expected\n");
    help(stderr);
    exit(1);
  }
  
  while (optind < argc) {
    if ((list = LP_VertexList_Read(argv[optind], scale)) == NULL)
      exit(1);
    if ((data = LP_VertexList_ListJoin(data, list)) == NULL)
      exit(1);
    optind++;
  }
  
  if (simplify > 0) {
    if (verbose)
      printf("\nSimplifying\n");
    for (count = 0, list = data; list != NULL; list = list->next, count++) {
      if ((list->vl = LP_Simplify(list->vl, simplify, 0)) == NULL)
	exit(1);
    }
  }
  
  if (convex) {
    if (verbose)
      printf("\nCalculating convex hulls\n");
    for (list = data; list != NULL; list = list->next) {
      if ((list->vl = LP_ConvexHull(list->vl)) == NULL)
	exit(1);
    }    
  }
  
  if (plane) {
    out = NULL;
    for (count = 0, list = data; list != NULL; list = list->next, count++) {
      if (verbose)
	printf("Cutting polyhedra %zu along plane\n", count);
      if ((list2 = LP_PlaneCut(list->vl, pval, pval[3])) == NULL)
	exit(1);
      if (verbose)
	printf("  -> Split into %zu polyhedra\n", LP_VertexList_ListLength(list2));
      if ((out = LP_VertexList_ListJoin(out, list2)) == NULL)
	exit(1);
    }
    LP_VertexList_ListFree(data);
    data = out;
  }
  
  if (decomp) {
    out = NULL;
    for (count = 0, list = data; list != NULL; list = list->next, count++) {
      if (verbose)
	printf("Decomposing polyhedra %zu\n", count);
      if ((list2 = LP_ConvexDecomp(list->vl, dval[0])) == NULL)
	exit(1);
      if (verbose)
	printf("  -> Split into %zu convex polyhedra\n", LP_VertexList_ListLength(list2));
      if ((out = LP_VertexList_ListJoin(out, list2)) == NULL)
	exit(1);
    }
    LP_VertexList_ListFree(data);
    data = out;
  }
  
  if (mass_prop) {
    if (verbose)
      printf("\nCalculating mass properies\n");
    for (count = 0, list = data; list != NULL; list = list->next, count++) {
      LP_MassProperties(list->vl, &mp);
      printf("Properties for polyhedra %zu:\n", count);
      printf("  Vertices: %u, Indices: %zu\n",
	     LP_VertexList_NumVert(list->vl),
	     LP_VertexList_NumInd(list->vl));
      printf("  Volume:         %g\n", mp.volume);
      printf("  Center of mass: (%g, %g, %g)\n", mp.center_of_mass[0], mp.center_of_mass[1], mp.center_of_mass[2]);
      printf("  Inertia Tensor:\n");
      printf("    [%20g, %20g, %20g]\n", mp.inertia_tensor[0], mp.inertia_tensor[1], mp.inertia_tensor[2]);
      printf("    [%20g, %20g, %20g]\n", mp.inertia_tensor[3], mp.inertia_tensor[4], mp.inertia_tensor[5]);
      printf("    [%20g, %20g, %20g]\n\n", mp.inertia_tensor[6], mp.inertia_tensor[7], mp.inertia_tensor[8]);
    }
  }
  
  if (*outfile != '\0') {
    LP_VertexList_Write(outfile, data, 1.0);
  }
  
  exit(0);
}
