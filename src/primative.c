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
  fprintf(out, "%s: Generate primative shapes\n", PACKAGE_STRING);
  fprintf(out, "  primative -t <type> [-n <number>]\n");
  fprintf(out, "    [-x <xsize>] [-y <ysize>] [-z <zsize>] <outfile>\n");
  fprintf(out, "  primative -h\n");
  fprintf(out, "  Generates a primative polyhedron.  The output is written to <outfile>.\n");
  fprintf(out, "  The generated polyhedron is controlled by the following options:\n");
  fprintf(out, "  -h\n");
  fprintf(out, "      Print this help and exit\n");
  fprintf(out, "  -n <number>\n");
  fprintf(out, "      Parameter that controls how the shape is generated.  See -t.\n");
  fprintf(out, "  -t <type>\n");
  fprintf(out, "      Type to generate. Supported types are:\n");
  fprintf(out, "       * cube: Generates a rectangular prism.  <number> is unused");
  fprintf(out, "       * cylinder: Generates a cylinder along the z-axis.  Diameter is\n");
  fprintf(out, "         <xsize>, height is <zsize>.  <number> is the number of points per\n");
  fprintf(out, "         revolution.  The default is 3.");
  fprintf(out, "       * uvsphere: Generates a sphere of diameter <xsize>.  <number> is the\n");
  fprintf(out, "         number of segments and the number of rings.  The default is 3.\n");
  fprintf(out, "       * icosphere: Generates a sphere of diameter <xsize>.  <number>\n");
  fprintf(out, "         represents the number of subdivisions.  The number of faces is\n");
  fprintf(out, "          20 * 4^<number>.  Default number of subdivisions is zero.");
  fprintf(out, "  -x <xsize>\n");
  fprintf(out, "      Size of primative in the x direction.  Default is 1.\n");
  fprintf(out, "  -y <xsize>\n");
  fprintf(out, "      Size of primative in the y direction.  Default is 1.\n");
  fprintf(out, "  -z <xsize>\n");
  fprintf(out, "      Size of primative in the z direction.  Default is 1.\n");
}

enum prim_type {
  pt_cube,
  pt_cylinder,
  pt_uvsphere,
  pt_icosphere
};

int main(int argc, char *argv[]) {
  int num_specified = 0;
  long number = 0;
  int type_specified = 0;
  enum prim_type type = 0;
  float xsize = 1.0f, ysize = 1.0f, zsize = 1.0f;
  char *end;
  size_t count;
  int opt;
  struct lp_vertex_list *vl = NULL;
  struct lp_vl_list list;

#ifdef HAVE_SETLOCALE
  setlocale(LC_NUMERIC, "C");
#endif
  
  while ((opt = getopt(argc, argv, "hn:t:x:y:z:")) >= 0) {
    switch (opt) {
    case 'h':
      help(stdout);
      exit(0);
      
    case 'n':
      num_specified = 1;
      number = strtol(optarg, &end, 0);
      if (*end != '\0') {
	fprintf(stderr, "Error: expected integer for -n argument: %s\n", optarg);
	help(stderr);
	exit(1);
      }
      break;

    case 't':
      type_specified = 1;
      if (strcasecmp(optarg, "cube") == 0)
	type = pt_cube;
      else if (strcasecmp(optarg, "cylinder") == 0)
	type = pt_cylinder;
      else if (strcasecmp(optarg, "uvsphere") == 0)
	type = pt_uvsphere;
      else if (strcasecmp(optarg, "icosphere") == 0)
	type = pt_icosphere;
      else {
	fprintf(stderr, "Error: Unknown type: %s\n", optarg);
	help(stderr);
	exit(1);
      }
      break;

    case 'x':
      xsize = strtof(optarg, &end);
      if (*end != '\0') {
	fprintf(stderr, "Error: expected floating point number for -x argument: %s\n", optarg);
	help(stderr);
	exit(1);
      }
      break;      
      
    case 'y':
      ysize = strtof(optarg, &end);
      if (*end != '\0') {
	fprintf(stderr, "Error: expected floating point number for -y argument: %s\n", optarg);
	help(stderr);
	exit(1);
      }
      break;      
      
    case 'z':
      zsize = strtof(optarg, &end);
      if (*end != '\0') {
	fprintf(stderr, "Error: expected floating point number for -z argument: %s\n", optarg);
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
  
  if (optind + 1 != argc) {
    fprintf(stderr, "Error: Exactly one output file expected\n");
    help(stderr);
    exit(1);
  }
  
  if (!type_specified) {
    fprintf(stderr, "Error: -t argument required\n");
    help(stderr);
    exit(1);
  }
  
  switch (type) {
  case pt_cube:
    vl = LP_Cube(xsize / 2, ysize / 2, zsize / 2);
    break;
  case pt_cylinder:
    vl = LP_Cylinder(xsize / 2, zsize, num_specified ? number : 3);
    break;
  case pt_uvsphere:
    vl = LP_UVSphere(xsize / 2, num_specified ? number : 3, num_specified ? number : 3);
    break;
  case pt_icosphere:
    vl = LP_IcoSphere(xsize / 2, num_specified ? number : 0);
    break;
  }
  
  if (vl == NULL) {
    fprintf(stderr, "Error: Unable to generate shape\n");
    help(stderr);
    exit(1);
  }
  
  list.vl   = vl;
  list.next = NULL;

  if (LP_VertexList_Write(argv[optind], &list, 1.0) < 0) {
    fprintf(stderr, "Error writing to file: %s\n", argv[optind]);
    help(stderr);
    exit(1);
  }
  
  exit(0);
}
