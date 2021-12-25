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

#ifndef LIBPOLYHEDRA_H
#define LIBPOLYHEDRA_H

#ifdef __cplusplus
extern "C" {
#endif

/* At present, only polyhedra with triangular faces are supported */

/*********************** Vertex List ********************************/
/* Basic polyhedra data structure */
enum primative_type {
  lp_pt_point,
  lp_pt_line,
  lp_pt_triangle,
  lp_pt_unspecified
};

struct lp_vertex_list;

struct lp_vertex_list *LP_VertexList_New(size_t floats_per_vert, enum primative_type pt);
void LP_VertexList_Free(struct lp_vertex_list *vl);
void LP_VertexList_Clear(struct lp_vertex_list *vl);
struct lp_vertex_list *LP_VertexList_Copy(const struct lp_vertex_list *vl, size_t new_floats_per_vert);

/* Return index, UINT_MAX from limits.h on error */
unsigned int LP_VertexList_Add(struct lp_vertex_list *vl, const float *vert);
unsigned int LP_VertexList_AddIndex(struct lp_vertex_list *vl, unsigned int index);

/* Optional: No more vertices will be added, free extra memory */
void LP_VertexList_Finalize(struct lp_vertex_list *vl);

size_t LP_VertexList_FloatsPerVert(const struct lp_vertex_list *vl);
enum primative_type LP_VertexList_PrimativeType(const struct lp_vertex_list *vl);
unsigned int LP_VertexList_NumVert(const struct lp_vertex_list *vl);
size_t LP_VertexList_NumInd(const struct lp_vertex_list *vl);
float *LP_VertexList_GetVert(const struct lp_vertex_list *vl);
unsigned int *LP_VertexList_GetInd(const struct lp_vertex_list *vl);

float *LP_VertexList_LookupVert(const struct lp_vertex_list *vl, unsigned int index);

struct lp_vl_list {
  struct lp_vertex_list *vl;
  struct lp_vl_list *next;
};

struct lp_vl_list *LP_VertexList_ListAppend(struct lp_vl_list *list, struct lp_vertex_list *vl);
struct lp_vl_list *LP_VertexList_ListJoin(struct lp_vl_list *list1, struct lp_vl_list *list2);
size_t LP_VertexList_ListLength(struct lp_vl_list *list);
void LP_VertexList_ListFree(struct lp_vl_list *list); /* Also frees the vertex lists */

/*         |      3D      |      2D      |
 * Format  | Read | Write | Read | Write |
 * --------|------|-------|------|-------|
 * .obj    |  x   |   x   |      |       |
 * .stl*   |  x   |   x   |      |       |
 * .svg    |      |       |      |   x   |
 * 
 * *binary stl only, ascii stl not supported
 */
struct lp_vl_list *LP_VertexList_Read(const char *filename, float scale);
int LP_VertexList_Write(const char *filename, struct lp_vl_list *list, float scale);

/****************** Triangulate 2D polygons w/ holes ****************/
/* Input:  a list of 2D line segments in any order. Polygons can share
 *         points and edges, but must not intersect or overlap.
 * Output: a list of 2D triangles
 */
struct lp_vertex_list *LP_Triangulate2D(const struct lp_vertex_list *in);

/********************** Mass properties ****************************/
/* Assumes closed polyhedra with uniform desnity */
struct lp_mass_properties {
  double volume;
  double center_of_mass[3];
  double inertia_tensor[9]; /* about center of mass, for unit density */
};

void LP_MassProperties(const struct lp_vertex_list *in, struct lp_mass_properties *properties);

/*********************** Transform *********************************/
struct lp_transform;

struct lp_transform *LP_Transform_New(void);
void LP_Transform_Free(struct lp_transform *trans);

void LP_Transform_Copy(struct lp_transform *dest, const struct lp_transform *src);
void LP_Transform_SetToIdentity(struct lp_transform *trans);
void LP_Transform_SetToMatrix4x4(struct lp_transform *trans, const float *m);

void LP_Transform_Translate(struct lp_transform *trans, float dx, float dy, float dz);
void LP_Transform_Rotate(struct lp_transform *trans, float angle_rad, float axis_x, float axis_y, float axis_z);
void LP_Transform_ApplyQauternion(struct lp_transform *trans, const float *wxyz);
void LP_Transform_Combine(struct lp_transform *dest, const struct lp_transform *a, const struct lp_transform *b);

void LP_Transform_Invert(struct lp_transform *trans);

#define LP_TRANSFORM_NONE      0
#define LP_TRANSFORM_NO_OFFSET 1
#define LP_TRANSFORM_INVERT    2
void LP_Transform_Point(const struct lp_transform *trans, float *dest, const float *src, int options);
struct lp_vertex_list *LP_Transform_VertexList(const struct lp_transform *trans, const struct lp_vertex_list *src, int options);

/*********************** Primatives ********************************/
struct lp_vertex_list *LP_Cube(float x, float y, float z);
struct lp_vertex_list *LP_Cylinder(float r, float h, int pts_per_rev);
struct lp_vertex_list *LP_UVSphere(float radius, int segs, int rings);
/* Num triangles = 20 * 4^num_subdiv */
struct lp_vertex_list *LP_IcoSphere(float radius, int num_subdiv);
  
/*********************** Simplify **********************************/
struct lp_vertex_list *LP_Simplify(const struct lp_vertex_list *in, size_t num_faces_out, float aggregation_thresh);

/*********************** Convex Hull *******************************/
struct lp_vertex_list *LP_ConvexHull(const struct lp_vertex_list *in);

/*********************** Plane Cut *********************************/
/* Cut a polyhedron into two pieces along a plane */
struct lp_vl_list *LP_PlaneCut(const struct lp_vertex_list *in, const float *norm, float dist);

/*************** Approx Convex Decompisition ***********************/
/* Decomposes a polyhedron into convex polyhedra */
struct lp_vl_list *LP_ConvexDecomp(const struct lp_vertex_list *in, float threshold);

#ifdef __cplusplus
}
#endif

#endif
