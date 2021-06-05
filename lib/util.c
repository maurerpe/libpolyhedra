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

#include <math.h>

#include "util.h"

float Norm2(const float *v) {
  return Dot(v, v);
}

float Norm(const float *v) {
  return sqrtf(Norm2(v));
}

float Norm2d2(const float *v) {
  return Dot2d(v, v);
}

float Norm2d(const float *v) {
  return sqrtf(Norm2d2(v));
}

float Normalize2d(float *v) {
  float norm, factor;

  norm = Norm2d(v);
  if (norm == 0)
    factor = 0;
  else
    factor = 1 / norm;

  v[0] *= factor;
  v[1] *= factor;
  
  return norm;
}

float Normalize(float *v) {
  float norm, factor;

  norm = Norm(v);
  if (norm == 0)
    factor = 0;
  else
    factor = 1 / norm;

  v[0] *= factor;
  v[1] *= factor;
  v[2] *= factor;
  
  return norm;
}

float Normalize4d(float *v) {
  float norm, factor;
  
  norm = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
  if (norm == 0)
    factor = 0;
  else
    factor = 1 / norm;

  v[0] *= factor;
  v[1] *= factor;
  v[2] *= factor;
  v[3] *= factor;
  
  return norm;
}

float Dist2d2(const float *p1, const float *p2) {
  float x, y;

  x = p1[0] - p2[0];
  y = p1[1] - p2[1];
  
  return x * x + y * y;
}

float Dist2d(const float *p1, const float *p2) {
  return sqrtf(Dist2d2(p1, p2));
}

float Dist2(const float *p1, const float *p2) {
  float x, y, z;

  x = p1[0] - p2[0];
  y = p1[1] - p2[1];
  z = p1[2] - p2[2];
  
  return x * x + y * y + z * z;
}

float Dist(const float *p1, const float *p2) {
  return sqrtf(Dist2(p1, p2));
}

float Dot2d(const float *a, const float *b) {
  return a[0] * b[0] + a[1] * b[1];
}

float Dot(const float *a, const float *b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

float Cross2d(const float *a, const float *b) {
  return a[0] * b[1] - a[1] * b[0];
}

static void CrossRaw(float *result, const float *a, const float *b) {
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];
}

void Cross(float *result, const float *a, const float *b) {
  float rr[3];
  
  if (result == a || result == b) {
    CrossRaw(rr, a, b);
    result[0] = rr[0];
    result[1] = rr[1];
    result[2] = rr[2];
  } else {
    CrossRaw(result, a, b);
  }
}

void BasisVectors(float *result_x, float *result_y, const float *norm) {
  float yy[3], min_val, max_val, abs;
  int min, max;

  yy[0] = norm[0];
  yy[1] = norm[1];
  yy[2] = norm[2];
  
  max = min = 0;
  max_val = min_val = fabsf(norm[0]);
  if ((abs = fabsf(norm[1])) < min_val) {
    min = 1;
    min_val = abs;
  }
  if (abs >= max_val) {
    max = 1;
    max_val = abs;
  }
  if ((abs = fabsf(norm[2])) < min_val) {
    min = 2;
    min_val = abs;
  }
  if (abs >= max_val) {
    max = 2;
    max_val = abs;
  }
  if (max == min) {
    min = (max + 1) % 3;
  }
  
  yy[min] = copysignf(norm[max], -norm[min]);
  yy[max] = copysignf(norm[min], -norm[max]);
  
  Cross(result_x, yy, norm);
  Normalize(result_x);
  Cross(result_y, norm, result_x);
  Normalize(result_y);
}

void PlaneNorm(float *norm, const float *p1, const float *p2, const float *p3) {
  float v1[3], v2[3];
  
  v1[0] = p2[0] - p1[0];
  v1[1] = p2[1] - p1[1];
  v1[2] = p2[2] - p1[2];
  
  v2[0] = p3[0] - p2[0];
  v2[1] = p3[1] - p2[1];
  v2[2] = p3[2] - p2[2];
  
  Cross(norm, v1, v2);
  Normalize(norm);
}

int Solve2x2(float *out, const float *m, const float *bb) {
  float a, b, c, d, det;
  
  a = m[2*0+0];
  b = m[2*0+1];
  c = m[2*1+0];
  d = m[2*1+1];
  det = a * d - b * c;
  if (det == 0)
    return -1;
  
  det = 1 / det;
  
  out[0] = det * (bb[0] * d - bb[1] * b);
  out[1] = det * (bb[1] * a - bb[0] * c);
  
  return 0;
}

/* https://github.com/willnode/N-Matrix-Programmer/blob/master/Info/Matrix_3x3.txt */
int Solve3x3(float *out, const float *m, const float *bb) {
  float a, b, c, det;
  
  a = m[3*1+1] * m[3*2+2] - m[3*1+2] * m[3*2+1];
  b = m[3*1+0] * m[3*2+2] - m[3*1+2] * m[3*2+0];
  c = m[3*1+0] * m[3*2+1] - m[3*1+1] * m[3*2+0];
  
  det = m[3*0+0] * a - m[3*0+1] * b + m[3*0+2] * c;
  if (det == 0)
    return -1;
  det = 1 / det;

  out[0] = det * (bb[0] *   a +
		  bb[1] * - (m[3*0+1] * m[3*2+2] - m[3*0+2] * m[3*2+1]) +
		  bb[2] *   (m[3*0+1] * m[3*1+2] - m[3*0+2] * m[3*1+1]));
  
  out[1] = det * (bb[0] * - b +
		  bb[1] *   (m[3*0+0] * m[3*2+2] - m[3*0+2] * m[3*2+0]) +
		  bb[2] * - (m[3*0+0] * m[3*1+2] - m[3*0+2] * m[3*1+0]));
  
  out[2] = det * (bb[0] *   c +
		  bb[1] * - (m[3*0+0] * m[3*2+1] - m[3*0+1] * m[3*2+0]) +
		  bb[2] *   (m[3*0+0] * m[3*1+1] - m[3*0+1] * m[3*1+0]));
  
  return 0;
}
