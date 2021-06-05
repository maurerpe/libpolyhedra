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

#ifndef LP_UTIL_H
#define LP_UTIL_H

float Norm2(const float *v);
float Norm(const float *v);
float Norm2d2(const float *v);
float Norm2d(const float *v);
float Normalize2d(float *v);
float Normalize(float *v);
float Normalize4d(float *v);

float Dist2d2(const float *p1, const float *p2);
float Dist2d(const float *p1, const float *p2);
float Dist2(const float *p1, const float *p2);
float Dist(const float *p1, const float *p2);

float Dot2d(const float *a, const float *b);
float Dot(const float *a, const float *b);
float Cross2d(const float *a, const float *b);
void Cross(float *result, const float *a, const float *b);
void BasisVectors(float *result_x, float *result_y, const float *norm);

/* For CCW winding order */
void PlaneNorm(float *norm, const float *p1, const float *p2, const float *p3);

int Solve2x2(float *out, const float *m, const float *bb);
int Solve3x3(float *out, const float *m, const float *bb);

#endif
