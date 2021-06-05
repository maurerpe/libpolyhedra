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

#include "array.h"

struct array {
  size_t len;
  size_t num_alloc;
  void **data;
  array_free_func_t free_func;
};

struct array *Array_New(size_t initial_alloc, array_free_func_t free_func) {
  struct array *arr;

  if ((arr = malloc(sizeof(*arr))) == NULL)
    goto err;
  memset(arr, 0, sizeof(*arr));
  arr->num_alloc = initial_alloc;
  arr->free_func = free_func;
  
  if ((arr->data = calloc(initial_alloc, sizeof(*arr->data))) == NULL)
    goto err2;

  return arr;

 err2:
  free(arr);
 err:
  return NULL;
}
    
void Array_Free(struct array *arr) {
  size_t count;
  
  if (arr == NULL)
    return;
  
  if (arr->free_func)
    for (count = 0; count < arr->len; count++)
      arr->free_func(arr->data[count]);
  
  free(arr->data);
  free(arr);
}

int Array_Add(struct array *arr, void *data) {
  void **new_data;
  size_t new_alloc;
  
  if (arr->len >= arr->num_alloc) {
    if (arr->num_alloc >= SIZE_MAX >> 1)
      return -1;

    new_alloc = arr->num_alloc << 1;
    if ((new_data = calloc(new_alloc, sizeof(*new_data))) == NULL)
      return -1;
    
    memcpy(new_data, arr->data, sizeof(*new_data) * arr->num_alloc);
    free(arr->data);
    arr->data = new_data;
    arr->num_alloc = new_alloc;
  }
  
  arr->data[arr->len++] = data;
  
  return 0;
}

void Array_Remove(struct array *arr, ssize_t idx) {
  if (idx < 0)
    idx += arr->len;
  
  if (idx < 0 || idx >= arr->len)
    return;

  if (arr->free_func)
    arr->free_func(arr->data[idx]);
  
  arr->data[idx] = arr->data[--arr->len];
}

size_t Array_Length(const struct array *arr) {
  return arr->len;
}

void **Array_Data(struct array *arr) {
  return arr->data;
}
