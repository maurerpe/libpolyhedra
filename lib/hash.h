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

#ifndef LP_HASH_H
#define LP_HASH_H

#include <stdint.h>

struct hash;

typedef uint64_t (*hash_hash_func_t)(const void *user, const unsigned char secret[16], const void *key);
typedef int (*hash_cmp_func_t)(const void *user, const void *key_a, const void *key_b); /* Returns 0 if objects are equal */
typedef void *(*hash_copy_func_t)(void *user, const void *);
typedef void (*hash_free_func_t)(void *user, void *);

struct hash *Hash_New(void *user, hash_hash_func_t hash_func, hash_cmp_func_t cmp, hash_copy_func_t copy_key, hash_free_func_t free_key, hash_copy_func_t copy_data, hash_free_func_t free_data, hash_free_func_t free_user);
struct hash *Hash_NewString(void *user, hash_copy_func_t copy_data, hash_free_func_t free_data, hash_free_func_t free_user);
struct hash *Hash_NewPtr(void *user, hash_free_func_t free_key, hash_copy_func_t copy_data, hash_free_func_t free_data, hash_free_func_t free_user);
struct hash *Hash_NewFixed(size_t size, void *user, hash_copy_func_t copy_data, hash_free_func_t free_data, hash_free_func_t free_user);
void Hash_Free(struct hash *hash);
void Hash_Clear(struct hash *hash);

size_t Hash_NumEntries(const struct hash *hash);

void *Hash_Lookup(const struct hash *hash, const void *key, int *was_found);
int Hash_Insert(struct hash *hash, const void *key, const void *data, void **key_out);
int Hash_Remove(struct hash *hash, const void *key);

struct hash_iterator;
struct hash_iterator *Hash_IteratorNew(struct hash *hash);
void Hash_IteratorFree(struct hash_iterator *hi);

int Hash_IteratorNext(struct hash_iterator *hi);
void *Hash_IteratorGetKey(struct hash_iterator *hi);
void *Hash_IteratorGetData(struct hash_iterator *hi);

void *Hash_GetFirstKey(struct hash *hash);

#endif
