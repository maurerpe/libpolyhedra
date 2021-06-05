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

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <limits.h>
#include <string.h>

#include "hash.h"
#include "random.h"
#include "SipHash/siphash.h"

struct element {
  uint64_t hash_val;
  void *key;
  void *data;
  struct element *next;
};

struct hash {
  size_t num_buckets;
  size_t num_items;

  void *user;
  hash_hash_func_t hash_func;
  hash_cmp_func_t cmp;
  hash_copy_func_t copy_key;
  hash_free_func_t free_key;
  hash_copy_func_t copy_data;
  hash_free_func_t free_data;
  hash_free_func_t free_user;
  
  struct element **buckets;
  
  unsigned char secret[16];
};

struct hash *Hash_New(void *user, hash_hash_func_t hash_func, hash_cmp_func_t cmp, hash_copy_func_t copy_key, hash_free_func_t free_key, hash_copy_func_t copy_data, hash_free_func_t free_data, hash_free_func_t free_user) {
  struct hash *hash;
  
  if (hash_func == NULL || cmp == NULL) {
    fprintf(stderr, "Error: Hash function and compare function cannot be NULL\n");
    goto err;
  }
  
  if ((hash = malloc(sizeof(*hash))) == NULL) {
    fprintf(stderr, "Error: Could not allocate space for hash table\n");
    goto err;
  }
  memset(hash, 0, sizeof(*hash));
  
  hash->user       = user;
  hash->hash_func  = hash_func;
  hash->cmp        = cmp;
  hash->copy_key   = copy_key;
  hash->free_key   = free_key;
  hash->copy_data  = copy_data;
  hash->free_data  = free_data;
  hash->free_user  = free_user;
  
  hash->num_buckets = 128;
  if ((hash->buckets = calloc(hash->num_buckets, sizeof(void *))) == NULL) {
    fprintf(stderr, "Error: Could not allocate space for hash table buckets\n");
    goto err2;
  }
  
  Random(hash->secret, sizeof(hash->secret));
  
  return hash;
  
 err2:
  free(hash);
 err:
  return NULL;
}

static uint64_t StringHash(const void *user, const unsigned char secret[16], const void *key) {
  return siphash(secret, (const unsigned char *) key, strlen((const char *) key));
}

static int StringCmp(const void *user, const void *key_a, const void *key_b) {
  return strcmp((const char *) key_a, (const char *) key_b);
}

static void *StringCopy(void *user, const void *key) {
  return strdup((const char *) key);
}

static void StringFree(void *user, void *key) {
  free(key);
}

struct hash *Hash_NewString(void *user, hash_copy_func_t copy_data, hash_free_func_t free_data, hash_free_func_t free_user) {
  return Hash_New(user, &StringHash, &StringCmp, &StringCopy, &StringFree, copy_data, free_data, free_user);
}

static uint64_t PtrHash(const void *user, const unsigned char secret[16], const void *key) {
  return siphash(secret, (unsigned char *) &key, sizeof(key));
}

static int PtrCmp(const void *user, const void *key_a, const void *key_b) {
  return !(key_a == key_b);
}

struct hash *Hash_NewPtr(void *user, hash_free_func_t free_key, hash_copy_func_t copy_data, hash_free_func_t free_data, hash_free_func_t free_user) {
  return Hash_New(user, &PtrHash, &PtrCmp, NULL, free_key, copy_data, free_data, free_user);
}

struct fixed_data {
  size_t size;
  void *user;
  hash_copy_func_t copy_data;
  hash_free_func_t free_data;
  hash_free_func_t free_user;
};

static uint64_t FixedHash(const void *user, const unsigned char secret[16], const void *key) {
  const struct fixed_data *fd = (const struct fixed_data *) user;
  
  return siphash(secret, (const unsigned char *) key, fd->size);
}

static int FixedCmp(const void *user, const void *key_a, const void *key_b) {
  const struct fixed_data *fd = (const struct fixed_data *) user;
  
  return memcmp((const char *) key_a, (const char *) key_b, fd->size);
}

static void *FixedCopyKey(void *user, const void *key) {
  struct fixed_data *fd = (struct fixed_data *) user;  
  void *mem;

  if ((mem = malloc(fd->size)) == NULL)
    return NULL;
  memcpy(mem, key, fd->size);
  
  return mem;
}

static void FixedFreeKey(void *user, void *key) {
  free(key);
}

static void *FixedCopyData(void *user, const void *data) {
  struct fixed_data *fd = (struct fixed_data *) user;
  
  if (fd->copy_data)
    return fd->copy_data(fd->user, data);
  
  return (void *) data;
}

static void FixedFreeData(void *user, void *data) {
  struct fixed_data *fd = (struct fixed_data *) user;
  
  if (fd->free_data)
    fd->free_data(fd->user, data);
}

static void FixedFreeUser(void *user, void *data) {
  struct fixed_data *fd = (struct fixed_data *) user;
  
  if (fd == NULL)
    return;

  if (fd->free_user)
    fd->free_user(fd->user, data);
  
  free(fd);
}

struct hash *Hash_NewFixed(size_t size, void *user, hash_copy_func_t copy_data, hash_free_func_t free_data, hash_free_func_t free_user) {
  struct fixed_data *fd;
  struct hash *hash;

  if ((fd = malloc(sizeof(*fd))) == NULL)
    goto err;
  memset(fd, 0, sizeof(*fd));

  fd->size      = size;
  fd->user      = user;
  fd->copy_data = copy_data;
  fd->free_data = free_data;
  fd->free_user = free_user;
  
  if ((hash = Hash_New(fd, &FixedHash, &FixedCmp, &FixedCopyKey, &FixedFreeKey, &FixedCopyData, &FixedFreeData, &FixedFreeUser)) == NULL)
    goto err2;
  
  return hash;

 err2:
  free(fd);
 err:
  return NULL;
}

static void FreeElement(struct hash *hash, struct element *elem) {
  if (hash->free_key)
    hash->free_key(hash->user, elem->key);
  if (hash->free_data)
    hash->free_data(hash->user, elem->data);
  free(elem);
}

void Hash_Free(struct hash *hash) {
  size_t count;
  struct element *cur, *next;

  if (hash == NULL)
    return;
  
  for (count = 0; count < hash->num_buckets; count++) {
    cur = hash->buckets[count];
    while (cur) {
      next = cur->next;
      FreeElement(hash, cur);
      cur = next;
    }
  }

  if (hash->free_user)
    hash->free_user(hash->user, NULL);
  free(hash->buckets);
  free(hash);
}

void Hash_Clear(struct hash *hash) {
  size_t count;
  struct element *cur, *next;
  
  for (count = 0; count < hash->num_buckets; count++) {
    cur = hash->buckets[count];
    hash->buckets[count] = NULL;
    while (cur) {
      next = cur->next;
      FreeElement(hash, cur);
      cur = next;
    }
  }
}

size_t Hash_NumEntries(const struct hash *hash) {
  return hash->num_items;
}

static void Rehash(struct hash *hash) {
  struct element **new_buckets, **old_buckets, *cur, *next, **dest;
  size_t new_num_buckets, count;
  
  if (hash->num_buckets > SIZE_MAX / 2)
    return;
  
  new_num_buckets = hash->num_buckets << 1;
  
  if ((new_buckets = calloc(new_num_buckets, sizeof(void *))) == NULL) {
    fprintf(stderr, "Could not allcoate memory for hash buckets");
    return;
  }
  
  old_buckets = hash->buckets;
  
  for (count = 0; count < hash->num_buckets; count++) {
    cur = old_buckets[count];
    while (cur) {
      next = cur->next;
      dest = new_buckets + (cur->hash_val % new_num_buckets);
      cur->next = *dest;
      *dest = cur;
      cur = next;
    }
  }

  hash->buckets = new_buckets;
  hash->num_buckets = new_num_buckets;
  free(old_buckets);
}

static struct element **Find(const struct hash *hash, const void *key, uint64_t hash_val) {
  size_t bucket;
  struct element **cur;
  
  bucket = hash_val % hash->num_buckets;

  cur = hash->buckets + bucket;
  while (*cur) {
    if ((*cur)->hash_val == hash_val) {
      if (hash->cmp(hash->user, (*cur)->key, key) == 0)
	return cur;
    }
    
    cur = &(*cur)->next;
  }
  
  return cur;
}

void *Hash_Lookup(const struct hash *hash, const void *key, int *was_found) {
  uint64_t hash_val;
  struct element **loc;
  
  hash_val = hash->hash_func(hash->user, hash->secret, key);
  loc = Find(hash, key, hash_val);
  if (!*loc) {
    if (was_found)
      *was_found = 0;
    return NULL;
  }
  
  if (was_found)
    *was_found = 1;
  
  return (*loc)->data;
}

int Hash_Insert(struct hash *hash, const void *key, const void *data, void **key_out) {
  uint64_t hash_val;
  struct element **loc, *new_elem;
  void *new_data, *new_key;
  
  hash_val = hash->hash_func(hash->user, hash->secret, key);
  loc = Find(hash, key, hash_val);
  if (*loc) {
    if (hash->free_data)
      hash->free_data(hash->user, (*loc)->data);
    if (hash->copy_data) {
      if ((new_data = hash->copy_data(hash->user, data)) == NULL) {
	fprintf(stderr, "Could not copy data into existing hash element\n");
	return -1;
      }
      (*loc)->data = new_data;
    } else {
      (*loc)->data = (void *) data;
    }
    
    if (key_out)
      *key_out = (*loc)->key;
    return 0;
  }
  
  if ((new_elem = malloc(sizeof(*new_elem))) == NULL) {
    fprintf(stderr, "Could not allocate memory for hash element\n");
    return -1;
  }
  memset(new_elem, 0, sizeof(*new_elem));
  
  new_elem->hash_val = hash_val;
  new_elem->next = NULL;
  
  if (hash->copy_key) {
    if ((new_key = hash->copy_key(hash->user, key)) == NULL) {
      fprintf(stderr, "Could not copy key into new hash element\n");
      return -1;
    }
    new_elem->key = new_key;
  } else {
    new_elem->key = (void *) key;
  }

  if (hash->copy_data) {
    if ((new_data = hash->copy_data(hash->user, data)) == NULL) {
      fprintf(stderr, "Could not copy data into new hash element\n");
      return -1;
    }
    new_elem->data = new_data;
  } else {
    new_elem->data = (void *) data;
  }

  *loc = new_elem;

  if ((hash->num_items >> 1) > hash->num_buckets)
    Rehash(hash);
  
  hash->num_items++;  
  if (key_out)
    *key_out = new_elem->key;
  return 1;
}

int Hash_Remove(struct hash *hash, const void *key) {
  uint64_t hash_val;
  struct element **loc, *rem;
  
  hash_val = hash->hash_func(hash->user, hash->secret, key);
  loc = Find(hash, key, hash_val);
  if (!*loc)
    return 0;
  
  rem = *loc;
  *loc = rem->next;
  
  FreeElement(hash, rem);
  hash->num_items--;
  
  return 1;
}

struct hash_iterator {
  struct hash *hash;
  size_t bucket;
  struct element *elem;
};

struct hash_iterator *Hash_IteratorNew(struct hash *hash) {
  struct hash_iterator *hi;

  if ((hi = malloc(sizeof(*hi))) == NULL) {
    perror("Error: Could not allocate memory for hash table iterator");
    goto err;
  }
  memset(hi, 0, sizeof(*hi));
  
  hi->hash = hash;
  hi->bucket = -1;
  
  return hi;
  
 err:
  return NULL;
}

void Hash_IteratorFree(struct hash_iterator *hi) {
  if (hi == NULL)
    return;
  
  free(hi);
}

int Hash_IteratorNext(struct hash_iterator *hi) {
  if (hi->elem && hi->elem->next) {
    hi->elem = hi->elem->next;
    return 1;
  }
  
  while (hi->bucket < hi->hash->num_buckets || hi->bucket == SIZE_MAX) {
    if ((hi->elem = hi->hash->buckets[++hi->bucket]))
      return 1;
  }
  
  return 0;
}

void *Hash_IteratorGetKey(struct hash_iterator *hi) {
  if (hi->elem == NULL)
    return NULL;

  return hi->elem->key;
}

void *Hash_IteratorGetData(struct hash_iterator *hi) {
  if (hi->elem == NULL)
    return NULL;

  return hi->elem->data;
}

void *Hash_GetFirstKey(struct hash *hash) {
  struct hash_iterator *hi;
  void *key;

  if ((hi = Hash_IteratorNew(hash)) == NULL)
    goto err;
  if (!Hash_IteratorNext(hi))
    goto err2;
  key = Hash_IteratorGetKey(hi);
  Hash_IteratorFree(hi);
  
  return key;

 err2:
  Hash_IteratorFree(hi);
 err:
  return NULL;
}
