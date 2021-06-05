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

#include <string.h>

#include "hash.h"
#include "unique_queue.h"

struct uq_list {
  void *value;
  struct uq_list *next;
};

struct unique_queue {
  struct hash *hash;
  struct uq_list *head;
  struct uq_list *tail;
};

#define PRESENT ((void *) 1)

struct unique_queue *UniqueQueue_New(void) {
  struct unique_queue *uq;
  
  if ((uq = malloc(sizeof(*uq))) == NULL)
    goto err;
  memset(uq, 0, sizeof(*uq));
  
  if ((uq->hash = Hash_NewPtr(NULL, NULL, NULL, NULL, NULL)) == NULL)
    goto err2;
  
  return uq;

 err2:
  free(uq);
 err:
  return NULL;
}

void UniqueQueue_Free(struct unique_queue *uq) {
  if (uq == NULL)
    return;
  
  UniqueQueue_Clear(uq);
  Hash_Free(uq->hash);
  free(uq);
}

void UniqueQueue_Clear(struct unique_queue *uq) {
  struct uq_list *cur, *next;
  
  cur = uq->head;
  while (cur) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  
  uq->head = NULL;
  uq->tail = NULL;
  Hash_Clear(uq->hash);
}

size_t UniqueQueue_NumEntries(const struct unique_queue *uq) {
  return Hash_NumEntries(uq->hash);
}

int UniqueQueue_Push(struct unique_queue *uq, void *value) {
  struct uq_list *elem;
  
  if (Hash_Lookup(uq->hash, value, NULL))
    return 0;
  
  if (Hash_Insert(uq->hash, value, PRESENT, NULL) < 0)
    goto err;
  
  if ((elem = malloc(sizeof(*elem))) == NULL)
    goto err2;
  memset(elem, 0, sizeof(*elem));
  
  elem->value = value;
  elem->next  = uq->head;
  uq->head = elem;
  if (uq->tail == NULL)
    uq->tail = elem;
  
  return 1;

 err2:
  Hash_Remove(uq->hash, value);
 err:
  return -1;
}

int UniqueQueue_PushBack(struct unique_queue *uq, void *value) {
  struct uq_list *elem;
  
  if (Hash_Lookup(uq->hash, value, NULL))
    return 0;
  
  if (Hash_Insert(uq->hash, value, PRESENT, NULL) < 0)
    goto err;
  
  if ((elem = malloc(sizeof(*elem))) == NULL)
    goto err2;
  memset(elem, 0, sizeof(*elem));
  
  elem->value = value;
  elem->next  = NULL;
  if (uq->tail) {
    uq->tail->next = elem;
  } else {
    uq->head = elem;
  }
  uq->tail = elem;
  
  return 1;

 err2:
  Hash_Remove(uq->hash, value);
 err:
  return -1;
}

void *UniqueQueue_Pop(struct unique_queue *uq) {
  struct uq_list *elem;
  void *value;
  
  if (uq->head == NULL)
    return NULL;
  
  elem = uq->head;
  uq->head = elem->next;
  value = elem->value;
  if (uq->head == NULL)
    uq->tail = NULL;
  
  free(elem);
  Hash_Remove(uq->hash, value);
  
  return value;
}
