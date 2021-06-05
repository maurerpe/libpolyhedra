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

#include "queue.h"

struct elem {
  void *value;
  struct elem *next;
};

struct queue {
  struct elem *head;
  struct elem *tail;
  size_t len;
};

struct queue *Queue_New(void) {
  struct queue *queue;

  if ((queue = malloc(sizeof(*queue))) == NULL)
    goto err;
  memset(queue, 0, sizeof(*queue));
  
  return queue;
  
 err:
  return NULL;
}

void Queue_Free(struct queue *queue) {
  if (queue == NULL)
    return;

  Queue_Clear(queue);
  free(queue);
}

void Queue_Clear(struct queue *queue) {
  struct elem *cur, *next;
  
  cur = queue->head;
  while (cur) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  
  memset(queue, 0, sizeof(*queue));
}

size_t Queue_Length(const struct queue *queue) {
  return queue->len;
}

int Queue_Push(struct queue *queue, void *value) {
  struct elem *elem;

  if ((elem = malloc(sizeof(*elem))) == NULL)
    goto err;
  memset(elem, 0, sizeof(*elem));
  
  elem->value = value;
  elem->next = queue->head;
  queue->head = elem;
  if (queue->len++ == 0)
    queue->tail = elem;
  
  return 0;
  
 err:
  return -1;
}

int Queue_PushBack(struct queue *queue, void *value) {
  struct elem *elem;
  
  if ((elem = malloc(sizeof(*elem))) == NULL)
    goto err;
  memset(elem, 0, sizeof(*elem));

  elem->value = value;
  if (queue->tail)
    queue->tail->next = elem;
  else
    queue->head = elem;
  queue->tail = elem;
  queue->len++;
  
  return 0;
  
 err:
  return -1;
}

void *Queue_Pop(struct queue *queue) {
  struct elem *elem;
  void *value;
  
  if ((elem = queue->head) == NULL)
    return NULL;
  
  value = elem->value;
  queue->head = elem->next;
  free(elem);
  if (--queue->len == 0)
    queue->tail = NULL;
  
  return value;
}

void *Queue_Peak(struct queue *queue) {
  if (queue->head == NULL)
    return NULL;
  
  return queue->head->value;
}

void *Queue_PeakBack(struct queue *queue) {
  if (queue->tail == NULL)
    return NULL;
  
  return queue->tail->value;
}
