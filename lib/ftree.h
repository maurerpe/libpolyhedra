/*
  Copyright (C) 2020-2021 Paul Maurer

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

#ifndef LP_FTREE_H
#define LP_FTREE_H

struct ftree;
struct ftree_node;

typedef void *(*ftree_copy_func_t)(const void *);
typedef void (*ftree_free_func_t)(void *);
typedef float (*ftree_key_func_t)(void *, void *);

struct ftree *FTree_New(ftree_copy_func_t copy_data, ftree_free_func_t free_data, ftree_key_func_t dyn_key);
void FTree_Free(struct ftree *ftree);
void FTree_Clear(struct ftree *ftree);

size_t FTree_ItemCount(struct ftree *ftree);

struct ftree_node *FTree_Insert(struct ftree *ftree, float key, void *data, void *user);
void FTree_Delete(struct ftree *ftree, struct ftree_node *node);
void FTree_Rekey(struct ftree *ftree, struct ftree_node *node, float new_key, void *user);

struct ftree_node *FTree_Lowest(struct ftree *tree);
struct ftree_node *FTree_Highest(struct ftree *tree);
struct ftree_node *FTree_Next(struct ftree *tree, struct ftree_node *node);
struct ftree_node *FTree_Prev(struct ftree *tree, struct ftree_node *node);
struct ftree_node *FTree_Median(struct ftree *tree);
float FTree_GetKey(struct ftree_node *node);
void *FTree_GetData(struct ftree_node *node);
float FTree_GetKeyDyn(struct ftree *tree, struct ftree_node *node, void *user);

typedef void (*ftree_map_func)(struct ftree_node *, void *);
void FTree_Map(struct ftree *ftree, ftree_map_func map_func, void *user);

/* Debug */
void FTree_Check(const struct ftree *ftree, void (*check_data)(void *));
int  FTree_HasData(const struct ftree *ftree, void *data);

#endif
