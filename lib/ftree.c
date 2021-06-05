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

#include "ftree.h"

struct ftree {
  struct ftree_node *root;
  ftree_copy_func_t copy_data;
  ftree_free_func_t free_data;
  ftree_key_func_t  dyn_key;
};

struct ftree_node {
  float key;
  void *data;
  size_t items;
  int height;
  struct ftree_node *parent;
  struct ftree_node *left;
  struct ftree_node *right;
};

struct ftree *FTree_New(ftree_copy_func_t copy_data, ftree_free_func_t free_data, ftree_key_func_t dyn_key) {
  struct ftree *ftree;

  if ((ftree = malloc(sizeof(*ftree))) == NULL)
    goto err;
  memset(ftree, 0, sizeof(*ftree));
  ftree->copy_data = copy_data;
  ftree->free_data = free_data;
  ftree->dyn_key   = dyn_key;

  if ((ftree->root = malloc(sizeof(*ftree->root))) == NULL)
    goto err2;
  memset(ftree->root, 0, sizeof(*ftree->root));
  ftree->root->items = 1;
  
  return ftree;

 err2:
  free(ftree);
 err:
  return NULL;
}

static void Node_Free(struct ftree_node *node, ftree_free_func_t free_data) {
  if (node == NULL)
    return;

  if (free_data)
    free_data(node->data);
  Node_Free(node->left, free_data);
  Node_Free(node->right, free_data);
  free(node);
}

void FTree_Free(struct ftree *ftree) {
  if (ftree == NULL)
    return;

  Node_Free(ftree->root, ftree->free_data);
  free(ftree);
}

void FTree_Clear(struct ftree *ftree) {
  Node_Free(ftree->root->left, ftree->free_data);
  memset(ftree->root, 0, sizeof(*ftree->root));
  ftree->root->items = 1;
}

#define HEIGHT(n) ((n) ? ((n)->height) : -1)
#define ITEMS(n)  ((n) ? ((n)->items) : 0)

size_t FTree_ItemCount(struct ftree *ftree) {
  return ITEMS(ftree->root->left);
}

static void FixHeight(struct ftree_node *n) {
  int lh, rh;

  if (n == NULL)
    return;
  
  lh = HEIGHT(n->left);
  rh = HEIGHT(n->right);
  n->height = (lh > rh ? lh : rh) + 1;
  
  n->items = ITEMS(n->left) + ITEMS(n->right) + 1;
}

/*    x            z            y
 *   / \          / \         /   \
 *      z   to   x      or   x     z
 *     / \      / \         / \   / \
 *    y            y          yl yr
 */
static void RotateLeft(struct ftree_node **n) {
  struct ftree_node *x, *y, *z;
  
  x = *n;
  z = x->right;
  y = z->left;
  
  if (HEIGHT(y) > HEIGHT(z->right)) {
    x->right = y->left;
    z->left  = y->right;
    y->left  = x;
    y->right = z;
    *n = y;
    y->parent = x->parent;
    x->parent = y;
    z->parent = y;
    if (x->right)
      x->right->parent = x;
    if (z->left)
      z->left->parent = z;
    FixHeight(x);
    FixHeight(z);
    FixHeight(y);
  } else {
    x->right = y;
    z->left  = x;
    *n = z;
    z->parent = x->parent;
    x->parent = z;
    if (y)
      y->parent = x;
    FixHeight(x);
    FixHeight(z);
  }
}

static void RotateRight(struct ftree_node **n) {
  struct ftree_node *x, *y, *z;
  
  x = *n;
  z = x->left;
  y = z->right;
  
  if (HEIGHT(y) > HEIGHT(z->left)) {
    x->left  = y->right;
    z->right = y->left;
    y->right = x;
    y->left  = z;
    *n = y;
    y->parent = x->parent;
    x->parent = y;
    z->parent = y;
    if (x->left)
      x->left->parent = x;
    if (z->right)
      z->right->parent = z;
    FixHeight(x);
    FixHeight(z);
    FixHeight(y);
  } else {
    x->left  = y;
    z->right = x;
    *n = z;
    z->parent = x->parent;
    x->parent = z;
    if (y)
      y->parent = x;
    FixHeight(x);
    FixHeight(z);
  }
}

static struct ftree_node **Node_Ptr(struct ftree_node *node) {
  struct ftree_node *parent = node->parent;
  
  if (parent == NULL)
    return NULL;
  
  if (parent->left == node)
    return &parent->left;
  
  return &parent->right;
}

static void Balance(struct ftree_node *node) {
  struct ftree_node **n = Node_Ptr(node);
  int lh, rh;
  
  while (n) {
    lh = HEIGHT((*n)->left);
    rh = HEIGHT((*n)->right);
    if (lh > rh + 1) {
      RotateRight(n);
    } else if (rh > lh + 1) {
      RotateLeft(n);
    } else {
      (*n)->height = (lh > rh ? lh : rh) + 1;
      (*n)->items = ITEMS((*n)->left) + ITEMS((*n)->right) + 1;
    }
    n = Node_Ptr((*n)->parent);
  }
}

static void Place(struct ftree *ftree, struct ftree_node *node, void *user) {
  struct ftree_node **cur, *parent;
  float key, cur_key;
  
  if (ftree->dyn_key)
    key = ftree->dyn_key(node->data, user);
  else
    key = node->key;
  
  parent = ftree->root;
  cur = &parent->left;
  while (*cur != NULL) {
    parent = *cur;
    
    if (ftree->dyn_key)
      cur_key = ftree->dyn_key((*cur)->data, user);
    else
      cur_key = (*cur)->key;
    
    if (key < cur_key)
      cur = &(*cur)->left;
    else
      cur = &(*cur)->right;
  }
  
  *cur = node;
  node->parent = parent;
  
  Balance(parent);
}

struct ftree_node *FTree_Insert(struct ftree *ftree, float key, void *data, void *user) {
  struct ftree_node *node;
  
  if ((node = malloc(sizeof(*node))) == NULL)
    goto err;
  memset(node, 0, sizeof(*node));

  node->key = key;
  node->items = 1;
  if (ftree->copy_data)
    node->data = ftree->copy_data(data);
  else
    node->data = data;
  
  Place(ftree, node, user);
  
  return node;
  
 err:
  return NULL;
}

static void RemoveFromTree(struct ftree *ftree, struct ftree_node *node) {
  struct ftree_node **n, *parent, *cur;
  
  if (node->right) {
    parent = node->right;
    n = &node->right;
    while ((*n)->left) {
      parent = *n;
      n = &parent->left;
    }
    cur = *n;
    if ((*n = cur->right))
      (*n)->parent = parent;
    cur->left   = node->left;
    cur->right  = node->right;
    cur->parent = node->parent;
    cur->height = node->height;
    cur->items  = node->items - 1;
    *Node_Ptr(node) = cur;
    if (cur->right)
      cur->right->parent = cur;
    if (cur->left)
      cur->left->parent = cur;
  } else {
    parent = node->parent;
    n = Node_Ptr(node);
    if ((*n = node->left))
      (*n)->parent = parent;
  }
  Balance(parent);
}

void FTree_Delete(struct ftree *ftree, struct ftree_node *node) {
  RemoveFromTree(ftree, node);
  
  if (ftree->free_data)
    ftree->free_data(node->data);
  free(node);  
}

void FTree_Rekey(struct ftree *ftree, struct ftree_node *node, float new_key, void *user) {
  RemoveFromTree(ftree, node);
  
  node->key = new_key;
  node->height = 0;
  node->items  = 1;
  node->left   = NULL;
  node->right  = NULL;
  
  Place(ftree, node, user);
}

struct ftree_node *FTree_Lowest(struct ftree *tree) {
  struct ftree_node *node = tree->root->left;
  
  if (node == NULL)
    return NULL;
  
  while (node->left) {
    node = node->left;
  }
  
  return node;
}

struct ftree_node *FTree_Highest(struct ftree *tree) {
  struct ftree_node *node = tree->root->left;
  
  if (node == NULL)
    return NULL;
  
  while (node->right) {
    node = node->right;
  }
  
  return node;
}

struct ftree_node *FTree_Next(struct ftree *tree, struct ftree_node *node) {
  struct ftree_node *prev;
  
  if (node->right) {
    node = node->right;
  
    while (node->left != NULL)
      node = node->left;
  
    return node;
  }
  
  prev = NULL;
  while (node->right == prev) {
    prev = node;
    node = node->parent;
    if (node->parent == NULL)
      return NULL;
  }
  
  return node;
}

struct ftree_node *FTree_Prev(struct ftree *tree, struct ftree_node *node) {
  struct ftree_node *prev;
  
  if (node->left) {
    node = node->left;
  
    while (node->right != NULL)
      node = node->right;
  
    return node;
  }
  
  prev = NULL;
  while (node->left == prev) {
    prev = node;
    node = node->parent;
    if (node->parent == NULL)
      return NULL;
  }
  
  return node;
}

struct ftree_node *FTree_Median(struct ftree *tree) {
  size_t target, pos;
  struct ftree_node *cur;
  
  target = FTree_ItemCount(tree) / 2;
  cur = tree->root->left;

  while (cur) {
    pos = ITEMS(cur->left);
    if (target == pos)
      return cur;

    if (target < pos) {
      cur = cur->left;
    } else {
      target -= pos + 1;
      cur = cur->right;
    }
  }
  
  return NULL;
}

float FTree_GetKey(struct ftree_node *node) {
  return node->key;
}

void *FTree_GetData(struct ftree_node *node) {
  return node->data;
}

float FTree_GetKeyDyn(struct ftree *tree, struct ftree_node *node, void *user) {
  if (tree->dyn_key)
    return tree->dyn_key(node->data, user);
  
  return node->key;
}

/* Debug */
static void CheckHeight(const struct ftree_node *node) {
  int lh, rh, hh, items;
  
  lh = HEIGHT(node->left);
  rh = HEIGHT(node->right);
  hh = (lh > rh ? lh : rh) + 1;
  if (node->height != hh)
    printf("Incorrect height %d -> (%d, %d)\n", node->height, lh, rh);
  
  if (lh > rh + 1 || rh > lh + 1)
    printf("Tree out of balance (%d, %d)\n", lh, rh);

  items = ITEMS(node->left) + ITEMS(node->right) + 1;
  if (items != node->items)
    printf("Incorrect count %zd -> (%zd, %zd)\n", node->items, ITEMS(node->left), ITEMS(node->right));
}

static void CheckNode(const struct ftree_node *node, const struct ftree_node *parent, void (*check_data)(void *)) {
  if (node == NULL)
    return;
  
  CheckHeight(node);
  
  if (node->parent != parent)
    printf("Incorrect parent %p != %p\n", node->parent, parent);
  
  if (check_data)
    check_data(node->data);
}

void FTree_Check(const struct ftree *ftree, void (*check_data)(void *)) {
  if (ftree->root == NULL) {
    printf("Invalid NULL root\n");
    return;
  }

  if (ftree->root->parent != NULL)
    printf("Invalid root->parent = %p\n", ftree->root->parent);
  
  if (ftree->root->key != 0)
    printf("Invalid root->key = %g\n", ftree->root->key);

  if (ftree->root->data != NULL)
    printf("Invalid root->data = %p\n", ftree->root->data);
  
  if (ftree->root->right != NULL)
    printf("Invalid root->right = %p\n", ftree->root->right);

  CheckNode(ftree->root->left, ftree->root, check_data);
}

static int HasData(const struct ftree_node *node, void *data) {
  if (node == NULL)
    return 0;

  if (node->data == data)
    return 1;

  if (HasData(node->left, data))
    return 1;

  return HasData(node->right, data);
}

int FTree_HasData(const struct ftree *ftree, void *data) {
  return HasData(ftree->root->left, data);
}
