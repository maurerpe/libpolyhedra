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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <string.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_BCRYPT_H
#include <Bcrypt.h>
#endif

#ifdef HAVE_PTHREADS
#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#endif
#define MUTEX pthread_mutex_t
#define MUTEX_INIT PTHREAD_MUTEX_INITIALIZER
#elif HAVE_WINDOWS_H
#include <Windows.h>
#define MUTEX HANDLE
#define MUTEX_INIT NULL
#endif

#include "random.h"
#include "SipHash/siphash.h"

static MUTEX mutex = MUTEX_INIT;
static unsigned char key[16];
static uint64_t counter = -1;

void Random_Init(void) {
#ifdef HAVE_PTHREADS
#else
#ifdef HAVE_CREATEMUTEXA
  if ((mutex = CreateMutexA(NULL, FALSE, NULL)) == NULL) {
    fprintf(stderr, "Could not create mutex for random number generator\n");
    exit(1);
  }
#endif
#endif
}

static void Random_Seed(void) {
#ifdef HAVE_GETENTROPY
  if (getentropy(key, sizeof(key)) < 0) {
    fprintf(stderr, "Could not get entropy to seed random number generator\n");
    exit(1);
  }
#else
#ifdef HAVE_BCRYPTGENRANDOM
  if (BCryptGenRandom(NULL, key, sizeof(key), BCRYPT_USE_SYSTEM_PREFERRED_RNG) != STATUS_SUCCESS) {
    fprintf(stderr, "Could not get bcrypt seed for random number generator\n");
    exit(1);
  }
#else
  FILE *in;
  
  if ((in = fopen("/dev/random", "rb")) == NULL) {
    fprintf(stderr, "Could not open /dev/random to seed random number generator\n");
    exit(1);
  }
  if (fread(key, sizeof(key), 1, in) <= 0) {
    fprintf(stderr, "Could not read from /dev/random to seed random number generator\n");
    exit(1);
  }
  fclose(in);
#endif
#endif
  
  counter = 0;
}

#ifdef HAVE_PTHREADS
static void Mutex_Lock(MUTEX mutex) {
  pthread_mutex_lock(&mutex);
}
static void Mutex_Unlock(MUTEX mutex) {
  pthread_mutex_unlock(&mutex);
}
#else
#ifdef HAVE_CREATEMUTEXA
static void Mutex_Lock(MUTEX mutex) {
  WaitForSingleObject(mutex, INFINITE);
}
static void Mutex_Unlock(MUTEX mutex) {
  ReleaseMutex(mutex);
}
#endif
#endif

uint64_t Random_Integer(void) {
  uint64_t data;

  Random(&data, sizeof(data));
  
  return data;
}

static inline uint8_t RandByte(void) {
  uint64_t hash;

  hash = siphash(key, (unsigned char *) &counter, sizeof(counter));
  counter++;
  
  hash ^= hash >> 32;
  hash ^= hash >> 16;
  hash ^= hash >> 8;
  
  return hash;
}

int Random(void *data, size_t len) {
  void *end = data + len;
  
#ifndef HAVE_PTHREADS
  if (mutex == MUTEX_INIT) {
    fprintf(stderr, "Warning: random number generator not initialzied properly\n");
    Random_Init();
  }
#endif
  
  Mutex_Lock(mutex);
  
  if (counter == -1)
    Random_Seed();
  
  while (data < end) {
    *((uint8_t *) data) = RandByte();
    data++;
  }
  Mutex_Unlock(mutex);
  
  return 0;
}
