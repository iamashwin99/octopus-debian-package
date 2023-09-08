/*
 Copyright (C) 2019,2023 S. Ohlmann

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/
#include "fortran_types.h"
#include "vectors.h"
#include <config.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_CUDA
#include <cuda.h>

#define CUDA_SAFE_CALL(x)                                                      \
  do {                                                                         \
    CUresult result = x;                                                       \
    if (result != CUDA_SUCCESS) {                                              \
      const char *msg;                                                         \
      cuGetErrorName(result, &msg);                                            \
      std::cerr << "\nerror: " #x " failed with error " << msg << '\n';        \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#endif

#ifdef __GLIBC__
#include <malloc.h>
static int initialized = 0;
static size_t threshold = 102400; // initial threshold at 100 kB
#endif

using namespace std;

void *allocate_aligned(fint8 size_bytes) {
#ifdef __GLIBC__
  // do this only for glibc (which is default on Linux)
  if(!initialized || (size_t)size_bytes < threshold) {
    // set threshold to use mmap to less than the size to be allocated
    // to make sure that all batches are allocated using mmap
    // this circumvents heap fragmentation problems with aligned memory
    // also make sure we set this at least once in the beginning
    threshold = min(threshold, (size_t)(0.9*size_bytes));
    int success = mallopt(M_MMAP_THRESHOLD, threshold);
    if(!success) {
      printf("Error setting mmap threshold option!\n");
      printf("You might run into an out-of-memory condition because of heap fragmentation.\n");
      printf("This can be caused by allocating aligned memory with posix_memalign in between other allocations.\n");
    }
    initialized = 1;
  }
#endif
#ifdef DEBUG_ALLOC
  printf("Allocating %d bytes, unpinned.\n", (size_t)size_bytes);
#endif
  void *aligned;
  int status;
  // align on vector size to improve vectorization
  status = posix_memalign(&aligned, (size_t)sizeof(double) * VEC_SIZE,
                          (size_t)size_bytes);
  if (status != 0) {
    printf("Error allocating aligned memory!\n");
    return NULL;
  }
  return aligned;
}

extern "C" void *dallocate_aligned(fint8 size) {
  return allocate_aligned(sizeof(double) * size);
}

extern "C" void *zallocate_aligned(fint8 size) {
  return allocate_aligned(sizeof(double) * 2 * size);
}

extern "C" void *sallocate_aligned(fint8 size) {
  return allocate_aligned(sizeof(float) * size);
}

extern "C" void *callocate_aligned(fint8 size) {
  return allocate_aligned(sizeof(float) * 2 * size);
}

extern "C" void deallocate_aligned(void *array) {
#ifdef DEBUG_ALLOC
  printf("Deallocating unpinned.\n");
#endif
  free(array);
}

void *allocate_pinned(fint8 size_bytes) {
#ifdef HAVE_CUDA
#ifdef DEBUG_ALLOC
  printf("Allocating %d bytes, pinned.\n", (size_t)size_bytes);
#endif
  void *pinned;
  CUDA_SAFE_CALL(cuMemAllocHost(&pinned, (size_t)size_bytes));
  return pinned;
#else
  printf("Error! Pinned memory requested, although CUDA not available. "
         "Returning aligned memory.\n");
  return allocate_aligned(size_bytes);
#endif
}

extern "C" void *dallocate_pinned(fint8 size) {
  return allocate_pinned(sizeof(double) * size);
}

extern "C" void *zallocate_pinned(fint8 size) {
  return allocate_pinned(sizeof(double) * 2 * size);
}

extern "C" void *sallocate_pinned(fint8 size) {
  return allocate_pinned(sizeof(float) * size);
}

extern "C" void *callocate_pinned(fint8 size) {
  return allocate_pinned(sizeof(float) * 2 * size);
}

extern "C" void deallocate_pinned(void *array) {
#ifdef HAVE_CUDA
#ifdef DEBUG_ALLOC
  printf("Deallocating pinned.\n");
#endif
  CUDA_SAFE_CALL(cuMemFreeHost(array));
#endif
}
