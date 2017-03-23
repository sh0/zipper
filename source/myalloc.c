/*

My memory allocation.

---------------------------------------------------------------

Copyright (c) 1994 The Board of Trustees of The Leland Stanford
Junior University.  All rights reserved.

Permission to use, copy, modify and distribute this software and its
documentation for any purpose is hereby granted without fee, provided
that the above copyright notice and this permission notice appear in
all copies of this software and that you do not sell the software.

THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

*/

#include <stdio.h>
#include <stdlib.h>
#include <myalloc.h>


/******************************************************************************
Allocate some memory.

Entry:
  size  - amount of memory requested (in bytes)
  lnum  - line number from which memory was requested
  fname - file name from which memory was requested
******************************************************************************/

char* my_alloc(size, lnum, fname)
int size;
int lnum;
char* fname;
{
    char* ptr;

    ptr = (char*) malloc(size);

    if (ptr == 0) {
        fprintf(stderr, "Memory allocation bombed on line %d in %s\n", lnum, fname);
    }

    return (ptr);
}


/******************************************************************************
Check to see if we can allocate anything.
******************************************************************************/

alloc_check(str)
char* str;
{
    char* ptr;

    ptr = (char*) malloc(1);

    if (ptr == 0) {
        fprintf(stderr, "can't allocate at '%s'\n", str);
    } else
        free(ptr);
}

