/*
 * Copyright (c) 1995-2017, Stanford University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of Stanford University nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY STANFORD UNIVERSITY ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL STANFORD UNIVERSITY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ZIPPER_PLY_H
#define ZIPPER_PLY_H

// External
#include <stdio.h>
#include <stddef.h>

#define PLY_ASCII      1        /* ascii PLY file */
#define PLY_BINARY_BE  2        /* binary PLY file, big endian */
#define PLY_BINARY_LE  3        /* binary PLY file, little endian */

#define PLY_OKAY    0           /* ply routine worked okay */
#define PLY_ERROR  -1           /* error in ply routine */

/* scalar data types supported by PLY format */
#define PLY_START_TYPE 0
#define PLY_CHAR       1
#define PLY_SHORT      2
#define PLY_INT        3
#define PLY_UCHAR      4
#define PLY_USHORT     5
#define PLY_UINT       6
#define PLY_FLOAT      7
#define PLY_DOUBLE     8
#define PLY_END_TYPE   9

#define  PLY_SCALAR  0
#define  PLY_LIST    1

// Description of a property
typedef struct PlyProperty {
    char* name;                           /* property name */
    int external_type;                    /* file's data type */
    int internal_type;                    /* program's data type */
    int offset;                           /* offset bytes of prop in a struct */

    int is_list;                          /* 1 = list, 0 = scalar */
    int count_external;                   /* file's count type */
    int count_internal;                   /* program's count type */
    int count_offset;                     /* offset byte for list count */
} PlyProperty;

// Description of an element
typedef struct PlyElement {
    char* name;                   /* element name */
    int num;                      /* number of elements in this object */
    int size;                     /* size of element (bytes) or -1 if variable */
    int nprops;                   /* number of properties for this element */
    PlyProperty** props;          /* list of properties in the file */
    char* store_prop;             /* flags: property wanted by user? */
    int other_offset;             /* offset to un-asked-for props, or -1 if none*/
    int other_size;               /* size of other_props structure */
} PlyElement;

// Describes other properties in an element
typedef struct PlyOtherProp {
    char* name;                   /* element name */
    int size;                     /* size of other_props */
    int nprops;                   /* number of properties in other_props */
    PlyProperty** props;          /* list of properties in other_props */
} PlyOtherProp;

// For storing other_props for an other element
typedef struct OtherData {
    void* other_props;
} OtherData;

// Data for one "other" element
typedef struct OtherElem {
    char* elem_name;             /* names of other elements */
    int elem_count;              /* count of instances of each element */
    OtherData** other_data;      /* actual property data for the elements */
    PlyOtherProp* other_props;   /* description of the property data */
} OtherElem;

// "Other" elements, not interpreted by user
typedef struct PlyOtherElems {
    int num_elems;                /* number of other elements */
    OtherElem* other_list;        /* list of data for other elements */
} PlyOtherElems;

// Description of PLY file
typedef struct PlyFile {
    FILE* fp;                     /* file pointer */
    int file_type;                /* ascii or binary */
    float version;                /* version number of file */
    int nelems;                   /* number of elements of object */
    PlyElement** elems;           /* list of elements */
    int num_comments;             /* number of comments */
    char** comments;              /* list of comments */
    int num_obj_info;             /* number of items of object information */
    char** obj_info;              /* list of object info items */
    PlyElement* which_elem;       /* which element we're currently writing */
    PlyOtherElems* other_elems;   /* "other" elements from a PLY file */
} PlyFile;

// Declarations
PlyFile* ply_write(FILE*, int, char**, int);
PlyFile* ply_open_for_writing(char*, int, char**, int, float*);
void ply_describe_element(PlyFile*, char*, int, int, PlyProperty*);
void ply_describe_property(PlyFile*, char*, PlyProperty*);
void ply_element_count(PlyFile*, char*, int);
void ply_header_complete(PlyFile*);
void ply_put_element_setup(PlyFile*, char*);
void ply_put_element(PlyFile*, void*);
void ply_put_comment(PlyFile*, char*);
void ply_put_obj_info(PlyFile*, char*);
PlyFile* ply_read(FILE*, int*, char***);
PlyFile* ply_open_for_reading(char*, int*, char***, int*, float*);
PlyProperty** ply_get_element_description(PlyFile*, char*, int*, int*);
void ply_get_element_setup(PlyFile*, char*, int, PlyProperty*);
void ply_get_property(PlyFile*, char*, PlyProperty*);
PlyOtherProp* ply_get_other_properties(PlyFile*, char*, int);
void ply_get_element(PlyFile*, void*);
char** ply_get_comments(PlyFile*, int*);
char** ply_get_obj_info(PlyFile*, int*);
void ply_close(PlyFile*);
void ply_get_info(PlyFile*, float*, int*);
PlyOtherElems* ply_get_other_element(PlyFile*, char*, int);
void ply_describe_other_elements(PlyFile*, PlyOtherElems*);
void ply_put_other_elements(PlyFile*);
void ply_free_other_elements(PlyOtherElems*);

int equal_strings(char*, char*);

#endif
