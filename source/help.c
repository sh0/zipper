/*

Help information about the program.

Greg Turk, August 1994

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
#include <math.h>
#include <string.h>
#include "strings.h"
#include "malloc.h"
#include "cyfile.h"
#include "zipper.h"

static char *help_strings[] = {
"quit - quit the program",
"help or ?   - print help about commands",
"help <word> - print any help line that contains a given word",
"",
"anchor <mesh> - specify the anchor mesh for alignment",
"next_align <mesh> - specify the next scan to align",
"align  [<mesh>] - align a mesh with the current anchor mesh",
"align_step <mesh> - align a mesh with the current anchor mesh",
"target <mesh> - specify the target mesh to merge into",
"next_merge <mesh> - specify the next scan to eat/merge",
"merge  [<mesh>] - merge a mesh into the target mesh",
"consensus <mesh> <res> - compute consensus geometry at a given resolution",
"eat_edges  <mesh>  - eat away edges of one mesh (first part of merge)",
"only_merge <mesh>  - perform zippering only (second part of merge)",
"fast_ops <off/on>  - use faster operations? (doesn't draw much)",
"parallel_max <num> - use up to this many processors to do eat & align",
"",
"undo - undo last movement operation",
"redo - redo last undo",
"beep - make a beep sound",
"mesh_resolution <mesh> <res> - set the mesh resolution (1, 2, 3, or 4)",
"remove_mesh <mesh> - remove a mesh from memory",
"",
"bmesh    <file> [angle]  - read in a binary mesh file",
"mesh     <file> [angle]  - read in Cyberware mesh file (not used much)",
"bpolygon <file> [angle]  - read in a PLY polygon file",
"bpwrite  <mesh> <file>   - write a mesh to a file",
"color_write <off/on>     - flag: write colors to file?",
"intensity_write <off/on> - flag: write intensity to file?",
"normals_write <off/on>   - flag: write normals & consensus tags to file?",
"",
"rename_mesh <old-mesh-name> <new-mesh-name> -  renames a mesh",
"remove_mesh <mesh> -  removes a mesh",
"",
"print         - print out positions of meshes",
"print <file>  - write out positions of meshes to a file",
"mprint        - print out matrices of meshes",
"mprint <file> - write out matrices of meshes to a file",
"source <file> - read in lines of a text file and execute them",
"",
"translate <mesh> <x> <y> <z> - translate a mesh",
"xrotate <mesh> <angle> - rotate a mesh by some degrees around the x-axis",
"yrotate <mesh> <angle> - rotate a mesh by some degrees around the y-axis",
"zrotate <mesh> <angle> - rotate a mesh by some degrees around the z-axis",
"qrotate <mesh> <q0> <q1> <q2> <q3> - rotate a mesh using a quaternion",
"speed <num>   - speed of object motion due to mouse (default = 1)",
"cull <off/on> - turn off/on backface rejection in line mode (default is on)",
"reset_view    - reset the view position",
"axes <off/on> - turn off or on drawing of axes",
"camera <cam_params> - set the camera position",
"background <r> <g> <b> - set the background color",
"singlebuffer <off/on> - choose between single- and double-buffering",
"orthographic <off/on> - choose between orthographic and perspective viewing",
};

static int nhelp = sizeof(help_strings) / sizeof(char*);


/******************************************************************************
Print out help information.
******************************************************************************/

help(word)
  char *word;
{
  int i;
  int found = 0;

  /* if command is null, print out all commands */

  if (word == NULL) {
    printf ("\n");
    for (i = 0; i < nhelp; i++)
      printf ("%s\n", help_strings[i]);
    printf ("\n");
    return;
  }

  /* find any mention of a given word in the help */

  printf ("\n");

  for (i = 0; i < nhelp; i++) {
    if (help_strings[i][0] != '\0' && strstr (help_strings[i], word) != NULL) {
      printf ("%s\n", help_strings[i]);
      found = 1;
    }
  }

  if (!found)
    printf ("'%s' not found in help\n", word);

  printf ("\n");
}

