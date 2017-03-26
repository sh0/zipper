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

#ifndef ZIPPER_PLY_WRAPPER_H
#define ZIPPER_PLY_WRAPPER_H

// Internal
#include "zipper.h"
#include "matrix.h"
#include "raw.h"

// Parameters
void set_range_data_sigma_factor(float factor);
float get_range_data_sigma_factor();
void set_range_data_min_intensity(float intensity);
float get_range_data_min_intensity();
void set_range_data_horizontal_erode(int erode);
int get_range_data_horizontal_erode();

// Declarations
void write_ply(Scan* sc, char* filename, int writeInfo);
int is_range_grid_file(char* filename);
int read_ply(char* filename);
RangeData* read_ply_geom(char* name);
int erode_forward(RangeData* plydata, char* continuity, int xx, int yy, int erodeMax);
int erode_backward(RangeData* plydata, char* continuity, int xx, int yy, int erodeMax);
int decide_continuity(RangeData* plydata, int xx, int yy);
void delete_ply_geom(RangeData* plydata);

#endif
