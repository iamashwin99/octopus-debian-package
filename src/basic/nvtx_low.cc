/*
 Copyright (C) 2019, 2021 S. Ohlmann

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

/* This is a wrapper around the NVTX (NVIDIA Tools Extension) profiling
 * functions */

#include <config.h>

#ifdef HAVE_NVTX
#include <nvToolsExt.h>
/* These are colors from the "light" qualitative color scheme
 * from https://personal.sron.nl/~pault/ */
const uint32_t colors[] = {0xff77aadd, 0xff99ddff, 0xff44bb99,
                           0xffbbcc33, 0xffaaaa00, 0xffeedd88,
                           0xffee8866, 0xffffaabb, 0xffdddddd};
const int num_colors = sizeof(colors) / sizeof(uint32_t);
#endif

#include "string_f.h" /* fortran <-> c string compatibility issues */

#include <fortran_types.h>

using namespace std;

extern "C" void FC_FUNC_(nvtx_range_push,
                         NVTX_RANGE_PUSH)(STR_F_TYPE range_name,
                                          const fint *idx STR_ARG1) {
#ifdef HAVE_NVTX
  char *range_name_c;
  TO_C_STR1(range_name, range_name_c);

  /* The code for the colored ranges is taken from a blog post by Jiri Kraus:
   * https://developer.nvidia.com/blog/cuda-pro-tip-generate-custom-application-profile-timelines-nvtx/
   */
  int color_id = *idx;
  color_id = color_id % num_colors;
  nvtxEventAttributes_t eventAttrib = {0};
  eventAttrib.version = NVTX_VERSION;
  eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib.colorType = NVTX_COLOR_ARGB;
  eventAttrib.color = colors[color_id];
  eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib.message.ascii = range_name_c;
  nvtxRangePushEx(&eventAttrib);

  free(range_name_c);
#endif
}

extern "C" void FC_FUNC_(nvtx_range_pop, NVTX_RANGE_POP)() {
#ifdef HAVE_NVTX
  nvtxRangePop();
#endif
}
