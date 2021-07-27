/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef SWIFT_LIGHTCONE_MAP_H
#define SWIFT_LIGHTCONE_MAP_H

/* Define this to write the expected sum of the map to the
   output file. This is to check that the SPH smoothing and
   communication code is conserving the quantity added to
   the map. */
#define LIGHTCONE_MAP_CHECK_TOTAL

/* Standard headers */
#include <math.h>
#include <limits.h>

/* Config parameters. */
#include "../config.h"

/* HDF5 */
#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

/* Local headers */
#include "units.h"

/**
 * @brief Struct to store a single lightcone healpix map
 */
struct lightcone_map {

  /*! Healpix nside parameter */
  int nside;

  /*! Total pixels in the map */
  size_t total_nr_pix;

  /*! Number of pixels stored on this node */
  size_t local_nr_pix;
  
  /*! Offset of the first pixel stored on this rank */
  size_t local_pix_offset;

  /*! Number of pixels per rank (last node has any extra) */
  size_t pix_per_rank;

  /*! Local healpix map data */
  double *data;

  /*! Inner radius */
  double r_min;

  /*! Outer radius */
  double r_max;

  /*! Units of this map */
  enum unit_conversion_factor units;

#ifdef LIGHTCONE_MAP_CHECK_TOTAL
  /*! Total quantity accumulated to this map, for consistency check */
  double total;
#endif

};


void lightcone_map_init(struct lightcone_map *map, const int nside, const size_t total_nr_pix,
                        const size_t pix_per_rank, const size_t local_nr_pix,
                        const size_t local_pix_offset, const double r_min, const double r_max,
                        enum unit_conversion_factor units);

void lightcone_map_clean(struct lightcone_map *map);

void lightcone_map_struct_dump(const struct lightcone_map *map, FILE *stream);

void lightcone_map_struct_restore(struct lightcone_map *map, FILE *stream);

void lightcone_map_allocate_pixels(struct lightcone_map *map, const int zero_pixels);

void lightcone_map_free_pixels(struct lightcone_map *map);

#ifdef HAVE_HDF5
void lightcone_map_write(struct lightcone_map *map, const hid_t loc_id, const char *name,
                         const struct unit_system *internal_units,
                         const struct unit_system *snapshot_units,
                         const int collective);
#endif

#endif /* #ifndef SWIFT_LIGHTCONE_MAP_H */
