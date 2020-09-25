/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 James Willis (james.s.willis@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

#ifdef WITH_FOF

/* Some standard headers. */
#include <errno.h>
#include <libgen.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "fof6d.h"

/* Local headers. */
#include "fof.h"
#include "black_holes.h"
#include "common_io.h"
#include "engine.h"
#include "hashmap.h"
#include "memuse.h"
#include "proxy.h"
#include "threadpool.h"

void fof6d_calc_vel_disp(struct fof_props *props, const size_t num_parts_in_groups, struct space *s) {

  const int num_groups = props->num_groups;
  struct gpart *gparts = s->gparts;
  const size_t nr_gparts = s->nr_gparts;
  //size_t *group_index = props->group_index;
  //size_t *group_size = props->group_size;
  double *group_mass = props->group_mass;
  double *v_disp = NULL;
  double *v_mean[3] = {NULL, NULL, NULL};
  size_t *part_index = NULL;

  /* Allocate and initialise a velocity dispersion array. */
  if (swift_memalign("6dfof_v_disp", (void **)&v_disp, 64,
                     num_groups * sizeof(double)) != 0)
    error("Failed to allocate list of group velocity dispersions for 6DFOF search.");

  if (swift_memalign("6dfof_v_mean", (void **)&v_mean, 64,
                     num_groups * 3 * sizeof(double)) != 0)
    error("Failed to allocate list of mean velocity for 6DFOF search.");

  if (swift_memalign("6dfof_part_index", (void **)&part_index, 64,
                     num_parts_in_groups * sizeof(size_t)) != 0)
    error("Failed to allocate list of particle indices in groups for 6DFOF search.");

  bzero(v_disp, num_groups * sizeof(double));
  bzero(v_mean, num_groups * sizeof(double));
  bzero(part_index, num_parts_in_groups * sizeof(size_t));
  
  size_t part_ctr = 0;

  /* Calculate the mean velocity for each group. */
  for (size_t i = 0; i < nr_gparts; i++) {
    
    const size_t group_id = gparts[i].fof_data.group_id;

    if(group_id != fof_props_default_group_id) {
       v_mean[group_id][0] += gparts[i].mass * gparts[i].v_full[0]; 
       v_mean[group_id][1] += gparts[i].mass * gparts[i].v_full[1]; 
       v_mean[group_id][2] += gparts[i].mass * gparts[i].v_full[2];

       /* JSW TODO: Could be calculated in fof_search_tree */
       part_index[part_ctr++] = i;
    } 
  }
  
  for (int i = 0; i < num_groups; i++) {
    const double one_over_mass = 1.0 / group_mass[i];

    v_mean[i][0] *= one_over_mass; 
    v_mean[i][1] *= one_over_mass; 
    v_mean[i][2] *= one_over_mass;
  }

  /* Calculate the velocity dispersion for each group. */
  for (size_t i = 0; i < num_parts_in_groups; i++) {
    
    const size_t index = part_index[i];
    const size_t group_id = gparts[index].fof_data.group_id;

    const double v_diff[3] = {gparts[index].v_full[0] - v_mean[group_id][0],
                              gparts[index].v_full[1] - v_mean[group_id][1],
                              gparts[index].v_full[2] - v_mean[group_id][2]};

    v_disp[group_id] += (v_diff[0] * v_diff[0] + 
                         v_diff[1] * v_diff[1] + 
                         v_diff[2] * v_diff[2]) * gparts[index].mass;
  }

  for (int i = 0; i < num_groups; i++) {
    v_disp[i] /= group_mass[i];
  }

}

#endif /* WITH_FOF6D */