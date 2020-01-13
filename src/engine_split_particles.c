/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* This object's header. */
#include "engine.h"

/* Local headers */
#include "active.h"
#include "black_holes_struct.h"
#include "chemistry.h"
#include "cooling.h"
#include "hydro.h"
#include "random.h"
#include "star_formation.h"
#include "tracers.h"

const int particle_split_factor = 2;
const double displacement_factor = 0.2;

/**
 * @brief Data structure used by the counter mapper function
 */
struct data_count {
  const struct engine *const e;
  const float mass_threshold;
  size_t counter;
};

/**
 * @brief Data structure used by the split mapper function
 */
struct data_split {
  const struct engine *const e;
  const float mass_threshold;
  size_t *const k_parts;
  size_t *const k_gparts;
  swift_lock_type lock;
};

/**
 * @brief Mapper function to count the number of #part above the mass threshold
 * for splitting.
 */
void engine_split_gas_particle_count_mapper(void *restrict map_data, int count,
                                            void *restrict extra_data) {

  /* Unpack the data */
  struct part *parts = (struct part *)map_data;
  struct data_count *data = (struct data_count *)extra_data;
  const struct engine *e = data->e;
  const float mass_threshold = data->mass_threshold;

  size_t counter = 0;

  for (int i = 0; i < count; ++i) {

    const struct part *p = &parts[i];

    /* Ignore inhibited particles */
    if (part_is_inhibited(p, e)) continue;

    /* Is the mass of this particle larger than the threshold? */
    const float gas_mass = hydro_get_mass(p);
    if (gas_mass > mass_threshold) ++counter;
  }

  /* Increment the global counter */
  atomic_add(&data->counter, counter);
}

/**
 * @brief Mapper function to split the #part above the mass threshold.
 */
void engine_split_gas_particle_split_mapper(void *restrict map_data, int count,
                                            void *restrict extra_data) {

  /* Unpack the data */
  struct part *parts = (struct part *)map_data;
  struct data_split *data = (struct data_split *)extra_data;
  const float mass_threshold = data->mass_threshold;
  const struct engine *e = data->e;
  const int with_gravity = (e->policy & engine_policy_self_gravity) ||
                           (e->policy & engine_policy_external_gravity);

  /* Additional thread-global particle arrays */
  const struct space *s = e->s;
  struct part *global_parts = s->parts;
  struct xpart *global_xparts = s->xparts;
  struct gpart *global_gparts = s->gparts;

  /* xpart array corresponding to the thread-local particle array */
  const ptrdiff_t offset = parts - global_parts;
  struct xpart *xparts = global_xparts + offset;

  /* Loop over the chunk of the part array assigned to this thread */
  for (int i = 0; i < count; ++i) {

    /* Is the mass of this particle larger than the threshold? */
    struct part *p = &parts[i];

    /* Ignore inhibited particles */
    if (part_is_inhibited(p, e)) continue;

    const float gas_mass = hydro_get_mass(p);
    const float h = p->h;

    /* Found a particle to split */
    if (gas_mass > mass_threshold) {

      /* Make sure only one thread is here */
      lock_lock(&data->lock);

      /* Where should we put the new particle in the array? */
      const size_t k_parts = *data->k_parts;
      const size_t k_gparts = *data->k_gparts;

      /* Make the next slot available for the next particle */
      (*data->k_parts)++;
      if (with_gravity) (*data->k_gparts)++;

      /* Release the lock and continue in parallel */
      if (lock_unlock(&data->lock) != 0)
        error("Impossible to unlock particle splitting");

      /* We now know where to put the new particle we are going to create. */

      /* Current other fields associated to this particle */
      struct xpart *xp = &xparts[i];
      struct gpart *gp = p->gpart;

      /* Start by copying over the particles */
      memcpy(&global_parts[k_parts], p, sizeof(struct part));
      memcpy(&global_xparts[k_parts], xp, sizeof(struct xpart));

      if (with_gravity) {
        memcpy(&global_gparts[k_gparts], gp, sizeof(struct gpart));
      }

      /* Update the IDs */
      global_parts[k_parts].id += 2;

      /* Re-link everything */
      if (with_gravity) {
        global_parts[k_parts].gpart = &global_gparts[k_gparts];
        global_gparts[k_gparts].id_or_neg_offset = -k_parts;
      }

      /* Displacement unit vector */
      const double delta_x = random_unit_interval(p->id, e->ti_current,
                                                  (enum random_number_type)0);
      const double delta_y = random_unit_interval(p->id, e->ti_current,
                                                  (enum random_number_type)1);
      const double delta_z = random_unit_interval(p->id, e->ti_current,
                                                  (enum random_number_type)2);

      /* Displace the old particle */
      p->x[0] += delta_x * displacement_factor * h;
      p->x[1] += delta_y * displacement_factor * h;
      p->x[2] += delta_z * displacement_factor * h;

      if (with_gravity) {
        gp->x[0] += delta_x * displacement_factor * h;
        gp->x[1] += delta_y * displacement_factor * h;
        gp->x[2] += delta_z * displacement_factor * h;
      }

      /* Displace the new particle */
      global_parts[k_parts].x[0] -= delta_x * displacement_factor * h;
      global_parts[k_parts].x[1] -= delta_y * displacement_factor * h;
      global_parts[k_parts].x[2] -= delta_z * displacement_factor * h;

      if (with_gravity) {
        global_gparts[k_gparts].x[0] -= delta_x * displacement_factor * h;
        global_gparts[k_gparts].x[1] -= delta_y * displacement_factor * h;
        global_gparts[k_gparts].x[2] -= delta_z * displacement_factor * h;
      }

      /* Divide the mass */
      const double new_mass = gas_mass * 0.5;
      hydro_set_mass(p, new_mass);
      hydro_set_mass(&global_parts[k_parts], new_mass);

      if (with_gravity) {
        gp->mass = new_mass;
        global_gparts[k_gparts].mass = new_mass;
      }

      /* Split the chemistry fields */
      chemistry_split_part(p, particle_split_factor);
      chemistry_split_part(&global_parts[k_parts], particle_split_factor);

      /* Split the cooling fields */
      cooling_split_part(p, xp, particle_split_factor);
      cooling_split_part(&global_parts[k_parts], &global_xparts[k_parts],
                         particle_split_factor);

      /* Split the star formation fields */
      star_formation_split_part(p, xp, particle_split_factor);
      star_formation_split_part(&global_parts[k_parts], &global_xparts[k_parts],
                                particle_split_factor);

      /* Split the tracer fields */
      tracers_split_part(p, xp, particle_split_factor);
      tracers_split_part(&global_parts[k_parts], &global_xparts[k_parts],
                         particle_split_factor);

      /* Mark the particles as not having been swallowed */
      black_holes_mark_part_as_not_swallowed(&p->black_holes_data);
      black_holes_mark_part_as_not_swallowed(
          &global_parts[k_parts].black_holes_data);
    }
  }
}

/**
 * @brief Identify all the gas particles above a given mass threshold
 * and split them into 2.
 *
 * This may reallocate the space arrays as new particles are created.
 * This is an expensive function as it has to loop multiple times over
 * the local array of particles. In case of reallocations, it may
 * also have to loop over the gravity and other arrays.
 *
 * This is done on a node-by-node basis. No MPI required here.
 *
 * @param e The #engine.
 */
void engine_split_gas_particles(struct engine *e) {

  /* Abort if we are not doing any splitting */
  if (!e->hydro_properties->particle_splitting) return;
  if (e->s->nr_parts == 0) return;

  /* Time this */
  const ticks tic = getticks();

  /* Get useful constants */
  struct space *s = e->s;
  const int with_gravity = (e->policy & engine_policy_self_gravity) ||
                           (e->policy & engine_policy_external_gravity);
  const float mass_threshold =
      e->hydro_properties->particle_splitting_mass_threshold;
  const size_t nr_parts_old = s->nr_parts;

  /* Quick check to avoid problems */
  if (particle_split_factor != 2) {
    error(
        "Invalid splitting factor. Can currently only split particles into 2!");
  }

  /* Start by counting how many particles are above the threshold
   * for splitting (this is done in parallel over the threads) */
  struct data_count data_count = {e, mass_threshold, 0};
  threadpool_map(&e->threadpool, engine_split_gas_particle_count_mapper,
                 s->parts, nr_parts_old, sizeof(struct part), 0, &data_count);
  const size_t counter = data_count.counter;

  /* Early abort? */
  if (counter == 0) {

    if (e->verbose)
      message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
              clocks_getunit());
    return;
  }

  /* Be verbose about this. This is an important event */
  message("Splitting %zd particles above the mass threshold", counter);

  /* Number of particles to create */
  const size_t count_new_gas = counter * particle_split_factor;

  /* Do we need to reallocate the gas array for the new particles? */
  if (s->nr_parts + count_new_gas > s->size_parts) {

    const size_t nr_parts_new = s->nr_parts + count_new_gas;
    s->size_parts = engine_parts_size_grow * nr_parts_new;

    if (e->verbose) message("Reallocating the part array!");

    /* Allocate a larger array and copy over */
    struct part *parts_new = NULL;
    if (swift_memalign("parts", (void **)&parts_new, part_align,
                       sizeof(struct part) * s->size_parts) != 0)
      error("Failed to allocate new part data.");

    /* Tell the compiler that the arrays are aligned */
    swift_align_information(struct part, s->parts, part_align);
    swift_align_information(struct part, parts_new, part_align);

    memcpy(parts_new, s->parts, sizeof(struct part) * s->nr_parts);
    swift_free("parts", s->parts);

    /* Allocate a larger array and copy over */
    struct xpart *xparts_new = NULL;
    if (swift_memalign("xparts", (void **)&xparts_new, xpart_align,
                       sizeof(struct xpart) * s->size_parts) != 0)
      error("Failed to allocate new part data.");

    /* Tell the compiler that the arrays are aligned */
    swift_align_information(struct xpart, s->xparts, xpart_align);
    swift_align_information(struct xpart, xparts_new, xpart_align);

    memcpy(xparts_new, s->xparts, sizeof(struct xpart) * s->nr_parts);
    swift_free("xparts", s->xparts);

    s->xparts = xparts_new;
    s->parts = parts_new;
  }

  /* Do we need to reallocate the gpart array for the new particles? */
  if (with_gravity && s->nr_gparts + count_new_gas > s->size_gparts) {

    const size_t nr_gparts_new = s->nr_gparts + count_new_gas;
    s->size_gparts = engine_parts_size_grow * nr_gparts_new;

    if (e->verbose) message("Reallocating the gpart array!");

    /* Allocate a larger array and copy over */
    struct gpart *gparts_new = NULL;
    if (swift_memalign("gparts", (void **)&gparts_new, gpart_align,
                       sizeof(struct gpart) * s->size_gparts) != 0)
      error("Failed to allocate new gpart data.");

    /* Tell the compiler that the arrays are aligned */
    swift_align_information(struct gpart, s->gparts, gpart_align);
    swift_align_information(struct gpart, gparts_new, gpart_align);

    /* Copy the particles */
    memcpy(gparts_new, s->gparts, sizeof(struct gpart) * s->nr_gparts);
    swift_free("gparts", s->gparts);

    /* We now need to correct all the pointers of the other particle arrays */
    part_relink_all_parts_to_gparts(gparts_new, s->nr_gparts, s->parts,
                                    s->sparts, s->bparts);
    s->gparts = gparts_new;
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that whatever reallocation happened we are still having correct
   * links */
  part_verify_links(s->parts, s->gparts, s->sparts, s->bparts, s->nr_parts,
                    s->nr_gparts, s->nr_sparts, s->nr_bparts, e->verbose);
#endif

  /* We now have enough memory in the part array to accomodate the new
   * particles. We can start the splitting procedure */

  size_t k_parts = s->nr_parts;
  size_t k_gparts = s->nr_gparts;

  /* Loop over the particles again to split them */
  struct data_split data_split = {e, mass_threshold, &k_parts, &k_gparts, 0};
  lock_init(&data_split.lock);
  threadpool_map(&e->threadpool, engine_split_gas_particle_split_mapper,
                 s->parts, nr_parts_old, sizeof(struct part), 0, &data_split);
  lock_destroy(&data_split.lock);

  /* Update the local counters */
  s->nr_parts = k_parts;
  s->nr_gparts = k_gparts;

#ifdef SWIFT_DEBUG_CHECKS
  if (s->nr_parts != nr_parts_old + (particle_split_factor - 1) * counter) {
    error("Incorrect number of particles created");
  }
#endif

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
