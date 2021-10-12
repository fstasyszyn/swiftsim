/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_SPHENIX_HYDRO_IACT_H
#define SWIFT_SPHENIX_HYDRO_IACT_H

/**
 * @file SPHENIX/hydro_iact.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added SPHENIX physics (Borrow 2020) (interaction routines)
 */

#include "adiabatic_index.h"
#include "hydro_parameters.h"
#include "minmax.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  float wi, wj, wi_dx, wj_dx;
  float dv[3], curlvr[3];

  const float r = sqrtf(r2);

#ifdef MHD_BASE
  double dB[3];
#ifdef MHD_EULER
  double dalpha, dbeta;
#endif
#endif


  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  /* Now we need to compute the div terms */
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->viscosity.div_v -= faci * dvdr;
  pj->viscosity.div_v -= facj * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

  /* Negative because of the change in sign of dx & dv. */
  pj->density.rot_v[0] += facj * curlvr[0];
  pj->density.rot_v[1] += facj * curlvr[1];
  pj->density.rot_v[2] += facj * curlvr[2];

#ifdef MHD_BASE
  for(int i=0;i<3;++i)
  	dB[i]= pi->bfld.B_pred[i] - pj->bfld.B_pred[i];
  const double dBdr = dB[0]*dx[0] + dB[1]*dx[1] + dB[2]*dx[2];
  pi->bfld.divB -= faci * dBdr;
  pj->bfld.divB -= facj * dBdr;
  pi->bfld.Bsm[0] += mj * wi * pj->bfld.B_pred[0];
  pi->bfld.Bsm[1] += mj * wi * pj->bfld.B_pred[1];
  pi->bfld.Bsm[2] += mj * wi * pj->bfld.B_pred[2];
  pj->bfld.Bsm[0] += mi * wj * pi->bfld.B_pred[0];
  pj->bfld.Bsm[1] += mi * wj * pi->bfld.B_pred[1];
  pj->bfld.Bsm[2] += mi * wj * pi->bfld.B_pred[2];

#ifdef MHD_DI
  pi->bfld.dBdt[0] += faci * ((pi->bfld.B_pred[0] * dv[1] - pi->bfld.B_pred[1] * dv[0]) * dx[1]
  		        +(pi->bfld.B_pred[0] * dv[2] - pi->bfld.B_pred[2] * dv[0]) * dx[2]);
  pi->bfld.dBdt[1] += faci * ((pi->bfld.B_pred[1] * dv[2] - pi->bfld.B_pred[2] * dv[1]) * dx[2]
  		        +(pi->bfld.B_pred[1] * dv[0] - pi->bfld.B_pred[0] * dv[1]) * dx[0]);
  pi->bfld.dBdt[2] += faci * ((pi->bfld.B_pred[2] * dv[0] - pi->bfld.B_pred[0] * dv[2]) * dx[0]
  		        +(pi->bfld.B_pred[2] * dv[1] - pi->bfld.B_pred[1] * dv[2]) * dx[1]);
  pj->bfld.dBdt[0] += facj * ((pj->bfld.B_pred[0] * dv[1] - pj->bfld.B_pred[1] * dv[0]) * dx[1]
  		        +(pj->bfld.B_pred[0] * dv[2] - pj->bfld.B_pred[2] * dv[0]) * dx[2]);
  pj->bfld.dBdt[1] += facj * ((pj->bfld.B_pred[1] * dv[2] - pj->bfld.B_pred[2] * dv[1]) * dx[2]
  		        +(pj->bfld.B_pred[1] * dv[0] - pj->bfld.B_pred[0] * dv[1]) * dx[0]);
  pj->bfld.dBdt[2] += facj * ((pj->bfld.B_pred[2] * dv[0] - pj->bfld.B_pred[0] * dv[2]) * dx[0]
  		        +(pj->bfld.B_pred[2] * dv[1] - pj->bfld.B_pred[1] * dv[2]) * dx[1]);
#endif
#ifdef MHD_EULER
  dalpha = pi->bfld.ep[0] - pj->bfld.ep[0];
  dbeta  = pi->bfld.ep[1] - pj->bfld.ep[1];

#if MHD_EULER_TEST == 1
// BrioWu
  const float LBOX=1.0;
  dalpha = (( dalpha > LBOX/2.0 ) ? dalpha-LBOX : ( ( dalpha < -LBOX/2.0 ) ? dalpha+LBOX: dalpha));
  dbeta  = (( dbeta > 0.75*LBOX/2.0 ) ? dbeta-0.75*LBOX : ( ( dbeta < -0.75*LBOX/2.0 ) ? dbeta+0.75*LBOX: dbeta));
#endif
#if MHD_EULER_TEST == 2
// VORTEX
  const float LBOX=1.0;
  dbeta  = (( dbeta  > LBOX/2.0 ) ? dbeta-LBOX  : ( ( dbeta  < -LBOX/2.0 ) ? dbeta+LBOX : dbeta));
#endif
  for(int i=0;i<3;++i)  
  pi->bfld.Grad_ep[0][i] += faci * dalpha*dx[i];

  for(int i=0;i<3;++i)  
  pi->bfld.Grad_ep[1][i] += faci * dbeta*dx[i];
  
  for(int i=0;i<3;++i)  
  pj->bfld.Grad_ep[0][i] += facj * dalpha*dx[i];

  for(int i=0;i<3;++i)  
  pj->bfld.Grad_ep[1][i] += facj * dbeta*dx[i];
#endif  /* MHD_EULER */
#endif  /* MHD_BASE */

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_density += wi;
  pj->n_density += wj;
  pi->N_density++;
  pj->N_density++;
#endif
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, const struct part* restrict pj, const float a,
    const float H) {

  float wi, wi_dx;
  float dv[3], curlvr[3];
#ifdef MHD_BASE
  double dB[3];
#ifdef MHD_EULER
  double dalpha, dbeta;
#endif
#endif

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->viscosity.div_v -= faci * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

#ifdef MHD_BASE
  for(int i=0;i<3;++i)
  	dB[i]= pi->bfld.B_pred[i] - pj->bfld.B_pred[i];
  const double dBdr = dB[0]*dx[0] + dB[1]*dx[1] + dB[2]*dx[2];
  pi->bfld.divB -= faci * dBdr;

  pi->bfld.Bsm[0] += mj * wi * pj->bfld.B_pred[0];
  pi->bfld.Bsm[1] += mj * wi * pj->bfld.B_pred[1];
  pi->bfld.Bsm[2] += mj * wi * pj->bfld.B_pred[2];
#ifdef MHD_DI
  pi->bfld.dBdt[0] += faci * ((pi->bfld.B_pred[0] * dv[1] - pi->bfld.B_pred[1] * dv[0]) * dx[1]
  		        +(pi->bfld.B_pred[0] * dv[2] - pi->bfld.B_pred[2] * dv[0]) * dx[2]);
  pi->bfld.dBdt[1] += faci * ((pi->bfld.B_pred[1] * dv[2] - pi->bfld.B_pred[2] * dv[1]) * dx[2]
  		        +(pi->bfld.B_pred[1] * dv[0] - pi->bfld.B_pred[0] * dv[1]) * dx[0]);
  pi->bfld.dBdt[2] += faci * ((pi->bfld.B_pred[2] * dv[0] - pi->bfld.B_pred[0] * dv[2]) * dx[0]
  		        +(pi->bfld.B_pred[2] * dv[1] - pi->bfld.B_pred[1] * dv[2]) * dx[1]);
#endif

#ifdef MHD_EULER
  dalpha = pi->bfld.ep[0] - pj->bfld.ep[0];
  dbeta  = pi->bfld.ep[1] - pj->bfld.ep[1];
  
#if MHD_EULER_TEST == 1
// BrioWu
  const float LBOX=1.0;
  dalpha = (( dalpha > LBOX/2.0 ) ? dalpha-LBOX : ( ( dalpha < -LBOX/2.0 ) ? dalpha+LBOX: dalpha));
  dbeta  = (( dbeta > 0.75*LBOX/2.0 ) ? dbeta-0.75*LBOX : ( ( dbeta < -0.75*LBOX/2.0 ) ? dbeta+0.75*LBOX: dbeta));
#endif
#if MHD_EULER_TEST == 2
// VORTEX
  const float LBOX=1.0;
  dbeta  = (( dbeta  > LBOX/2.0 ) ? dbeta-LBOX  : ( ( dbeta  < -LBOX/2.0 ) ? dbeta+LBOX : dbeta));
#endif
  for(int i=0;i<3;++i)  
     pi->bfld.Grad_ep[0][i] += faci * dalpha * dx[i];
  
  for(int i=0;i<3;++i)  
     pi->bfld.Grad_ep[1][i] += faci * dbeta  * dx[i];
  
#endif  /* MHD_EULER */
#endif  /* MHD_BASE */

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_density += wi;
  pi->N_density++;
#endif
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Cosmology terms for the signal velocity */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */

  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float new_v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig);
  pj->viscosity.v_sig = max(pj->viscosity.v_sig, new_v_sig);

  /* Calculate Del^2 u for the thermal diffusion coefficient. */
  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx, wj, wj_dx;

  const float ui = r / hi;
  const float uj = r / hj;

  kernel_deval(ui, &wi, &wi_dx);
  kernel_deval(uj, &wj, &wj_dx);

  const float delta_u_factor = (pi->u - pj->u) * r_inv;
  pi->diffusion.laplace_u += pj->mass * delta_u_factor * wi_dx / pj->rho;
  pj->diffusion.laplace_u -= pi->mass * delta_u_factor * wj_dx / pi->rho;

  /* Set the maximal alpha from the previous step over the neighbours
   * (this is used to limit the diffusion in hydro_prepare_force) */
  const float alpha_i = pi->viscosity.alpha;
  const float alpha_j = pj->viscosity.alpha;
  pi->force.alpha_visc_max_ngb = max(pi->force.alpha_visc_max_ngb, alpha_j);
  pj->force.alpha_visc_max_ngb = max(pj->force.alpha_visc_max_ngb, alpha_i);

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_gradient += wi;
  pj->n_gradient += wj;
  pi->N_gradient++;
  pj->N_gradient++;
#endif
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Cosmology terms for the signal velocity */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */

  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float new_v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig);

  /* Calculate Del^2 u for the thermal diffusion coefficient. */
  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx;

  const float ui = r / hi;

  kernel_deval(ui, &wi, &wi_dx);

  const float delta_u_factor = (pi->u - pj->u) * r_inv;
  pi->diffusion.laplace_u += pj->mass * delta_u_factor * wi_dx / pj->rho;

  /* Set the maximal alpha from the previous step over the neighbours
   * (this is used to limit the diffusion in hydro_prepare_force) */
  const float alpha_j = pj->viscosity.alpha;
  pi->force.alpha_visc_max_ngb = max(pi->force.alpha_visc_max_ngb, alpha_j);

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_gradient += wi;
  pi->N_gradient++;
#endif
}

/**
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

#ifdef MHD_BASE
  const float MU0_1 = 1.0/(4.0*M_PI);
#endif
  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;
  const float mi = pi->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Includes the hubble flow term; not used for du/dt */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
#ifndef MHD_BASE
  const float v_sig = pi->force.soundspeed + pj->force.soundspeed -
                      const_viscosity_beta * mu_ij;
#else
  // CHECK MU0
  const float ci   = pi->force.soundspeed;
  const float cj   = pj->force.soundspeed;
  const float b2_i = (pi->bfld.B_pred[0]*pi->bfld.B_pred[0] + pi->bfld.B_pred[1]*pi->bfld.B_pred[1] + pi->bfld.B_pred[2]*pi->bfld.B_pred[2] );
  const float b2_j = (pj->bfld.B_pred[0]*pj->bfld.B_pred[0] + pj->bfld.B_pred[1]*pj->bfld.B_pred[1] + pj->bfld.B_pred[2]*pj->bfld.B_pred[2] ); 
  //float vcsa2_i = ci * ci + min(MU0_1 * b2_i/rhoi,10.0*ci*ci); 
  //float vcsa2_j = cj * cj + min(MU0_1 * b2_j/rhoj,10.0*cj*cj); 
  const float vcsa2_i = ci * ci + MU0_1 * b2_i/rhoi; 
  const float vcsa2_j = cj * cj + MU0_1 * b2_j/rhoj; 
  float Bpro2_i = (pi->bfld.B_pred[0]*dx[0]+ pi->bfld.B_pred[1]*dx[1]+ pi->bfld.B_pred[2]*dx[2]) * r_inv;
        Bpro2_i *= Bpro2_i;
  float mag_speed_i = sqrtf(0.5 * (vcsa2_i + 
  		      sqrtf(max(  (vcsa2_i * vcsa2_i - 4. * ci * ci * Bpro2_i * MU0_1 / rhoi),0.0))));
  float Bpro2_j = (pj->bfld.B_pred[0]*dx[0]+ pj->bfld.B_pred[1]*dx[1]+ pj->bfld.B_pred[2]*dx[2]) * r_inv;
        Bpro2_j *= Bpro2_j;
  float mag_speed_j = sqrtf(0.5 * (vcsa2_j + 
  		      sqrtf(max(  (vcsa2_j * vcsa2_j - 4. * cj * cj * Bpro2_j * MU0_1 / rhoj),0.0))));

  const float v_sig = 3.0*(mag_speed_i + mag_speed_j - const_viscosity_beta/2.0 * mu_ij);
#endif

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Construct the full viscosity term */
  const float rho_ij = rhoi + rhoj;
  const float alpha = pi->viscosity.alpha + pj->viscosity.alpha;
  const float visc =
      -0.25f * alpha * v_sig * mu_ij * (balsara_i + balsara_j) / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  pj->a_hydro[0] += mi * acc * dx[0];
  pj->a_hydro[1] += mi * acc * dx[1];
  pj->a_hydro[2] += mi * acc * dx[2];

  /* Eventually got the MHD accel */ 
#ifdef MHD_BASE
  const float mag_faci = MU0_1 * f_ij * wi_dr * r_inv /(rhoi*rhoi);
  const float mag_facj = MU0_1 * f_ji * wj_dr * r_inv /(rhoj*rhoj);
//  float mm_i[3][3],mm_j[3][3];
//  
//  for(int i=0;i<3;i++)
//    for(int j=0;j<3;j++)
//  	{
//	mm_i[i][j] = pi->Bfld[i]*pi->Bfld[j];
 // 	
//	mm_j[i][j] = pj->Bfld[i]*pj->Bfld[j];
//	}
//  for(int j=0;j<3;j++)
//  	{mm_i[j][j] -= 0.5*pi->Bfld[j]*pi->Bfld[j];
//  	 mm_j[j][j] -= 0.5*pj->Bfld[j]*pj->Bfld[j];}
///////////////////////////////
/*  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
    {  
       //pi->a_hydro[i] -= mj * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj) * dx[j];
       //pj->a_hydro[i] += mi * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj) * dx[j];
          pi->a_hydro[i] -= mj * ((pi->Bfld[i]*pi->Bfld[j])*mag_faci
	  + ( pj->Bfld[i]*pj->Bfld[j])*mag_facj)*dx[j];
	 if(i==j)
          pi->a_hydro[i] += mj * 0.5 * ((pi->Bfld[i]*pi->Bfld[j])*mag_faci
	  + ( pj->Bfld[i]*pj->Bfld[j])*mag_facj)*dx[j];

         pj->a_hydro[i] += mi * ((pi->Bfld[i]*pi->Bfld[j])*mag_faci
	  + ( pj->Bfld[i]*pj->Bfld[j])*mag_facj)*dx[j];
	 if(i==j)
         pj->a_hydro[i] -= mi * 0.5 * ((pi->Bfld[i]*pi->Bfld[j])*mag_faci
	  + ( pj->Bfld[i]*pj->Bfld[j])*mag_facj)*dx[j];
     }*/
  pi->a_hydro[0] -= mj * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[0]) * mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pi->a_hydro[0] -= mj * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[1]) * mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pi->a_hydro[0] -= mj * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[2]) * mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
  pi->a_hydro[1] -= mj * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[0]) * mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pi->a_hydro[1] -= mj * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[1]) * mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pi->a_hydro[1] -= mj * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[2]) * mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
  pi->a_hydro[2] -= mj * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[0]) * mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pi->a_hydro[2] -= mj * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[1]) * mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pi->a_hydro[2] -= mj * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[2]) * mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
  pi->a_hydro[0] += mj * 0.5 * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[0])*mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pi->a_hydro[1] += mj * 0.5 * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[1])*mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pi->a_hydro[2] += mj * 0.5 * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[2])*mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
  
  pj->a_hydro[0] += mi * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[0]) * mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pj->a_hydro[0] += mi * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[1]) * mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pj->a_hydro[0] += mi * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[2]) * mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
  pj->a_hydro[1] += mi * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[0]) * mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pj->a_hydro[1] += mi * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[1]) * mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pj->a_hydro[1] += mi * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[2]) * mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
  pj->a_hydro[2] += mi * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[0]) * mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pj->a_hydro[2] += mi * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[1]) * mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pj->a_hydro[2] += mi * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[2]) * mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
  pj->a_hydro[0] -= mi * 0.5 * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[0])*mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pj->a_hydro[1] -= mi * 0.5 * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[1])*mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pj->a_hydro[2] -= mi * 0.5 * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[2])*mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
//Take out the divergence term  
/*  for(int i=0;i<3;++i)
    for(int j=0;j<3;j++)
    {
     pi->a_hydro[i] += mj * pi->Bfld[i] * (pi->Bfld[j]*mag_faci+pj->Bfld[j]*mag_facj)*dx[j];
     pj->a_hydro[i] -= mi * pj->Bfld[i] * (pi->Bfld[j]*mag_faci+pj->Bfld[j]*mag_facj)*dx[j];
    }*/
  pi->a_hydro[0] += mj * pi->bfld.B_pred[0] * (pi->bfld.B_pred[0]*mag_faci+pj->bfld.B_pred[0]*mag_facj)*dx[0];
  pi->a_hydro[0] += mj * pi->bfld.B_pred[0] * (pi->bfld.B_pred[1]*mag_faci+pj->bfld.B_pred[1]*mag_facj)*dx[1];
  pi->a_hydro[0] += mj * pi->bfld.B_pred[0] * (pi->bfld.B_pred[2]*mag_faci+pj->bfld.B_pred[2]*mag_facj)*dx[2];
  pi->a_hydro[1] += mj * pi->bfld.B_pred[1] * (pi->bfld.B_pred[0]*mag_faci+pj->bfld.B_pred[0]*mag_facj)*dx[0];
  pi->a_hydro[1] += mj * pi->bfld.B_pred[1] * (pi->bfld.B_pred[1]*mag_faci+pj->bfld.B_pred[1]*mag_facj)*dx[1];
  pi->a_hydro[1] += mj * pi->bfld.B_pred[1] * (pi->bfld.B_pred[2]*mag_faci+pj->bfld.B_pred[2]*mag_facj)*dx[2];
  pi->a_hydro[2] += mj * pi->bfld.B_pred[2] * (pi->bfld.B_pred[0]*mag_faci+pj->bfld.B_pred[0]*mag_facj)*dx[0];
  pi->a_hydro[2] += mj * pi->bfld.B_pred[2] * (pi->bfld.B_pred[1]*mag_faci+pj->bfld.B_pred[1]*mag_facj)*dx[1];
  pi->a_hydro[2] += mj * pi->bfld.B_pred[2] * (pi->bfld.B_pred[2]*mag_faci+pj->bfld.B_pred[2]*mag_facj)*dx[2];
 
  pj->a_hydro[0] -= mi * pj->bfld.B_pred[0] * (pi->bfld.B_pred[0]*mag_faci+pj->bfld.B_pred[0]*mag_facj)*dx[0];
  pj->a_hydro[0] -= mi * pj->bfld.B_pred[0] * (pi->bfld.B_pred[1]*mag_faci+pj->bfld.B_pred[1]*mag_facj)*dx[1];
  pj->a_hydro[0] -= mi * pj->bfld.B_pred[0] * (pi->bfld.B_pred[2]*mag_faci+pj->bfld.B_pred[2]*mag_facj)*dx[2];
  pj->a_hydro[1] -= mi * pj->bfld.B_pred[1] * (pi->bfld.B_pred[0]*mag_faci+pj->bfld.B_pred[0]*mag_facj)*dx[0];
  pj->a_hydro[1] -= mi * pj->bfld.B_pred[1] * (pi->bfld.B_pred[1]*mag_faci+pj->bfld.B_pred[1]*mag_facj)*dx[1];
  pj->a_hydro[1] -= mi * pj->bfld.B_pred[1] * (pi->bfld.B_pred[2]*mag_faci+pj->bfld.B_pred[2]*mag_facj)*dx[2];
  pj->a_hydro[2] -= mi * pj->bfld.B_pred[2] * (pi->bfld.B_pred[0]*mag_faci+pj->bfld.B_pred[0]*mag_facj)*dx[0];
  pj->a_hydro[2] -= mi * pj->bfld.B_pred[2] * (pi->bfld.B_pred[1]*mag_faci+pj->bfld.B_pred[1]*mag_facj)*dx[1];
  pj->a_hydro[2] -= mi * pj->bfld.B_pred[2] * (pi->bfld.B_pred[2]*mag_faci+pj->bfld.B_pred[2]*mag_facj)*dx[2];
#endif

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;
  const float sph_du_term_j = P_over_rho2_j * dvdr * r_inv * wj_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Diffusion term */
  /* Combine the alpha_diff into a pressure-based switch -- this allows the
   * alpha from the highest pressure particle to dominate, so that the
   * diffusion limited particles always take precedence - another trick to
   * allow the scheme to work with thermal feedback. */
  const float alpha_diff =
      (pressurei * pi->diffusion.alpha + pressurej * pj->diffusion.alpha) /
      (pressurei + pressurej);
  const float v_diff = alpha_diff * 0.5f *
                       (sqrtf(2.f * fabsf(pressurei - pressurej) / rho_ij) +
                        fabsf(fac_mu * r_inv * dvdr_Hubble));
  /* wi_dx + wj_dx / 2 is F_ij */
  const float diff_du_term =
      v_diff * (pi->u - pj->u) * (f_ij * wi_dr / rhoi + f_ji * wj_dr / rhoj);

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;
  const float du_dt_j = sph_du_term_j + visc_du_term - diff_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

  /* Update if we need to; this should be guaranteed by the gradient loop but
   * due to some possible synchronisation problems this is here as a _quick
   * fix_. Added: 14th August 2019. To be removed by 1st Jan 2020. (JB) */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, v_sig);
  pj->viscosity.v_sig = max(pj->viscosity.v_sig, v_sig);

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pj->n_force += wi + wj;
  pi->N_force++;
  pj->N_force++;
#endif
}

/**
 * @brief Force interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, const struct part* restrict pj, const float a,
    const float H) {

#ifdef MHD_BASE
  const float MU0_1 = 1.0/(4.0*M_PI);
#endif
  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Includes the hubble flow term; not used for du/dt */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
#ifndef MHD_BASE
  const float v_sig = pi->force.soundspeed + pj->force.soundspeed -
                      const_viscosity_beta * mu_ij;
#else
  const float ci   = pi->force.soundspeed;
  const float cj   = pj->force.soundspeed;
  const float b2_i = (pi->bfld.B_pred[0]*pi->bfld.B_pred[0] + pi->bfld.B_pred[1]*pi->bfld.B_pred[1] + pi->bfld.B_pred[2]*pi->bfld.B_pred[2] );
  const float b2_j = (pj->bfld.B_pred[0]*pj->bfld.B_pred[0] + pj->bfld.B_pred[1]*pj->bfld.B_pred[1] + pj->bfld.B_pred[2]*pj->bfld.B_pred[2] ); 
  //float vcsa2_i = ci * ci + min(MU0_1 * b2_i/rhoi,10.0*ci*ci); 
  //float vcsa2_j = cj * cj + min(MU0_1 * b2_j/rhoj,10.0*cj*cj); 
  const float vcsa2_i = ci * ci + MU0_1 * b2_i/rhoi; 
  const float vcsa2_j = cj * cj + MU0_1 * b2_j/rhoj; 
  float Bpro2_i = (pi->bfld.B_pred[0]*dx[0]+ pi->bfld.B_pred[1]*dx[1]+ pi->bfld.B_pred[2]*dx[2]) * r_inv;
        Bpro2_i *= Bpro2_i;
  float mag_speed_i = sqrtf(0.5 * (vcsa2_i + 
  		      sqrtf(max(  (vcsa2_i * vcsa2_i - 4. * ci * ci * Bpro2_i * MU0_1 / rhoi),0.0))));
  float Bpro2_j = (pj->bfld.B_pred[0]*dx[0]+ pj->bfld.B_pred[1]*dx[1]+ pj->bfld.B_pred[2]*dx[2]) * r_inv;
        Bpro2_j *= Bpro2_j;
  float mag_speed_j = sqrtf(0.5 * (vcsa2_j + 
  		      sqrtf(max(  (vcsa2_j * vcsa2_j - 4. * cj * cj * Bpro2_j * MU0_1 / rhoj),0.0))));

  const float v_sig = 3.0*(mag_speed_i + mag_speed_j - const_viscosity_beta/2.0 * mu_ij);
#endif


  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Construct the full viscosity term */
  const float rho_ij = rhoi + rhoj;
  const float alpha = pi->viscosity.alpha + pj->viscosity.alpha;
  const float visc =
      -0.25f * alpha * v_sig * mu_ij * (balsara_i + balsara_j) / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Eventually got the MHD accel */ 
#ifdef MHD_BASE
  const float mag_faci = MU0_1 * f_ij * wi_dr * r_inv /(rhoi*rhoi);
  const float mag_facj = MU0_1 * f_ji * wj_dr * r_inv /(rhoj*rhoj);
  //float mm_i[3][3],mm_j[3][3];
  
 // for(int i=0;i<3;i++)
 //   for(int j=0;j<3;j++)
//  	{
//	mm_i[i][j] = pi->Bfld[i]*pi->Bfld[j];
//	mm_j[i][j] = pj->Bfld[i]*pj->Bfld[j];
//	 }
//  for(int j=0;j<3;j++)
//  	{mm_i[j][j]-=0.5*pi->Bfld[j]*pi->Bfld[j];
//  	 mm_j[j][j]-=0.5*pj->Bfld[j]*pj->Bfld[j];}
//////////////////////////////////
/*  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
//          pi->a_hydro[i] -= mj * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj)*dx[j];
          pi->a_hydro[i] -= mj * ((pi->Bfld[i]*pi->Bfld[j])*mag_faci
	  + ( pj->Bfld[i]*pj->Bfld[j])*mag_facj)*dx[j];
	 if(i==j)
          pi->a_hydro[i] += mj * 0.5 * ((pi->Bfld[i]*pi->Bfld[j])*mag_faci
	  + ( pj->Bfld[i]*pj->Bfld[j])*mag_facj)*dx[j];
	 }
*/
  pi->a_hydro[0] -= mj * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[0]) * mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pi->a_hydro[0] -= mj * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[1]) * mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pi->a_hydro[0] -= mj * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[2]) * mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
  pi->a_hydro[1] -= mj * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[0]) * mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pi->a_hydro[1] -= mj * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[1]) * mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pi->a_hydro[1] -= mj * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[2]) * mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
  pi->a_hydro[2] -= mj * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[0]) * mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pi->a_hydro[2] -= mj * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[1]) * mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pi->a_hydro[2] -= mj * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[2]) * mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
  pi->a_hydro[0] += mj * 0.5 * ((pi->bfld.B_pred[0]*pi->bfld.B_pred[0])*mag_faci + ( pj->bfld.B_pred[0]*pj->bfld.B_pred[0])*mag_facj)*dx[0];
  pi->a_hydro[1] += mj * 0.5 * ((pi->bfld.B_pred[1]*pi->bfld.B_pred[1])*mag_faci + ( pj->bfld.B_pred[1]*pj->bfld.B_pred[1])*mag_facj)*dx[1];
  pi->a_hydro[2] += mj * 0.5 * ((pi->bfld.B_pred[2]*pi->bfld.B_pred[2])*mag_faci + ( pj->bfld.B_pred[2]*pj->bfld.B_pred[2])*mag_facj)*dx[2];
//Take out the divergence term  
//  for(int i=0;i<3;++i)
//    for(int j=0;j<3;++j)
//	  pi->a_hydro[i] += mj * pi->Bfld[i] * (pi->Bfld[j]*mag_faci+pj->Bfld[j]*mag_facj)*dx[j];
  pi->a_hydro[0] += mj * pi->bfld.B_pred[0] * (pi->bfld.B_pred[0]*mag_faci+pj->bfld.B_pred[0]*mag_facj)*dx[0];
  pi->a_hydro[0] += mj * pi->bfld.B_pred[0] * (pi->bfld.B_pred[1]*mag_faci+pj->bfld.B_pred[1]*mag_facj)*dx[1];
  pi->a_hydro[0] += mj * pi->bfld.B_pred[0] * (pi->bfld.B_pred[2]*mag_faci+pj->bfld.B_pred[2]*mag_facj)*dx[2];
  pi->a_hydro[1] += mj * pi->bfld.B_pred[1] * (pi->bfld.B_pred[0]*mag_faci+pj->bfld.B_pred[0]*mag_facj)*dx[0];
  pi->a_hydro[1] += mj * pi->bfld.B_pred[1] * (pi->bfld.B_pred[1]*mag_faci+pj->bfld.B_pred[1]*mag_facj)*dx[1];
  pi->a_hydro[1] += mj * pi->bfld.B_pred[1] * (pi->bfld.B_pred[2]*mag_faci+pj->bfld.B_pred[2]*mag_facj)*dx[2];
  pi->a_hydro[2] += mj * pi->bfld.B_pred[2] * (pi->bfld.B_pred[0]*mag_faci+pj->bfld.B_pred[0]*mag_facj)*dx[0];
  pi->a_hydro[2] += mj * pi->bfld.B_pred[2] * (pi->bfld.B_pred[1]*mag_faci+pj->bfld.B_pred[1]*mag_facj)*dx[1];
  pi->a_hydro[2] += mj * pi->bfld.B_pred[2] * (pi->bfld.B_pred[2]*mag_faci+pj->bfld.B_pred[2]*mag_facj)*dx[2];
#endif

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Diffusion term */
  /* Combine the alpha_diff into a pressure-based switch -- this allows the
   * alpha from the highest pressure particle to dominate, so that the
   * diffusion limited particles always take precedence - another trick to
   * allow the scheme to work with thermal feedback. */
  const float alpha_diff =
      (pressurei * pi->diffusion.alpha + pressurej * pj->diffusion.alpha) /
      (pressurei + pressurej);
  const float v_diff = alpha_diff * 0.5f *
                       (sqrtf(2.f * fabsf(pressurei - pressurej) / rho_ij) +
                        fabsf(fac_mu * r_inv * dvdr_Hubble));
  /* wi_dx + wj_dx / 2 is F_ij */
  const float diff_du_term =
      v_diff * (pi->u - pj->u) * (f_ij * wi_dr / rhoi + f_ji * wj_dr / rhoj);

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

  /* Update if we need to; this should be guaranteed by the gradient loop but
   * due to some possible synchronisation problems this is here as a _quick
   * fix_. Added: 14th August 2019. To be removed by 1st Jan 2020. (JB) */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, v_sig);

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pi->N_force++;
#endif
}
#endif /* SWIFT_SPHENIX_HYDRO_IACT_H */

