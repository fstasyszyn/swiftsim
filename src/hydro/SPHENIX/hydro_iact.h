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
  double dB[3];
  for(int i=0;i<3;++i)
  	dB[i]= pi->BPred[i] - pj->BPred[i];
  const double dBdr = dB[0]*dx[0] + dB[1]*dx[1] + dB[2]*dx[2];
  pi->divB -= faci * dBdr;
  pj->divB -= facj * dBdr;
#ifdef MHD_EULER
  double dalpha, dbeta;
  dalpha = pi->ep[0] - pj->ep[0];
  dbeta  = pi->ep[1] - pj->ep[1];

#if MHD_EULER_TEST == 1
// BrioWu
  dalpha = nearest(dalpha, engine_extra_dims[2]);
  dbeta  = nearest(dbeta, 0.75*engine_extra_dims[1]);
#endif
#if MHD_EULER_TEST == 2
// VORTEX
  dbeta = nearest(dbeta, engine_extra_dims[2]);
#endif
  for(int i=0;i<3;++i) 
  { 
  pi->Grad_ep[0][i] += faci * dalpha*dx[i];

  pi->Grad_ep[1][i] += faci * dbeta*dx[i];
  
  pj->Grad_ep[0][i] += facj * dalpha*dx[i];

  pj->Grad_ep[1][i] += facj * dbeta*dx[i];
  }
#endif  /* MHD_EULER */
#ifdef MHD_VECPOT
  double dA[3];
  for(int i=0;i<3;++i)
  	//dA[i]= pi->APred[i] - pj->APred[i];
  	dA[i]= pi->APot[i] - pj->APot[i];
  const double dAdr = dA[0]*dx[0] + dA[1]*dx[1] + dA[2]*dx[2];
  pi->divA -= faci * dAdr;
  pj->divA -= facj * dAdr;
  for(int i=0;i<3;++i)
  { 
     pi->BPred[i] += faci * (dA[(i+1)%3]*dx[(i+2)%3] - dA[(i+2)%3]*dx[(i+1)%3] ) ;
     pj->BPred[i] += facj * (dA[(i+1)%3]*dx[(i+2)%3] - dA[(i+2)%3]*dx[(i+1)%3] ) ;
  }
#endif
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
  double dB[3];
  for(int i=0;i<3;++i)
  	dB[i]= pi->BPred[i] - pj->BPred[i];
  const double dBdr = dB[0]*dx[0] + dB[1]*dx[1] + dB[2]*dx[2];
  pi->divB -= faci * dBdr;
#ifdef MHD_EULER
  double dalpha, dbeta;
  dalpha = pi->ep[0] - pj->ep[0];
  dbeta  = pi->ep[1] - pj->ep[1];
  
#if MHD_EULER_TEST == 1
// BrioWu
  dalpha = nearest(dalpha, engine_extra_dims[2]);
  dbeta  = nearest(dbeta, 0.75*engine_extra_dims[1]);
#endif
#if MHD_EULER_TEST == 2
// VORTEX
  dbeta = nearest(dbeta, engine_extra_dims[2]);
#endif
  for(int i=0;i<3;++i){  
  pi->Grad_ep[0][i] += faci * dalpha * dx[i];
  
  pi->Grad_ep[1][i] += faci * dbeta  * dx[i];
  }
#endif  /* MHD_EULER */
#ifdef MHD_VECPOT
  double dA[3];
  for(int i=0;i<3;++i)
  	dA[i]= pi->APot[i] - pj->APot[i];
  	//dA[i]= pi->APred[i] - pj->APred[i];
  const double dAdr = dA[0]*dx[0] + dA[1]*dx[1] + dA[2]*dx[2];
  pi->divA -= faci * dAdr;
  for(int i=0;i<3;++i)
     pi->BPred[i] += faci * (dA[(i+1)%3]*dx[(i+2)%3] - dA[(i+2)%3]*dx[(i+1)%3] ) ;
#endif 
#endif  /* MHD */

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
#ifdef MHD_VECPOT
  for(int i=0;i<3;i++){
     pi->Bfld[i] +=  pj->mass * wi * pj->BPred[i];
     pj->Bfld[i] +=  pi->mass * wj * pi->BPred[i];
     }
  pi->Q0 += pj->mass * wi;
  pj->Q0 += pi->mass * wj;
#endif
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
#ifdef MHD_VECPOT  
  for(int i=0;i<3;i++)
     pi->Bfld[i] += pj->mass * wi * pj->BPred[i];
  pi->Q0 += pj->mass * wi;
#endif
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
#ifdef MHD_ORESTIS
  float Bi[3];
  float Bj[3];
  Bi[0] = pi->BPred[0] * rhoi;
  Bi[1] = pi->BPred[1] * rhoi;
  Bi[2] = pi->BPred[2] * rhoi;
  Bj[0] = pj->BPred[0] * rhoj;
  Bj[1] = pj->BPred[1] * rhoj;
  Bj[2] = pj->BPred[2] * rhoj;
#endif
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
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
#ifndef MHD_BASE
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;
#else
  const float b2_i = (pi->BPred[0]*pi->BPred[0] + pi->BPred[1]*pi->BPred[1] + pi->BPred[2]*pi->BPred[2] );
  const float b2_j = (pj->BPred[0]*pj->BPred[0] + pj->BPred[1]*pj->BPred[1] + pj->BPred[2]*pj->BPred[2] ); 
  const float vcsa2_i = ci * ci + MU0_1 * b2_i/rhoi; 
  const float vcsa2_j = cj * cj + MU0_1 * b2_j/rhoj; 
  float Bpro2_i = (pi->BPred[0]*dx[0]+ pi->BPred[1]*dx[1]+ pi->BPred[2]*dx[2]) * r_inv;
        Bpro2_i *= Bpro2_i;
  float mag_speed_i = sqrtf(0.5 * (vcsa2_i + 
  		      sqrtf(max(  (vcsa2_i * vcsa2_i - 4. * ci * ci * Bpro2_i * MU0_1 / rhoi),0.0))));
  float Bpro2_j = (pj->BPred[0]*dx[0]+ pj->BPred[1]*dx[1]+ pj->BPred[2]*dx[2]) * r_inv;
        Bpro2_j *= Bpro2_j;
  float mag_speed_j = sqrtf(0.5 * (vcsa2_j + 
  		      sqrtf(max(  (vcsa2_j * vcsa2_j - 4. * cj * cj * Bpro2_j * MU0_1 / rhoj),0.0))));

  const float v_sig = (mag_speed_i + mag_speed_j - const_viscosity_beta * mu_ij);
  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, v_sig);
  pj->viscosity.v_sig = max(pj->viscosity.v_sig, v_sig);
  /* Grab balsara switches */
  //const float balsara_i = 1.f;
  //const float balsara_j = 1.f;
#endif
  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;
  /* Construct the full viscosity term */
  const float rho_ij = rhoi + rhoj;
#ifdef MHD_BASE
  const float alpha = 2.f * (pi->viscosity.alpha + pj->viscosity.alpha);
#else
  const float alpha = pi->viscosity.alpha + pj->viscosity.alpha;
#endif
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
  float mm_i[3][3],mm_j[3][3];
//  
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
  	{
	mm_i[i][j] = pi->BPred[i]*pi->BPred[j];
	mm_j[i][j] = pj->BPred[i]*pj->BPred[j];
	}
  for(int j=0;j<3;j++)
  	{
	 mm_i[j][j] -= 0.5 * (pi->BPred[0]*pi->BPred[0]+pi->BPred[1]*pi->BPred[1]+pi->BPred[2]*pi->BPred[2]);
  	 mm_j[j][j] -= 0.5 * (pj->BPred[0]*pj->BPred[0]+pj->BPred[1]*pj->BPred[1]+pj->BPred[2]*pj->BPred[2]);
	 }
///////////////////////////////
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
    {  
       pi->a_hydro[i] += mj * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj) * dx[j];
       pj->a_hydro[i] -= mi * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj) * dx[j];
       pi->a_hydro[i] -= pi->Q0 * mj * pi->BPred[i] * (pi->BPred[j]*mag_faci+pj->BPred[j]*mag_facj)*dx[j];
       pj->a_hydro[i] += pj->Q0 * mi * pj->BPred[i] * (pi->BPred[j]*mag_faci+pj->BPred[j]*mag_facj)*dx[j];
     }
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

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pj->n_force += wi + wj;
  pi->N_force++;
  pj->N_force++;
#endif
#ifdef MHD_DI
  const float mag_Indi = wi_dr * r_inv / rhoi;
  const float mag_Indj = wj_dr * r_inv / rhoj;
  float dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  for(int i=0;i<3;i++) 
  {
  pi->dBdt[i] += mj * mag_Indi * ((pi->BPred[i] * dv[(i+1)%3] - pi->BPred[(i+1)%3] * dv[i]) * dx[(i+1)%3]
			        + (pi->BPred[i] * dv[(i+2)%3] - pi->BPred[(i+2)%3] * dv[i]) * dx[(i+2)%3]);
  pj->dBdt[i] += mi * mag_Indj * ((pj->BPred[i] * dv[(i+1)%3] - pj->BPred[(i+1)%3] * dv[i]) * dx[(i+1)%3]
			        + (pj->BPred[i] * dv[(i+2)%3] - pj->BPred[(i+2)%3] * dv[i]) * dx[(i+2)%3]);
  pi->dBdt[i] += pi->Q1 * mj * mag_Indi * (pi->phi - pj->phi) * dx[i];
  pj->dBdt[i] += pj->Q1 * mi * mag_Indj * (pi->phi - pj->phi) * dx[i];
  }
#endif
#ifdef MHD_VECPOT
  const float mag_Indi = wi_dr * r_inv / rhoi;
  const float mag_Indj = wj_dr * r_inv / rhoj;
  // ADVECTIVE GAUGE
  //float dv[3];
  //dv[0] = pi->v[0] - pj->v[0];
  //dv[1] = pi->v[1] - pj->v[1];
  //dv[2] = pi->v[2] - pj->v[2];
  //const float SourceAi = dv[0]*pi->APred[0] + dv[1]*pi->APred[1] + dv[2]*pi->APred[2];
  //const float SourceAj = dv[0]*pj->APred[0] + dv[1]*pj->APred[1] + dv[2]*pj->APred[2];
  // Normal Gauge
  float dA[3];
  for(int i=0;i<3;i++)
     dA[i] = pi->APred[i] - pj->APred[i];
     //dA[i] = pi->APred[i]/(rhoi*rhoi) + pj->APred[i]/(rhoj*rhoj);
  const float SourceAi = -(dA[0]*pi->v[0] + dA[1]*pi->v[1] + dA[2]*pi->v[2]);
  const float SourceAj = -(dA[0]*pj->v[0] + dA[1]*pj->v[1] + dA[2]*pj->v[2]);
  float SAi = SourceAi + (pi->GauPred - pj->GauPred);
  float SAj = SourceAj + (pi->GauPred - pj->GauPred);
///  float SAi = (pi->GauPred - pj->GauPred);
//  float SAj = (pi->GauPred - pj->GauPred);
     
  for(int i=0;i<3;i++)
  {
//     pi->dAdt[i] += mj *mag_Indi* (
//         pi->v[(i+1)%3]*(dA[(i+1)%3]*dx[i] - dA[i]*dx[(i+1)%3])
//       - pi->v[(i+2)%3]*(dA[(i+2)%3]*dx[i] - dA[i]*dx[(i+2)%3]) );
//     pj->dAdt[i] += mi *mag_Indj* (
//         pj->v[(i+1)%3]*(dA[(i+1)%3]*dx[i] - dA[i]*dx[(i+1)%3])
//       - pj->v[(i+2)%3]*(dA[(i+2)%3]*dx[i] - dA[i]*dx[(i+2)%3]) );
     
     pi->dAdt[i] += mj *mag_Indi* SAi *dx[i];
     pj->dAdt[i] += mi *mag_Indj* SAj *dx[i];
  }
  //Dissipation
  const float Deta = 0.0005f;
  const float mag_Disi = wi_dr * r_inv * rhoi / (rho_ij * rho_ij);
  const float mag_Disj = wj_dr * r_inv * rhoj / (rho_ij * rho_ij);
  for(int i=0;i<3;i++)
  {
     pi->dAdt[i] += mj * 2.0 * Deta * mag_Disi* dA[i];
     pj->dAdt[i] += mi * 2.0 * Deta * mag_Disj* dA[i];
  }
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
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
#ifndef MHD_BASE
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;
#else
  const float b2_i = (pi->BPred[0]*pi->BPred[0] + pi->BPred[1]*pi->BPred[1] + pi->BPred[2]*pi->BPred[2] );
  const float b2_j = (pj->BPred[0]*pj->BPred[0] + pj->BPred[1]*pj->BPred[1] + pj->BPred[2]*pj->BPred[2] ); 
  const float vcsa2_i = ci * ci + MU0_1 * b2_i/rhoi; 
  const float vcsa2_j = cj * cj + MU0_1 * b2_j/rhoj; 
  float Bpro2_i = (pi->BPred[0]*dx[0]+ pi->BPred[1]*dx[1]+ pi->BPred[2]*dx[2]) * r_inv;
        Bpro2_i *= Bpro2_i;
  float mag_speed_i = sqrtf(0.5 * (vcsa2_i + 
  		      sqrtf(max(  (vcsa2_i * vcsa2_i - 4. * ci * ci * Bpro2_i * MU0_1 / rhoi),0.0))));
  float Bpro2_j = (pj->BPred[0]*dx[0]+ pj->BPred[1]*dx[1]+ pj->BPred[2]*dx[2]) * r_inv;
        Bpro2_j *= Bpro2_j;
  float mag_speed_j = sqrtf(0.5 * (vcsa2_j + 
  		      sqrtf(max(  (vcsa2_j * vcsa2_j - 4. * cj * cj * Bpro2_j * MU0_1 / rhoj),0.0))));

  const float v_sig = (mag_speed_i + mag_speed_j - const_viscosity_beta * mu_ij);
//  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;
  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, v_sig);
  /* Grab balsara switches */
  //const float balsara_i = 1.f;
  //const float balsara_j = 1.f;
#endif
  /* Grab balsara switches */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;


  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  /* Construct the full viscosity term */
  const float rho_ij = rhoi + rhoj;
#ifdef MHD_BASE
  const float alpha = 2.f * (pi->viscosity.alpha + pj->viscosity.alpha);
#else
  const float alpha = pi->viscosity.alpha + pj->viscosity.alpha;
#endif
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
  /* Eventually got the MHD accel */ 
//_FORCE
  const float mag_faci = MU0_1 * f_ij * wi_dr * r_inv /(rhoi*rhoi);
  const float mag_facj = MU0_1 * f_ji * wj_dr * r_inv /(rhoj*rhoj);
  float mm_i[3][3],mm_j[3][3];
  
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
  	{
	mm_i[i][j] = pi->BPred[i]*pi->BPred[j];
	mm_j[i][j] = pj->BPred[i]*pj->BPred[j];
	 }
  for(int j=0;j<3;j++){
	 mm_i[j][j] -= 0.5 * (pi->BPred[0]*pi->BPred[0]+pi->BPred[1]*pi->BPred[1]+pi->BPred[2]*pi->BPred[2]);
  	 mm_j[j][j] -= 0.5 * (pj->BPred[0]*pj->BPred[0]+pj->BPred[1]*pj->BPred[1]+pj->BPred[2]*pj->BPred[2]);}
//////////////////////////////////
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
          pi->a_hydro[i] += mj * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj)*dx[j];
	  pi->a_hydro[i] -= pi->Q0 * mj * pi->BPred[i] * (pi->BPred[j]*mag_faci+pj->BPred[j]*mag_facj)*dx[j];
	 }
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

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pi->N_force++;
#endif
#ifdef MHD_DI
  const float mag_Indi = wi_dr * r_inv / rhoi;
  //const float mag_faci = MU0_1 * f_ij * wi_dr * r_inv /(rhoi*rhoi);
  float dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  for(int i=0;i<3;i++) 
  pi->dBdt[i] += mj * mag_Indi * ((pi->BPred[i] * dv[(i+1)%3] - pi->BPred[(i+1)%3] * dv[i]) * dx[(i+1)%3]
			        + (pi->BPred[i] * dv[(i+2)%3] - pi->BPred[(i+2)%3] * dv[i]) * dx[(i+2)%3]);
  for(int i=0;i<3;i++) 
  pi->dBdt[i] += pi->Q1 * mj * mag_Indi * (pi->phi - pj->phi) * dx[i];
#endif
#ifdef MHD_VECPOT
  const float mag_Indi = wi_dr * r_inv / rhoi;
  // ADVECTIVE GAUGE
  //float dv[3];
  //dv[0] = pi->v[0] - pj->v[0];
  //dv[1] = pi->v[1] - pj->v[1];
  //dv[2] = pi->v[2] - pj->v[2];
  //const float SourceAi = dv[0]*pi->APred[0] + dv[1]*pi->APred[1] + dv[2]*pi->APred[2];
  // Normal Gauge
  float dA[3];
  for(int i=0;i<3;i++)
     dA[i] = pi->APred[i] - pj->APred[i];
     //dA[i] = pi->APred[i]/(rhoi*rhoi) + pj->APred[i]/(rhoj*rhoj);
  //const float SourceAi = -(rhoi*rhoi)*(dA[0]*pi->v[0] + dA[1]*pi->v[1] + dA[2]*pi->v[2]);
  const float SourceAi = -(dA[0]*pi->v[0] + dA[1]*pi->v[1] + dA[2]*pi->v[2]);
  float SAi = SourceAi + (pi->GauPred - pj->GauPred);
  //for(int i=0;i<3;i++)
  //   pi->dAdt[i] += mj *mag_Indi* (
  //      pi->v[(i+1)%3]*(dA[(i+1)%3]*dx[i] - dA[i]*dx[(i+1)%3])
  //    - pi->v[(i+2)%3]*(dA[(i+2)%3]*dx[i] - dA[i]*dx[(i+2)%3]));
  //float SAi = (pi->GauPred - pj->GauPred);
  for(int i=0;i<3;i++)
     pi->dAdt[i] += mj *mag_Indi* SAi *dx[i];
  //Dissipation
  const float Deta = 0.0005f;
  const float mag_Disi = wi_dr * r_inv * rhoi / (rho_ij * rho_ij);
  for(int i=0;i<3;i++)
     pi->dAdt[i] += mj * 2.0 * Deta * mag_Disi* dA[i];
#endif
}
#endif /* SWIFT_SPHENIX_HYDRO_IACT_H */

