/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_MINIMAL_HYDRO_IACT_H
#define SWIFT_MINIMAL_HYDRO_IACT_H

/**
 * @file Minimal/hydro_iact.h
 * @brief Minimal conservative implementation of SPH (Neighbour loop equations)
 *
 * The thermal variable is the internal energy (u). Simple constant
 * viscosity term with the Balsara (1995) switch. No thermal conduction
 * term is implemented.
 *
 * This corresponds to equations (43), (44), (45), (101), (103)  and (104) with
 * \f$\beta=3\f$ and \f$\alpha_u=0\f$ of Price, D., Journal of Computational
 * Physics, 2012, Volume 231, Issue 3, pp. 759-794.
 */

#include "adiabatic_index.h"
#include "hydro_parameters.h"
#include "minmax.h"

#ifdef MHD_EULER_TEST
#include "periodic.h"
#include "space.h"
#endif

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  float wi, wj, wi_dx, wj_dx;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

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

  /* Compute dv dot r */
  float dv[3], curlvr[3];

  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->density.div_v -= faci * dvdr;
  pj->density.div_v -= facj * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

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
  dalpha = nearest(dalpha, engine_extra_dims[2]);
  dbeta  = nearest(dbeta, 0.75*engine_extra_dims[1]);
#endif
#if MHD_EULER_TEST == 2
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
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

  float wi, wi_dx;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute dv dot r */
  float dv[3], curlvr[3];

  const float faci = mj * wi_dx * r_inv;

  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->density.div_v -= faci * dvdr;

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
  dalpha = nearest(dalpha, engine_extra_dims[2]);
  dbeta  = nearest(dbeta, 0.75*engine_extra_dims[1]);
#endif
#if MHD_EULER_TEST == 2
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
}

/**
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
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

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;
#ifdef MHD_ORESTIS
  /* Isotropic pressure */
  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];
  // const float isoPi = pressurei + 0.5f * B2i / const_vacuum_permeability;
  // const float isoPj = pressurej + 0.5f * B2j / const_vacuum_permeability;

  /* B dot r. */
  const float Bri = (Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2]);
  const float Brj = (Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2]);

  /* Compute gradient terms */
  const float over_rho2_i = 1.0f / (rhoi * rhoi) * f_ij;
  const float over_rho2_j = 1.0f / (rhoj * rhoj) * f_ji;
#else
  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;
#endif // MHD_ORESTIS

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
#ifndef MHD_BASE
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;
  /* Grab balsara switches */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;
#else
  const float b2_i = (pi->BPred[0]*pi->BPred[0] + pi->BPred[1]*pi->BPred[1] + pi->BPred[2]*pi->BPred[2] );
  const float b2_j = (pj->BPred[0]*pj->BPred[0] + pj->BPred[1]*pj->BPred[1] + pj->BPred[2]*pj->BPred[2] ); 
  float vcsa2_i = ci * ci + min(MU0_1 * b2_i/rhoi,10.0*ci*ci); 
  float vcsa2_j = cj * cj + min(MU0_1 * b2_j/rhoj,10.0*cj*cj); 
  //const float vcsa2_i = ci * ci + MU0_1 * b2_i/rhoi; 
  //const float vcsa2_j = cj * cj + MU0_1 * b2_j/rhoj; 
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
  /* Grab balsara switches */
  const float balsara_i = 1.f;
  const float balsara_j = 1.f;
#endif

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
#ifdef MHD_ORESTIS
  const float v_sig_hydro = ci + cj - const_viscosity_beta * mu_ij;
  const float visc =
      -0.25f * v_sig_hydro * (balsara_i + balsara_j) * mu_ij / rho_ij;
#else
  const float visc = -0.25f * v_sig * (balsara_i + balsara_j) * mu_ij / rho_ij;
#endif
  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;
#ifdef MHD_ORESTIS

  const float permeability_inv = MU0_1;
  const float monopole_beta = 1.f;

  /* SPH acceleration term in x direction, i_th particle */
  float sph_acc_term_i[3] = {0.f, 0.f, 0.f};

  /* Accelerations along X */

  /* Normal hydro SPH term */
  sph_acc_term_i[0] += pressurei * over_rho2_i * wi_dr * r_inv * dx[0];
  sph_acc_term_i[0] += pressurej * over_rho2_j * wj_dr * r_inv * dx[0];

  /* Isotropic MHD pressure term */
  sph_acc_term_i[0] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[0];
  sph_acc_term_i[0] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[0];

  /* Anisotropic MHD term */
  sph_acc_term_i[0] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[0];
  sph_acc_term_i[0] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[0];

  /* Accelerations along Y */

  /* Normal hydro SPH term */
  sph_acc_term_i[1] += pressurei * over_rho2_i * wi_dr * r_inv * dx[1];
  sph_acc_term_i[1] += pressurej * over_rho2_j * wj_dr * r_inv * dx[1];

  /* Isotropic MHD pressure term */
  sph_acc_term_i[1] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[1];
  sph_acc_term_i[1] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[1];

  /* Anisotropic MHD term */
  sph_acc_term_i[1] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[1];
  sph_acc_term_i[1] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[1];

  /* Accelerations along Z */

  /* Normal hydro SPH term */
  sph_acc_term_i[2] += pressurei * over_rho2_i * wi_dr * r_inv * dx[2];
  sph_acc_term_i[2] += pressurej * over_rho2_j * wj_dr * r_inv * dx[2];

  /* Isotropic MHD pressure term */
  sph_acc_term_i[2] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[2];
  sph_acc_term_i[2] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[2];

  /* Anisotropic MHD term */
  sph_acc_term_i[2] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[2];
  sph_acc_term_i[2] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[2];

  /* SPH acceleration term in x direction, j_th particle */
  float sph_acc_term_j[3];
  sph_acc_term_j[0] = -sph_acc_term_i[0];
  sph_acc_term_j[1] = -sph_acc_term_i[1];
  sph_acc_term_j[2] = -sph_acc_term_i[2];

  /* Divergence cleaning term */
  /* Manifestly *NOT* symmetric in i <-> j */

  sph_acc_term_i[0] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[0];
  sph_acc_term_i[0] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[0];

  sph_acc_term_i[1] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[1];
  sph_acc_term_i[1] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[1];

  sph_acc_term_i[2] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[2];
  sph_acc_term_i[2] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[2];

  sph_acc_term_j[0] -= monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bj[0];
  sph_acc_term_j[0] -= monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bj[0];

  sph_acc_term_j[1] -= monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bj[1];
  sph_acc_term_j[1] -= monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bj[1];

  sph_acc_term_j[2] -= monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bj[2];
  sph_acc_term_j[2] -= monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bj[2];

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * sph_acc_term_i[0] + mj * visc_acc_term * dx[0];
  pi->a_hydro[1] -= mj * sph_acc_term_i[1] + mj * visc_acc_term * dx[1];
  pi->a_hydro[2] -= mj * sph_acc_term_i[2] + mj * visc_acc_term * dx[2];

  pj->a_hydro[0] -= mi * sph_acc_term_j[0] - mi * visc_acc_term * dx[0];
  pj->a_hydro[1] -= mi * sph_acc_term_j[1] - mi * visc_acc_term * dx[1];
  pj->a_hydro[2] -= mi * sph_acc_term_j[2] - mi * visc_acc_term * dx[2];


  /* Get the time derivative for u. */
  const float sph_du_term_i = pressurei * over_rho2_i * dvdr * r_inv * wi_dr;
  const float sph_du_term_j = pressurej * over_rho2_j * dvdr * r_inv * wj_dr;
#else  // Basically MHD_BASE && NO MHD
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
  
#ifdef MHD_BASE 
  /* Eventually got the MHD accel */ 
//#_FORCE
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

#endif //MHD_ORESTIS
  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;
  const float du_dt_j = sph_du_term_j + visc_du_term;

  /* Internal energy time derivatibe */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr * f_ij;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr * f_ji;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);
#ifdef MHD_ORESTIS
  /* */
  const float dB_dt_pref_i = over_rho2_i * wi_dr * r_inv;
  const float dB_dt_pref_j = over_rho2_j * wj_dr * r_inv;

  /* */
  float dB_dt_i[3];
  dB_dt_i[0] = -Bri * (pi->v[0] - pj->v[0]);
  dB_dt_i[1] = -Bri * (pi->v[1] - pj->v[1]);
  dB_dt_i[2] = -Bri * (pi->v[2] - pj->v[2]);

  float dB_dt_j[3];
  dB_dt_j[0] = -Brj * (pi->v[0] - pj->v[0]);
  dB_dt_j[1] = -Brj * (pi->v[1] - pj->v[1]);
  dB_dt_j[2] = -Brj * (pi->v[2] - pj->v[2]);

  /* */
  pi->dBdt[0] += mj * dB_dt_pref_i * dB_dt_i[0];
  pi->dBdt[1] += mj * dB_dt_pref_i * dB_dt_i[1];
  pi->dBdt[2] += mj * dB_dt_pref_i * dB_dt_i[2];

  pj->dBdt[0] += mi * dB_dt_pref_j * dB_dt_j[0];
  pj->dBdt[1] += mi * dB_dt_pref_j * dB_dt_j[1];
  pj->dBdt[2] += mi * dB_dt_pref_j * dB_dt_j[2];
#endif //MHD_ORESTIS
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
  const float SourceAi = -(dA[0]*pi->v[0] + dA[1]*pi->v[1] + dA[2]*pi->v[2]);
  const float SourceAj = -(dA[0]*pj->v[0] + dA[1]*pj->v[1] + dA[2]*pj->v[2]);
  float SAi = SourceAi + (pi->GauPred - pj->GauPred);
  float SAj = SourceAj + (pi->GauPred - pj->GauPred);
  for(int i=0;i<3;i++)
  {
     pi->dAdt[i] += mj *mag_Indi* SAi *dx[i];
     pj->dAdt[i] += mi *mag_Indj* SAj *dx[i];
  }
  //Dissipation
/*  const float deta = 0.0;
  const float mag_Disi = wi_dr * r_inv * rhoi / (rho_ij * rho_ij);
  const float mag_Disj = wj_dr * r_inv * rhoj / (rho_ij * rho_ij);
  for(int i=0;i<3;i++)
  {
     pi->dAdt[i] += mj * deta * mag_Disi* dA[i];
     pj->dAdt[i] += mi * deta * mag_Disj* dA[i];
  }*/
#endif
}

/**
 * @brief Force interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {


#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
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

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;
#ifdef MHD_ORESTIS
  /* Isotropic pressure */
  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];
  // const float isoPi = pressurei + 0.5f * B2i / const_vacuum_permeability;
  // const float isoPj = pressurej + 0.5f * B2j / const_vacuum_permeability;

  /* B dot r. */
  const float Bri = (Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2]);
  const float Brj = (Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2]);

  /* Compute gradient terms */
  const float over_rho2_i = 1.0f / (rhoi * rhoi) * f_ij;
  const float over_rho2_j = 1.0f / (rhoj * rhoj) * f_ji;

#else
/* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;
#endif // MHD_ORESTIS
  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
#ifndef MHD_BASE
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;
  /* Grab balsara switches */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;
#else
  const float b2_i = (pi->BPred[0]*pi->BPred[0] + pi->BPred[1]*pi->BPred[1] + pi->BPred[2]*pi->BPred[2] );
  const float b2_j = (pj->BPred[0]*pj->BPred[0] + pj->BPred[1]*pj->BPred[1] + pj->BPred[2]*pj->BPred[2] ); 
  float vcsa2_i = ci * ci + min(MU0_1 * b2_i/rhoi,10.0*ci*ci); 
  float vcsa2_j = cj * cj + min(MU0_1 * b2_j/rhoj,10.0*cj*cj); 
  //const float vcsa2_i = ci * ci + MU0_1 * b2_i/rhoi; 
  //const float vcsa2_j = cj * cj + MU0_1 * b2_j/rhoj; 
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
  /* Grab balsara switches */
  const float balsara_i = 1.f;
  const float balsara_j = 1.f;
#endif

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
#ifdef MHD_ORESTIS
  const float v_sig_hydro = ci + cj - const_viscosity_beta * mu_ij;
  const float visc =
      -0.25f * v_sig_hydro * (balsara_i + balsara_j) * mu_ij / rho_ij;
#else
  const float visc = -0.25f * v_sig * (balsara_i + balsara_j) * mu_ij / rho_ij;
#endif

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;
#ifdef MHD_ORESTIS
  const float permeability_inv = MU0_1;
  const float monopole_beta = 1.f;

  /* SPH acceleration term in x direction, i_th particle */
  float sph_acc_term_i[3] = {0.f, 0.f, 0.f};

  /* Accelerations along X */

  /* Normal hydro SPH term */
  sph_acc_term_i[0] += pressurei * over_rho2_i * wi_dr * r_inv * dx[0];
  sph_acc_term_i[0] += pressurej * over_rho2_j * wj_dr * r_inv * dx[0];

  /* Isotropic MHD pressure term */
  sph_acc_term_i[0] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[0];
  sph_acc_term_i[0] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[0];

  /* Anisotropic MHD term */
  sph_acc_term_i[0] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[0];
  sph_acc_term_i[0] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[0];

  /* Accelerations along Y */

  /* Normal hydro SPH term */
  sph_acc_term_i[1] += pressurei * over_rho2_i * wi_dr * r_inv * dx[1];
  sph_acc_term_i[1] += pressurej * over_rho2_j * wj_dr * r_inv * dx[1];

  /* Isotropic MHD pressure term */
  sph_acc_term_i[1] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[1];
  sph_acc_term_i[1] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[1];

  /* Anisotropic MHD term */
  sph_acc_term_i[1] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[1];
  sph_acc_term_i[1] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[1];

  /* Accelerations along Z */

  /* Normal hydro SPH term */
  sph_acc_term_i[2] += pressurei * over_rho2_i * wi_dr * r_inv * dx[2];
  sph_acc_term_i[2] += pressurej * over_rho2_j * wj_dr * r_inv * dx[2];

  /* Isotropic MHD pressure term */
  sph_acc_term_i[2] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[2];
  sph_acc_term_i[2] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[2];

  /* Anisotropic MHD term */
  sph_acc_term_i[2] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[2];
  sph_acc_term_i[2] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[2];

  /* Divergence cleaning term */
  /* Manifestly *NOT* symmetric in i <-> j */

  sph_acc_term_i[0] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[0];
  sph_acc_term_i[0] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[0];

  sph_acc_term_i[1] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[1];
  sph_acc_term_i[1] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[1];

  sph_acc_term_i[2] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[2];
  sph_acc_term_i[2] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[2];

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * sph_acc_term_i[0] + mj * visc_acc_term * dx[0];
  pi->a_hydro[1] -= mj * sph_acc_term_i[1] + mj * visc_acc_term * dx[1];
  pi->a_hydro[2] -= mj * sph_acc_term_i[2] + mj * visc_acc_term * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i = pressurei * over_rho2_i * dvdr * r_inv * wi_dr;

#else
  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];
  
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
#endif // MHD_ORESTIS
  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;

  /* Internal energy time derivatibe */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr * f_ij;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);

#ifdef MHD_ORESTIS
  /* */
  const float dB_dt_pref_i = over_rho2_i * wi_dr * r_inv;

  /* */
  float dB_dt_i[3];
  dB_dt_i[0] = -Bri * (pi->v[0] - pj->v[0]);
  dB_dt_i[1] = -Bri * (pi->v[1] - pj->v[1]);
  dB_dt_i[2] = -Bri * (pi->v[2] - pj->v[2]);

  /* */
  pi->dBdt[0] += mj * dB_dt_pref_i * dB_dt_i[0];
  pi->dBdt[1] += mj * dB_dt_pref_i * dB_dt_i[1];
  pi->dBdt[2] += mj * dB_dt_pref_i * dB_dt_i[2];
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
  const float SourceAi = -(dA[0]*pi->v[0] + dA[1]*pi->v[1] + dA[2]*pi->v[2]);
  float SAi = SourceAi + (pi->GauPred - pj->GauPred);
  for(int i=0;i<3;i++)
     pi->dAdt[i] += mj *mag_Indi* SAi *dx[i];
  //Dissipation
/*  const float deta = 0.002;
  const float mag_Disi = wi_dr * r_inv * rhoi / (rho_ij * rho_ij);
  for(int i=0;i<3;i++)
     pi->dAdt[i] += mj * deta * mag_Disi* dA[i];
*/
#endif
}

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
