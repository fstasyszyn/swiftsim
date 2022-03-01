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
#ifdef MHD_BASE
  double dB[3];
#ifdef MHD_EULER
  double dalpha, dbeta;
#endif
#endif

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
  for(int i=0;i<3;++i)
  	dB[i]= pi->BPred[i] - pj->BPred[i];
  const double dBdr = dB[0]*dx[0] + dB[1]*dx[1] + dB[2]*dx[2];
  pi->divB -= faci * dBdr;
  pj->divB -= facj * dBdr;
#ifdef MHD_EULER
  dalpha = pi->ep[0] - pj->ep[0];
  dbeta  = pi->ep[1] - pj->ep[1];

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
  { 
  pi->Grad_ep[0][i] += faci * dalpha*dx[i];

  pi->Grad_ep[1][i] += faci * dbeta*dx[i];
  
  pj->Grad_ep[0][i] += facj * dalpha*dx[i];

  pj->Grad_ep[1][i] += facj * dbeta*dx[i];
  }
#endif  /* MHD_EULER */
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
#ifdef MHD_BASE
  double dB[3];
#ifdef MHD_EULER
  double dalpha, dbeta;
#endif
#endif

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
  for(int i=0;i<3;++i)
  	dB[i]= pi->BPred[i] - pj->BPred[i];
  const double dBdr = dB[0]*dx[0] + dB[1]*dx[1] + dB[2]*dx[2];
  pi->divB -= faci * dBdr;
#ifdef MHD_EULER
  dalpha = pi->ep[0] - pj->ep[0];
  dbeta  = pi->ep[1] - pj->ep[1];
  
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
  for(int i=0;i<3;++i){  
  pi->Grad_ep[0][i] += faci * dalpha * dx[i];
  
  pi->Grad_ep[1][i] += faci * dbeta  * dx[i];
  }
#endif  /* MHD_EULER */
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

#ifdef MHD_BASE
#ifdef MHD_EULER_TEST
  const float MU0_1 = 1.0;
#else
  const float MU0_1 = 1.0/(4.0*M_PI);
#endif
#endif

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

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;

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

  const float v_sig = (mag_speed_i + mag_speed_j - const_viscosity_beta/2.0 * mu_ij);
  /* Grab balsara switches */
  const float balsara_i = 1.f;
  const float balsara_j = 1.f;
#endif

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * (balsara_i + balsara_j) * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;

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
//#_FORCE
  const float mag_faci = MU0_1 * f_ij * wi_dr * r_inv /(rhoi*rhoi);
  const float mag_facj = MU0_1 * f_ji * wj_dr * r_inv /(rhoj*rhoj);
  float mm_i[3][3],mm_j[3][3];
//  
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
  	{
	mm_i[i][j] = pi->Bfld[i]*pi->Bfld[j];
	mm_j[i][j] = pj->Bfld[i]*pj->Bfld[j];
	}
  for(int j=0;j<3;j++)
  	{
	 mm_i[j][j] -= 0.5 * (pi->Bfld[0]*pi->Bfld[0]+pi->Bfld[1]*pi->Bfld[1]+pi->Bfld[2]*pi->Bfld[2]);
  	 mm_j[j][j] -= 0.5 * (pj->Bfld[0]*pj->Bfld[0]+pj->Bfld[1]*pj->Bfld[1]+pj->Bfld[2]*pj->Bfld[2]);
	 }
///////////////////////////////
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
    {  
       pi->a_hydro[i] += mj * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj) * dx[j];
       pj->a_hydro[i] -= mi * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj) * dx[j];
       pi->a_hydro[i] -= mj * pi->Bfld[i] * (pi->Bfld[j]*mag_faci+pj->Bfld[j]*mag_facj)*dx[j];
       pj->a_hydro[i] += mi * pj->Bfld[i] * (pi->Bfld[j]*mag_faci+pj->Bfld[j]*mag_facj)*dx[j];
     }
#endif
  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;
  const float sph_du_term_j = P_over_rho2_j * dvdr * r_inv * wj_dr;

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

#ifdef MHD_BASE
#ifdef MHD_EULER_TEST
  const float MU0_1 = 1.0;
#else
  const float MU0_1 = 1.0/(4.0*M_PI);
#endif
#endif

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

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;

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

  const float v_sig = (mag_speed_i + mag_speed_j - const_viscosity_beta/2.0 * mu_ij);
  /* Grab balsara switches */
  const float balsara_i = 1.f;
  const float balsara_j = 1.f;
#endif

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * (balsara_i + balsara_j) * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;

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
//_FORCE
  const float mag_faci = MU0_1 * f_ij * wi_dr * r_inv /(rhoi*rhoi);
  const float mag_facj = MU0_1 * f_ji * wj_dr * r_inv /(rhoj*rhoj);
  float mm_i[3][3],mm_j[3][3];
  
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
  	{
	mm_i[i][j] = pi->Bfld[i]*pi->Bfld[j];
	mm_j[i][j] = pj->Bfld[i]*pj->Bfld[j];
	 }
  for(int j=0;j<3;j++){
	 mm_i[j][j] -= 0.5 * (pi->Bfld[0]*pi->Bfld[0]+pi->Bfld[1]*pi->Bfld[1]+pi->Bfld[2]*pi->Bfld[2]);
  	 mm_j[j][j] -= 0.5 * (pj->Bfld[0]*pj->Bfld[0]+pj->Bfld[1]*pj->Bfld[1]+pj->Bfld[2]*pj->Bfld[2]);}
//////////////////////////////////
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
          pi->a_hydro[i] += mj * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj)*dx[j];
	  pi->a_hydro[i] -= mj * pi->Bfld[i] * (pi->Bfld[j]*mag_faci+pj->Bfld[j]*mag_facj)*dx[j];
	 }
#endif

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

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
}

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
