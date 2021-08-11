/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_GADGET2_HYDRO_IACT_H
#define SWIFT_GADGET2_HYDRO_IACT_H

/**
 * @file Gadget2/hydro_iact.h
 * @brief SPH interaction functions following the Gadget-2 version of SPH.
 *
 * The interactions computed here are the ones presented in the Gadget-2 paper
 * Springel, V., MNRAS, Volume 364, Issue 4, pp. 1105-1134.
 * We use the same numerical coefficients as the Gadget-2 code. When used with
 * the Spline-3 kernel, the results should be equivalent to the ones obtained
 * with Gadget-2 up to the rounding errors and interactions missed by the
 * Gadget-2 tree-code neighbours search.
 */

#include "cache.h"
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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float wi, wi_dx;
  float wj, wj_dx;
  float dv[3], curlvr[3];
#ifdef GADGET_MHD
  double dB[3];
#ifdef GADGET_MHD_EULER
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
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Compute the kernel function for pi */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the density */
  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute the kernel function for pj */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Compute contribution to the density */
  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);

  /* Compute contribution to the number of neighbours */
  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Compute dv dot r */
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

#ifdef GADGET_MHD
  for(int i=0;i<3;++i)
  	dB[i]= pi->Bfld[i] - pj->Bfld[i];
  const double dBdr = dB[0]*dx[0] + dB[1]*dx[1] + dB[2]*dx[2];
  pi->divB -= faci * dBdr;
  pj->divB -= facj * dBdr;
#endif
#ifdef GADGET_MHD_EULER
  dalpha = pi->ep[0] - pj->ep[0];
  dbeta  = pi->ep[1] - pj->ep[1];

#ifdef GADGET_MHD_EULER_TEST
 // BRIO_WU
#if GADGET_MHD_EULER_TEST==1
  float LBOX=1.0;
  dalpha = (( dalpha > LBOX/2.0 ) ? dalpha-LBOX : ( ( dalpha < -LBOX/2.0 ) ? dalpha+LBOX: dalpha));
  dbeta  = (( dbeta > 0.75*LBOX/2.0 ) ? dbeta-0.75*LBOX : ( ( dbeta < -0.75*LBOX/2.0 ) ? dbeta+0.75*LBOX: dbeta));
#else
// VORTEX
  const float LBOX=1.0;
  dbeta  = (( dbeta  > LBOX/2.0 ) ? dbeta-LBOX  : ( ( dbeta  < -LBOX/2.0 ) ? dbeta+LBOX : dbeta));
#endif
#endif
  for(int i=0;i<3;++i)  
  pi->Grad_ep[0][i] += faci * dalpha*dx[i];

  for(int i=0;i<3;++i)  
  pi->Grad_ep[1][i] += faci * dbeta*dx[i];
  
  for(int i=0;i<3;++i)  
  pj->Grad_ep[0][i] += facj * dalpha*dx[i];

  for(int i=0;i<3;++i)  
  pj->Grad_ep[1][i] += facj * dbeta*dx[i];
#endif

#ifdef DEBUG_INTERACTIONS_SPH
  /* Update ngb counters */
  if (pi->num_ngb_density < MAX_NUM_OF_NEIGHBOURS)
    pi->ids_ngbs_density[pi->num_ngb_density] = pj->id;
  ++pi->num_ngb_density;

  if (pj->num_ngb_density < MAX_NUM_OF_NEIGHBOURS)
    pj->ids_ngbs_density[pj->num_ngb_density] = pi->id;
  ++pj->num_ngb_density;
#endif
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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    const struct part *restrict pj, float a, float H) {

  float wi, wi_dx;
  float dv[3], curlvr[3];
#ifdef GADGET_MHD
  double dB[3];
#ifdef GADGET_MHD_EULER
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
  const float r_inv = 1.0f / r;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the density */
  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  const float fac = mj * wi_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  pi->density.div_v -= fac * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += fac * curlvr[0];
  pi->density.rot_v[1] += fac * curlvr[1];
  pi->density.rot_v[2] += fac * curlvr[2];

#ifdef GADGET_MHD
  for(int i=0;i<3;++i)
  	dB[i]= pi->Bfld[i] - pj->Bfld[i];
  const double dBdr = dB[0]*dx[0] + dB[1]*dx[1] + dB[2]*dx[2];
  pi->divB -= fac * dBdr;
#endif
#ifdef GADGET_MHD_EULER
  dalpha = pi->ep[0] - pj->ep[0];
  dbeta  = pi->ep[1] - pj->ep[1];
  
#ifdef GADGET_MHD_EULER_TEST
  //const float LBOX=0.5;
  //BRIO_WU
  //float LBOX=1.0;
  //dalpha = (( dalpha > LBOX/2.0 ) ? dalpha-LBOX : ( ( dalpha < -LBOX/2.0 ) ? dalpha+LBOX: dalpha));
  //LBOX=1.0;
  //dbeta  = (( dbeta > 0.75*LBOX/2.0 ) ? dbeta-0.75*LBOX : ( ( dbeta < -0.75*LBOX/2.0 ) ? dbeta+0.75*LBOX: dbeta));
  //VORTEX
  const float LBOX=1.0;
  dbeta  = (( dbeta  > LBOX/2.0 ) ? dbeta-LBOX  : ( ( dbeta  < -LBOX/2.0 ) ? dbeta+LBOX : dbeta));
#endif
  
  for(int i=0;i<3;++i)  
  pi->Grad_ep[0][i] += fac * dalpha * dx[i];
  
  for(int i=0;i<3;++i)  
  pi->Grad_ep[1][i] += fac * dbeta  * dx[i];
  
#endif


#ifdef DEBUG_INTERACTIONS_SPH
  /* Update ngb counters */
  if (pi->num_ngb_density < MAX_NUM_OF_NEIGHBOURS)
    pi->ids_ngbs_density[pi->num_ngb_density] = pj->id;
  ++pi->num_ngb_density;
#endif
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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float wi, wj, wi_dx, wj_dx;

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
  const float r_inv = 1.0f / r;

  /* Get some values in local variables. */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute h-gradient terms */
  const float f_i = pi->force.f;
  const float f_j = pj->force.f;

  /* Compute pressure terms */
  const float P_over_rho2_i = pi->force.P_over_rho2;
  const float P_over_rho2_j = pj->force.P_over_rho2;

  /* Compute sound speeds */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float MU0_1 = 1.0/(4.0*M_PI);
#ifndef GADGET_MHD
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;
#else
  // CHECK MU0
  const float b2_i = (pi->Bfld[0]*pi->Bfld[0] + pi->Bfld[1]*pi->Bfld[1] + pi->Bfld[2]*pi->Bfld[2] );
  const float b2_j = (pj->Bfld[0]*pj->Bfld[0] + pj->Bfld[1]*pj->Bfld[1] + pj->Bfld[2]*pj->Bfld[2] ); 
  float vcsa2_i = ci * ci + min(MU0_1 * b2_i/rhoi,10.0*ci*ci); 
  float vcsa2_j = cj * cj + min(MU0_1 * b2_j/rhoj,10.0*cj*cj); 
  //const float vcsa2_i = ci * ci + MU0_1 * b2_i/rhoi; 
  //const float vcsa2_j = cj * cj + MU0_1 * b2_j/rhoj; 
  float Bpro2_i = (pi->Bfld[0]*dx[0]+ pi->Bfld[1]*dx[1]+ pi->Bfld[2]*dx[2]) * r_inv;
        Bpro2_i *= Bpro2_i;
  float mag_speed_i = sqrt(0.5 * (vcsa2_i + 
  		      sqrt(max(  (vcsa2_i * vcsa2_i - 4. * ci * ci * Bpro2_i * MU0_1 / rhoi),0.0))));
  float Bpro2_j = (pj->Bfld[0]*dx[0]+ pj->Bfld[1]*dx[1]+ pj->Bfld[2]*dx[2]) * r_inv;
        Bpro2_j *= Bpro2_j;
  float mag_speed_j = sqrt(0.5 * (vcsa2_j + 
  		      sqrt(max(  (vcsa2_j * vcsa2_j - 4. * cj * cj * Bpro2_j * MU0_1 / rhoj),0.0))));

  const float v_sig = 3.0*(mag_speed_i + mag_speed_j - const_viscosity_beta/2.0 * mu_ij);
//  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;
#endif

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * mu_ij *10* (balsara_i + balsara_j) / rho_ij;

  /* Now, convolve with the kernel */
  const float visc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;
  const float sph_term =
      (f_i * P_over_rho2_i * wi_dr + f_j * P_over_rho2_j * wj_dr) * r_inv;

  /* Eventually got the acceleration */
  const float acc = visc_term + sph_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  pj->a_hydro[0] += mi * acc * dx[0];
  pj->a_hydro[1] += mi * acc * dx[1];
  pj->a_hydro[2] += mi * acc * dx[2];
  
  /* Eventually got the MHD accel */ 
#ifdef GADGET_MHD
  float mm_i[3][3],mm_j[3][3];
  const float mag_faci = MU0_1 * f_i * wi_dr * r_inv /(rhoi*rhoi);
  const float mag_facj = MU0_1 * f_j * wj_dr * r_inv /(rhoj*rhoj);
  //float magacc[3],magcorr[3];
  
  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
  	mm_i[i][j]=pi->Bfld[i]*pi->Bfld[j];
  for(int i=0;i<3;++i)
  	mm_i[i][i]-=0.5*b2_i;
 //   for(int j=0;j<3;++j)
  //	mm_i[i][i]-=0.5*pi->Bfld[j]*pi->Bfld[j];
  
  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
  	mm_j[i][j]=pj->Bfld[i]*pj->Bfld[j];
  for(int i=0;i<3;++i)
  	mm_j[i][i]-=0.5*b2_j;
  //  for(int j=0;j<3;++j)
  //	mm_j[i][i]-=0.5*pj->Bfld[j]*pj->Bfld[j];
  
  
  for(int i=0;i<2;++i)
    for(int j=0;j<3;++j)
     pi->a_hydro[i] += mj * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj)*dx[j];
 for(int i=0;i<2;++i)
    for(int j=0;j<3;++j)
     pi->a_hydro[i] -= mj * pi->Bfld[i] * (pi->Bfld[j]*mag_faci+pj->Bfld[j]*mag_facj)*dx[j];
  for(int i=0;i<2;++i)
    for(int j=0;j<3;++j)
     pj->a_hydro[i] += mi * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj)*dx[j];
  //Take out the divergence term  
  for(int i=0;i<2;++i)
    for(int j=0;j<3;++j)
     pj->a_hydro[i] -= mi * pj->Bfld[i] * (pi->Bfld[j]*mag_faci+pj->Bfld[j]*mag_facj)*dx[j];

#endif

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);

  /* Change in entropy */
  pi->entropy_dt += mj * visc_term * dvdr_Hubble;
  pj->entropy_dt += mi * visc_term * dvdr_Hubble;

#ifdef DEBUG_INTERACTIONS_SPH
  /* Update ngb counters */
  if (pi->num_ngb_force < MAX_NUM_OF_NEIGHBOURS)
    pi->ids_ngbs_force[pi->num_ngb_force] = pj->id;
  ++pi->num_ngb_force;

  if (pj->num_ngb_force < MAX_NUM_OF_NEIGHBOURS)
    pj->ids_ngbs_force[pj->num_ngb_force] = pi->id;
  ++pj->num_ngb_force;
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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    const struct part *restrict pj, float a, float H) {

  float wi, wj, wi_dx, wj_dx;

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
  const float r_inv = 1.0f / r;

  /* Get some values in local variables. */
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute h-gradient terms */
  const float f_i = pi->force.f;
  const float f_j = pj->force.f;

  /* Compute pressure terms */
  const float P_over_rho2_i = pi->force.P_over_rho2;
  const float P_over_rho2_j = pj->force.P_over_rho2;

  /* Compute sound speeds */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float MU0_1 = 1.0/(4.0*M_PI);
#ifndef GADGET_MHD
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;
#else
  // CHECK MU0
  const float b2_i = (pi->Bfld[0]*pi->Bfld[0] + pi->Bfld[1]*pi->Bfld[1] + pi->Bfld[2]*pi->Bfld[2] );
  const float b2_j = (pj->Bfld[0]*pj->Bfld[0] + pj->Bfld[1]*pj->Bfld[1] + pj->Bfld[2]*pj->Bfld[2] ); 
  float vcsa2_i = ci * ci + min(MU0_1 * b2_i/rhoi,10.0*ci*ci); 
  float vcsa2_j = cj * cj + min(MU0_1 * b2_j/rhoj,10.0*cj*cj); 
  //const float vcsa2_i = ci * ci + MU0_1 * b2_i/rhoi; 
  //const float vcsa2_j = cj * cj + MU0_1 * b2_j/rhoj; 
  float Bpro2_i = (pi->Bfld[0]*dx[0]+ pi->Bfld[1]*dx[1]+ pi->Bfld[2]*dx[2]) * r_inv;
        Bpro2_i *= Bpro2_i;
  float mag_speed_i = sqrt(0.5 * (vcsa2_i + 
  		      sqrt(max(  (vcsa2_i * vcsa2_i - 4. * ci * ci * Bpro2_i * MU0_1 / rhoi),0.0))));
  float Bpro2_j = (pj->Bfld[0]*dx[0]+ pj->Bfld[1]*dx[1]+ pj->Bfld[2]*dx[2]) * r_inv;
        Bpro2_j *= Bpro2_j;
  float mag_speed_j = sqrt(0.5 * (vcsa2_j + 
  		      sqrt(max(  (vcsa2_j * vcsa2_j - 4. * cj * cj * Bpro2_j * MU0_1 / rhoj),0.0))));

  const float v_sig = 3.0*(mag_speed_i + mag_speed_j - const_viscosity_beta/2.0 * mu_ij);
  //const float v_sig = ci + cj - const_viscosity_beta * mu_ij;
#endif

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * mu_ij *10* (balsara_i + balsara_j) / rho_ij;

  /* Now, convolve with the kernel */
  const float visc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;
  const float sph_term =
      (f_i * P_over_rho2_i * wi_dr + f_j * P_over_rho2_j * wj_dr) * r_inv;

  /* Eventually got the acceleration */
  const float acc = visc_term + sph_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];
  
  /* Eventually got the MHD accel */ 
#ifdef GADGET_MHD
  float mm_i[3][3],mm_j[3][3];
  const float mag_faci = MU0_1 * f_i * wi_dr * r_inv /(rhoi*rhoi);
  const float mag_facj = MU0_1 * f_j * wj_dr * r_inv /(rhoj*rhoj);
  
  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
  	mm_i[i][j]=pi->Bfld[i]*pi->Bfld[j];
  for(int i=0;i<3;++i)
  	mm_i[i][i]-=0.5*b2_i;
 //   for(int j=0;j<3;++j)
 // 	mm_i[i][i]-=0.5*pi->Bfld[j]*pi->Bfld[j];
  
  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
  	mm_j[i][j]=pj->Bfld[i]*pj->Bfld[j];
  for(int i=0;i<3;++i)
  	mm_j[i][i]-=0.5*b2_j;
//    for(int j=0;j<3;++j)
//  	mm_j[i][i]-=0.5*pj->Bfld[j]*pj->Bfld[j];
  
  
  for(int i=0;i<2;++i)
    for(int j=0;j<3;++j)
          pi->a_hydro[i] += mj * (mm_i[i][j]*mag_faci+mm_j[i][j]*mag_facj)*dx[j];
  //Take out the divergence term  
  for(int i=0;i<2;++i)
    for(int j=0;j<3;++j)
	  pi->a_hydro[i] -= mj * pi->Bfld[i] * (pi->Bfld[j]*mag_faci+pj->Bfld[j]*mag_facj)*dx[j];
#endif

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);

  /* Change in entropy */
  pi->entropy_dt += mj * visc_term * dvdr_Hubble;

#ifdef DEBUG_INTERACTIONS_SPH
  /* Update ngb counters */
  if (pi->num_ngb_force < MAX_NUM_OF_NEIGHBOURS)
    pi->ids_ngbs_force[pi->num_ngb_force] = pj->id;
  ++pi->num_ngb_force;
#endif
}

#endif /* SWIFT_GADGET2_HYDRO_IACT_H */
