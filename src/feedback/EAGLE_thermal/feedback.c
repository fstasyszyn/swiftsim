/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* This file's header */
#include "feedback.h"

/* Local includes. */
#include "../EAGLE/enrichment.h"
#include "../EAGLE/yield_tables.h"
#include "hydro_properties.h"
#include "inline.h"
#include "random.h"
#include "timers.h"

/**
 * @brief Return the change in temperature (in internal units) to apply to a
 * gas particle affected by SNe feedback.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 */
double eagle_feedback_temperature_change(const struct spart* sp,
                                         const struct feedback_props* props) {

  /* In the EAGLE REF model, the change of temperature is constant */
  return props->SNe_deltaT_desired;
}

/**
 * @brief Compute the critical heating temperature increase, dT_crit.
 *
 * @param ngb_n The ambient density of the star.
 * @param ngb_Z The ambient metallicity of the star.
 * @param props The feedback model properties.
 */
double compute_SNII_dT_crit(const double ngb_n,
                            const double ngb_Z,
                            const double mean_ngb_mass,
                            const struct feedback_props* props) {

  /* Base line: analytic value from Dalla Vecchia & Schaye (2012) */
  double dT_crit = 3.162e7 * pow(ngb_n * 0.1, 0.6666667) *
      pow(mean_ngb_mass * props->mass_to_solar_mass * 1e-6, 0.33333333);

  /* Limit dT_crit to a specified minimum */
  dT_crit = max(dT_crit, props->SNII_dTcrit_floor);

  /* Finally, apply a density and metallicity dependent additional floor */
  const double dT_crit_plane =
      pow(ngb_n, props->SNII_dTcrit_exp_nH) *
      pow(ngb_Z / props->Z_Sol, props->SNII_dTcrit_exp_Z) *
      props->SNII_dTcrit_norm;
  dT_crit = max(dT_crit, dT_crit_plane);

  return dT_crit;
}

/**
 * @brief Compute the limiting mimimum heating temperature increase, dT_limit.
 *
 * @param ngb_n The ambient density of the star.
 * @param ngb_Z The ambient metallicity of the star.
 * @param props The feedback model properties.
 */
double compute_SNII_dT_limit(const double ngb_n,
                             const double ngb_Z,
                             const struct feedback_props* props) {

  const double Z = max(ngb_Z / props->Z_Sol, props->SNII_dT_Zmin);
  const double dT_limit =
      pow(ngb_n, props->SNII_dTlimit_exp_nH) *
      pow(Z, props->SNII_dTlimit_exp_Z) * props->SNII_dTlimit_norm;

  return dT_limit;
}

/**
 * @brief Compute the target sampling temperature increase, dT_sample.
 *
 * @param ngb_n The ambient density of the star.
 * @param ngb_Z The ambient metallicity of the star.
 * @param props The feedback model properties.
 */
double compute_SNII_dT_sample(
  const double ngb_n, const double ngb_Z, const double gamma_star,
  const double frac_SNII, const double ngb_SFR, const double ngb_rho_phys,
  const double G_Newton, const double SNe_energy, const double mean_ngb_mass,
  const double dt, struct spart* sp, const struct feedback_props* props) {

  /* Extract properties from SN model */
  const double num_to_heat = props->SNII_delta_T_num_ngb_to_heat;

  /* -- Increase base line sampling factor if desired -- */

  double n_single = -1.;
  switch (props->SNII_with_oversampling_timescale) {

    case eagle_SNII_timescale_none:
      n_single = frac_SNII * num_to_heat;
      break;

    case eagle_SNII_timescale_gasconsum:
      {
        /* Time scale for gas consumption by star formation.
         * If SFR = 0, set it to very large, in which case it has no effect. */
        const double dt_consum = (ngb_SFR > 0) ?
            ngb_rho_phys / ngb_SFR : FLT_MAX;

        /* Maximally, ask for num_to_heat heating events, even for long dt */
        const double f_dt = min(dt / dt_consum, 1.);
        n_single = num_to_heat * max(frac_SNII, f_dt);
      }
      break;

    case eagle_SNII_timescale_freefall:
      {
        const double dt_freefall =
            sqrt(3. * M_PI /
                 (32. * G_Newton * ngb_rho_phys));

        /* Maximally, ask for num_to_heat heating events, even for long dt */
        const double f_dt = min(dt / dt_freefall, 1.);
        n_single = num_to_heat * max(frac_SNII, f_dt);
      }
      break;

    default:
      error("Invalid SNII oversampling requirement!!!");
  }

   /* Sanity check to make sure that we have set n_single... */
    if (n_single < 0)
    error("Tiny problem in adaptive SN-dT: N_single=%g.", n_single);

  /* -- Compute sampling reduction factor -- */

  const double rho_birth_phys = sp->birth_density;
  const double sfr_birth_phys = sp->birth_star_formation_rate;
  const double m_initial = sp->mass_init;
  double nu = rho_birth_phys * sfr_birth_phys / (m_initial * m_initial) /
      gamma_star;

  /* Apply specified limit to nu */
  nu = min(nu, props->SNII_maximum_nu);

  /* Normally, we want nu >= 1, but this limit can be disabled */
  if (!props->SNII_with_nu_below_one)
    nu = max(nu, 1.0);

  /* Record current sampling reduction factor */
  sp->nu = (float) nu;

  /* -- Calculate sampling temperature -- */
  const double dT_sample = SNe_energy /
      (mean_ngb_mass * props->temp_to_u_factor * n_single / nu);

  return dT_sample;
}

/* Returns the modelled SN efficiency at sub-critical dT.
 *
 * This is the fraction of the injected energy that is assumed to go into
 * doing actual work, rather than being lost to numerical overcooling.
 * 
 * @param dT The temperature increase for which to calculate the efficiency.
 * @param dT_crit The temperature above which the efficiency is 1.
 * @param dT_limit The temperature below which the efficiency is 0.
 * @param zeta The power-law index of eta between dT_limit and dT_crit.
 */
INLINE static double SN_eta(
  const double dT, const double dT_crit, const double dT_limit,
  const double zeta) {

  /* If dT < dT_limit, then all is lost */
  if (dT <= dT_limit)
    return 0.;

  /* If dT > dT_crit, then all is well */
  if (dT >= dT_crit)
    return 1.;

  /* Otherwise, use the appropriate power-law model (index zeta) */
  const double tau = (dT - dT_limit) / (dT_crit - dT_limit);
  return pow(tau, zeta);
}

/**
 * @brief Compute the "compromise" temperature increase that results in the
 *        right amount of energy received by the gas if dT must be < dT_crit.
 *
 * @param dT_sample The ideal dT for satisfying the sampling criterion.
 * @param dT_crit The temperature at which feedback would be numerically fully
 *                efficient.
 * @param props The feedback model properties.
 */
double compute_compromise_dT(
  const double dT_sample, const double dT_crit, const double dT_limit,
  const struct feedback_props* props) {

  /* Extract relevant model parameters. */
  const double zeta = props->SNII_efficiency_zeta;

  /* Compute minimum allowed dT */
  const double eta_min = props->SNII_efficiency_eta_min;
  const double dT_min = pow(eta_min, 1./zeta) * (dT_crit - dT_limit) + dT_limit;

  /* If energy compensation is off, stay as close to dT_sample as possible. */
  if (!props->SNII_with_energy_compensation)
    return max(dT_min, dT_sample);

  /* Check whether dT_min gives an effective dT that is above
   * dT_sample. If this is the case, we would have to use dT_comp < dT_min
   * to get an effective dT_sample, which is not allowed. Use dT_min. */
  if (dT_min * SN_eta(dT_min, dT_crit, dT_limit, zeta) >= dT_sample)
    return dT_min;

  /* If we get here, we can find a compromise dT between dT_min and dT_crit
   * that gives an effective dT_sample. Start in the middle of the two, and
   * iterate with a Newton-Raphson scheme (should always work...) */
  double dT = (dT_min + dT_crit) / 2.;
  for (int ii = 0; ii < 1000; ii++) {

    const double eta = SN_eta(dT, dT_crit, dT_limit, zeta);
    const double tau = (dT - dT_limit) / (dT_crit - dT_limit);

    /* Effective temperature offset from dT_sample and its derivative */
    const double f = dT * eta - dT_sample;
    const double f_dash =
        eta + dT * zeta * pow(tau, zeta - 1.) / (dT_crit - dT_limit);

    /* Updated estimate of dT according to Newton-Raphson scheme.
     * The value should never get below dT_min, but better be explicit... */
    dT = max(dT - f / f_dash, dT_min);
    
    /* Check whether we are sufficiently close to the true answer */
    const double dT_eff = SN_eta(dT, dT_crit, dT_limit, zeta) * dT;
    if ((fabs(dT_eff - dT_sample) / dT_sample) < 1e-4)
      break;
  }
  
  /* Make sure we have not gone too high somehow... */
  if (dT > dT_crit * 1.001)
    error("dT_comp=%g > dT_crit=%g!", dT, dT_crit);

  /* Make sure we have not somehow remained unconverged... */
  const double dT_eff = SN_eta(dT, dT_crit, dT_limit, zeta) * dT;
  if ((fabs(dT_eff - dT_sample) / dT_sample) >= 1e-4)
    error("Unconverged dT_comp=%g: dT_sample=%g, dT_eff=%g "
          "(dT_crit=%g, dT_limit=%g, zeta=%g)",
         dT, dT_sample, dT_eff, dT_crit, dT_limit, zeta);

  return dT;
}

/**
 * @brief Return the variable change in temperature (in internal units) to
 * apply to a gas particle affected by SNe feedback.
 *
 * @param sp The #spart.
 * @param ngb_gas_mass The total mass of the star's gas neighbours
 * @param num_gas_ngbs The total number of gas neighbours of the star
 * @param ngb_nH_cgs The SPH-smoothed gas density within the kernel [cgs]
 * @param p_SNe_energy Pointer to the SN energy fraction of this star
 * @param props The feedback model
 * @param frac_SNII The fraction of the star's SNII going off in this step
 * @param h_pkpc_inv Inverse of the star's smoothing length [kpc^-1]
 * @param dt The length of the star's current time step
 * @param G_Newton Newton's gravitational constant
 * @param ngb_SFR The SPH-smoothed SFR of the neighbour gas particles
 * @param ngb_rho_phys The SPH-smoothed gas density within the kernel [code]
 * @param ngb_Z The SPH-smoothed gas metallicity within the kernel
 */
double eagle_variable_feedback_temperature_change(
  struct spart* sp, const double ngb_gas_mass, const int num_gas_ngbs,
  const double ngb_nH_cgs, double* p_SNe_energy,
  const struct feedback_props* props, const double frac_SNII,
  const double h_pkpc_inv,
  const double dt, const double G_Newton, const double ngb_SFR,
  const double ngb_rho_phys, const double ngb_Z) {

  /* Safety check: should only get here if we run with the adaptive-dT model,
   * version 3 */
  if (props->SNII_use_variable_delta_T != 3)
    error("Attempting to compute variable SNII heating temperature "
          "without activating this model, v3. Cease and desist.");

  /* Model parameters */
  const double f_crit = props->SNII_T_crit_factor;
  const double dT_max = props->SNII_delta_T_max;  
  const double zeta = props->SNII_efficiency_zeta;
  const double omega_max = (props->SNII_with_energy_compensation) ?
      props->SNII_delta_T_omega_max : 1.0;
  const double gamma_star = props->SNII_gamma_star *
      ((props->SNII_sampling_reduction_within_smoothing_length) ?
       (h_pkpc_inv * h_pkpc_inv * h_pkpc_inv) : 1.0);

  /* Relevant star properties */
  const double mean_ngb_mass = ngb_gas_mass / ((double)num_gas_ngbs);
  const double n_phys = (props->SNII_use_instantaneous_density_for_dT ?
      ngb_nH_cgs : sp->birth_density * props->rho_to_n_cgs);
  const double Z_star = (props->SNII_use_instantaneous_Z_for_dT) ?
      ngb_Z : chemistry_get_star_total_metal_mass_fraction_for_feedback(sp);

  /* Calculate critical temperature */
  const double dT_crit = compute_SNII_dT_crit(
      n_phys, Z_star, mean_ngb_mass, props);

  /* Calculate minimum effective temperature */
  const double dT_limit = compute_SNII_dT_limit(n_phys, Z_star, props);

  /* Calculate target sampling temperature */
  const double dT_sample =
      compute_SNII_dT_sample(n_phys, Z_star, gamma_star, frac_SNII, ngb_SFR,
        ngb_rho_phys, G_Newton, *p_SNe_energy, mean_ngb_mass, dt, sp, props);

  /* Begin decision logic... */
  double dT = -1;
  double omega = 1.0;

  if (f_crit * dT_crit < dT_sample) {
    /* Case A: critical temperature below sampling temperature --> fine
     * (the critical temperature is always above dT_limit) */
    dT = dT_crit * f_crit;
  } else {
    /* Case B: critical temperature above sampling --> need to compromise */
    dT = compute_compromise_dT_v3(dT_sample, dT_crit, dT_limit, props);

    /* Compromise dT is guaranteed to be within partly efficient region,
     * and may be much higher than dT_sample. To avoid boosting the energy
     * too much, limit omega to the specified maximum. */
    omega = (props->SNII_omega_by_sampling) ? dT / dT_sample :
        1. / SN_eta(dT, dT_crit, dT_limit, zeta);
    omega = min(omega, omega_max); 
  }

  /* Apply maximum dT ceiling, correcting for expected cooling losses */
  if (dT > dT_max) {
    dT = dT_max;
    if (dT > dT_crit * 1.01)
      error("Internal logic error: forced to heat at dT_max=%g, "
            "but dT_crit=%g", dT_max, dT_crit);
    if (props->SNII_with_energy_compensation) {
      omega = 1. / SN_eta(dT, dT_crit, dT_limit, zeta);
 
      /* Limit omega energy boost to maximum (should rarely apply) */
      omega = min(omega, omega_max);
    }
  }

  /* Apply (potential) SNe energy increase to achieve desired sampling */
  if (omega > 1.01 && !props->SNII_with_energy_compensation)
    error("Running without energy compensation, but omega=%g!", omega);
  *p_SNe_energy *= omega;

  if (dT < 0)
    error("Ahm... I'm not going to set a negative heating temperature "
          "(star ID=%lld, dT=%g).",
          sp->id, dT);

  /* Record how delta_T compares to targets */
  const double critical_fraction = dT / dT_crit;
  const double sampling_fraction = dT_sample / dT * omega;

  sp->delta_T_min = (float) min(sp->delta_T_min, dT);
  sp->delta_T_max = (float) max(sp->delta_T_max, dT);
  sp->T_critical_fraction_min = (float) min(sp->T_critical_fraction_min,
                                    critical_fraction);
  sp->T_critical_fraction_max = (float) max(sp->T_critical_fraction_max,
                                    critical_fraction);
  sp->T_sampling_fraction_min = (float) min(sp->T_sampling_fraction_min,
                                    sampling_fraction);
  sp->T_sampling_fraction_max = (float) max(sp->T_sampling_fraction_max,
                                    sampling_fraction);

  sp->delta_T = (float) dT;
  sp->dT_crit = (float) dT_crit;
  sp->dT_sample = (float) dT_sample;
  sp->omega = (float) omega;
  sp->dT_limit = (float) dT_limit;
  sp->eta = (float) SN_eta(dT, dT_crit, dT_limit, zeta);  

  sp->ngb_rho = (float) ngb_rho_phys;
  sp->ngb_Z = (float) ngb_Z;

  sp->omega_min = (float) min(sp->omega_min, omega);
  sp->omega_max = (float) max(sp->omega_max, omega);

  return dT;
}

/**
 * @brief Computes the number of supernovae of type II exploding for a given
 * star particle assuming that all the SNII stars go off at once.
 *
 * @param sp The #spart.
 * @param props The properties of the stellar model.
 */
double eagle_feedback_number_of_SNII(const struct spart* sp,
                                     const struct feedback_props* props) {

  /* Note: For a Chabrier 2003 IMF and SNII going off
   * - between 6 and 100 M_sun, the first term is 0.017362 M_sun^-1 (EAGLE)
   * - between 8 and 100 M_sun, the first term is 0.011801 M_sun^-1 (EAGLE-XL)
   */
  return props->num_SNII_per_msun * sp->mass_init * props->mass_to_solar_mass;
}

/**
 * @brief Computes the number of supernovae of type II exploding for a given
 * star particle between two mass limits
 *
 * @param sp The #spart.
 * @param props The properties of the stellar model.
 * @param min_dying_mass_Msun Minimal star mass dying this step (in solar
 * masses).
 * @param max_dying_mass_Msun Maximal star mass dying this step (in solar
 * masses).
 */
double eagle_feedback_number_of_sampled_SNII(const struct spart* sp,
                                             const struct feedback_props* props,
                                             const double min_dying_mass_Msun,
                                             const double max_dying_mass_Msun) {

  /* The max dying star mass is below the SNII mass window
   * --> No SNII */
  if (max_dying_mass_Msun < props->SNII_min_mass_msun) return 0.;

  /* The min dying star mass is above the SNII mass window
   * --> No SNII */
  if (min_dying_mass_Msun > props->SNII_max_mass_msun) return 0.;

  /* Ok, we have some overlap with the SNII mass window. */

  double log10_min_mass_Msun = -10.;
  double log10_max_mass_Msun = -10.;

  /* The min dying star mass dies inside the SNII mass window */
  if (min_dying_mass_Msun <= props->SNII_max_mass_msun &&
      min_dying_mass_Msun > props->SNII_min_mass_msun) {

    /* Now, check the max dying star mass */

    /* The max dying star mass is also inside the SNII mass window */
    if (max_dying_mass_Msun <= props->SNII_max_mass_msun) {
      log10_min_mass_Msun = log10(min_dying_mass_Msun);
      log10_max_mass_Msun = log10(max_dying_mass_Msun);
    }

    /* The max dying star is above the SNII mass window */
    else {
      log10_min_mass_Msun = log10(min_dying_mass_Msun);
      log10_max_mass_Msun = props->log10_SNII_max_mass_msun;
    }

  }

  /* The min dying star mass dies below the SNII mass window */
  else if (min_dying_mass_Msun <= props->SNII_min_mass_msun) {

    /* Now, check the max dying star mass */

    /* The max dying star mass is inside the SNII mass window */
    if (max_dying_mass_Msun > props->SNII_min_mass_msun &&
        max_dying_mass_Msun <= props->SNII_max_mass_msun) {
      log10_min_mass_Msun = props->log10_SNII_min_mass_msun;
      log10_max_mass_Msun = log10(max_dying_mass_Msun);
    }

    /* The max dying star is above the SNII mass window */
    else if (max_dying_mass_Msun > props->SNII_max_mass_msun) {
      log10_min_mass_Msun = props->log10_SNII_min_mass_msun;
      log10_max_mass_Msun = props->log10_SNII_max_mass_msun;
    }

    /* The max dying star mass is also below the SNII mass window */
    else {

      /* We already excluded this at the star of the function */
#ifdef SWIFT_DEBUG_CHECKS
      error("Error in the logic");
#endif
    }
  }

  /* The min dying star mass dies above the SNII mass window */
  else {

    /* We already excluded this at the star of the function */
#ifdef SWIFT_DEBUG_CHECKS
    error("Error in the logic");
#endif
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (log10_min_mass_Msun == -10. || log10_max_mass_Msun == -10.)
    error("Something went wrong in the calculation of the number of SNII.");
#endif

  /* Calculate how many supernovae have exploded in this timestep
   * by integrating the IMF between the bounds we chose */
  const double num_SNII_per_msun =
      integrate_imf(log10_min_mass_Msun, log10_max_mass_Msun,
                    eagle_imf_integration_no_weight, NULL, props);

  return num_SNII_per_msun * sp->mass_init * props->mass_to_solar_mass;
}

/**
 * @brief Computes the number of supernovae of type Ia exploding for a given
 * star particle between time t0 and t1
 *
 * @param M_init The inital mass of the star particle in internal units.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param props The properties of the stellar model.
 */
double eagle_feedback_number_of_SNIa(const double M_init, const double t0,
                                     const double t1,
                                     const struct feedback_props* props) {

  double num_SNIa_per_Msun = 0.;

#ifdef SWIFT_DEBUG_CHECKS
  if (t1 < t0) error("Negative time range!");
  if (t0 < props->SNIa_DTD_delay_Gyr)
    error("Initial time smaller than the delay time!");
#endif

  switch (props->SNIa_DTD) {

    case eagle_feedback_SNIa_DTD_exponential: {

      /* We follow Foerster et al. 2006, MNRAS, 368 */

      /* The calculation is written as the integral between t0 and t1 of
       * eq. 3 of Schaye 2015 paper. */
      const double tau = props->SNIa_DTD_exp_timescale_Gyr_inv;
      const double nu = props->SNIa_DTD_exp_norm;
      num_SNIa_per_Msun = nu * (exp(-t0 * tau) - exp(-t1 * tau));
      break;
    }

    case eagle_feedback_SNIa_DTD_power_law: {

      /* We follow Graur et al. 2011, MNRAS, 417 */

      const double norm = props->SNIa_DTD_power_law_norm;
      num_SNIa_per_Msun = norm * log(t1 / t0);
      break;
    }

    default: {
      num_SNIa_per_Msun = 0.;
      error("Invalid choice of SNIa delay time distribution!");
    }
  }

  return num_SNIa_per_Msun * M_init * props->mass_to_solar_mass;
}

/**
 * @brief Computes the fraction of the available super-novae energy to
 * inject for a given event.
 *
 * Note that the fraction can be > 1.
 *
 * We use equation 7 of Schaye et al. 2015.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 * @param ngb_nH_cgs Hydrogen number density of the gas surrounding the star
 * (physical cgs units).
 * @param ngb_Z Metallicity (metal mass fraction) of the gas surrounding the
 * star.
 */
double eagle_feedback_energy_fraction(const struct spart* sp,
                                      const struct feedback_props* props,
                                      const double ngb_nH_cgs,
                                      const double ngb_Z) {

  /* Model parameters */
  const double f_E_max = props->f_E_max;
  const double f_E_min = props->f_E_min;
  const double Z_0 = props->Z_0;
  const double n_0 = props->n_0_cgs;
  const double n_Z = props->n_Z;
  const double n_n = props->n_n;

  /* Metallicity (metal mass fraction) at birth time of the star */
  const double Z_birth =
      chemistry_get_star_total_metal_mass_fraction_for_feedback(sp);

  /* Physical density of the gas at the star's birth time */
  const double rho_birth = sp->birth_density;
  const double n_birth_cgs = rho_birth * props->rho_to_n_cgs;

  /* Choose either the birth properties or current properties */
  const double nH =
      props->use_birth_density_for_f_th ? n_birth_cgs : ngb_nH_cgs;
  const double Z = props->use_birth_Z_for_f_th ? Z_birth : ngb_Z;

  /* Calculate f_E */
  const double Z_term = pow(max(Z, 1e-6) / Z_0, n_Z);
  const double n_term = pow(nH / n_0, -n_n);
  const double denonimator = 1. + Z_term * n_term;

  return f_E_min + (f_E_max - f_E_min) / denonimator;
}

/**
 * @brief Compute the properties of the SNII stochastic feedback energy
 * injection.
 *
 * Only does something if the particle reached the SNII age during this time
 * step.
 *
 * @param sp The star particle.
 * @param star_age Age of star at the beginning of the step in internal units.
 * @param dt Length of time-step in internal units.
 * @param ngb_gas_N Total (integer) number of gas neighbours within the
 * star's kernel.
 * @param ngb_gas_mass Total un-weighted mass in the star's kernel (internal
 * units)
 * @param ngb_nH_cgs Hydrogen number density of the gas surrounding the star
 * (physical cgs units).
 * @param ngb_Z Metallicity (metal mass fraction) of the gas surrounding the
 * star.
 * @param feedback_props The properties of the feedback model.
 * @param min_dying_mass_Msun Minimal star mass dying this step (in solar
 * masses).
 * @param max_dying_mass_Msun Maximal star mass dying this step (in solar
 * masses).
 * @param ti_begin ???
 * @param h_pkpc_inv The inverse of the smoothing length in kpc
 * @param G_Newton Newton's gravitational constant
 * @param ngb_SFR The SPH-smoothed SFR of the gas neighbour particles
 * @param ngb_rho_phys The SPH-smoothed density within the kernel
 */
INLINE static void compute_SNII_feedback(
    struct spart* sp, const double star_age, const double dt,
    const int ngb_gas_N, const float ngb_gas_mass, const double ngb_nH_cgs,
    const double ngb_Z, const struct feedback_props* feedback_props,
    const double min_dying_mass_Msun, const double max_dying_mass_Msun,
    const integertime_t ti_begin, const double h_pkpc_inv,
    const double G_Newton, const double ngb_SFR, const double ngb_rho_phys) {

  /* Are we sampling the delay function or using a fixed delay? */
  const int SNII_sampled_delay = feedback_props->SNII_sampled_delay;

  /* Time after birth considered for SNII feedback (internal units)
   * when using a fixed delay */
  const double SNII_wind_delay = feedback_props->SNII_wind_delay;

  /* Are we doing feedback this step?
   * Note that since the ages are calculated using an interpolation table we
   * must allow some tolerance here*/
  if ((SNII_sampled_delay) || (star_age <= SNII_wind_delay &&
                               (star_age + 1.001 * dt) > SNII_wind_delay)) {

    /* Make sure a star does not do feedback twice
     * when using a fixed delay! */
    if (!SNII_sampled_delay && sp->f_E != -1.f) {
#ifdef SWIFT_DEBUG_CHECKS
      message("Star has already done feedback! sp->id=%lld age=%e d=%e", sp->id,
              star_age, dt);
#endif
      return;
    }

    /* Properties of the model (all in internal units) */
    const double delta_T =
        eagle_feedback_temperature_change(sp, feedback_props);
    const double E_SNe = feedback_props->E_SNII;
    const double f_E =
        eagle_feedback_energy_fraction(sp, feedback_props, ngb_nH_cgs, ngb_Z);

    /* Number of SNe at this time-step */
    double N_SNe;
    double frac_SNII;
    if (SNII_sampled_delay) {
      N_SNe = eagle_feedback_number_of_sampled_SNII(
          sp, feedback_props, min_dying_mass_Msun, max_dying_mass_Msun);
      frac_SNII = N_SNe / eagle_feedback_number_of_SNII(sp, feedback_props);
    } else {
      N_SNe = eagle_feedback_number_of_SNII(sp, feedback_props);
      frac_SNII = 1.0;
    }

    /* Abort if there are no SNe exploding this step */
    if (N_SNe <= 0.) return;

    /* Conversion factor from T to internal energy */
    const double conv_factor = feedback_props->temp_to_u_factor;

    /* Calculate the heating temperature increment */
    double SNe_energy = f_E * E_SNe * N_SNe;
    const double delta_T = (feedback_props->SNII_use_variable_delta_T) ?
      eagle_variable_feedback_temperature_change(
          sp, ngb_gas_mass, ngb_gas_N, ngb_nH_cgs, &SNe_energy,
          feedback_props, frac_SNII, h_pkpc_inv, dt,
          G_Newton, ngb_SFR, ngb_rho_phys, ngb_Z) :
      eagle_feedback_temperature_change(sp, feedback_props);

    /* Calculate the default heating probability (accounting for round-off) */
    double prob = SNe_energy / (conv_factor * delta_T * ngb_gas_mass);
    prob = max(prob, 0.0);

    /* Calculate the change in internal energy of the gas particles that get
     * heated */
    double delta_u;

    /* Number of SNII events for this stellar particle */
    int number_of_SN_events = 0;

    if (prob <= 1.) {

      /* Normal case */
      delta_u = delta_T * conv_factor;

      for (int i = 0; i < ngb_gas_N; i++) {
        const double rand_thermal = random_unit_interval_part_ID_and_ray_idx(
            sp->id, i, ti_begin, random_number_stellar_feedback_3);
        if (rand_thermal < prob) number_of_SN_events++;
      }

    } else {

      /* Special case: we need to adjust the energy irrespective of the
         desired deltaT to ensure we inject all the available energy. */
      delta_u = f_E * E_SNe * N_SNe / ngb_gas_mass;

      /* Number of SNII events is equal to the number of Ngbs */
      number_of_SN_events = ngb_gas_N;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (f_E < feedback_props->f_E_min || f_E > feedback_props->f_E_max)
      error("f_E is not in the valid range! f_E=%f sp->id=%lld", f_E, sp->id);
#endif

    /* If we have more heating events than the maximum number of
     * rays (eagle_feedback_number_of_rays), then we cannot
     * distribute all of the heating events (since 1 event = 1 ray), so we need
     * to increase the thermal energy per ray and make the number of events
     * equal to the number of rays */
    if (number_of_SN_events > eagle_SNII_feedback_num_of_rays) {
      const double alpha_thermal =
          (double)number_of_SN_events / (double)eagle_SNII_feedback_num_of_rays;
      delta_u *= alpha_thermal;
      number_of_SN_events = eagle_SNII_feedback_num_of_rays;
    }

    /* Current total f_E for this star */
    double star_f_E = sp->f_E * sp->number_of_SNII_events;

    /* New total */
    star_f_E = (star_f_E + f_E) / (sp->number_of_SNII_events + 1.);

    /* Store all of this in the star for delivery onto the gas and recording */
    sp->f_E = star_f_E;
    sp->number_of_SNII_events++;
    sp->feedback_data.to_distribute.SNII_delta_u = delta_u;
    sp->feedback_data.to_distribute.SNII_num_of_thermal_energy_inj =
        number_of_SN_events;
  }
}

/**
 * @brief calculates stellar mass in spart that died over the timestep, calls
 * functions to calculate feedback due to SNIa, SNII and AGB
 *
 * @param feedback_props feedback_props data structure
 * @param phys_const The physical constants in internal units.
 * @param cosmo The cosmological model.
 * @param sp spart that we're evolving
 * @param us unit_system data structure
 * @param age age of spart at beginning of step
 * @param dt length of current timestep
 * @param ti_begin The current integer time (for random number hashing).
 */
void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct hydro_props* hydro_props,
                               const struct phys_const* phys_const,
                               const struct cosmology* cosmo, struct spart* sp,
                               const struct unit_system* us, const double age,
                               const double dt, const integertime_t ti_begin) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (age < 0.f) error("Negative age for a star.");

  if (sp->count_since_last_enrichment != 0)
    error("Computing feedback on a star that should not");
#endif

  /* Convert dt and stellar age from internal units to Gyr. */
  const double Gyr_inv = 1. / (phys_const->const_year * 1e9);
  const double dt_Gyr = dt * Gyr_inv;
  const double star_age_Gyr = age * Gyr_inv;

  /* Get the birth mass of the star */
  const double M_init = sp->mass_init;

  /* Get the total metallicity (metal mass fraction) at birth time and impose a
   * minimum */
  const double Z =
      max(chemistry_get_star_total_metal_mass_fraction_for_feedback(sp),
          exp10(log10_min_metallicity));

  /* Get the individual abundances (mass fractions at birth time) */
  const float* const abundances =
      chemistry_get_star_metal_mass_fraction_for_feedback(sp);

  /* Properties collected in the stellar density loop. */
  const int ngb_gas_N = sp->feedback_data.to_collect.ngb_N;
  const float ngb_gas_mass = sp->feedback_data.to_collect.ngb_mass;
  const float ngb_gas_Z = sp->feedback_data.to_collect.ngb_Z;
  const float ngb_gas_SFR = sp->feedback_data.to_collect.ngb_SFR;
  const float ngb_gas_rho_phys =
      sp->feedback_data.to_collect.ngb_rho * cosmo->a3_inv;
  const float ngb_gas_phys_nH_cgs =
      ngb_gas_rho_phys * feedback_props->rho_to_n_cgs;
  const float enrichment_weight_sum =
      sp->feedback_data.to_collect.enrichment_weight_sum;

  /* Check if there are neighbours, otherwise exit */
  if (ngb_gas_mass == 0.f || sp->density.wcount * pow_dimension(sp->h) < 1e-4) {
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (enrichment_weight_sum < 0.)
    error("Enrichment weight sum must be positive, not %g!",
          enrichment_weight_sum);
#endif

  /* Now we start filling the data structure for information to apply to the
   * particles. Do _NOT_ read from the to_collect substructure any more. */

  /* Zero all the output fields */
  feedback_reset_feedback(sp, feedback_props);

  /* Store the inverse sum of all enrichment weights as normalisation factor */
  sp->feedback_data.to_distribute.enrichment_normalisation =
      (enrichment_weight_sum > 0) ? 1.f / enrichment_weight_sum : 0.f;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.enrichment_normalisation < 0.)
    error("Enrichment normalisation must be positive, not %g",
          sp->feedback_data.to_distribute.enrichment_normalisation);
#endif

  /* Calculate mass of stars that has died from the star's birth up to the
   * beginning and end of timestep */
  const double max_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr, Z, feedback_props);
  const double min_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr + dt_Gyr, Z, feedback_props);

#ifdef SWIFT_DEBUG_CHECKS
  /* Sanity check. Worth investigating if necessary as functions for evaluating
   * mass of stars dying might be strictly decreasing.  */
  if (min_dying_mass_Msun > max_dying_mass_Msun)
    error("min dying mass is greater than max dying mass");
#endif

  /* Compute properties of the stochastic SNII feedback model. */
  if (feedback_props->with_SNII_feedback) {
    const double h_pkpc_inv = phys_const->const_parsec * 1e3 * cosmo->a_inv /
        sp->h;
    const double G_Newton = phys_const->const_newton_G;
    compute_SNII_feedback(sp, age, dt, ngb_gas_N, ngb_gas_mass,
                          ngb_gas_phys_nH_cgs, ngb_gas_Z, feedback_props,
                          min_dying_mass_Msun, max_dying_mass_Msun,
                          ti_begin, h_pkpc_inv, G_Newton,
                          ngb_gas_SFR, ngb_gas_rho_phys);
  }

  /* Integration interval is zero - this can happen if minimum and maximum
   * dying masses are above imf_max_mass_Msun. Return without doing any
   * enrichment. */
  if (min_dying_mass_Msun == max_dying_mass_Msun) return;

  /* Life is better in log */
  const double log10_max_dying_mass_Msun = log10(max_dying_mass_Msun);
  const double log10_min_dying_mass_Msun = log10(min_dying_mass_Msun);

  /* Compute elements, energy and momentum to distribute from the
   *  three channels SNIa, SNII, AGB */
  if (feedback_props->with_SNIa_enrichment) {
    evolve_SNIa(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
                feedback_props, star_age_Gyr, dt_Gyr, &sp->feedback_data);
  }
  if (feedback_props->with_SNII_enrichment) {
    evolve_SNII(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
                abundances, feedback_props, &sp->feedback_data);
  }
  if (feedback_props->with_AGB_enrichment) {
    evolve_AGB(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
               abundances, feedback_props, &sp->feedback_data);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.mass != 0.f)
    error("Injected mass will be lost");
#endif

  /* Compute the total mass to distribute (H + He  metals) */
  sp->feedback_data.to_distribute.mass =
      sp->feedback_data.to_distribute.total_metal_mass +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_H] +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_He];

  /* Compute energy change due to kinetic energy of ejectas */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass *
      feedback_props->AGB_ejecta_specific_kinetic_energy;

  /* Compute energy change due to kinetic energy of the star */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass * 0.5f *
      (sp->v[0] * sp->v[0] + sp->v[1] * sp->v[1] + sp->v[2] * sp->v[2]) *
      cosmo->a2_inv;

  TIMER_TOC(timer_do_star_evol);
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
void feedback_props_init(struct feedback_props* fp,
                         const struct phys_const* phys_const,
                         const struct unit_system* us,
                         struct swift_params* params,
                         const struct hydro_props* hydro_props,
                         const struct cosmology* cosmo) {

  /* Main operation modes ------------------------------------------------- */

  fp->with_SNII_feedback =
      parser_get_param_int(params, "EAGLEFeedback:use_SNII_feedback");

  fp->with_SNIa_feedback =
      parser_get_param_int(params, "EAGLEFeedback:use_SNIa_feedback");

  fp->with_AGB_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_AGB_enrichment");

  fp->with_SNII_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_SNII_enrichment");

  fp->with_SNIa_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_SNIa_enrichment");

  if (fp->with_SNIa_feedback && !fp->with_SNIa_enrichment) {
    error("Cannot run with SNIa feedback without SNIa enrichment.");
  }

  /* Properties of the IMF model ------------------------------------------ */

  /* Minimal and maximal mass considered */
  fp->imf_max_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:IMF_max_mass_Msun");
  fp->imf_min_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:IMF_min_mass_Msun");

  /* Check that it makes sense. */
  if (fp->imf_max_mass_msun < fp->imf_min_mass_msun) {
    error("Can't have the max IMF mass smaller than the min IMF mass!");
  }

  fp->log10_imf_max_mass_msun = log10(fp->imf_max_mass_msun);
  fp->log10_imf_min_mass_msun = log10(fp->imf_min_mass_msun);

  /* Properties of the SNII energy feedback model ------------------------- */

  char model[64];
  parser_get_param_string(params, "EAGLEFeedback:SNII_feedback_model", model);
  if (strcmp(model, "Random") == 0)
    fp->feedback_model = SNII_random_ngb_model;
  else if (strcmp(model, "Isotropic") == 0)
    fp->feedback_model = SNII_isotropic_model;
  else if (strcmp(model, "MinimumDistance") == 0)
    fp->feedback_model = SNII_minimum_distance_model;
  else if (strcmp(model, "MinimumDensity") == 0)
    fp->feedback_model = SNII_minimum_density_model;
  else
    error(
        "The SNII feedback model must be either 'Random', 'MinimumDistance', "
        "'MinimumDensity' or 'Isotropic', not %s",
        model);

  /* Are we sampling the SNII lifetimes for feedback or using a fixed delay? */
  fp->SNII_sampled_delay =
      parser_get_param_int(params, "EAGLEFeedback:SNII_sampled_delay");

  if (!fp->SNII_sampled_delay) {

    /* Set the delay time before SNII occur */
    fp->SNII_wind_delay =
        parser_get_param_double(params, "EAGLEFeedback:SNII_wind_delay_Gyr") *
        phys_const->const_year * 1e9;
  }

  /* Energy released by supernova type II */
  fp->E_SNII_cgs =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_erg");
  fp->E_SNII =
      fp->E_SNII_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Stellar mass limits for SNII feedback */
  const double SNII_min_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:SNII_min_mass_Msun");
  const double SNII_max_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:SNII_max_mass_Msun");

  /* Check that it makes sense. */
  if (SNII_max_mass_msun < SNII_min_mass_msun) {
    error("Can't have the max SNII mass smaller than the min SNII mass!");
  }

  fp->SNII_min_mass_msun = SNII_min_mass_msun;
  fp->SNII_max_mass_msun = SNII_max_mass_msun;
  fp->log10_SNII_min_mass_msun = log10(SNII_min_mass_msun);
  fp->log10_SNII_max_mass_msun = log10(SNII_max_mass_msun);

  /* Heating temperature parameters */
  fp->SNII_use_variable_delta_T =
      parser_get_param_int(params, "EAGLEFeedback:SNII_use_variable_delta_T");
  if (fp->SNII_use_variable_delta_T) {
    fp->SNII_use_instantaneous_density_for_dT =
        parser_get_param_int(
            params, "EAGLEFeedback:SNII_use_instantaneous_density_for_dT");
    fp->SNII_use_instantaneous_Z_for_dT =
        parser_get_param_int(
          params, "EAGLEFeedback:SNII_use_instantaneous_metallicity_for_dT");
    fp->SNII_delta_T_num_ngb_to_heat =
        parser_get_param_double(
          params, "EAGLEFeedback:SNII_delta_T_num_ngb_to_heat");
    fp->SNII_T_crit_factor =
        parser_get_param_double(params, "EAGLEFeedback:SNII_T_crit_factor");

    fp->SNII_delta_T_max =
        parser_get_param_double(params, "EAGLEFeedback:SNII_delta_T_max") /
        units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

    fp->SNII_with_energy_compensation =
        parser_get_param_int(
            params, "EAGLEFeedback:SNII_with_energy_compensation");
    if (fp->SNII_with_energy_compensation) {
      fp->SNII_delta_T_omega_max = parser_get_param_double(
        params, "EAGLEFeedback:SNII_delta_T_omega_max");
      fp->SNII_efficiency_zeta =
          parser_get_param_double(params, "EAGLEFeedback:SNII_efficiency_zeta");
    }

    fp->SNII_gamma_star =
        parser_get_param_double(params, "EAGLEFeedback:SNII_gamma") /
            (phys_const->const_year * 1e7 * 4.0 / 3.0 * M_PI *
            phys_const->const_parsec * phys_const->const_parsec *
            phys_const->const_parsec * 1e9);
    fp->SNII_sampling_reduction_within_smoothing_length =
        parser_get_param_int(params,
          "EAGLEFeedback:SNII_sampling_reduction_within_smoothing_length");
    fp->SNII_with_nu_below_one =
        parser_get_param_int(params, "EAGLEFeedback:SNII_with_nu_below_one");

    char temp[32];
    parser_get_param_string(
        params, "EAGLEFeedback:SNII_with_oversampling_timescale", temp);
    if (strcmp(temp, "None") == 0)
      fp->SNII_with_oversampling_timescale = eagle_SNII_timescale_none;
    else if (strcmp(temp, "GasConsumption") == 0)
      fp->SNII_with_oversampling_timescale = eagle_SNII_timescale_gasconsum;
    else if (strcmp(temp, "FreeFall") == 0)
      fp->SNII_with_oversampling_timescale = eagle_SNII_timescale_freefall;
    else
      error("Invalid value of "
            "EAGLEFeedback:SNII_with_oversampling_timescale: '%s'",
            temp);

    fp->SNII_dTcrit_floor =
        parser_get_param_float(params, "EAGLEFeedback:SNII_dTcrit_floor");
    fp->SNII_dTcrit_exp_nH =
        parser_get_param_float(params, "EAGLEFeedback:SNII_dTcrit_exp_nH");
    fp->SNII_dTcrit_exp_Z =
        parser_get_param_float(params, "EAGLEFeedback:SNII_dTcrit_exp_Z");
    fp->SNII_dTcrit_norm =
        parser_get_param_float(params, "EAGLEFeedback:SNII_dTcrit_norm");
    fp->SNII_dTlimit_exp_nH =
        parser_get_param_float(params, "EAGLEFeedback:SNII_dTlimit_exp_nH");
    fp->SNII_dTlimit_exp_Z =
        parser_get_param_float(params, "EAGLEFeedback:SNII_dTlimit_exp_Z");
    fp->SNII_dT_Zmin =
        parser_get_param_float(params, "EAGLEFeedback:SNII_dT_Zmin");
    fp->SNII_dTlimit_norm =
        parser_get_param_float(params, "EAGLEFeedback:SNII_dTlimit_norm");
    fp->Z_Sol = parser_get_param_float(params, "EAGLEFeedback:Z_Sol");

    fp->SNII_efficiency_eta_min = parser_get_param_float(
        params, "EAGLEFeedback:SNII_efficiency_eta_min");
    fp->SNII_maximum_nu = parser_get_param_float(
        params, "EAGLEFeedback:SNII_maximum_nu");
    fp->SNII_omega_by_sampling = parser_get_param_int(
       params, "EAGLEFeedback:SNII_omega_by_sampling");

  } else {    /* Constant dT model */
    fp->SNe_deltaT_desired =
        parser_get_param_float(params, "EAGLEFeedback:SNII_delta_T_K") /
        units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  }

  /* Properties of the energy fraction model */
  fp->f_E_min =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_min");
  fp->f_E_max =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_max");
  fp->Z_0 =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_Z_0");
  fp->n_0_cgs = parser_get_param_double(
      params, "EAGLEFeedback:SNII_energy_fraction_n_0_H_p_cm3");
  fp->n_n =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_n_n");
  fp->n_Z =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_n_Z");

  /* Check that it makes sense. */
  if (fp->f_E_max < fp->f_E_min) {
    error("Can't have the maximal energy fraction smaller than the minimal!");
  }

  /* Are we using the stars' birth properties or at feedback time? */
  fp->use_birth_density_for_f_th = parser_get_param_int(
      params, "EAGLEFeedback:SNII_energy_fraction_use_birth_density");
  fp->use_birth_Z_for_f_th = parser_get_param_int(
      params, "EAGLEFeedback:SNII_energy_fraction_use_birth_metallicity");

  /* Properties of the SNII enrichment model -------------------------------- */

  /* Set factors for each element adjusting SNII yield */
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    char buffer[50];
    sprintf(buffer, "EAGLEFeedback:SNII_yield_factor_%s",
            chemistry_get_element_name((enum chemistry_element)elem));

    fp->SNII_yield_factor[elem] =
        parser_get_opt_param_float(params, buffer, 1.f);
  }

  /* Properties of the SNIa enrichment model -------------------------------- */

  fp->SNIa_DTD_delay_Gyr =
      parser_get_param_double(params, "EAGLEFeedback:SNIa_DTD_delay_Gyr");

  char temp[32] = {0};
  parser_get_param_string(params, "EAGLEFeedback:SNIa_DTD", temp);

  if (strcmp(temp, "Exponential") == 0) {

    fp->SNIa_DTD = eagle_feedback_SNIa_DTD_exponential;

    /* Read SNIa exponential DTD model parameters */
    fp->SNIa_DTD_exp_norm = parser_get_param_float(
        params, "EAGLEFeedback:SNIa_DTD_exp_norm_p_Msun");
    fp->SNIa_DTD_exp_timescale_Gyr = parser_get_param_float(
        params, "EAGLEFeedback:SNIa_DTD_exp_timescale_Gyr");
    fp->SNIa_DTD_exp_timescale_Gyr_inv = 1.f / fp->SNIa_DTD_exp_timescale_Gyr;

  } else if (strcmp(temp, "PowerLaw") == 0) {

    fp->SNIa_DTD = eagle_feedback_SNIa_DTD_power_law;

    /* Read SNIa power-law DTD model parameters */
    fp->SNIa_DTD_power_law_norm = parser_get_param_float(
        params, "EAGLEFeedback:SNIa_DTD_power_law_norm_p_Msun");

    /* Renormalize everything such that the integral converges to
       'SNIa_DTD_power_law_norm' over 13.6 Gyr. */
    fp->SNIa_DTD_power_law_norm /= log(13.6 / fp->SNIa_DTD_delay_Gyr);

  } else {
    error("Invalid SNIa DTD model: '%s'", temp);
  }

  /* Energy released by supernova type Ia */
  fp->E_SNIa_cgs =
      parser_get_param_double(params, "EAGLEFeedback:SNIa_energy_erg");
  fp->E_SNIa =
      fp->E_SNIa_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Properties of the SNIa enrichment model -------------------------------- */

  /* Read AGB ejecta velocity */
  const float ejecta_velocity_km_p_s = parser_get_param_float(
      params, "EAGLEFeedback:AGB_ejecta_velocity_km_p_s");

  /* Convert to internal units */
  const float ejecta_velocity_cgs = ejecta_velocity_km_p_s * 1e5;
  const float ejecta_velocity =
      ejecta_velocity_cgs / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  /* Convert to specific thermal energy */
  fp->AGB_ejecta_specific_kinetic_energy =
      0.5f * ejecta_velocity * ejecta_velocity;

  /* Properties of the enrichment down-sampling ----------------------------- */

  fp->stellar_evolution_age_cut =
      parser_get_param_double(params,
                              "EAGLEFeedback:stellar_evolution_age_cut_Gyr") *
      phys_const->const_year * 1e9;

  fp->stellar_evolution_sampling_rate = parser_get_param_double(
      params, "EAGLEFeedback:stellar_evolution_sampling_rate");

  if (fp->stellar_evolution_sampling_rate < 1 ||
      fp->stellar_evolution_sampling_rate >= (1 << (8 * sizeof(char) - 1)))
    error("Stellar evolution sampling rate too large. Must be >0 and <%d",
          (1 << (8 * sizeof(char) - 1)));

  /* Check that we are not downsampling before reaching the SNII delay */
  if (!fp->SNII_sampled_delay &&
      fp->stellar_evolution_age_cut < fp->SNII_wind_delay)
    error(
        "Time at which the enrichment downsampling starts is lower than the "
        "SNII wind delay!");

  /* Gather common conversion factors --------------------------------------- */

  /* Calculate internal mass to solar mass conversion factor */
  const double Msun_cgs = phys_const->const_solar_mass *
                          units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  fp->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
  fp->solar_mass_to_mass = 1. / fp->mass_to_solar_mass;

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  fp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* Calculate conversion factor from rho to n_H
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  fp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Initialise the IMF ------------------------------------------------- */

  init_imf(fp);

  /* Calculate number of type II SN per unit solar mass based on our choice
   * of IMF and integration limits for type II SNe.
   * Note: No weighting by yields here. */
  fp->num_SNII_per_msun =
      integrate_imf(fp->log10_SNII_min_mass_msun, fp->log10_SNII_max_mass_msun,
                    eagle_imf_integration_no_weight,
                    /*(stellar_yields=)*/ NULL, fp);

  /* Initialise the yields ---------------------------------------------- */

  /* Read yield table filepath  */
  parser_get_param_string(params, "EAGLEFeedback:filename",
                          fp->yield_table_path);

  /* Allocate yield tables  */
  allocate_yield_tables(fp);

  /* Read the tables  */
  read_yield_tables(fp);

  /* Set yield_mass_bins array */
  const float imf_log10_mass_bin_size =
      (fp->log10_imf_max_mass_msun - fp->log10_imf_min_mass_msun) /
      (eagle_feedback_N_imf_bins - 1);

  for (int i = 0; i < eagle_feedback_N_imf_bins; i++)
    fp->yield_mass_bins[i] =
        imf_log10_mass_bin_size * i + fp->log10_imf_min_mass_msun;

  /* Resample yields from mass bins used in tables to mass bins used in IMF  */
  compute_yields(fp);

  /* Resample ejecta contribution to enrichment from mass bins used in tables to
   * mass bins used in IMF  */
  compute_ejecta(fp);

  message("initialized stellar feedback");
}

/**
 * @brief Clean-up the memory allocated for the feedback routines
 *
 * We simply free all the arrays.
 *
 * @param fp the feedback data structure.
 */
void feedback_clean(struct feedback_props* fp) {

  swift_free("imf-tables", fp->imf);
  swift_free("imf-tables", fp->imf_mass_bin);
  swift_free("imf-tables", fp->imf_mass_bin_log10);
  swift_free("feedback-tables", fp->yields_SNIa);
  swift_free("feedback-tables", fp->yield_SNIa_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.mass);
  swift_free("feedback-tables", fp->yield_AGB.metallicity);
  swift_free("feedback-tables", fp->yield_AGB.yield);
  swift_free("feedback-tables", fp->yield_AGB.yield_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.ejecta);
  swift_free("feedback-tables", fp->yield_AGB.ejecta_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.total_metals);
  swift_free("feedback-tables", fp->yield_AGB.total_metals_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.mass);
  swift_free("feedback-tables", fp->yield_SNII.metallicity);
  swift_free("feedback-tables", fp->yield_SNII.yield);
  swift_free("feedback-tables", fp->yield_SNII.yield_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.ejecta);
  swift_free("feedback-tables", fp->yield_SNII.ejecta_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.total_metals);
  swift_free("feedback-tables", fp->yield_SNII.total_metals_IMF_resampled);
  swift_free("feedback-tables", fp->lifetimes.mass);
  swift_free("feedback-tables", fp->lifetimes.metallicity);
  swift_free("feedback-tables", fp->yield_mass_bins);
  for (int i = 0; i < eagle_feedback_lifetime_N_metals; i++) {
    free(fp->lifetimes.dyingtime[i]);
  }
  free(fp->lifetimes.dyingtime);
  for (int i = 0; i < eagle_feedback_SNIa_N_elements; i++) {
    free(fp->SNIa_element_names[i]);
  }
  free(fp->SNIa_element_names);
  for (int i = 0; i < eagle_feedback_SNII_N_elements; i++) {
    free(fp->SNII_element_names[i]);
  }
  free(fp->SNII_element_names);
  for (int i = 0; i < eagle_feedback_AGB_N_elements; i++) {
    free(fp->AGB_element_names[i]);
  }
  free(fp->AGB_element_names);
}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream) {

  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the feedback routines. Helps debugging. */
  struct feedback_props feedback_copy = *feedback;

  /* zero AGB and SNII table pointers */
  zero_yield_table_pointers(&feedback_copy.yield_AGB);
  zero_yield_table_pointers(&feedback_copy.yield_SNII);

  /* zero SNIa table pointers */
  feedback_copy.yield_SNIa_IMF_resampled = NULL;
  feedback_copy.yields_SNIa = NULL;
  feedback_copy.yield_SNIa_total_metals_IMF_resampled = 0;

  /* zero element name tables */
  feedback_copy.SNIa_element_names = NULL;
  feedback_copy.SNII_element_names = NULL;
  feedback_copy.AGB_element_names = NULL;

  /* zero mass bins table */
  feedback_copy.yield_mass_bins = NULL;

  /* zero lifetime tracks */
  feedback_copy.lifetimes.mass = NULL;
  feedback_copy.lifetimes.metallicity = NULL;
  feedback_copy.lifetimes.dyingtime = NULL;

  /* zero IMF tables */
  feedback_copy.imf = NULL;
  feedback_copy.imf_mass_bin = NULL;
  feedback_copy.imf_mass_bin_log10 = NULL;

  restart_write_blocks((void*)&feedback_copy, sizeof(struct feedback_props), 1,
                       stream, "feedback", "feedback function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Read the structure from the stream and restore the feedback tables by
 * re-reading them.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream) {
  restart_read_blocks((void*)feedback, sizeof(struct feedback_props), 1, stream,
                      NULL, "feedback function");

  if (strlen(feedback->yield_table_path) != 0)
    feedback_restore_tables(feedback);
}
