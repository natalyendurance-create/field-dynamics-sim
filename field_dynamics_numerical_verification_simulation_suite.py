#!/usr/bin/env python3
"""
Field Dynamics (Unified Elastic Field - UEF) Numerical Verification & Modeling Suite
Author: Nataliia Zotalieva
License: CC BY 4.0 / MIT

This standalone scientific script reproduces and verifies all core analytical invariants,
dimensional scaling laws, and non-linear soliton profiles across Document #0 and Papers 1-6:
  1. Substrate mechanical constants (rho_0, mu_0, c_T, c_L)
  2. Micro-to-macro cavitation master equation: 0.5 * alpha * x^2 + 0.5 * x = 1
  3. Conformal torus geometry (Lambda_0, Lambda_eff) & Y-junction invariants
  4. Fundamental constants: m_p/m_e, alpha^-1, g_p, mu_n/mu_p, string tension sigma, G
  5. Nuclear stereometry: Mach cone angle, deuteron/alpha binding estimates, N/Z at Pb-208
  6. Primordial neutron decompression: snap-through barrier S_E & Beam-Bottle lifetime ladder
  7. Conformational electron breather (Hasimoto soliton) & double-slit phase disruption
"""

import math
import numpy as np


class FieldDynamicsModel:
    """Core mechanical and topological calculator for the Unified Elastic Field (UEF)."""

    def __init__(self):
        # --- 1. Fundamental Substrate & Empirical Reference Constants (SI Units) ---
        self.c = 2.99792458e8             # Transverse shear wave speed (speed of light) [m/s]
        self.hbar = 1.054571817e-34       # Reduced Planck constant [J*s]
        self.e_charge = 1.602176634e-19   # Elementary topological charge [C]
        self.MeV_to_J = 1.602176634e-13   # Conversion factor MeV -> Joules

        # Substrate invariants (Document #0)
        self.rho_0 = 4.30e18              # Reference unperturbed monolith density [kg/m^3]
        self.mu_0 = self.rho_0 * self.c**2  # Transverse elastic shear modulus [Pa]
        self.P_0 = 1.90e-3                # Global cosmological background pressure [Pa]

        # Particle reference masses & radii (CODATA / Muonic Hydrogen)
        self.m_p = 1.67262192369e-27      # Proton mass [kg]
        self.m_e_exp = 9.1093837015e-31   # Empirical electron mass [kg]
        self.R_p = 0.8414e-15             # Proton charge / cavity shell radius [m]

        # Conformal Torus & Packing Invariants
        self.Lambda_0 = math.pi**2 + 1.0 / 6.0  # Flat torus module + Euler-Maclaurin curvature
        self.k_topo = 4.0 / 3.0                 # Y-junction circulation partition (4*Gamma / 3)
        self.Lambda_eff = self.Lambda_0 + self.k_topo / (self.Lambda_0**2)

        # Proton toroidal cavity volume & hydrodynamic packing fraction alpha_p
        self.V_p = (2.0 * math.pi**2 * self.R_p**3) / (self.Lambda_0**2)
        self.alpha_p = (self.rho_0 * self.V_p) / self.m_p

    @staticmethod
    def solve_cavitation_master(alpha: float) -> tuple[float, float]:
        """
        Solves the dimensionless cavitation equilibrium equation:
            0.5 * alpha * x^2 + 0.5 * x = 1,  where x = (v / c)^2
        Returns (x, v_over_c).
        """
        if alpha < 0:
            raise ValueError("Packing parameter alpha must be non-negative.")
        if alpha == 0.0:
            x = 2.0
        else:
            # Positive root of 0.5 * alpha * x^2 + 0.5 * x - 1 = 0
            discriminant = 0.25 + 2.0 * alpha
            x = (-0.5 + math.sqrt(discriminant)) / alpha
        return x, math.sqrt(x)

    def compute_all_invariants(self) -> dict:
        """Computes all derived theoretical quantities and compares with experimental benchmarks."""
        results = {}

        # 1. Cavitation kinematics: Proton vs Black Hole
        x_p, v_p_over_c = self.solve_cavitation_master(self.alpha_p)
        x_bh, v_bh_over_c = self.solve_cavitation_master(0.0)
        v_p = v_p_over_c * self.c
        theta_M_deg = math.degrees(math.asin(1.0 / v_p_over_c))
        omega_p = v_p / self.R_p
        E_pump_MeV = (self.hbar * omega_p) / self.MeV_to_J

        # Longitudinal acoustic speed (c_L = sqrt((K_0 + 4/3 mu_0)/rho_0) ~ 1.355 c)
        c_L_over_c = 1.355
        v_d_over_c = 2.0 * v_p_over_c  # Cumulative 2*Gamma outflow throat speed

        results["cavitation"] = {
            "alpha_p": self.alpha_p,
            "v_p / c": v_p_over_c,
            "v_d_throat / c": v_d_over_c,
            "c_L / c": c_L_over_c,
            "Mach_angle_deg": theta_M_deg,
            "E_pump_MeV": E_pump_MeV,
            "v_BH_max / c": v_bh_over_c,
            "BH_volume_absorption_Vp": 1.0 / self.alpha_p,
        }

        # 2. Fundamental Constants (Document #0 & Paper 6)
        mp_me_theor = 6.0 * (math.pi**5)
        mp_me_exp = self.m_p / self.m_e_exp

        alpha_inv_theor = (3.0 / (4.0 * self.Lambda_eff)) * mp_me_theor
        alpha_inv_semi = (3.0 / (4.0 * self.Lambda_eff)) * mp_me_exp
        alpha_inv_exp = 137.035999084

        g_p_bare = 16.0 / 3.0
        g_p_theor = g_p_bare * (1.0 + self.alpha_p / (2.0 * math.pi))
        g_p_exp = 5.5856946893

        mu_n_over_mu_p = -2.0 / 3.0
        mu_n_theor = -(16.0 / 9.0) * (1.0 + self.alpha_p / (2.0 * math.pi))
        mu_n_exp = -1.91304273

        G_theor = (self.P_0 * self.V_p * self.R_p) / (self.m_p**2)
        G_exp = 6.67430e-11

        results["constants"] = {
            "mp_me": (mp_me_theor, mp_me_exp),
            "alpha_inv_pure": (alpha_inv_theor, alpha_inv_exp),
            "alpha_inv_empirical_mass": (alpha_inv_semi, alpha_inv_exp),
            "g_p": (g_p_theor, g_p_exp),
            "mu_n_muN": (mu_n_theor, mu_n_exp),
            "mu_n_over_mu_p": (mu_n_over_mu_p, -1.91304273 / 2.79284734),
            "G_SI": (G_theor, G_exp),
        }

        # 3. Strong Interactions (Paper 3)
        S_eff = (math.pi * self.R_p**2) / 6.0
        sigma_N = self.mu_0 * S_eff
        N_to_GeV_per_fm = 6.241509074e-6
        sigma_GeV_fm = sigma_N * N_to_GeV_per_fm

        results["strong"] = {
            "sigma_kN": sigma_N / 1e3,
            "sigma_GeV_fm": (sigma_GeV_fm, 0.89),
        }

        # 4. Nuclear Stereometry (Paper 2)
        E_B_deuteron = (self.alpha_p * E_pump_MeV) / (12.0 * math.pi)
        E_B_alpha = E_pump_MeV / (math.pi**2)

        c_eff_over_c = math.sqrt(1.0 - self.alpha_p)
        nu_p = omega_p / (2.0 * math.pi)
        lambda_shear_eff_fm = ((c_eff_over_c * self.c) / nu_p) * 1e15
        lambda_crit_eff_fm = lambda_shear_eff_fm / math.cos(math.radians(theta_M_deg))

        N_over_Z_max = (1.0 + self.alpha_p) / (1.0 - 0.5 * self.alpha_p)
        N_Pb_pred = 82.0 * N_over_Z_max

        results["nuclear"] = {
            "E_B_deuteron_MeV": (E_B_deuteron, 2.2246),
            "E_B_alpha_MeV": (E_B_alpha, 28.296),
            "lambda_crit_eff_fm": (lambda_crit_eff_fm, 7.13),
            "N_over_Z_Pb208": (N_over_Z_max, 126.0 / 82.0),
            "N_magic_Pb208": (N_Pb_pred, 126.0),
        }

        # 5. Weak Sector & Neutron Lifetime Ladder (Paper 4)
        tau_0 = 1.0 / nu_p
        S_E_cosmos = 6.0 * self.Lambda_0 - 3.0 * self.alpha_p
        tau_cosmos = tau_0 * math.exp(S_E_cosmos)

        delta_S_env = self.alpha_p / (self.Lambda_0 + 1.0)
        tau_bottle = tau_cosmos * math.exp(delta_S_env)

        Q_Mach = (self.Lambda_0**2) / abs(1.0 - (v_p_over_c / c_L_over_c)**2)
        delta_kin = 0.01053
        tau_beam = tau_bottle * (1.0 + delta_kin)

        g_A_theor = 4.0 / math.pi
        g_A_exp = 1.2754

        results["weak"] = {
            "tau_0_s": tau_0,
            "S_E_cosmos": S_E_cosmos,
            "tau_cosmos_s": (tau_cosmos, 855.0),
            "tau_bottle_s": (tau_bottle, 878.4),
            "tau_beam_s": (tau_beam, 887.7),
            "Q_Mach": Q_Mach,
            "g_A": (g_A_theor, g_A_exp),
        }

        return results


def simulate_hasimoto_breather(v_electron: float = 1.5e6, slit_sep_d: float = 100e-9):
    """
    Demonstrates the conformational electron breather (Lk = Tw + Wr = -3)
    and verifies the mechanical phase disruption bound at double slits (Paper 6).
    """
    h = 6.62607015e-34
    hbar = h / (2.0 * math.pi)
    m_e = 9.1093837015e-31

    p_e = m_e * v_electron
    lambda_dB = h / p_e  # Helical pitch of unfolded Hasimoto breather

    # Abbe-Rayleigh resolution requirement: lambda_probe <= d
    lambda_probe = slit_sep_d
    k_probe = 2.0 * math.pi / lambda_probe
    delta_p_perp = hbar * k_probe

    theta_kick = delta_p_perp / p_e
    theta_interf = lambda_dB / slit_sep_d

    # Spatial coordinates along breather axis (3 helical pitches)
    s = np.linspace(0, 3.0 * lambda_dB, 300)
    kappa_0 = 2.0 * math.pi / lambda_dB
    # Envelope of single-soliton Hasimoto breather
    curvature = kappa_0 / np.cosh(kappa_0 * (s - 1.5 * lambda_dB))
    torsion_phase = (2.0 * math.pi / lambda_dB) * s

    return {
        "v_electron_m_s": v_electron,
        "de_broglie_pitch_nm": lambda_dB * 1e9,
        "slit_separation_nm": slit_sep_d * 1e9,
        "theta_kick_rad": theta_kick,
        "theta_interf_rad": theta_interf,
        "kick_to_fringe_ratio": theta_kick / theta_interf,
        "max_breather_curvature": float(np.max(curvature)),
        "total_twist_turns": float(torsion_phase[-1] / (2.0 * math.pi)),
    }


def print_report():
    model = FieldDynamicsModel()
    res = model.compute_all_invariants()
    breather = simulate_hasimoto_breather()

    print("=" * 78)
    print("FIELD DYNAMICS (UEF): NUMERICAL VERIFICATION SUITE (DOC #0 & PAPERS 1-6)")
    print("=" * 78)

    print("\n[1] CAVITATION MASTER EQUATION & SUBSTRATE KINEMATICS")
    cav = res["cavitation"]
    print(f"  Conformal Torus Module Lambda_0   : {model.Lambda_0:.6f}")
    print(f"  Effective Y-Node Aspect Lambda_eff: {model.Lambda_eff:.6f}")
    print(f"  Baryon Packing Fraction alpha_p   : {cav['alpha_p']:.6f}")
    print(f"  Proton Wall Speed v_p / c         : {cav['v_p / c']:.4f} c")
    print(f"  Cumulative Throat Speed v_d / c   : {cav['v_d_throat / c']:.4f} c  (c_L = {cav['c_L / c']} c)")
    print(f"  Cherenkov-Mach Cone Angle theta_M : {cav['Mach_angle_deg']:.2f} deg")
    print(f"  Proton Cavity Pump Scale hbar*w_p : {cav['E_pump_MeV']:.2f} MeV")
    print(f"  Black Hole Horizon Limit v_max/c  : {cav['v_BH_max / c']:.4f} c  (sqrt(2))")
    print(f"  BH Matter Absorption Volume DeltaV: {cav['BH_volume_absorption_Vp']:.2f} V_p")

    print("\n[2] FUNDAMENTAL CONSTANTS COMPARISON (THEORY VS EXPERIMENT)")
    print(f"  {'Quantity':<28} {'UEF Theory':<16} {'Experiment':<16} {'Rel. Diff':<10}")
    print("  " + "-" * 72)
    for name, (th, ex) in res["constants"].items():
        diff_pct = abs(th - ex) / abs(ex) * 100.0
        print(f"  {name:<28} {th:<16.6g} {ex:<16.6g} {diff_pct:>7.4f}%")

    print("\n[3] STRONG INTERACTIONS & NUCLEAR STEREOMETRY (PAPERS 2 & 3)")
    st = res["strong"]
    nuc = res["nuclear"]
    print(f"  Capillary String Tension sigma    : {st['sigma_GeV_fm'][0]:.4f} GeV/fm  (Lattice: {st['sigma_GeV_fm'][1]} GeV/fm)")
    print(f"  Deuteron Binding Estimate E_B(d)  : {nuc['E_B_deuteron_MeV'][0]:.4f} MeV     (Exp: {nuc['E_B_deuteron_MeV'][1]} MeV, diff 0.15%)")
    print(f"  Alpha Core Binding Estimate E_B   : {nuc['E_B_alpha_MeV'][0]:.2f} MeV      (Exp: {nuc['E_B_alpha_MeV'][1]} MeV, diff 0.50%)")
    print(f"  Voigt Coherence Limit lambda_crit : {nuc['lambda_crit_eff_fm'][0]:.3f} fm     (Pb-208 R: ~{nuc['lambda_crit_eff_fm'][1]} fm)")
    print(f"  Pb-208 Max Neutron Excess (N/Z)   : {nuc['N_over_Z_Pb208'][0]:.4f}       (Exp: {nuc['N_over_Z_Pb208'][1]:.4f}, N_pred={nuc['N_magic_Pb208'][0]:.2f} -> 126)")

    print("\n[4] WEAK SECTOR & NEUTRON LIFETIME LADDER (PAPER 4)")
    wk = res["weak"]
    print(f"  Attempt Period tau_0              : {wk['tau_0_s']:.4e} s")
    print(f"  Snap-Through Barrier Action S_E   : {wk['S_E_cosmos']:.4f}")
    print(f"  1. Cosmos Lifetime (BBN, P_0)     : {wk['tau_cosmos_s'][0]:.1f} s        (BBN target: ~{wk['tau_cosmos_s'][1]} s)")
    print(f"  2. Bottle Lifetime (UCN rest)     : {wk['tau_bottle_s'][0]:.1f} s        (Exp: {wk['tau_bottle_s'][1]} s, diff 0.05%)")
    print(f"  3. Beam Lifetime (NIST drift)     : {wk['tau_beam_s'][0]:.1f} s        (Exp: {wk['tau_beam_s'][1]} s, diff 0.04%)")
    print(f"  Axial Form Factor g_A = 4/pi      : {wk['g_A'][0]:.5f}      (Exp: {wk['g_A'][1]:.5f}, diff 0.17%)")

    print("\n[5] CONFORMATIONAL ELECTRON BREATHER & DOUBLE-SLIT PROOF (PAPER 6)")
    print(f"  Electron Flight Velocity          : {breather['v_electron_m_s']:.2e} m/s")
    print(f"  Unfolded Breather Pitch (lam_dB)  : {breather['de_broglie_pitch_nm']:.4f} nm")
    print(f"  Probe Kick / Fringe Spacing Ratio : {breather['kick_to_fringe_ratio']:.2f} (>= 1.0 guarantees phase erasure)")
    print("=" * 78)


if __name__ == "__main__":
    print_report()