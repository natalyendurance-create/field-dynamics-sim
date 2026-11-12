[repository_documentation_github_zenodo.md](https://github.com/user-attachments/files/34763149/repository_documentation_github_zenodo.md)
# Field Dynamics (Unified Elastic Field - UEF): Numerical Verification Suite

**Author:** Nataliia Zotalieva (`nataly.endurance@gmail.com`)  
**Parent Monograph (Document #0):** [10.5281/zenodo.23857436](https://doi.org/10.5281/zenodo.23857436)

## Overview

This repository provides a self-contained, reproducible Python verification suite (`field_dynamics_sim.py`) accompanying the **Field Dynamics (Unified Elastic Field)** preprint series (**Document #0** and **Papers 1–6**).

It computes all dimensionless geometric invariants, cavitation boundary limits, and topological mass/coupling ratios directly from two primary mechanical parameters of the three-dimensional Euclidean ($\mathbb{R}^3$) micropolar Cosserat substrate:
* Reference unperturbed density: $\rho_0 \approx 4.30 \times 10^{18}\text{ kg/m}^3$
* Transverse elastic shear modulus: $\mu_0 = \rho_0 c^2 \approx 3.865 \times 10^{35}\text{ Pa}$

## Quick Start

### Requirements
* Python 3.9+
* `numpy` (`pip install numpy`)

### Running the Verification Suite
```bash
python field_dynamics_sim.py
```

## Summary of Verified Predictions

| Observable / Constant | UEF Analytical Expression | UEF Value | Empirical / CODATA | Precision |
| :--- | :--- | :--- | :--- | :--- |
| **Proton-Electron Mass Ratio ($m_p/m_e$)** | $6\pi^5$ | `1836.118` | `1836.153` | `0.002%` |
| **Fine-Structure Constant ($\alpha^{-1}$)** | $\frac{3}{4\Lambda_{\text{eff}}}\frac{m_p}{m_e}$ | `137.0304–137.0330` | `137.0360` | `0.002%` |
| **Proton Landé Factor ($g_p$)** | $\frac{16}{3}\left(1 + \frac{\alpha_p}{2\pi}\right)$ | `5.5899` | `5.5857` | `0.075%` |
| **Neutron Magnetic Ratio ($\mu_n/\mu_p$)** | $Y$-node phase inversion | `-0.6667` | `-0.6850` | `2.67%` |
| **QCD String Tension ($\sigma$)** | $\mu_0 \frac{\pi R_p^2}{6}$ | `0.8942 GeV/fm` | `~0.89 GeV/fm` | `0.5%` |
| **Deuteron Binding ($E_B(\mathrm{d})$)** | $\frac{\alpha_p \hbar\omega_p}{12\pi}$ | `2.228 MeV` | `2.2246 MeV` | `0.15%` |
| **Alpha Binding ($E_B(^4\text{He})$)** | $\frac{\hbar\omega_p}{\pi^2}$ | `28.15 MeV` | `28.296 MeV` | `0.50%` |
| **Lead-208 Neutron Excess ($N/Z$)** | $\frac{1 + \alpha_p}{1 - \alpha_p/2}$ | `1.5343` ($N=126$) | `1.5366` | `0.16%` |
| **Neutron Lifetime (Bottle UCN)** | $\tau_0 e^{S_E^{(0)} + \Delta S_{\text{env}}}$ | `878.8 s` | `878.4 ± 0.5 s` | `0.05%` |
| **Neutron Lifetime (Beam NIST)** | $\tau_{\text{bottle}}(1 + \delta_{\text{kin}})$ | `888.1 s` | `887.7 ± 2.2 s` | `0.04%` |
| **Axial Vector Coupling ($g_A$)** | $4/\pi$ | `1.27324` | `1.2754` | `0.17%` |
| **Newton's Constant ($G$)** | $\frac{P_0 V_p R_p}{m_p^2}$ | `6.669e-11` | `6.6743e-11` | `0.07%` |

## License

This software and documentation are released under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.
