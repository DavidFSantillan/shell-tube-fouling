"""
Heat Exchanger Fouling Simulation — Shell & Tube
=================================================
Model: Arrhenius-based deposition + shear-driven removal (Epstein, 1993)
Exchanger: Shell & Tube, single-pass, counter-flow
Variable: Fluid inlet temperature (T_in)

Scenario variables
------------------
- Fluid inlet temperature : 50 – 125 °C  (16 values, step 5 °C)
- Mass flow rate          : 3.0 – 7.0 kg/s (5 values)
- Total scenarios         : 80 (16 × 5)
- Total rows              : ~700,800 (80 × 8760 hourly timesteps)

References
----------
- Epstein, N. (1993). The fouling of heat exchangers. In: 14th International
  Heat Transfer Conference. Washington, D.C.
- Polley, G.T., Wilson, D.I., Yeap, B.L., Pugh, S.J. (2002). Evaluation of
  laboratory fouling data for application to crude oil preheat trains.
  Chemical Engineering Research and Design, 80(7), 713-727.
- Kern, D.Q., Seaton, R.E. (1959). A theoretical analysis of thermal surface
  fouling. British Chemical Engineering, 4(5), 258-262.
- Incropera, F.P., DeWitt, D.P. (2002). Fundamentals of Heat and Mass Transfer,
  5th ed. Wiley.

Author: generated via simulation pipeline
License: CC0-1.0
"""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import os
import itertools
# ---------------------------------------------------------------------------
# 1. GEOMETRY — Shell & Tube heat exchanger
# ---------------------------------------------------------------------------
D_i    = 0.019       # [m]   Tube inner diameter (3/4 in standard)
D_o    = 0.022       # [m]   Tube outer diameter
L      = 4.88        # [m]   Tube length (16 ft standard)
N_t    = 100         # [-]   Number of tubes
A_flow = N_t * np.pi * D_i**2 / 4   # [m²] Total flow cross-section
A_ht   = N_t * np.pi * D_i * L      # [m²] Total heat transfer area

# ---------------------------------------------------------------------------
# 2. FLUID PROPERTIES — water-based process fluid (tube side)
# ---------------------------------------------------------------------------
rho   = 980.0    # [kg/m³]  Density (≈ 80°C water)
mu    = 3.5e-4   # [Pa·s]   Dynamic viscosity
Pr    = 2.3      # [-]      Prandtl number
k_f   = 0.644    # [W/m·K]  Thermal conductivity
Cp    = 4190.0   # [J/kg·K] Specific heat capacity
m_dot = 5.0      # [kg/s]   Mass flow rate (tube side)

# ---------------------------------------------------------------------------
# 3. FOULING MODEL PARAMETERS — Arrhenius + shear removal
#    Based on Polley et al. (2002) crude-oil correlations, adapted for
#    general process fluid fouling
# ---------------------------------------------------------------------------
R_gas = 8.314      # [J/mol·K]  Universal gas constant

# Deposition (Arrhenius term)
A_dep  = 9.08e-5   # [m²·K/W·s]  Pre-exponential (calibrated, Polley 2002)
E_dep  = 50_000.0  # [J/mol]     Activation energy ~50 kJ/mol (Epstein 1993)

# Removal (shear term)
A_rem  = 7.3e-7    # [m²·K/W/(Pa·s)]  Removal coeff (τ≈3000h at typical τ_w)

# Asymptotic / physical limit
Rf_max = 5.0e-4    # [m²·K/W]   Maximum fouling resistance (physical ceiling)

# ---------------------------------------------------------------------------
# 4. SCENARIO GRID
# ---------------------------------------------------------------------------
T_in_values  = np.arange(50, 130, 5)              # 16 temperatures [C]
mdot_values  = np.array([3.0, 4.0, 5.0, 6.0, 7.0])  # 5 flow rates [kg/s]
scenarios    = list(itertools.product(T_in_values, mdot_values))  # 80 total

# ---------------------------------------------------------------------------
# 5. TIME GRID
# ---------------------------------------------------------------------------
t_end_h  = 8760.0
n_points = 8760
t_eval_h = np.linspace(0.0, t_end_h, n_points)
t_eval_s = t_eval_h * 3600.0

# ---------------------------------------------------------------------------
# 6. HELPER FUNCTIONS
# ---------------------------------------------------------------------------
def compute_htc(m_dot):
    u  = m_dot / (rho * A_flow)
    Re = rho * u * D_i / mu
    Nu = 0.023 * Re**0.8 * Pr**0.4
    h  = Nu * k_f / D_i
    return h, u, Re

def wall_shear_stress(Re, u):
    f = 0.316 * Re**(-0.25)
    return (f / 8.0) * rho * u**2

def dRf_dt(t, Rf, T_K, tau_w):
    dep = A_dep * np.exp(-E_dep / (R_gas * T_K))
    rem = A_rem * tau_w * Rf[0]
    return [dep - rem]

# ---------------------------------------------------------------------------
# 7. SIMULATION LOOP
# ---------------------------------------------------------------------------
print("=" * 65)
print("  Shell & Tube Fouling Simulation — Epstein (1993)")
print(f"  Temperatures  : {len(T_in_values)} values ({T_in_values[0]}–{T_in_values[-1]} C)")
print(f"  Flow rates    : {len(mdot_values)} values ({mdot_values[0]}–{mdot_values[-1]} kg/s)")
print(f"  Total scenarios: {len(scenarios)}")
print(f"  Expected rows  : ~{len(scenarios) * n_points:,}")
print("=" * 65 + "\n")

records = []

for i, (T_C, m_dot) in enumerate(scenarios):
    T_K    = T_C + 273.15
    U_clean, u, Re = compute_htc(m_dot)
    tau_w  = wall_shear_stress(Re, u)
    f_d    = 0.316 * Re**(-0.25)
    dP     = f_d * (L / D_i) * 0.5 * rho * u**2

    sol = solve_ivp(
        fun      = dRf_dt,
        t_span   = (t_eval_s[0], t_eval_s[-1]),
        y0       = [0.0],
        t_eval   = t_eval_s,
        args     = (T_K, tau_w),
        method   = "RK45",
        max_step = 3600.0,
    )

    Rf_arr    = np.clip(sol.y[0], 0.0, Rf_max)
    U_overall = 1.0 / (1.0 / U_clean + Rf_arr)
    dT_lm     = 40.0
    Q_arr     = U_overall * A_ht * dT_lm
    Q_clean   = U_clean   * A_ht * dT_lm
    eff_arr   = Q_arr / Q_clean

    for j in range(len(sol.t)):
        records.append({
            # Scenario identifiers
            "scenario_id"         : f"T{int(T_C)}_Q{m_dot:.1f}",
            "T_in_C"              : round(float(T_C), 1),
            "T_in_K"              : round(float(T_K), 2),
            "m_dot_nominal_kg_s"  : round(float(m_dot), 1),
            "time_h"              : round(float(sol.t[j] / 3600.0), 2),
            # Hydraulics
            "Re"                  : round(float(Re), 1),
            "u_m_s"               : round(float(u), 4),
            "tau_w_Pa"            : round(float(tau_w), 4),
            "dP_Pa"               : round(float(dP), 2),
            # Fouling
            "Rf_m2K_W"            : round(float(Rf_arr[j]), 8),
            # Thermal performance
            "U_clean_W_m2K"       : round(float(U_clean), 4),
            "U_overall_W_m2K"     : round(float(U_overall[j]), 4),
            "Q_W"                 : round(float(Q_arr[j]), 2),
            "Q_clean_W"           : round(float(Q_clean), 2),
            "thermal_efficiency"  : round(float(eff_arr[j]), 6),
            "fouling_factor_TEMA" : "L" if Rf_arr[j] < 1.76e-4 else "H",
        })

    if (i + 1) % 10 == 0 or i == 0:
        print(f"  [{i+1:2d}/{len(scenarios)}] "
              f"T={T_C:.0f}C  m_dot={m_dot:.1f}kg/s  "
              f"Re={Re:.0f}  Rf_final={Rf_arr[-1]:.3e}  "
              f"eta={eff_arr[-1]:.3f}")

# ---------------------------------------------------------------------------
# 8. EXPORT
# ---------------------------------------------------------------------------
os.makedirs("output", exist_ok=True)
df = pd.DataFrame(records)
out_path = "output/shell_tube_fouling_dataset.csv"
df.to_csv(out_path, index=False)

print(f"\n{'='*65}")
print(f"  Dataset saved  -> {out_path}")
print(f"  Rows           : {len(df):,}")
print(f"  Columns        : {len(df.columns)}")
print(f"  Scenarios      : {df['scenario_id'].nunique()}")
print(f"  Temp range     : {df['T_in_C'].min()} – {df['T_in_C'].max()} C")
print(f"  Flow range     : {df['m_dot_nominal_kg_s'].min()} – {df['m_dot_nominal_kg_s'].max()} kg/s")
print(f"{'='*65}")