"""
Heat Exchanger Fouling Simulation — Shell & Tube
=================================================
Model: Arrhenius-based deposition + shear-driven removal (Epstein, 1993)
Exchanger: Shell & Tube, single-pass, counter-flow
Variable: Fluid inlet temperature (T_in)

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
# 4. CLEAN HEAT TRANSFER COEFFICIENT — Dittus-Boelter correlation
# ---------------------------------------------------------------------------
u_tube = m_dot / (rho * A_flow)          # [m/s] flow velocity
Re     = rho * u_tube * D_i / mu         # [-]   Reynolds number
Nu     = 0.023 * Re**0.8 * Pr**0.4       # [-]   Nusselt (heating, n=0.4)
h_i    = Nu * k_f / D_i                  # [W/m²·K] inner HTC
U_clean = h_i                            # Simplified: tube-side controlling


def wall_shear_stress(Re, rho, u):
    """
    Blasius friction factor for turbulent flow in smooth tubes.
    f = 0.316 * Re^(-0.25)  valid for 4000 < Re < 10^5
    τ_w = f/8 * ρ * u²
    """
    f = 0.316 * Re**(-0.25)
    return (f / 8.0) * rho * u**2


# ---------------------------------------------------------------------------
# 5. ODE — dRf/dt = deposition_rate - removal_rate
# ---------------------------------------------------------------------------
def dRf_dt(t, Rf, T_K, tau_w):
    """
    Fouling resistance ODE (Epstein model):
        dRf/dt = A_dep * exp(-E_dep / RT) - A_rem * tau_w * Rf
    """
    deposition = A_dep * np.exp(-E_dep / (R_gas * T_K))
    removal    = A_rem * tau_w * Rf[0]
    return [deposition - removal]


# ---------------------------------------------------------------------------
# 6. SIMULATION — sweep over T_in values
# ---------------------------------------------------------------------------
T_in_values_C = np.arange(50, 130, 5)      # [°C] 50 → 125 °C (16 scenarios)
t_start = 0.0
t_end   = 8760.0 * 3600.0                  # [s]  1 year
n_points = 8760                             # hourly resolution
t_eval  = np.linspace(t_start, t_end, n_points)

records = []

print(f"{'='*60}")
print(f"  Shell & Tube Fouling Simulation — Epstein (1993) model")
print(f"  Re = {Re:.0f}  |  u = {u_tube:.2f} m/s  |  U_clean = {U_clean:.1f} W/m²K")
print(f"  Scenarios: {len(T_in_values_C)} temperature values")
print(f"{'='*60}\n")

for T_C in T_in_values_C:
    T_K   = T_C + 273.15
    tau_w = wall_shear_stress(Re, rho, u_tube)

    sol = solve_ivp(
        fun      = dRf_dt,
        t_span   = (t_start, t_end),
        y0       = [0.0],
        t_eval   = t_eval,
        args     = (T_K, tau_w),
        method   = "RK45",
        max_step = 3600.0
    )

    Rf_arr = np.clip(sol.y[0], 0.0, Rf_max)    # physical ceiling

    # Derived quantities
    U_arr      = 1.0 / (1.0 / U_clean + Rf_arr)          # [W/m²K]  overall HTC
    dT_lm      = 40.0                                      # [K]  assumed constant LMTD
    Q_arr      = U_arr * A_ht * dT_lm                     # [W]  heat duty
    Q_clean    = U_clean * A_ht * dT_lm
    eff_arr    = Q_arr / Q_clean                           # [-]  thermal efficiency
    # Pressure drop (Darcy-Weisbach + fouling narrows effective diameter)
    f_darcy    = 0.316 * Re**(-0.25)
    dP_arr     = f_darcy * (L / D_i) * 0.5 * rho * u_tube**2  # Pa (constant D approx)

    hours_arr  = sol.t / 3600.0

    for i in range(len(hours_arr)):
        records.append({
            "T_in_C"              : round(float(T_C), 1),
            "T_in_K"              : round(float(T_K), 2),
            "time_h"              : round(float(hours_arr[i]), 2),
            "Re"                  : round(float(Re), 1),
            "u_m_s"               : round(float(u_tube), 4),
            "tau_w_Pa"            : round(float(tau_w), 4),
            "Rf_m2K_W"            : round(float(Rf_arr[i]), 8),
            "U_overall_W_m2K"     : round(float(U_arr[i]), 4),
            "U_clean_W_m2K"       : round(float(U_clean), 4),
            "Q_W"                 : round(float(Q_arr[i]), 2),
            "Q_clean_W"           : round(float(Q_clean), 2),
            "thermal_efficiency"  : round(float(eff_arr[i]), 6),
            "dP_Pa"               : round(float(dP_arr), 2),
            "fouling_factor_TEMA" : "L" if Rf_arr[i] < 1.76e-4 else "H",
        })

    Rf_final = Rf_arr[-1]
    print(f"  T_in={T_C:5.1f}°C | Rf_final={Rf_final:.4e} m²K/W "
          f"| U_final={1/(1/U_clean+Rf_final):.1f} W/m²K "
          f"| η={1/(1+U_clean*Rf_final):.3f}")

# ---------------------------------------------------------------------------
# 7. EXPORT
# ---------------------------------------------------------------------------
os.makedirs("output", exist_ok=True)
df = pd.DataFrame(records)
out_path = "output/shell_tube_fouling_dataset.csv"
df.to_csv(out_path, index=False)

print(f"\n{'='*60}")
print(f"  Dataset saved → {out_path}")
print(f"  Rows: {len(df):,}  |  Columns: {len(df.columns)}")
print(f"  Temperature scenarios: {df['T_in_C'].nunique()}")
print(f"  Time span: {df['time_h'].min():.0f} – {df['time_h'].max():.0f} h")
print(f"{'='*60}")
