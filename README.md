# Aero Solver

Physics-based aircraft sizing and performance estimation tool for fixed-wing UAV conceptual design.

## Overview

This project sizes a full fixed-wing UAV configuration from first principles by combining:
- aerodynamic theory,
- component-level empirical drag models, and
- iterative numerical convergence for cruise-point estimation.

The solver generates a flyable configuration by coupling wing, fuselage, and tail geometry with aerodynamic and stability constraints.

## What the model includes

- **Airfoil modeling**  
  Built-in wing airfoils:
  - `S1223` (default)
  - `Clark-Y`
  
  Built-in tail airfoils:
  - `NACA0012` for horizontal tail
  - `NACA0012` for vertical tail

  The model uses embedded parameters such as:
  - maximum lift coefficient (`C_Lmax`)
  - stall angle
  - zero-lift angle of attack
  - thickness and max-thickness location terms used in drag form factors

- **Geometry synthesis**  
  Computes wing, fuselage, horizontal tail, and vertical tail geometry using:
  - span, root/tip chord, taper ratio, sweep
  - mean aerodynamic chord (MAC)
  - aspect ratio
  - tail volume coefficient methods (`Vol_ht`, `Vol_vt`)
  - tail moment arms and areas

- **Aerodynamic force model**  
  Lift and drag are computed using standard aerodynamic relations.  
  Total drag is decomposed into:
  - **parasitic drag (`C_D0`)** from component-level wetted area, skin friction, form factor, and interference correction
  - **induced drag (`C_Dind`)** via aspect ratio and Oswald efficiency

- **Iterative cruise solver**  
  Cruise velocity is solved iteratively to converge to a minimum-drag operating condition for the specified aircraft setup.
  
  The solver includes:
  - Oswald efficiency modeling
  - Reynolds number effects
  - Mach-dependent corrections

- **Second-order effects and controllability context**  
  Includes downwash/effective angle adjustments and tail-volume-based sizing to maintain realistic stability/control geometry relationships while minimizing drag.

## Outputs

Running the script prints:
- Airfoil parameters and effective angles
- Full geometry for wing/HT/VT/fuselage
- Cruise speed, stall speed, Mach number
- Reynolds number and wing loading
- `C_D0`, `C_Dind`, total `C_D`, cruise `C_L`, Oswald efficiency
- Lift, drag, thrust margin, and lift-to-drag ratio

## Run

### 1) Install dependencies

```bash
pip install numpy matplotlib
```

### 2) Execute

```bash
python main.py
```

## Notes

- All calculations are in SI units unless otherwise stated by output labels.
- Current configuration values are hardcoded in `main.py` for rapid iteration during conceptual design.
- You can switch wing airfoils by changing the selected function call near the top of `main.py` (e.g., `s1223()` vs `clarky()`).
