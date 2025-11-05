# Maraak2026Imperatives
Supplementary software to "Imperatives for anaerobic high cell density culturing by denitrification: technical and physiological perspectives" by Maråk et al.

## Modeling pH, NO3-, and CO2 kinetics
A dynamic model, Modeling_pHstat.py, was developed in Python to simulate the kinetics of CO₂ production, pH variation, and nitrate (NO₃⁻) consumption in the bioreactor operated under pH-stat control. The model was designed to describe the coupled biochemical and physicochemical processes governing the denitrification-driven cultivation of Paracoccus pantotrophus.

The model was implemented in Python (v3.7) using the pandas library. Model parameters (e.g., specific growth rate, yield coefficients, equilibrium constants) were estimated from experimental data obtained in this study and in Maråk et. al 2024.

Simulations were performed using experimental conditions (temperature, gas flow rate, agitation speed, and initial concentrations) identical to those in the bioreactor experiments. The model output included the time courses of CO₂ concentration, pH, and NO₃⁻ concentration in the liquid phase.

## Medium Calculator
A medium calculator, Rx_MediumCalculator.xlsx, was designed in Excel (Excel version 2508) to calculate the steady-state biomass for a high cell density culture and the mass balance for all elements. The parameters used were estimated from experimental data obtained in this study and in Maråk et. al 2024.

## References
Maråk, M. M., Kellermann, R., Bergaust, L. L., & Bakken, L. R. (2024). High cell density cultivation by anaerobic respiration. Microbial Cell Factories, 23(1), 320.
