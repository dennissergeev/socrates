"""
Calculate approximate pathlengths for a given .raw file.

Sample output:

Global Average H2-H2 Path: 9.4950e+05 kg^2/m^5
Global Average H2-He Path: 3.3274e+05 kg^2/m^5

Table of absorbers with computed annual global mmr values and total column paths

Gas    Mean MMR     Total column path (kg/m2)
H2     7.4042e-01   4.1737e+05
He     2.5947e-01   1.4626e+05
CH4    1.8531e-01   1.0446e+05
H2O    9.2701e-02   5.2255e+04
C2H6   2.3059e-09   1.2998e-03
C2H2   1.1778e-09   6.6390e-04
CO     1.2139e-09   6.8428e-04
"""

import pandas as pd

# Physical constants
g = 8.87  # m/s2
R = 3615.0  # J/kg/K

# Read data from file
file_path = "uranus_atm.raw"
df = pd.read_csv(file_path, sep=r"\s+", comment="*", nrows=29)
print(df.head())
print()

# Pressure differences between layers
df["delta_p"] = df["PRESS(PA)"].diff()
df["delta_p"] = df["delta_p"].fillna(df["PRESS(PA)"])

# Calculate layer mass and the total mass
df["LAYER_MASS(KG/M2)"] = df["delta_p"] / g

total_atmosphere_mass = df["LAYER_MASS(KG/M2)"].sum()
total_atmosphere_mass.item()

# Calculate density
df["DENSITY(KG/M3)"] = df["PRESS(PA)"] / (R * df["TEMP(K)"])

# Calculate continuum pathlengths
df["CONT_H2_H2"] = (
    df["H2(NONE)"] * df["H2(NONE)"] * df["LAYER_MASS(KG/M2)"] * df["DENSITY(KG/M3)"]
)
df["CONT_H2_He"] = (
    df["He(NONE)"] * df["H2(NONE)"] * df["LAYER_MASS(KG/M2)"] * df["DENSITY(KG/M3)"]
)

total_h2h2 = df["CONT_H2_H2"].sum()
total_h2he = df["CONT_H2_He"].sum()

print(f"Global Average H2-H2 Path: {total_h2h2:.4e} kg^2/m^5")
print(f"Global Average H2-He Path: {total_h2he:.4e} kg^2/m^5")

# Calculate pathlengths for normal absorbers (kg/m2) and their mean MMR
gases = [
    "H2",
    "He",
    "CH4",
    "H2O",
    "NH3",
    "H2S",
    "C2H6",
    "C2H2",
    "C2H4",
    "CO",
    "CO2",
    "PH3",
]

print(
    "\nTable of absorbers with computed annual global mmr values and total column paths"
)
print("\nGas    Mean MMR     Total column path (kg/m2)")
for gas in gases:
    column = f"{gas}(NONE)"
    if column in df.columns:
        total_path = (df[column] * df["LAYER_MASS(KG/M2)"]).sum()
        mean_mmr = total_path / total_atmosphere_mass
        print(f"{gas:6} {mean_mmr:10.4e}   {total_path:10.4e}")
