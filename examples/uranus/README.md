# Scripts for creating spectral files for Uranus

## Contributors

- Denis Sergeev (University of Bristol)
- Evie Cushing (University of Bristol)
- James Manners (Met Office)

## Key features
- Clear-sky
- 5 gases: CH4 (6), H2 (23), He (24), C2H2 (34), C2H6 (46)
- 2 CIA: H2-H2, H2-He
- High-resolution files: 
  - SW: 293 bands from 0.3 to 15 microns
  - LW: 319 bands from 3.25 to 1000 microns
- Low-resolution files:
  - SW: 26 bands
  - LW: 20 bands

## Files

| File | Description |
| --- | --- |
| `env_uranus` | Shared environment variables sourced by the `mk_*` scripts: data paths (HITRAN, LBL directories), band counts, column mass/path scaling factors, and P-T table filenames. |
| `calc_pathlengths.py` | Computes approximate gas column paths (kg/m2) and CIA continuum paths (kg2/m5) from `uranus_atm.raw`, used to set the `COL_MASS_K_*` / `COL_H2H2_C` / `COL_H2HE_C` constants in `env_uranus`. |
| `mk_*_lbl_*_highres_ec_uranus` | Per-gas, per-band line-by-line (LBL) cross-section generators (e.g. `mk_ch4_lbl_lw_highres_ec_uranus`, `mk_c2h2_lbl_sw_highres_ec_uranus`). Each runs `Ccorr_k` over every high-res band and concatenates the results with `lblcat` into a single LBL `.nc` file for that gas/region. |
| `mk_sp_*_skel` | Builds the skeleton spectral files (band structure, gas/CIA list) using SOCRATES's `prep_spec`. |
| `mk_sp_*_uranus` | Top-level driver scripts (e.g. `mk_sp_lw_highres_ec_uranus`, `mk_sp_sw_lowres_ec_uranus`) that build the skeleton, run `Ccorr_k` for each gas and CIA continuum using the pre-built LBL files, assemble everything into a spectral file with `prep_spec`, and tidy it with `tidy_90`. |
| `planet.nml` | SOCRATES namelist with Uranus physical constants (radius, gravity, gas constant, specific heat) and radiance/PL control flags used by `Cl_run_cdf`. |
| `mk_uranus_atm` | Converts a raw atmosphere profile (`uranus_atm.raw`) plus solar zenith angle/flux into the CDF/NetCDF input files (`uranus.h2`, `uranus.ch4`, `uranus.t`, `uranus.szen`, `uranus.stoa`, `uranus.surf`, etc.) needed to run SOCRATES. |
| `pt_uranus` | P-T table centred on a representative temperature profile +/-10,20K either side, used by `Ccorr_k -F` when generating gas absorption coefficients. |
| `pt_cont_uranus` | Temperature-only table used for generating CIA continuum absorption coefficients. |
| `ref_pt_uranus` | Reference pressure/temperature points per gas index, used to normalise absorption coefficients during spectral file construction. |
| `run_uranus` | End-to-end example run: builds the atmosphere with `mk_uranus_atm` and runs LW/SW radiative transfer with `Cl_run_cdf`, outputing fluxes/heating rates per level. |
