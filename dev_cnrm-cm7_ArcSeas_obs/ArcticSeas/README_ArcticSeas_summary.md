# ArcticSeas component: summary of work so far

Context file to start a new discussion on **adding observational references** to the ArcticSeas atlas.
Working directory: `dev_cnrm-cm7_ArcSeas_obs` (duplicated from `dev_cnrm-cm7_ArcSeas_grids`, branch `myversion`).
Status: tested with REF1, REF3 and CNRM-CM6.2 (runoffs corr); everything works as expected.
Files of the component: `ArcticSeas/params_ArcticSeas.py` (user parameters) and `ArcticSeas/diagnostics_ArcticSeas.py` (logic).

## 1. Goal

A C-ESM-EP (CliMAF) component producing an html page with, for each Arctic sea, time series of **sea ice area** (from `sic`)
and **sea ice volume** (from `sit`), one curve per simulation of `models` (see `datasets_setup.py`). Written in the house
style of the other components (MyOwnDiag, MainTimeSeries, NEMO_main).

## 2. How it works

External scripts from the `sea_ice_diag_tools` repo (`ArcticSeas_tools_dir`), called through `subprocess.run(check=True)`:

- `masks/create_mask_regions.py --maskfile --gridfile --out`: builds the sea-mask netcdf (one 2D variable per sea, `long_name` = sea name). Run only if the mask file is missing.
- `masks/check_masks.py --gridfile --mask --label --latcutoff`: plots the Arctic seas on a polar map (always writes `check_masks_arctic.png`, which is renamed to `<mask name without .nc>.png` in `workdir`).
- `comp_seaiceindex.py --data --var --mask --grid --dxvar --dyvar --meanORsum --out`: computes per-sea area-weighted sum (or mean) of a variable.

Processing flow in `diagnostics_ArcticSeas.py` (inside `if do_ArcticSeas_timeseries:`):

1. Per simulation, resolve grid + mask files (`resolve_grid_and_mask`); build the mask if absent (`ensure_mask`).
2. Top of the page: one "Arctic seas" check-mask figure **per distinct grid** (not per simulation); its label lists every simulation sharing that grid.
3. For each concept in `ARCTIC_SEAS_VARIABLES = ['sic', 'sit']` and each simulation: resolve the real netcdf variable name, build the period with `get_period_manager(wmodel, diag='ts')`, check the cache, otherwise `ds(**wmodel)` (optionally `ccdo yearmean`), `cfile()`, then `comp_seaiceindex.py`.
4. Plot per sea x variable: read the cached netcdf of every simulation, build a matplotlib figure, embed with `cell(...)`.

## 3. Parameters (`params_ArcticSeas.py`)

| Parameter | Role |
|---|---|
| `ArcticSeas_tools_dir` | clone of `sea_ice_diag_tools` |
| `ArcticSeas_gridfile` / `ArcticSeas_maskfile` | **dicts keyed by simulation `customname`** (or `experiment` if none): one grid file and one mask file per simulation. Simulations sharing a grid repeat the same paths (mask built once). A simulation absent from `ArcticSeas_gridfile` is skipped with a log message. If a simulation has a grid but no mask entry, the mask path is derived from the grid file name (`mesh_mask`/`meshmask` stripped) under `ArcticSeas_cache_dir/masks/` and built automatically |
| `ArcticSeas_dxvar`, `ArcticSeas_dyvar` | cell size variables (`e1t`, `e2t`) |
| `ArcticSeas_cache_dir` | cache of computed per-sea series: one sub-directory per simulation, one netcdf per variable (`sia`, `siv`, `sic`, `sit` depending on `meanORsum`) |
| `ArcticSeas_variable_names` | dict-of-dicts per simulation, e.g. `{'REF1': {'sic': 'siconc', 'sit': 'sithic'}}` (the real netcdf variable names can differ between simulations, e.g. N4CPL `sithic` vs N3CPL `sivolu`; the values currently in the params file were validated on REF1, REF3 and CNRM-CM6.2 and must be kept as they are) |
| `ArcticSeas_meanORsum` | `'sum'` (area/volume) or `'mean'` (mean concentration/thickness) |
| `ArcticSeas_seas_list` | subset of seas to plot (21 codes: `arcticoc`, 16 IHO subdivisions, `margseas`, `centrarc`, `wcentarc`, `ecentarc`) |
| `ArcticSeas_annual_mean` | annual mean before `comp_seaiceindex.py` |
| `ArcticSeas_on_missing_simulations` | `'error'` (error message instead of the figure if some simulations lack data) or `'partial'` (plot the available ones). If none has data: always an error message |
| `ArcticSeas_thumbnail_size` | thumbnail size |

`ARCTIC_SEAS_VARIABLES` (`'sic'`, `'sit'`) is a fixed internal constant of `diagnostics_ArcticSeas.py`, deliberately **not** a user parameter.
The `reference = 'default'` and `season` / `latcutoff` parameters in the params file are inherited from other examples and are **not used** by ArcticSeas yet.

## 4. Key design decisions / behaviours

- **Cache**: a simulation is recomputed only if no cache exists or if the simulation now extends further in time than the cached file (`latest_year_of_wmodel` vs `latest_year_in_ncfile`).
- **Robustness**: missing data / grid / mask / variable name for a simulation is logged and skipped, not fatal. A guard `if 'period' not in wmodel: continue` avoids `build_plot_title` crashing on simulations without data.
- **Sea ice concentration units**: if `sic` max > 1.5 (percent), a corrected copy divided by 100 is written in `workdir` before calling `comp_seaiceindex.py` (CliMAF's cache file is never modified).
- **Display units** (local to the plot, cache file untouched): in `'sum'` mode, `sic` and `sit` values are divided by `1e12` (m2 -> millions km2; m3 -> thousand km3) to match the y labels "Sea ice area (Millions km2)" / "Sea ice volume (Thousand km3)". The two cases share a single `if variable in ('sic', 'sit') and ArcticSeas_meanORsum == 'sum'`.
- Naming: loop variable `variable` is the concept (`sic`/`sit`), `real_variable` the actual netcdf name, `plot_var` an entry of `variable_plot_specs`, `model_label_of(wmodel)` returns `customname` or `experiment`.

## 5. Pitfalls found (useful for the next step)

- `get_period_manager(wmodel, diag=...)` only understands `None`, `'clim'`, `'ts'`; `'ts'` is the one resolving `ts_period` into a concrete `period`.
- `dict.get(k, default)` evaluates `default` eagerly (it crashed through `build_plot_title` when `period` was unresolved).
- Framework function is `build_period_str` (not `build_str_period`).
- `check_masks.py` always writes the same relative PNG name, hence the rename to the mask name to avoid collisions between grids.
- Careful with blind `replace_all` renames (`plot_spec` is a substring of `variable_plot_specs`).
- Project aliasing in `project_N4cpl.py`: `calias("N4CPL", 'sit', 'sithic')`, `calias("N3CPL", 'sit', 'sivolu')`, `calias(..., 'sic', 'siconc', scale=100.)`.
- Unrelated issue diagnosed earlier: SST missing in `MainTimeSeries` came from a simulation removed from the `models` list in `datasets_setup.py`.

## 6. Current simulations (`datasets_setup.py`)

`REF1` (project N4CPL, `AOGCM_drv632_NEMO422_LR_t7`), `REF3` (project N3CPL, `CM67_N4_3_REF3`), plus the CNRM-CM6.2 "runoffs corr" run, using the grid/mask files in `params_ArcticSeas.py`
(cnrmcm6 vs cnrmcm7 meshmasks and Arctic-seas masks).

## 7. Next topic: adding observational references

Open questions to settle at the start of the new discussion:

- Which observational products (sea ice area/extent/volume: e.g. NSIDC, OSI SAF, PIOMAS for volume...), which variables and time coverage.
- Observations may be provided as **concentration/thickness fields on their own grid** (then need their own grid file, mask, `ArcticSeas_variable_names` entry, unit handling and possibly the 0-100 correction), or as **already-integrated per-sea series** (then they bypass `comp_seaiceindex.py` and can only be plotted).
- How to declare them: an entry in `models`, or a dedicated parameter (e.g. `ArcticSeas_obs`) with its own grid/mask/variable/period, and how it interacts with `reference = 'default'`.
- Plot style (distinct color/linestyle, label), behaviour with the `ArcticSeas_on_missing_simulations` policy (should a missing obs count as a missing simulation?), and cache layout for obs.
- Period handling: obs period vs simulations' `ts_period` (and annual mean consistency).
