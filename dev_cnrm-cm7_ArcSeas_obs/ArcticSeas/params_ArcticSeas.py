# ------------------------------------------------------------------------------------------ \
# --                                                                                        - \
# --                                                                                         - \
# --      User Interface for:                                                                 - \
# --          CliMAF Earth System Model Evaluation Platform                                    - \
# --             - component: ArcticSeas                                                       - |
# --                                                                                           - |
# --      Contact : virginie.guemas@meteo.fr                                                   - |
# --      History : Created September 2026   -    virginie.guemas@meteo.fr                    - |
# --                                                                                           - /
# --------------------------------------------------------------------------------------------- /

# -- Preliminary settings: import module, set the verbosity and the 'safe mode'
# ---------------------------------------------------------------------------- >
from os import getcwd

# -- Set the verbosity of CliMAF (minimum is 'critical', maximum is 'debug', intermediate -> 'warning')
verbose = 'debug'
# -- Safe Mode (set to False and verbose='debug' if you want to debug)
safe_mode = False
# -- Set to True to clean the CliMAF cache
clean_cache = False
# -- Patterns to clean the cache at the end of the execution of the atlas
routine_cache_cleaning = [dict(age='+20')]
# -- Parallel and memory instructions
do_parallel = False
#
# WARNING : The arguments above are copied from other example, not sure they actually
# have an impact

# -- Set the reference against which we plot the diagnostics 
# ---------------------------------------------------------------------------- >
# --    -> 'default' uses variable2reference to point to a default
# --       reference dataset (obs and reanalyses)
# --    -> you can set reference to a dictionary that will point any other
# --       climaf dataset
# --       For instance, you can set it to models[0] if you want to see the
# --       differences relative to the first simulation of the list 'models'
reference = 'default'

# WARNING : I put default in case there had been some sea ice data in C-ESM-EP
# for sea ice but there isn't - need to add something

# -- Head title of the atlas
# ---------------------------------------------------------------------------- >
atlas_head_title = "Arctic Seas - Sea Ice Area and Volume"


# -- Set the overall season
# ---------------------------------------------------------------------------- >
# -> Choose among all the possible values taken by clim_average (see help(clim_average)) like JFM, December,...
season    = 'ANM'
# -- Minimum latitude for stereographic projection
latcutoff =  50. 
#
# WARNING : For now the arguments season and latcutoff are not used, to be improved

# ---------------------------------------------------------------------------- >
# -- Arctic Seas : per-sea sea ice area / volume time series
# -- Relies on the external script comp_seaiceindex.py from the sea_ice_diag_tools
# -- repository, which computes, for a variable, either the area-weighted mean or
# -- the area-weighted sum (area-integral) over each region described in a mask
# -- netcdf file (one region = one data variable in the mask file), weighting by
# -- the grid cell area (dxvar*dyvar) taken from a grid file.
# --   -> applied to 'sic' with meanORsum='sum'  => actual sea ice area  per sea
# --   -> applied to 'sit' with meanORsum='sum'  => actual sea ice volume per sea
# -- (meanORsum='mean' instead gives the area-weighted mean concentration/thickness)
#
# -- The sea mask file can be computed by create_mask_regions.py from the same
# -- sea_ice_diag_tools repository, given the a grid description file holding
# -- the latitudes and longitudes for each grid point and a land-sea mask to know
# -- where the coastline is for a particular grid
# ---------------------------------------------------------------------------- >
do_ArcticSeas_timeseries = True
# WARNING : put there to do as other examples do, but why is that useful ?

# -- Path to the sea_ice_diag_tools clone holding :
# -- comp_seaiceindex.py, create_mask_regions.py and check_masks.py
# -- https://github.com/virginieguemas/sea_ice_diag_tools.git
# -- git@github.com:virginieguemas/sea_ice_diag_tools.git
ArcticSeas_tools_dir = '/cnrm/ioga/Users/guemas/diag_tools/sea_ice_diag_tools'

# -- Observational references, declared here rather than in `models` (datasets_setup.py)
# -- so that the other components do not try to process them as simulations.
# -- Same syntax as an entry of `models`: a CliMAF dataset description (the project is
# -- declared in my_projects.py), with 'frequency' and 'ts_period' ('full' = every
# -- available date). Each observation is then looked up by its 'customname' in
# -- ArcticSeas_gridfile, ArcticSeas_maskfile and ArcticSeas_variable_names exactly like
# -- a simulation; it only needs an ArcticSeas_variable_names entry for the variables
# -- it provides (e.g. no 'sit' for a sea ice concentration product).
# -- Observations are plotted in black/grey with distinct line styles, and a missing
# -- observation never counts as a missing simulation (see ArcticSeas_on_missing_simulations).
ArcticSeas_obs = [
    dict(project='NSIDC_CDR_V6', customname='nsidc_cdr_v6', frequency='monthly', ts_period='full'),
]

# -- How the observations are plotted, since their dates generally do not match those of
# -- the simulations:
# --   'climato' -> climatology of the observations over ArcticSeas_obs_clim_period,
# --                repeated every year along the dates of the simulations (a single mean
# --                value if ArcticSeas_annual_mean, otherwise the mean of each calendar month)
# --   'real'    -> observations plotted at their actual dates
# -- ArcticSeas_obs_clim_period : 'full' (every available year) or 'YYYY-YYYY' (only used
# -- with 'climato').
# -- Both can be set for one observation only with the keys 'obs_dates' and 'obs_clim_period'
# -- in its ArcticSeas_obs entry, e.g. dict(..., obs_dates='real') or
# -- dict(..., obs_clim_period='1991-2020').
ArcticSeas_obs_dates = 'climato'
ArcticSeas_obs_clim_period = 'full'


# -- Grid description file (NEMO mesh_mask-like: land-sea mask, lon/lat, cell sizes
# -- e1t/e2t) and matching Arctic-seas mask file (one data variable (2D) per sea/region,
# -- with a 'long_name' attribute giving the sea name ; built on demand by
# -- ArcticSeas_tools_dir/masks/create_mask_regions.py if missing), one entry each for
# -- every simulation in `models` (see datasets_setup.py), keyed by its 'customname'.
# -- Simulations that share the same grid should simply repeat the same paths: the mask
# -- is then only built once and reused.
# -- A simulation without a 'customname' is instead keyed by its 'experiment'.
# -- A simulation missing from these dictionaries is skipped, with a message in the log.
ArcticSeas_gridfile = {
    'CNRM-CM6.2 (runoffs corr)': '/cnrm/ioga/Users/guemas/gridfiles/meshmask/mesh_mask.cnrmcm6.nc',
    'REF1': '/cnrm/ioga/Users/guemas/gridfiles/meshmask/mesh_mask.cnrmcm7.nc',
    'REF3': '/cnrm/ioga/Users/guemas/gridfiles/meshmask/mesh_mask.cnrmcm7.nc',
    # -- NSIDC polar stereographic 25 km grid (psn25, shared with CryoSat-2 RDEFT4),
    # -- built by sea_ice_diag_tools/masks/create_grid_psn25.py
    'nsidc_cdr_v6': '/cnrm/ioga/Users/guemas/gridfiles/meshmask/mesh_mask.psn25.nc',
}

ArcticSeas_maskfile = {
    'CNRM-CM6.2 (runoffs corr)': '/cnrm/ioga/Users/guemas/gridfiles/seas/mask.ArcticSeas.cnrmcm6.nc',
    'REF1': '/cnrm/ioga/Users/guemas/gridfiles/seas/mask.ArcticSeas.cnrmcm7.nc',
    'REF3': '/cnrm/ioga/Users/guemas/gridfiles/seas/mask.ArcticSeas.cnrmcm7.nc',
    'nsidc_cdr_v6': '/cnrm/ioga/Users/guemas/gridfiles/seas/mask.ArcticSeas.psn25.nc',
}
ArcticSeas_dxvar = 'e1t'
ArcticSeas_dyvar = 'e2t'

# -- Where the per-sea sea ice area/volume time series computed by comp_seaiceindex.py
# -- are cached, one sub-directory per simulation (based on its customname), one
# -- netcdf file per variable inside. Reused across atlas runs: a simulation is only
# -- (re-)processed by comp_seaiceindex.py if no cache exists yet, or if the
# -- simulation now extends further in time than what is cached (i.e. it was
# -- lengthened since the cache was last built).
ArcticSeas_cache_dir = '/cnrm/ioga/Users/guemas/NO_SAVE/ArcticSeas_cache'

# -- Variable name to find the sea ice concentration (sic) and the sea ice tickness (sit)
# -- in the output files, one entry per simulation (keyed by its 'customname',
# -- or its 'experiment' if it has none -- same key as ArcticSeas_gridfile/
# -- ArcticSeas_maskfile). A simulation missing from this dictionary is skipped, 
# -- with a message in the log. Ex : N4CPL's volume variable is 'sithic', but
# -- N3CPL's is 'sivolu' 
ArcticSeas_variable_names = {
    'CNRM-CM6.2 (runoffs corr)': {'sic': 'siconc', 'sit': 'sivolu'},  
    'REF1': {'sic': 'siconc', 'sit': 'sithic'}, 
    'REF3': {'sic': 'siconc', 'sit': 'sithic'},
    'nsidc_cdr_v6': {'sic': 'cdr_seaice_conc_monthly'},
}

# -- 'sum'  -> area-weighted sum (area-integral): sic gives the actual sea ice
# --           area, sit gives the actual sea ice volume, in each sea
# -- 'mean' -> area-weighted mean: sic gives the mean sea ice concentration,
# --           sit gives the mean sea ice thickness, in each sea
ArcticSeas_meanORsum = 'sum'

# -- Restrict to a subset of seas (list of the netcdf variable names found in
# -- ArcticSeas_maskfile); set to None (or an empty list) to use every region found
# -- in the mask file.
ArcticSeas_seas_list = [
    'arcticoc',  # Arctic Ocean (overall basin)
    'eastsibe',  # East Siberian Sea
    'laptevse',  # Laptev Sea
    'karaseax',  # Kara Sea
    'barentse',  # Barents Sea
    'whitesea',  # White Sea
    'greenlds',  # Greenland Sea
    'norwegia',  # Norwegian Sea
    'icelands',  # Iceland Sea
    'davisstr',  # Davis Strait
    'hudsonst',  # Hudson Strait
    'hudsonba',  # Hudson Bay
    'baffinba',  # Baffin Bay
    'lincolns',  # Lincoln Sea
    'nwpassag',  # Northwestern Passages
    'beaufort',  # Beaufort Sea
    'chukchis',  # Chukchi Sea
    'margseas',  # Arctic Marginal Seas (union of the 16 seas above)
    'centrarc',  # Central Arctic
    'wcentarc',  # Western Central Arctic
    'ecentarc',  # Eastern Central Arctic
]

# -- Compute annual means before calling comp_seaiceindex.py (recommended: shorter,
# -- more readable multi-decadal time series). Set to False to keep the native
# -- (monthly) frequency.
ArcticSeas_annual_mean = True

# -- What to do, for a given sea/variable plot, when data is available for only
# -- some of the simulations (the others failed: no file found, missing grid file,
# -- comp_seaiceindex.py error...):
# --   'error'   -> replace the whole figure by an error message listing the
# --                simulations that are missing
# --   'partial' -> plot the simulations that do have data anyway, silently
# --                leaving out the missing ones
# -- Regardless of this setting, if NO simulation has data at all for a given
# -- sea/variable, an error message is always shown instead of an (empty) figure.
# -- Observations (ArcticSeas_obs) are not concerned: a missing observation is simply
# -- not plotted.
ArcticSeas_on_missing_simulations = 'error'

# -- Thumbnail size for the time series plots
ArcticSeas_thumbnail_size = '450*300'

# -- Some settings -- customization
# ---------------------------------------------------------------------------- >
# -- Name of the html file
# -- if index_name is set to None, it will be build as user_comparisonname_season
# -- with comparisonname being the name of the parameter file without 'params_'
# -- (and '.py' of course)
# ---------------------------------------------------------------------------- >
index_name = None

# ---------------------------------------------------------------------------------------- #
# -- END                                                                                -- #
# ---------------------------------------------------------------------------------------------- #
