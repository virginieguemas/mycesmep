# ------------------------------------------------------------------------------------------ \
# --                                                                                        - \
# --                                                                                         - \
# --      User Interface for:                                                                 - \
# --          CliMAF Earth System Model Evaluation Platform                                     - \
# --             - component: ArcticSeas                                                    - |
# --                                                                                             - |
# --      Developed within the ANR Convergence Project                                           - |
# --      CNRM GAME, IPSL, CERFACS                                                               - |
# --                                                                                             - |
# --                                                                                            - /
# --------------------------------------------------------------------------------------------- /

# -- Preliminary settings: import module, set the verbosity and the 'safe mode'
# ---------------------------------------------------------------------------- >
from os import getcwd

# -- Set the verbosity of CliMAF (minimum is 'critical', maximum is 'debug', intermediate -> 'warning')
verbose = 'debug'
# -- Safe Mode (set to False and verbose='debug' if you want to debug)
safe_mode = True
# -- Set to True to clean the CliMAF cache
clean_cache = False
# -- Patterns to clean the cache at the end of the execution of the atlas
routine_cache_cleaning = [dict(age='+20')]
# -- Parallel and memory instructions
do_parallel = False

# -- Head title of the atlas
# ---------------------------------------------------------------------------- >
atlas_head_title = "Arctic Seas - Sea Ice Area and Volume"


# -- Set the overall season, region and geographical domain
# ---------------------------------------------------------------------------- >
season = 'ANM'
proj = 'NH'
domain = dict()


# ---------------------------------------------------------------------------- >
# -- Arctic Seas : per-sea sea ice area / volume time series
# -- Relies on the external script comp_seaiceindex.py from the sea_ice_diag_tools
# -- repository, which averages a variable over each region described in a mask
# -- netcdf file (one region = one data variable in the mask file).
# --   -> applied to 'siconc'  => (mean) sea ice area  per sea
# --   -> applied to 'sithick' => (mean) sea ice volume per sea
# ---------------------------------------------------------------------------- >
do_ArcticSeas_timeseries = True

# -- Path to the sea_ice_diag_tools clone holding comp_seaiceindex.py
# -- (https://git.meteo.fr/cnrm-cerfacs-esm/seaice/sea_ice_diag_tools)
# -- /!\ TO BE ADAPTED to the actual location on the machine running the atlas
ArcticSeas_tools_dir = '/home/guemas/sea_ice_diag_tools'

# -- Mask netcdf file: one data variable (2D, on the model t-grid) per sea/region,
# -- with a 'long_name' attribute giving the sea name; built by
# -- ArcticSeas_tools_dir/masks/create_mask_regions.py
# -- /!\ As of 2026-08, create_mask_regions.py stops (sys.exit()) before writing out
# --     the individual Arctic seas (Barents, Kara, Laptev...): that part of the
# --     script is still legacy/unreachable code. Point this to a mask file that
# --     actually contains the per-sea breakdown once that script is completed.
ArcticSeas_maskfile = '/home/guemas/tmp/test_regions/mask.ArcticSeas.cnrmcm7.nc'

# -- Variables used to compute the per-sea indices
# --   siconc  -> sea ice area  (mean sea ice concentration over the sea)
# --   sithick -> sea ice volume (mean sea ice thickness over the sea)
ArcticSeas_variables = ['siconc', 'sithick']

# -- Restrict to a subset of seas (list of the 'long_name'/variable names found in
# -- ArcticSeas_maskfile); set to None (or empty list) to use every region found in
# -- the mask file
ArcticSeas_seas_list = [
    'Central Arctic',
    'Barents Sea',
    'Kara Sea',
    'Laptev Sea',
    'East Siberian Sea',
    'Chukchi Sea',
    'Beaufort Sea',
    'Baffin Bay',
    'Hudson',
    'Labrador Sea',
    'Bering',
    'Okhotsk',
    'Nordic Seas',
]

# -- Compute annual means before calling comp_seaiceindex.py (recommended: shorter,
# -- more readable multi-decadal time series). Set to False to keep the native
# -- (monthly) frequency.
ArcticSeas_annual_mean = True

# -- Thumbnail size for the time series plots
ArcticSeas_thumbnail_size = '450*300'


# -- Some settings -- customization
# ---------------------------------------------------------------------------- >
# -- Name of the html file
# -- if index_name is set to None, it will be build as user_comparisonname_season
index_name = None


# ---------------------------------------------------------------------------------------- #
# -- END                                                                                -- #
# ---------------------------------------------------------------------------------------------- #
