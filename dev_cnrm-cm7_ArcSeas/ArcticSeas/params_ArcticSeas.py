# ------------------------------------------------------------------------------------------ \
# --                                                                                        - \
# --                                                                                         - \
# --      User Interface for:                                                                 - \
# --          CliMAF Earth System Model Evaluation Platform                                    - \
# --             - component: ArcticSeas                                                       - |
# --                                                                                           - |
# --      Contact : virginie.guemas@meteo.fr                                                   - |
# --      History : Creating September 2026   -    virginie.guemas@meteo.fr                    - |
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

# -- Set the reference against which we plot the diagnostics 
# ---------------------------------------------------------------------------- >
# --    -> 'default' uses variable2reference to point to a default
# --       reference dataset (obs and reanalyses)
# --    -> you can set reference to a dictionary that will point any other
# --       climaf dataset
# --       For instance, you can set it to models[0] if you want to see the
# --       differences relative to the first simulation of the list 'models'
reference = 'default'

# -- Head title of the atlas
# ---------------------------------------------------------------------------- >
atlas_head_title = "Arctic Seas - Sea Ice Area and Volume"


# -- Set the overall season
# ---------------------------------------------------------------------------- >
# -> Choose among all the possible values taken by clim_average (see help(clim_average)) like JFM, December,...
season    = 'ANM'
# -- Minimum latitude for stereographic projection
latcutoff =  50. 


# ---------------------------------------------------------------------------- >
# -- Arctic Seas : per-sea sea ice area / volume time series
# -- Relies on the external script comp_seaiceindex.py from the sea_ice_diag_tools
# -- repository, which computes, for a variable, either the area-weighted mean or
# -- the area-weighted sum (area-integral) over each region described in a mask
# -- netcdf file (one region = one data variable in the mask file), weighting by
# -- the grid cell area (dxvar*dyvar) taken from a grid file.
# --   -> applied to 'siconc'  with meanORsum='sum'  => actual sea ice area  per sea
# --   -> applied to 'sithick' with meanORsum='sum'  => actual sea ice volume per sea
# -- (meanORsum='mean' instead gives the area-weighted mean concentration/thickness)
#
# -- The sea mask file can be computed by create_mask_regions.py from the same
# -- sea_ice_diag_tools repository, given the a grid description file holding
# -- the latitudes and longitudes for each grid point and a land-sea mask to know
# -- where the coastline is for a particular grid
# ---------------------------------------------------------------------------- >
do_ArcticSeas_timeseries = True

# -- Path to the sea_ice_diag_tools clone holding :
# -- comp_seaiceindex.py, create_mask_regions.py and check_masks.py
# -- https://github.com/virginieguemas/sea_ice_diag_tools.git
# -- git@github.com:virginieguemas/sea_ice_diag_tools.git
ArcticSeas_tools_dir = '/cnrm/ioga/Users/guemas/diag_tools/sea_ice_diag_tools'

# -- Mask netcdf file: one data variable (2D, on the model t-grid) per sea/region,
# -- with a 'long_name' attribute giving the sea name; built by
# -- ArcticSeas_tools_dir/masks/create_mask_regions.py (grid='cnrmcm7' case).
# -- If the mask netcdf file is missing, diagnostics_ArcticSeas.py runs 
# -- create_mask_regions.py itself
# -- (from this same directory, so its relative 'mask.ArcticSeas.cnrmcm7.nc' output
# -- lands here) to build it before going any further.
ArcticSeas_maskfile = '/cnrm/ioga/Users/guemas/gridfiles/seas/mask.ArcticSeas.cnrmcm7.nc'
#
# WARNING : Script create_mask_regions to be made generic enough to send grid description
# file as input arguments (right now hardcoded and sea masks built by hand)
# WARNING : As coded now, ArcticSeas can only compare simulations with the same grid
# Need generalization

# -- Grid description netcdf file (NEMO mesh_mask/mesh_hgr-like file) giving the
# -- size of the grid cells (used to area-weight the mean computed over each sea).
# -- Used as a fallback only: if the model dictionary (see datasets_setup.py) already
# -- defines a 'mesh_hgr' (or 'gridfile') key, that per-model file is used instead,
# -- since different simulations may not share the same grid/resolution.
ArcticSeas_gridfile = '/cnrm/ioga/Users/guemas/gridfiles/meshmask/mesh_mask.cnrmcm7.nc'
ArcticSeas_dxvar = 'e1t'
ArcticSeas_dyvar = 'e2t'
#
# WARNING : Tested only on comparison of simulations with the same grids

# -- Where the per-sea sea ice area/volume time series computed by comp_seaiceindex.py
# -- are cached, one sub-directory per simulation (based on its customname), one
# -- netcdf file per variable inside. Reused across atlas runs: a simulation is only
# -- (re-)processed by comp_seaiceindex.py if no cache exists yet, or if the
# -- simulation now extends further in time than what is cached (i.e. it was
# -- lengthened since the cache was last built).
ArcticSeas_cache_dir = '/cnrm/ioga/Users/guemas/ArcticSeas_cache'

# -- Variables used to compute the per-sea indices
# --   siconc  -> sea ice area
# --   sithick -> sea ice volume
ArcticSeas_variables = ['siconc', 'sithick']

# -- 'sum'  -> area-weighted sum (area-integral): siconc gives the actual sea ice
# --           area, sithick gives the actual sea ice volume, in each sea
# -- 'mean' -> area-weighted mean: siconc gives the mean sea ice concentration,
# --           sithick gives the mean sea ice thickness, in each sea
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
