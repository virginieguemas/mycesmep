# ------------------------------------------------------------------------------------------------------ \
# --                                                                                                    - \
# --      Scientific diagnostics for the                                                                 - \
# --          CliMAF Earth System Model Evaluation Platform                                               - |
# --      diagnostics_ArcticSeas.py                                                                       - |
# --        ==> add html code to 'index' (initialized with 'header')                                      - |
# --            using the CliMAF html toolbox (start_line, cell, close_table... )                         - |
# --            to create the Arctic Seas atlas page                                                      - |
# --                                                                                                      - |
# --      Time series of sea ice area (sic) and sea ice volume (sit) computed as an                       - |
# --      area-weighted sum over each Arctic sea, using the external script comp_seaiceindex.py           - |
#         (sea_ice_diag_tools repository) applied to a mask file describing the individual seas and       - |
#         a grid file giving the cell areas. If not provided by the user, the mask file for the           - |
#         individual seas is computed by create_mask_regions.py from the sea_ice_diag_tools repository    - |
# --                                                                                                      - /
# --                                                                                                     - /
# --      Contact : virginie.guemas@meteo.fr                                                            - /
# --      History : Created September 2026   -    virginie.guemas@meteo.fr                             - /
# --                                                                                                  - /
# ---------------------------------------------------------------------------------------------------- /

from os import getcwd
import os
import re
import copy
import subprocess
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# -- Init html index - Example from MyOwnDiag - fonction from CLIMAF library 
# -----------------------------------------------------------------------------------
index = header(atlas_head_title, style_file=style_file)


def sanitize(name):
    """Turn a sea/model name into a safe token for file names."""
    # Replace every character other than letter or number by _ (/-:) and join
    # together in a safe name all the initial and replaced characters
    cleaned = "".join(c if c.isalnum() else "_" for c in name)
    return re.sub(r"_+", "_", cleaned).strip("_")


def model_label_of(wmodel):
    """Identifier for a simulation, used to look it up in the ArcticSeas_*
    per-simulation dictionaries (params_ArcticSeas.py) and in log/error messages: its
    customname if defined, else experiment if defined, else unknown simulation ."""
    return wmodel.get('customname', wmodel.get('experiment', 'unknown simulation'))


# -- Generic name for the variables on which the average or sum should be
# -- computed : 'sic' -> sea ice area, 'sit' -> sea ice volume.
# -- Fixed names (not a user parameter). A correspondance is set in 
# -- params_ArcticSeas.py to find actual netcdf variable for a given
# -- simulation through ArcticSeas_variable_names.
ARCTIC_SEAS_VARIABLES = ['sic', 'sit']

# -- File names for the cache files (stored in one directory per experiment)
# -- sia/siv = sea ice area/volume (area-weighted sum)
# -- sic/sit = sea ice concentration/thickness (area-weighted mean)
CACHE_VARIABLE_LABELS = {
    ('sic', 'sum'): 'sia',
    ('sit', 'sum'): 'siv',
    ('sic', 'mean'): 'sic',
    ('sit', 'mean'): 'sit',
}

# -- (color, linestyle) for the observations (ArcticSeas_obs), so that
# -- they stand out from the simulations (default matplotlib color cycle, solid lines)
OBS_STYLES = [('black', '-'), ('black', '--'), ('black', ':'), ('black', '-.'),
              ('grey', '-'), ('grey', '--'), ('grey', ':'), ('grey', '-.')]


def latest_year_of_wmodel(wmodel):
    """Latest year covered by a period-managed model dict, from build_period_str()."""
    # build_period_str return a chain of characters holding the first and last years
    # then re.findall(r'\d{4}' finds the groups of 4 consecutive digits, months are dropped
    # at this stage because there are only 2 digits.
    # the first and last simulation years are stored in years
    years = re.findall(r'\d{4}', str(build_period_str(wmodel)))
    return int(years[-1]) if years else None


def latest_year_in_ncfile(ncfile):
    """Latest year found in the (single) time-like dimension of a netcdf file."""
    if not os.path.exists(ncfile):
        return None
    # If there is no netcdf file yet, this function returns None
    try:
        cached = xr.open_dataset(ncfile)
        first_var = next(iter(cached.data_vars), None)
        # iter iterates on the keys of cached.data_vars dictionary and next takes the
        # first, if there is none, fist_var becomes None (very unlikely but in case ...)
        if first_var is None:
            cached.close()
            return None
        # If the netcdf file is empty, this function returns None
        time_dim = cached[first_var].dims[0]
        times = cached[time_dim].values
        cached.close()
        return int(str(times[-1])[:4]) if len(times) else None
        # times[-1] is the last date of the file, its first four digits are the year
    except Exception:
        return None
    # If for any reason the last year of the file can not be computed, this function returns None


# -- Message issued by comp_seaiceindex.py when the field has missing values within a sea
SEA_WARNING_PATTERN = re.compile(r"Missing values found in field for sea '([^']+)' where mask == 1")
# [^']+ means one or more character that are not '
# WARNING : if the comp_seaiceindex message were to change, the change should be reflected in the
# line above otherwise the ArcticSeas option could not detect it

def warnings_file_of(cache_file):
    """Text file, next to a per-sea cache file, holding the 'missing values' warnings
    issued by comp_seaiceindex.py when it computed that cache file (one per line, may be empty)."""
    return cache_file + '.warnings'


def write_sea_warnings(cache_file, stderr):
    """Keep the 'missing values' warnings found in the stderr of comp_seaiceindex.py."""
    with open(warnings_file_of(cache_file), 'w') as f:
    # A cach_file_name.warnings text file is opened, overwritten if there was one before
        for match in SEA_WARNING_PATTERN.finditer(stderr):
        # We find all occurrences (finditer) of the warning message about missing values 
        # (SEA_WARNING_PATTERN) in a sea in the standard error output of comp_seaiceindex.py
            f.write(match.group(0) + '\n')
            # match.group(0) is the full error message output by comp_seaiceindex.py without
            # the path of the script nor the warnings.warn added by python


def read_sea_warnings(cache_file):
    """Builds a dictionary holding for each sea name as a key : the 'missing values' warning 
    message, as kept by write_sea_warnings."""
    sea_warnings = dict()
    if os.path.exists(warnings_file_of(cache_file)):
    # If the cache_file_name.warnings file exist, we read it
        with open(warnings_file_of(cache_file)) as f:
            for line in f:
                # Each line in this file should correspond with a sea with missing values
                match = SEA_WARNING_PATTERN.search(line)
                # SEA_WARNING_PATTERN is the generic error message from comp_seaiceindex
                # match is a Match object holding the full error message with the actual sea 
                # name in match.group(0) and the sea name in match.group(1) because a set of 
                # parenthesis was included in SEA_WARNING_PATTERN to be able to isolate that
                # sea name
                if match:
                    sea_warnings[match.group(1)] = match.group(0)
                    # If we find an occurrence with add an element to the dictionary
    return sea_warnings


def sea_areas_text(sea, labels, seaindex_files, area_attr, area_divisor):
    """ From an input sea (ex: centrarc), and for a list of experiment and observation
    labels (typically customname), this function will print a line listing the area
    of that sea in each of the experiment and observation dataset (it is not necessarily the
    same area because of the land/sea mask). It needs a dictionary containing for each label
    and each variable, its cache file : seaindex_files.

    One line listing, for each simulation/observation of labels, the area attribute
    area_attr (written by comp_seaiceindex.py) of sea in its cached per-sea file, divided by
    area_divisor. 'n/a' if not found (no data, or cache computed by an older version of
    comp_seaiceindex.py)."""
    parts = []
    # We define a list which will hold the parts of the final line to print
    for label in labels:
        # We loop on the experiment and observation labels
        area = None
        for variable in ARCTIC_SEAS_VARIABLES:
            # We loop on the Arctic variables just in case some are missing because it is the
            # same land/sea mask for all variables of a given experiment but observation often
            # have only of the variables
            outfile = seaindex_files.get(variable, dict()).get(label)
            # We get the cache file corresponding to that label and that variable, NONE if there
            # isn't any cache. We need to use get to avoid a crash when a variable isn't there.
            if outfile is None or not os.path.exists(outfile):
                continue # If we don't have any cache file, we go to next variable
            with xr.open_dataset(outfile) as out_ds:
                if sea in out_ds.data_vars and area_attr in out_ds[sea].attrs:
                # We look for the area_attr attribute in the cache file we found
                # It should contain the area of the sea in that experiment/observation
                    area = float(out_ds[sea].attrs[area_attr]) / area_divisor
                # area_divisor is to convert to Thousand km3 or Millons km2
            if area is not None:
                break
            # If we found the area of the sea for one variable, we skip the next variable
            # and go directly to next label
        parts.append('%s %s' % (label, 'n/a' if area is None else '%.0f' % area))
        # We add to parts a label followed by the area if found None if not
    return ', '.join(parts)


def obs_climatology_curve(da, clim_period, annual, sim_times):
    """Climatology of an observed per-sea series da (time as first dimension) over
    clim_period ('full' or 'YYYY-YYYY'), repeated along sim_times (the time values of the
    simulations): a single mean value if the series are annual means, otherwise the mean
    of each calendar month. Returns (values along sim_times, first year, last year of the
    climatology), or None if no observation falls within clim_period."""
    time_dim = da.dims[0]
    years = da[time_dim].dt.year
    if clim_period != 'full':
        first, last = [int(year) for year in clim_period.split('-')]
        da = da.where((years >= first) & (years <= last), drop=True)
        if da.sizes[time_dim] == 0:
            return None
        years = da[time_dim].dt.year
    if annual:
        values = np.full(len(sim_times), float(da.mean()))
    else:
        # -- Mean of each calendar month (NaN for a month never observed, e.g. summer
        # -- for a winter-only product), then picked for the month of each simulation date
        monthly_clim = da.groupby(time_dim + '.month').mean().reindex(month=range(1, 13))
        sim_months = xr.DataArray(sim_times, dims='simtime').dt.month
        values = monthly_clim.sel(month=sim_months).values
    return values, int(years.min()), int(years.max())


if do_ArcticSeas_timeseries:
    # 
    # ==> -- Open the section and an html file - Function from CLIMAF library
    # -----------------------------------------------------------------------------------------
    index += section("Sea ice area and volume per Arctic sea", level=4)
    #
    # ==> -- Control the size of the thumbnail -> thumbN_size
    # Note: ArcticSeas_thumbnail_size defined in params_ArcticSeas.py
    # Different from what is done in other diagnostics to avoid overwriting plot parameters
    # for all diagnostics running in parallel
    # -----------------------------------------------------------------------------------------
    if thumbnail_size:
        thumbN_size = thumbnail_size
    elif 'ArcticSeas_thumbnail_size' in dir() and ArcticSeas_thumbnail_size:
        thumbN_size = ArcticSeas_thumbnail_size
    else:
        thumbN_size = thumbnail_size_global
    # 
    # Create a working directory under the execution directory if non-existant already
    # -----------------------------------------------------------------------------------------
    workdir = os.path.join(getcwd(), 'ArcticSeas_workfiles')
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    #
    # Complete path to the scripts required for the diagnostics under ArcticSeas
    # -----------------------------------------------------------------------------------------
    comp_seaiceindex_script = os.path.join(ArcticSeas_tools_dir, 'comp_seaiceindex.py')
    create_mask_regions_script = os.path.join(ArcticSeas_tools_dir, 'masks', 'create_mask_regions.py')
    check_masks_script = os.path.join(ArcticSeas_tools_dir, 'masks', 'check_masks.py')
    #
    # Writing error messages in the html page - functions used from CLIMAF library
    # -----------------------------------------------------------------------------------------
    if not os.path.exists(comp_seaiceindex_script):
        index += open_table()
        index += start_line('Error')
        index += "comp_seaiceindex.py not found at %s ; check ArcticSeas_tools_dir in params_ArcticSeas.py" \
                  % comp_seaiceindex_script
        index += close_line() + close_table()

    else:
        # -- Simulations may not all share the same grid/land-sea mask (different NEMO
        # -- configurations, resolutions...). ArcticSeas_gridfile and ArcticSeas_maskfile
        # -- (params_ArcticSeas.py) are dictionaries giving, for each simulation (keyed by
        # -- its 'customname', or its 'experiment' if it has none), the grid description file
        # -- to use and the matching Arctic-seas mask file (built on demand if missing). A
        # -- simulation with a grid file but no mask file entry gets one derived automatically
        # -- from its grid file's name, under ArcticSeas_cache_dir.
        # -- Simulations sharing the same grid should repeat the same paths: the mask itself is
        # -- then only built once.
        # -----------------------------------------------------------------------------------------
        def resolve_grid_and_mask(wmodel):
            """Returns the complete path for the gridfile if defined in params_ArcticSeas.py (None
               if not), the matching maskfile (declared, or derived from the gridfile name if not),
               and the model_label to use for the plots"""
            model_label = model_label_of(wmodel)
            if model_label not in ArcticSeas_gridfile:
                return None, None, model_label
            gridfile = ArcticSeas_gridfile[model_label]
            if model_label in ArcticSeas_maskfile:
                maskfile = ArcticSeas_maskfile[model_label]
            else:
                # -- No mask file declared for this simulation in ArcticSeas_maskfile: derive
                # -- one automatically from its grid file's name, under ArcticSeas_cache_dir.
                # -- It does not need to exist yet -- ensure_mask() below builds it on demand,
                # -- exactly like any explicitly declared mask file.
                grid_name = os.path.splitext(os.path.basename(gridfile))[0]
                # -- Drop 'mesh_mask'/'meshmask' if present: the grid name is typically
                # -- provided before or after it (e.g. 'mesh_mask.cnrmcm7')
                for token in ('mesh_mask', 'meshmask'):
                    grid_name = grid_name.replace(token, '')
                grid_name = grid_name.strip('._-')
                maskfile = os.path.join(ArcticSeas_cache_dir, 'masks', 'mask.ArcticSeas.%s.nc' % grid_name)
            return gridfile, maskfile, model_label

        def ensure_mask(gridfile, maskfile):
            """Build maskfile from gridfile with create_mask_regions.py if it is missing."""
            if os.path.exists(maskfile):
                return True
            if not os.path.exists(create_mask_regions_script):
                print("ArcticSeas: create_mask_regions.py not found at %s ; cannot build %s"
                      % (create_mask_regions_script, maskfile))
                return False
            mask_dir = os.path.dirname(maskfile)
            if not os.path.isdir(mask_dir):
                os.makedirs(mask_dir)
            try:
                # We need to use a subprocess here because 1. the script uses variable names that
                # could overwrite local variables otherwise, 2. it ends with a sys.exit() which
                # would stop the diagnostics otherwise
                subprocess.run(['python3', create_mask_regions_script,
                                 '--maskfile', gridfile, '--gridfile', gridfile,
                                 '--out', maskfile], check=True)
            except Exception as e:
                print("ArcticSeas: create_mask_regions.py failed for %s -> %s" % (maskfile, e))
            return os.path.exists(maskfile)

        def ensure_check_masks_plot(gridfile, maskfile, label):
            """(Re-)run check_masks.py for this grid if its plot is missing or predates maskfile.
            check_masks.py always writes fixed relative filenames ('check_masks_arctic.png',
            'check_masks_antarctic.png') in its current directory -- run it in workdir and
            immediately rename the Arctic one to <mask file name without .nc>.png, so every
            distinct grid ends up with its own plot name. MPLBACKEND=Agg avoids the script's
            closing plt.show() blocking/failing headless."""
            png = os.path.join(workdir, os.path.splitext(os.path.basename(maskfile))[0] + '.png')
            need = not os.path.exists(png) or os.path.getmtime(png) < os.path.getmtime(maskfile)
            if need and os.path.exists(check_masks_script):
                try:
                    subprocess.run(['python3', check_masks_script,
                                     '--gridfile', gridfile, '--mask', maskfile, '--label', label],
                                    cwd=workdir, check=True, env=dict(os.environ, MPLBACKEND='Agg'))
                    generated = os.path.join(workdir, 'check_masks_arctic.png')
                    if os.path.exists(generated):
                        os.replace(generated, png)
                except Exception as e:
                    print("ArcticSeas: check_masks.py failed for %s -> %s" % (label, e))
            return png if os.path.exists(png) else None
            # WARNING : Remove plot for missing seas

        # -- Find every distinct grid/mask actually needed by the configured simulations,
        # -- then build/locate each of them (deduplicated by mask file)
        # -----------------------------------------------------------------------------------------
        # From MyOwnDiagnostics, not sure we need to copy those
        Wmodels = copy.deepcopy(models)
        # -- Observational references (params_ArcticSeas.py) go through exactly the same
        # -- processing as the simulations; they are only told apart (obs_labels) when
        # -- checking for missing simulations and when plotting
        Wobs = copy.deepcopy(ArcticSeas_obs) if 'ArcticSeas_obs' in dir() and ArcticSeas_obs else []
        # -- How each observation is plotted (see ArcticSeas_obs_dates/ArcticSeas_obs_clim_period
        # -- in params_ArcticSeas.py, overridden by 'obs_dates'/'obs_clim_period' in its own
        # -- entry). These keys are removed from the entry so that they never reach ds()
        default_obs_dates = ArcticSeas_obs_dates if 'ArcticSeas_obs_dates' in dir() else 'climato'
        default_obs_clim_period = ArcticSeas_obs_clim_period if 'ArcticSeas_obs_clim_period' in dir() else 'full'
        obs_settings = dict()  # obs label -> (obs_dates, obs_clim_period)
        for obs in Wobs:
            obs_dates = obs.pop('obs_dates', default_obs_dates)
            obs_clim_period = str(obs.pop('obs_clim_period', default_obs_clim_period))
            if obs_dates not in ('climato', 'real'):
                print("ArcticSeas: obs_dates '%s' for %s should be 'climato' or 'real', using 'climato'"
                      % (obs_dates, model_label_of(obs)))
                obs_dates = 'climato'
            if obs_clim_period != 'full' and not re.match(r'^\d{4}-\d{4}$', obs_clim_period):
                print("ArcticSeas: obs_clim_period '%s' for %s should be 'full' or 'YYYY-YYYY', using 'full'"
                      % (obs_clim_period, model_label_of(obs)))
                obs_clim_period = 'full'
            obs_settings[model_label_of(obs)] = (obs_dates, obs_clim_period)
        obs_labels = [model_label_of(obs) for obs in Wobs]
        for label in obs_labels:
            if label in [model_label_of(model) for model in Wmodels]:
                print("ArcticSeas: observation %s has the same customname as a simulation; "
                      "their per-simulation settings and cache will be mixed up" % label)
        Wdatasets = Wmodels + Wobs
        grids_needed = dict()  # maskfile -> (gridfile, [label, ...]) -- one label per simulation
        # sharing this maskfile, in the order they are encountered in Wdatasets
        for model in Wdatasets:
            gridfile, maskfile, label = resolve_grid_and_mask(model)
            if gridfile is None:
                print("ArcticSeas: no ArcticSeas_gridfile/ArcticSeas_maskfile entry for %s, skipping"
                      % label)
                continue
            if maskfile in grids_needed:
                grids_needed[maskfile][1].append(label)
            else:
                grids_needed[maskfile] = (gridfile, [label])

        available_masks = dict()  # maskfile -> (gridfile, [label, ...]), only those built/found
        for maskfile, (gridfile, labels) in grids_needed.items():
            if ensure_mask(gridfile, maskfile):
                available_masks[maskfile] = (gridfile, labels)
        # Build all required maskfiles only once thanks to the information stored in grids_needed
        # If maskfile properly built fill in available_masks

        if not available_masks:
            index += open_table()
            index += start_line('Error')
            index += "No Arctic seas mask file could be found or built for any of the configured " \
                      "simulations (see the job log for details)."
            index += close_line() + close_table()

        # Finally !!! We have everything we need to run the actual diagnostics
        # --------------------------------------------------------------------------------------------
        else:
            # -- Sanity-check plot of the mask geometry, one per distinct grid (not one per
            # -- simulation: the label already lists every simulation sharing that grid) 
            # -----------------------------------------------------------------------------------------
            index += section("Arctic seas", level=5)
            index += open_table()
            index += start_line('Arctic seas')
            for maskfile, (gridfile, shared_labels) in available_masks.items():
                combined_label = ', '.join(shared_labels)
                png = ensure_check_masks_plot(gridfile, maskfile, combined_label)
                if png:
                    index += cell(combined_label, png, thumbnail=thumbN_size, hover=hover,
                                   **alternative_dir)
                # Options thumbnail, hover, alternative_dir set globally in C-ESM-EP, check those in
                # case of issues with the plotting on the html
            index += close_line() + close_table()

            # -- Get the list of seas from any one of the mask files: they all describe the same
            # -- set of named seas (netcdf variable names, e.g. 'barentse'), just discretized on
            # -- different grids; their 'long_name' attribute (e.g. 'Barents Sea') labels the page
            # -----------------------------------------------------------------------------------------
            reference_maskfile = next(iter(available_masks))
            mask_ds = xr.open_dataset(reference_maskfile)
            available_seas = list(mask_ds.data_vars)
            sea_display_names = {sea: mask_ds[sea].attrs.get('long_name', sea) for sea in available_seas}
            mask_ds.close()
            #
            # -- Determine which seas to plot according to params_ArcticSeas.py
            # -----------------------------------------------------------------------------------------
            if ArcticSeas_seas_list:
                seas = [sea for sea in ArcticSeas_seas_list if sea in available_seas]
                missing_seas = [sea for sea in ArcticSeas_seas_list if sea not in available_seas]
            else:
                seas = available_seas
                missing_seas = []

            # ==> -- For each model/variable, reuse the cached
            # ==> -- per-sea time series from ArcticSeas_cache_dir if it already covers the
            # ==> -- simulation's latest available year, otherwise (re-)compute it with
            # ==> -- comp_seaiceindex.py and update the cache
            # -----------------------------------------------------------------------------------------
            seaindex_files = dict()
            # seaindex_files is set to contain the file names holding the sea ice area / volume  
            # organized by simulation and variable
            model_labels = []
            #
            # ARCTIC_SEAS_VARIABLES = ['sic', 'sit'] -- these are not necessarily the actual
            # netcdf variable name in the netcdf files that varies by simulation (e.g. N3CPL's volume
            # variable is 'sivolu', not 'sit') and is looked up below in ArcticSeas_variable_names
            # (params_ArcticSeas.py), a dict-of-dicts keyed by simulation then by ('sit'/'sic').
            for variable in ARCTIC_SEAS_VARIABLES:
                seaindex_files[variable] = dict()
                # Wdatasets is a list of dictionaries holding information set in datasetsetup.py
                # (simulations) followed by those set in ArcticSeas_obs (observations)
                for model in Wdatasets:
                    wmodel = model.copy()
                    model_label = model_label_of(wmodel)
                    if model_label not in model_labels:
                        model_labels.append(model_label)

                    # -- Netcdf variable name to find in output files for this simulation/variable
                    # -----------------------------------------------------------------------------------------
                    real_variable = ArcticSeas_variable_names.get(model_label, dict()).get(variable)
                    if real_variable is None:
                        print("ArcticSeas: no ArcticSeas_variable_names entry for %s / %s, skipping"
                              % (model_label, variable))
                        continue
                        # We skip to next model here if variable not declared in param_ArcticSeas.py

                    # Climaf needs the variable to find the input netcdf files
                    wmodel.update(dict(variable=real_variable))
                    # -- get_period_manager() with diag='ts' is what actually
                    # -- resolves ts_period into a concrete 'period' (date1-date2) usable by ds()
                    # -- period included in the wmodel dictionary
                    wmodel = get_period_manager(wmodel, diag='ts')

                    # -- when it finds no file for a variable/simulation, get_period_manager()
                    # -- prints 'Error in get_period_manager => No File found for ...'
                    # -- and leaves 'period' unset)
                    if 'period' not in wmodel:
                        print("ArcticSeas: no data found for %s / %s (variable %s), skipping"
                              % (model_label, variable, real_variable))
                        continue
                        # We skip to next model here if period is not available = no input data found

                    # Build a name for the output files (sea ice volume and sea ice area)
                    cache_label = CACHE_VARIABLE_LABELS[(variable, ArcticSeas_meanORsum)]
                    cache_file = os.path.join(ArcticSeas_cache_dir, sanitize(model_label), cache_label + '.nc')
                    # Finding the last simulation year both in the input and output files
                    simulation_latest_year = latest_year_of_wmodel(wmodel)
                    cache_latest_year = latest_year_in_ncfile(cache_file)
                    up_to_date = (cache_latest_year is not None and simulation_latest_year is not None
                                  and cache_latest_year >= simulation_latest_year)
                    # The output sea ice index file does not need to be updated.

                    if up_to_date:
                        seaindex_files[variable][model_label] = cache_file
                        continue
                        # We do not need to compute this diagnostic so we skip to next model

                    # -- Grid and matching sea mask for this simulation specifically (see
                    # -- resolve_grid_and_mask above: may differ from the defaults)
                    gridfile, maskfile, label = resolve_grid_and_mask(wmodel)

                    if gridfile is None:
                        print("ArcticSeas: no ArcticSeas_gridfile/ArcticSeas_maskfile entry for %s, "
                              "skipping" % model_label)
                        continue
                    if not os.path.exists(gridfile):
                        print("ArcticSeas: grid file not found for %s -> %s" % (model_label, gridfile))
                        continue
                        # Here we skip to the next model because we do not need to update the diagnostics
                    if maskfile not in available_masks:
                        print("ArcticSeas: no sea mask available for %s (mask %s), skipping"
                              % (model_label, maskfile))
                        continue

                    try:
                        # Gather in dataset dat all informations about model read in datasetup.py and
                        # the period of the netcdf files computed by get_period_manager
                        dat = ds(**wmodel)
                        if ArcticSeas_annual_mean:
                            dat = ccdo(dat, operator='yearmean')
                            # Determine which command line to run if needed and store it in dat together
                            # with all the information about input files for this model
                        datafile = cfile(dat)
                        # Actual computation of annual means if needed

                        # -- in case sea ice concentration is published as a percentage (0-100) 
                        # -- turn it into a fraction (0-1). Detect and write the corrected copy
                        # -- to our own workdir rather than overwriting CliMAF's cache file.
                        if variable == 'sic':
                            with xr.open_dataset(datafile) as check_ds:
                            # datafile is opened only in this loop, closed automatically when loop closes
                                is_percent = float(check_ds[real_variable].max()) > 1.5
                            if is_percent:
                                fraction_file = os.path.join(
                                    workdir, 'sic_fraction_%s.nc' % sanitize(model_label))
                                with xr.open_dataset(datafile) as check_ds:
                                # idem datafile opened only in this loop
                                    fraction_ds = check_ds.copy()
                                    fraction_ds[real_variable] = fraction_ds[real_variable] / 100.
                                    fraction_ds.to_netcdf(fraction_file)
                                datafile = fraction_file

                        # Create a cache directory to hold the sea ice index output files if not there
                        if not os.path.isdir(os.path.dirname(cache_file)):
                            os.makedirs(os.path.dirname(cache_file))

                        cmd = ['python3', comp_seaiceindex_script,
                               '--data', datafile,
                               '--var', real_variable,
                               '--mask', maskfile,
                               '--grid', gridfile,
                               '--dxvar', ArcticSeas_dxvar,
                               '--dyvar', ArcticSeas_dyvar,
                               '--meanORsum', ArcticSeas_meanORsum,
                               '--out', cache_file]
                        # -- Output captured (and echoed to the log) to keep its 'missing values'
                        # -- warnings next to the cache file, so that they are shown on the atlas
                        # -- page also when the cache is reused later
                        result = subprocess.run(cmd, capture_output=True, text=True)
                        print(result.stdout + result.stderr)
                        result.check_returncode()
                        write_sea_warnings(cache_file, result.stderr)
                        seaindex_files[variable][model_label] = cache_file
                    except Exception as e:
                        print("ArcticSeas: failed to compute %s (variable %s) for %s -> %s"
                              % (variable, real_variable, model_label, e))

            # ==> -- Build one table row per sea; for each sea, one plot for the ice area (sic) and
            # ==> -- one for the ice volume (sit), overlaying all the simulations
            # -----------------------------------------------------------------------------------------
            index += open_table()
            if ArcticSeas_meanORsum == 'sum':
                variable_plot_specs = [
                    dict(variable='sic', title='Sea ice area',
                         ylabel='Sea ice area (Millions km2)'),
                    dict(variable='sit', title='Sea ice volume',
                         ylabel='Sea ice volume (Thousand km3)'),
                ]
            else:
                variable_plot_specs = [
                    dict(variable='sic', title='Sea ice concentration',
                         ylabel='Area-weighted mean sea ice concentration (%)'),
                    dict(variable='sit', title='Sea ice thickness',
                         ylabel='Area-weighted mean sea ice thickness (m)'),
                ]
            # -- comp_seaiceindex.py stores, as an attribute of each sea in its output file, the
            # -- total area of the sea ('sum') or the mean grid cell area within the sea ('mean')
            # -- on the grid of that simulation/observation
            if ArcticSeas_meanORsum == 'sum':
                area_attr, area_divisor, area_title = 'total_area_m2', 1e9, 'Sea area (thousand km2)'
            else:
                area_attr, area_divisor, area_title = 'mean_cell_area_m2', 1e6, 'Mean grid cell area (km2)'
            # -- 'Missing values' warnings of comp_seaiceindex.py, per sea, for every
            # -- simulation/observation and variable plotted
            variable_titles = dict((spec['variable'], spec['title']) for spec in variable_plot_specs)
            warnings_per_sea = dict()  # sea -> ['label (variable title): message', ...]
            for variable, files in seaindex_files.items():
                for model_label, cache_file in files.items():
                    for warned_sea, message in read_sea_warnings(cache_file).items():
                        warnings_per_sea.setdefault(warned_sea, []).append(
                            '%s (%s): %s' % (model_label, variable_titles.get(variable, variable), message))
            for sea in seas:
                sea_name = sea_display_names.get(sea, sea)
                # If sea_display_names is not defined for sea, sea_name is sea
                title = sea_name + '<br><small>' + area_title + ': ' + \
                    sea_areas_text(sea, model_labels, seaindex_files, area_attr, area_divisor) + '</small>'
                if sea in warnings_per_sea:
                    title += '<br><small><font color="red">Warning - ' + \
                             '<br>Warning - '.join(warnings_per_sea[sea]) + '</font></small>'
                index += start_line(title)
                for plot_var in variable_plot_specs:
                    # variable_plot_specs is a list of dictionaries holding plot titles and ylabels
                    variable = plot_var['variable']

                    # -- Gather the data available for this sea/variable first, so we know
                    # -- which simulations (if any) are missing before deciding whether to
                    # -- plot the available ones or to report an error instead
                    # -----------------------------------------------------------------------------------------
                    curves = []  # (model_label, time_values, data_values) for simulations with data
                    obs_series = dict()  # obs label -> per-sea DataArray, for observations with data
                    present_labels = []
                    for model_label, outfile in seaindex_files.get(variable, dict()).items():
                        if not os.path.exists(outfile):
                            continue
                        out_ds = xr.open_dataset(outfile)
                        if sea in out_ds.data_vars:
                            da = out_ds[sea].load()
                            time_dim = da.dims[0]
                            if variable in ('sic', 'sit') and ArcticSeas_meanORsum == 'sum':
                                # -- sea ice area (m^2 -> millions km^2) or volume (m^3 ->
                                # -- thousand km^3), to match the plot's ylabel; purely local
                                # -- to this plot, the cache file on disk is left untouched
                                da = da / 1e12
                            if model_label in obs_labels:
                                obs_series[model_label] = da
                            else:
                                curves.append((model_label, da[time_dim].values, da.values))
                            present_labels.append(model_label)
                        out_ds.close()

                    # -- Observations: either their actual dates ('real'), or their climatology
                    # -- over obs_clim_period repeated along the simulation dates ('climato').
                    # -- Without any simulation to follow, a climatology falls back to real dates
                    # -----------------------------------------------------------------------------------------
                    obs_curves = []  # (obs label, legend label, time_values, data_values)
                    sim_times = None
                    if curves:
                        sim_times = np.unique(np.concatenate([time_values for _, time_values, _ in curves]))
                    for obs_label, da in obs_series.items():
                        obs_dates, obs_clim_period = obs_settings[obs_label]
                        time_dim = da.dims[0]
                        if obs_dates == 'real' or sim_times is None:
                            obs_curves.append((obs_label, obs_label, da[time_dim].values, da.values))
                            continue
                        clim = obs_climatology_curve(da, obs_clim_period, ArcticSeas_annual_mean, sim_times)
                        if clim is None:
                            print("ArcticSeas: no %s data within obs_clim_period %s for %s, not plotted"
                                  % (obs_label, obs_clim_period, sea))
                            continue
                        values, first_year, last_year = clim
                        obs_curves.append((obs_label, '%s (clim %d-%d)' % (obs_label, first_year, last_year),
                                           sim_times, values))
                    # model_labels holds every simulation ArcticSeas attempted (see the compute
                    # loop above); any of them missing from present_labels has no data here.
                    # Observations are left out: a missing observation is just not plotted
                    missing_labels = [m for m in model_labels
                                      if m not in present_labels and m not in obs_labels]

                    fig, ax = plt.subplots(figsize=(6, 4))
                    # New figure created for each sea, each variable and associated diagnostics
                    if not curves and not obs_curves:
                        # -- Always an error message here, regardless of ArcticSeas_on_missing_simulations
                        ax.text(0.5, 0.5, 'No data available for any simulation',
                                ha='center', va='center', wrap=True, transform=ax.transAxes)
                    elif missing_labels and ArcticSeas_on_missing_simulations == 'error':
                        ax.text(0.5, 0.5, 'Data missing for: %s' % ', '.join(missing_labels),
                                ha='center', va='center', wrap=True, transform=ax.transAxes)
                    else:
                        # -- ArcticSeas_on_missing_simulations == 'partial' (or nothing missing):
                        # -- plot the simulations that do have data, silently skipping the rest
                        # -- Simulations first (default color cycle), then observations on top
                        # -- with their own fixed style (same style for a given observation
                        # -- on every plot, from its rank in ArcticSeas_obs)
                        for model_label, time_values, data_values in curves:
                            ax.plot(time_values, data_values, lw=1.5, label=model_label)
                        for obs_label, legend_label, time_values, data_values in obs_curves:
                            color, linestyle = OBS_STYLES[obs_labels.index(obs_label) % len(OBS_STYLES)]
                            ax.plot(time_values, data_values, lw=2, color=color, linestyle=linestyle,
                                    label=legend_label, zorder=3)
                        ax.legend(fontsize=8)

                    ax.set_title("%s - %s" % (sea_name, plot_var['title']))
                    ax.set_xlabel('Time')
                    ax.set_ylabel(plot_var['ylabel'])
                    fig.tight_layout()

                    png_path = os.path.join(workdir, 'ts_%s_%s.png' % (sanitize(sea), variable))
                    fig.savefig(png_path, dpi=100)
                    plt.close(fig)

                    index += cell(plot_var['title'], png_path, thumbnail=thumbN_size, hover=hover, **alternative_dir)
                index += close_line()
            index += close_table()

            if missing_seas:
                index += open_table()
                index += start_line('Note')
                index += "Seas requested in ArcticSeas_seas_list but not found in the mask file " \
                          "%s: %s" % (reference_maskfile, ", ".join(missing_seas))
                index += close_line() + close_table()

# -----------------------------------------------------------------------------------
# --   End
# -----------------------------------------------------------------------------------
