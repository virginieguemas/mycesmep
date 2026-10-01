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
    return "".join(c if c.isalnum() else "_" for c in name)


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

# -- sia/siv = sea ice area/volume (area-weighted sum)
# -- sic/sit = sea ice concentration/thickness (area-weighted mean)
CACHE_VARIABLE_LABELS = {
    ('sic', 'sum'): 'sia',
    ('sit', 'sum'): 'siv',
    ('sic', 'mean'): 'sic',
    ('sit', 'mean'): 'sit',
}


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
        grids_needed = dict()  # maskfile -> (gridfile, [label, ...]) -- one label per simulation
        # sharing this maskfile, in the order they are encountered in Wmodels
        for model in Wmodels:
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
                # Wmodels is a list of dictionaries holding information set in datasetsetup.py
                for model in Wmodels:
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
                        subprocess.run(cmd, check=True)
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
            for sea in seas:
                sea_name = sea_display_names.get(sea, sea)
                # If sea_display_names is not defined for sea, sea_name is sea
                index += start_line(sea_name)
                for plot_var in variable_plot_specs:
                    # variable_plot_specs is a list of dictionaries holding plot titles and ylabels
                    variable = plot_var['variable']

                    # -- Gather the data available for this sea/variable first, so we know
                    # -- which simulations (if any) are missing before deciding whether to
                    # -- plot the available ones or to report an error instead
                    # -----------------------------------------------------------------------------------------
                    curves = []  # (model_label, time_values, data_values) for simulations with data
                    present_labels = []
                    for model_label, outfile in seaindex_files.get(variable, dict()).items():
                        if not os.path.exists(outfile):
                            continue
                        out_ds = xr.open_dataset(outfile)
                        if sea in out_ds.data_vars:
                            da = out_ds[sea]
                            time_dim = da.dims[0]
                            data_values = da.values
                            if variable in ('sic', 'sit') and ArcticSeas_meanORsum == 'sum':
                                # -- sea ice area (m^2 -> millions km^2) or volume (m^3 ->
                                # -- thousand km^3), to match the plot's ylabel; purely local
                                # -- to this plot, the cache file on disk is left untouched
                                data_values = data_values / 1e12
                            curves.append((model_label, out_ds[time_dim].values, data_values))
                            present_labels.append(model_label)
                        out_ds.close()
                    # model_labels holds every simulation ArcticSeas attempted (see the compute
                    # loop above); any of them missing from present_labels has no data here
                    missing_labels = [m for m in model_labels if m not in present_labels]

                    fig, ax = plt.subplots(figsize=(6, 4))
                    # New figure created for each sea, each variable and associated diagnostics
                    if not curves:
                        # -- Always an error message here, regardless of ArcticSeas_on_missing_simulations
                        ax.text(0.5, 0.5, 'No data available for any simulation',
                                ha='center', va='center', wrap=True, transform=ax.transAxes)
                    elif missing_labels and ArcticSeas_on_missing_simulations == 'error':
                        ax.text(0.5, 0.5, 'Data missing for: %s' % ', '.join(missing_labels),
                                ha='center', va='center', wrap=True, transform=ax.transAxes)
                    else:
                        # -- ArcticSeas_on_missing_simulations == 'partial' (or nothing missing):
                        # -- plot the simulations that do have data, silently skipping the rest
                        for model_label, time_values, data_values in curves:
                            ax.plot(time_values, data_values, lw=1.5, label=model_label)
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
