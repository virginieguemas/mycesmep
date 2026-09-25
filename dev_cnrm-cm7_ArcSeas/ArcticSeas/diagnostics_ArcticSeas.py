# ------------------------------------------------------------------------------------------------------ \
# --                                                                                                    - \
# --      Scientific diagnostics for the                                                                  - \
# --          CliMAF Earth System Model Evaluation Platform                                               - |
# --      diagnostics_ArcticSeas.py                                                                       - |
# --        ==> add html code to 'index' (initialized with 'header')                                      - |
# --            using the CliMAF html toolbox (start_line, cell, close_table... )                         - |
# --            to create the Arctic Seas atlas page                                                      - |
# --                                                                                                      - |
# --      Time series of sea ice area (siconc) and sea ice volume (sithick) computed as an                - |
# --      area-weighted sum over each Arctic sea, using the external script comp_seaiceindex.py           - |
#         (sea_ice_diag_tools repository) applied to a mask file describing the individual seas and       - |
#         a grid file giving the cell areas. If not provided by the user, the mask file for the           - |
#         individual seas is computed by create_mask_regions.py from the sea_ice_diag_tools repository    - |
# --                                                                                                      - /
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

# -- Init html index
# -----------------------------------------------------------------------------------
index = header(atlas_head_title, style_file=style_file)


def sanitize(name):
    """Turn a sea/model name into a safe token for file names."""
    return "".join(c if c.isalnum() else "_" for c in name)


def latest_year_of_wmodel(wmodel):
    """Latest year covered by a period-managed model dict, from build_str_period()."""
    years = re.findall(r'\d{4}', str(build_str_period(wmodel)))
    return int(years[-1]) if years else None


def latest_year_in_ncfile(ncfile):
    """Latest year found in the (single) time-like dimension of a netcdf file."""
    if not os.path.exists(ncfile):
        return None
    try:
        cached = xr.open_dataset(ncfile)
        first_var = next(iter(cached.data_vars), None)
        if first_var is None:
            cached.close()
            return None
        time_dim = cached[first_var].dims[0]
        times = cached[time_dim].values
        cached.close()
        return int(str(times[-1])[:4]) if len(times) else None
    except Exception:
        return None


if do_ArcticSeas_timeseries:
    # 
    # ==> -- Open the section and an html file
    # -----------------------------------------------------------------------------------------
    # WARNING : Different from MyOwnDiagnostics
    index += section("Sea ice area and volume per Arctic sea", level=4)
    #
    # ==> -- Control the size of the thumbnail -> thumbN_size
    # Note: ArcticSeas_thumbnail_size defined in params_ArcticSeas.py
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
    #
    # Ecriture des messages d'erreur dans la page html
    # -----------------------------------------------------------------------------------------
    if not os.path.exists(comp_seaiceindex_script):
        index += open_table()
        index += start_line('Error')
        index += "comp_seaiceindex.py not found at %s ; check ArcticSeas_tools_dir in params_ArcticSeas.py" \
                  % comp_seaiceindex_script
        index += close_line() + close_table()

    elif not os.path.exists(ArcticSeas_maskfile) and not os.path.exists(create_mask_regions_script):
        index += open_table()
        index += start_line('Error')
        index += "Mask file not found: %s, and create_mask_regions.py not found at %s to build it ; " \
                  "check ArcticSeas_maskfile and ArcticSeas_tools_dir in params_ArcticSeas.py" \
                  % (ArcticSeas_maskfile, create_mask_regions_script)
        index += close_line() + close_table()

    else:
        mask_dir = os.path.dirname(ArcticSeas_maskfile)

        # -- Build the mask file if it is not there yet: create_mask_regions.py takes no
        # -- argument and writes its output ('mask.ArcticSeas.cnrmcm7.nc') as a relative
        # -- path in its current directory, so it is run from ArcticSeas_maskfile's directory
        # -----------------------------------------------------------------------------------------
        if not os.path.exists(ArcticSeas_maskfile):
            if not os.path.isdir(mask_dir):
                os.makedirs(mask_dir)
            try:
                subprocess.run(['python3', create_mask_regions_script], cwd=mask_dir, check=True)
            except Exception as e:
                print("ArcticSeas: create_mask_regions.py failed -> %s" % e)
        #
        # WARNING : code to be generalized so that it takes command line arguments for the different 
        # potential grids

        if not os.path.exists(ArcticSeas_maskfile):
            index += open_table()
            index += start_line('Error')
            index += "Mask file still missing after running create_mask_regions.py: %s" % ArcticSeas_maskfile
            index += close_line() + close_table()
 
        # Finally !!! We have everything we need to run the actual diagnostics
        # --------------------------------------------------------------------------------------------
        else:
            # -- Use check_masks.py from the sea_ice_diag_tools repository to draw
            # -- a map colouring each named Arctic sea. check_masks.py also produces an Antarctic
            # -- map (same script, unrelated hemisphere), which is left out of this Arctic-only
            # -- page. Only (re-)run it when there is no plot yet or it predates the mask file
            # -- (e.g. the mask was just rebuilt above); MPLBACKEND=Agg avoids the script's
            # -- closing plt.show() blocking/failing headless.
            # -----------------------------------------------------------------------------------------
            check_masks_script = os.path.join(ArcticSeas_tools_dir, 'masks', 'check_masks.py')
            check_masks_arctic_png = os.path.join(mask_dir, 'check_masks_arctic.png')
            need_check_plot = (not os.path.exists(check_masks_arctic_png)
                                or os.path.getmtime(check_masks_arctic_png) < os.path.getmtime(ArcticSeas_maskfile))
            # WARNING : Need to simplify by removing the plt.show() in check_masks.py and the
            #           MPLPACKEND=Agg here
            # WARNING : Same generalization as create_mask_regions.py needed
            # WARNING : Penser à ne pas plotter les mers manquantes

            if need_check_plot and os.path.exists(check_masks_script):
                try:
                    subprocess.run(['python3', check_masks_script], cwd=mask_dir, check=True,
                                    env=dict(os.environ, MPLBACKEND='Agg'))
                except Exception as e:
                    print("ArcticSeas: check_masks.py failed -> %s" % e)
          
            #
            # ==> -- Add the Arctic seas plot to the html page
            # -----------------------------------------------------------------------------------------
            if os.path.exists(check_masks_arctic_png):
                index += section("Arctic seas", level=5)
                index += open_table()
                index += start_line('Arctic seas')
                index += cell('Arctic seas', check_masks_arctic_png, thumbnail=thumbN_size, hover=hover,
                               **alternative_dir)
                index += close_line() + close_table()

            # -- Get the list of seas from the mask file: netcdf variable names (e.g. 'barentse')
            # -- are used to index the data, their 'long_name' attribute (e.g. 'Barents Sea') to
            # -- label the html page and the plots
            # -----------------------------------------------------------------------------------------
            mask_ds = xr.open_dataset(ArcticSeas_maskfile)
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
            #
            # ==> -- Apply the period_for_diag_manager 
            # -----------------------------------------------------------------------------------------
            # WARNING : From MyOwnDiagnostics, check whether we need that
            Wmodels = copy.deepcopy(models)

            # ==> -- For each model/variable, reuse the cached
            # ==> -- per-sea time series from ArcticSeas_cache_dir if it already covers the
            # ==> -- simulation's latest available year, otherwise (re-)compute it with
            # ==> -- comp_seaiceindex.py and update the cache
            # -----------------------------------------------------------------------------------------
            seaindex_files = dict()
            model_labels = []
            for variable in ArcticSeas_variables:
                seaindex_files[variable] = dict()
                for model in Wmodels:
                    wmodel = model.copy()
                    wmodel.update(dict(variable=variable))
                    wmodel = get_period_manager(wmodel, diag='ArcticSeas')

                    model_label = wmodel.get('customname', build_plot_title(wmodel, None))
                    if model_label not in model_labels:
                        model_labels.append(model_label)

                    cache_file = os.path.join(ArcticSeas_cache_dir, sanitize(model_label), variable + '.nc')

                    simulation_latest_year = latest_year_of_wmodel(wmodel)
                    cache_latest_year = latest_year_in_ncfile(cache_file)
                    up_to_date = (cache_latest_year is not None and simulation_latest_year is not None
                                  and cache_latest_year >= simulation_latest_year)

                    if up_to_date:
                        seaindex_files[variable][model_label] = cache_file
                        continue

                    # -- Prefer a per-model grid description file (models may not share the
                    # -- same grid/resolution) over the global fallback from params_ArcticSeas.py
                    gridfile = wmodel.get('mesh_hgr', wmodel.get('gridfile', ArcticSeas_gridfile))

                    if not os.path.exists(gridfile):
                        print("ArcticSeas: grid file not found for %s -> %s" % (model_label, gridfile))
                        continue

                    try:
                        dat = ds(**wmodel)
                        if ArcticSeas_annual_mean:
                            dat = ccdo(dat, operator='yearmean')
                        datafile = cfile(dat)

                        if not os.path.isdir(os.path.dirname(cache_file)):
                            os.makedirs(os.path.dirname(cache_file))

                        cmd = ['python3', comp_seaiceindex_script,
                               '--data', datafile,
                               '--var', variable,
                               '--mask', ArcticSeas_maskfile,
                               '--grid', gridfile,
                               '--dxvar', ArcticSeas_dxvar,
                               '--dyvar', ArcticSeas_dyvar,
                               '--meanORsum', ArcticSeas_meanORsum,
                               '--out', cache_file]
                        subprocess.run(cmd, check=True)
                        seaindex_files[variable][model_label] = cache_file
                    except Exception as e:
                        print("ArcticSeas: failed to compute %s for %s -> %s" % (variable, model_label, e))

            # ==> -- Build one table row per sea; for each sea, one plot for the ice area (siconc) and
            # ==> -- one for the ice volume (sithick), overlaying all the simulations
            # -----------------------------------------------------------------------------------------
            index += open_table()
            if ArcticSeas_meanORsum == 'sum':
                variable_plot_specs = [
                    dict(variable='siconc', title='Sea ice area',
                         ylabel='Sea ice area (millions km2)'),
                    dict(variable='sithick', title='Sea ice volume',
                         ylabel='Sea ice volume (Thousand km3)'),
                ]
            else:
                variable_plot_specs = [
                    dict(variable='siconc', title='Sea ice concentration',
                         ylabel='Area-weighted mean sea ice concentration'),
                    dict(variable='sithick', title='Sea ice thickness',
                         ylabel='Area-weighted mean sea ice thickness (m)'),
                ]
            for sea in seas:
                sea_name = sea_display_names.get(sea, sea)
                index += start_line(sea_name)
                for spec in variable_plot_specs:
                    variable = spec['variable']
                    fig, ax = plt.subplots(figsize=(6, 4))
                    has_curve = False
                    for model_label, outfile in seaindex_files.get(variable, dict()).items():
                        if not os.path.exists(outfile):
                            continue
                        out_ds = xr.open_dataset(outfile)
                        if sea in out_ds.data_vars:
                            da = out_ds[sea]
                            time_dim = da.dims[0]
                            ax.plot(out_ds[time_dim].values, da.values, lw=1.5, label=model_label)
                            has_curve = True
                        out_ds.close()

                    ax.set_title("%s - %s" % (sea_name, spec['title']))
                    ax.set_xlabel('Time')
                    ax.set_ylabel(spec['ylabel'])
                    if has_curve:
                        ax.legend(fontsize=8)
                    else:
                        ax.text(0.5, 0.5, 'no data', ha='center', va='center', transform=ax.transAxes)
                    fig.tight_layout()

                    png_path = os.path.join(workdir, 'ts_%s_%s.png' % (sanitize(sea), variable))
                    fig.savefig(png_path, dpi=100)
                    plt.close(fig)

                    index += cell(spec['title'], png_path, thumbnail=thumbN_size, hover=hover, **alternative_dir)
                index += close_line()
            index += close_table()

            if missing_seas:
                index += open_table()
                index += start_line('Note')
                index += "Seas requested in ArcticSeas_seas_list but not found in the mask file " \
                          "%s: %s" % (ArcticSeas_maskfile, ", ".join(missing_seas))
                index += close_line() + close_table()

# -----------------------------------------------------------------------------------
# --   End
# -----------------------------------------------------------------------------------
