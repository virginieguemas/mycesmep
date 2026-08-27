# ------------------------------------------------------------------------------------------------------ \
# --                                                                                                    - \
# --      Scientific diagnostics for the                                                                  - \
# --          CliMAF Earth System Model Evaluation Platform                                               - |
# --      diagnostics_ArcticSeas.py                                                                       - |
# --                                                                                                      - |
# --      Time series of sea ice area (siconc) and sea ice volume (sithick) averaged                     - |
# --      over each Arctic sea, computed with the external script                                        - |
# --      comp_seaiceindex.py (sea_ice_diag_tools repository) applied to a mask file                      - |
# --      describing the individual seas.                                                                - |
# --                                                                                                      - /
# ---------------------------------------------------------------------------------------------------- /

from os import getcwd
import os
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


if do_ArcticSeas_timeseries:
    # ==> -- Open the section
    # -----------------------------------------------------------------------------------------
    index += section("Sea ice area and volume time series per Arctic sea", level=4)

    if thumbnail_size:
        thumbN_size = thumbnail_size
    elif 'ArcticSeas_thumbnail_size' in dir() and ArcticSeas_thumbnail_size:
        thumbN_size = ArcticSeas_thumbnail_size
    else:
        thumbN_size = thumbnail_size_global

    workdir = os.path.join(getcwd(), 'ArcticSeas_workfiles')
    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    comp_seaiceindex_script = os.path.join(ArcticSeas_tools_dir, 'comp_seaiceindex.py')

    if not os.path.exists(comp_seaiceindex_script):
        index += open_table()
        index += start_line('Error')
        index += "comp_seaiceindex.py not found at %s ; check ArcticSeas_tools_dir in params_ArcticSeas.py" \
                  % comp_seaiceindex_script
        index += close_line() + close_table()

    elif not os.path.exists(ArcticSeas_maskfile):
        index += open_table()
        index += start_line('Error')
        index += "Mask file not found: %s ; check ArcticSeas_maskfile in params_ArcticSeas.py" % ArcticSeas_maskfile
        index += close_line() + close_table()

    else:
        # -- Get the list of seas from the mask file (long_name if available, else var name)
        # -----------------------------------------------------------------------------------------
        mask_ds = xr.open_dataset(ArcticSeas_maskfile)
        available_seas = list(mask_ds.data_vars)
        mask_ds.close()

        if ArcticSeas_seas_list:
            seas = [sea for sea in ArcticSeas_seas_list if sea in available_seas]
            missing_seas = [sea for sea in ArcticSeas_seas_list if sea not in available_seas]
        else:
            seas = available_seas
            missing_seas = []

        Wmodels = copy.deepcopy(models)

        # ==> -- Apply the period_manager and materialize a netcdf file per model/variable, then
        # ==> -- run comp_seaiceindex.py on it to get the per-sea time series
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

                try:
                    dat = ds(**wmodel)
                    if ArcticSeas_annual_mean:
                        dat = ccdo(dat, operator='yearmean')
                    datafile = cfile(dat)

                    outfile = os.path.join(
                        workdir, 'seaindex_%s_%s.nc' % (variable, sanitize(model_label)))

                    cmd = ['python3', comp_seaiceindex_script,
                           '--maskfile', ArcticSeas_maskfile,
                           '--datafile', datafile,
                           '--variable', variable,
                           '--outfile', outfile]
                    subprocess.run(cmd, check=True)
                    seaindex_files[variable][model_label] = outfile
                except Exception as e:
                    print("ArcticSeas: failed to compute %s for %s -> %s" % (variable, model_label, e))

        # ==> -- Build one table row per sea; for each sea, one plot for the ice area (siconc) and
        # ==> -- one for the ice volume (sithick), overlaying all the simulations
        # -----------------------------------------------------------------------------------------
        index += open_table()
        variable_plot_specs = [
            dict(variable='siconc', title='Sea ice area', ylabel='Mean sea ice concentration'),
            dict(variable='sithick', title='Sea ice volume', ylabel='Mean sea ice thickness (m)'),
        ]
        for sea in seas:
            index += start_line(sea)
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

                ax.set_title("%s - %s" % (sea, spec['title']))
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
