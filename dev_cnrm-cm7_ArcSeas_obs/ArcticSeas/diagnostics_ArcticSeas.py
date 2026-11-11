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
# -- The cached series are always monthly (seasonal means are computed when plotting,
# -- see seasonal_mean_complete_years)
CACHE_VARIABLE_LABELS = {
    ('sic', 'sum'): 'sia',
    ('sit', 'sum'): 'siv',
    ('sic', 'mean'): 'sic',
    ('sit', 'mean'): 'sit',
}

# -- (color, linestyle) for the observations (ArcticSeas_obs), so that
# -- they stand out from the simulations (colors of SIM_COLORS, solid lines)
OBS_STYLES = [('black', '-'), ('black', '--'), ('black', ':'), ('black', '-.'),
              ('grey', '-'), ('grey', '--'), ('grey', ':'), ('grey', '-.')]

# -- Colors of the simulations (solid lines): the default matplotlib colors without its grey
# -- (kept for the observations). A simulation gets the color of its rank among all the
# -- simulations of datasets_setup.py, so that it keeps the same color on every plot, even
# -- when other simulations are missing from some of them (colors repeat beyond 9 simulations)
SIM_COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
              'tab:brown', 'tab:pink', 'tab:olive', 'tab:cyan']


def latest_year_of_wmodel(wmodel):
    """Latest year covered by an experiment or observation described in wmodel dictionary after its 
    period has been detected by get_period_manager. This function relies on build_period_str()."""
    # build_period_str return a chain of characters holding the first and last dates
    # then re.findall(r'\d{4}' finds the groups of 4 consecutive digits, months are dropped
    # at this stage because there are only 2 digits.
    # WARNING : this function would not work if CLIMAF starts handling daily data
    years = re.findall(r'\d{4}', str(build_period_str(wmodel)))
    return int(years[-1]) if years else None


def latest_year_in_ncfile(ncfile):
    """Latest year of a netcdf file whose variables have time as first dimension (e.g. the
    per-sea cache files written by comp_seaiceindex.py): last value of the first dimension
    of its first variable."""
    if not os.path.exists(ncfile):
        return None
    # If there is no netcdf file yet, this function returns None
    try:
        cached = xr.open_dataset(ncfile)
        first_var = next(iter(cached.data_vars), None)
        # iter iterates on the keys of cached.data_vars dictionary and next takes the
        # first, if there is none, first_var becomes None (very unlikely but in case ...)
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


def report_failed_comput(message_failed_comput, message):
    """Prints a message (why a diagnostic (sea ice area/volume) could not be computed from
    a simulation or observation) in the log and stores it in the list of messages named 
    message_failed_comput, to be shown on the atlas page. A message already
    in the list (e.g. same missing grid file for sic and sit) is not added twice."""
    print("ArcticSeas: " + message)
    if message not in message_failed_comput:
        message_failed_comput.append(message)


def warnings_file_of(cache_file):
    """Defines a filename for the text file, next to a per-sea cache file, holding the 
    'missing values' warnings issued by comp_seaiceindex.py when it computed that cache 
    file (one per line, may be empty).
    This file allows to store information that will be printed on the Atlas html page later."""
    return cache_file + '.warnings'


def write_sea_warnings(cache_file, stderr):
    """From the standard error of comp_seaiceindex.py, all the 'missing values' warnings 
    related to a cache_file are stored into a text file named cache_file.warnings."""
    with open(warnings_file_of(cache_file), 'w') as f:
    # A cache_file_name.warnings text file is opened, overwritten if there was one before
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
    # If the cache_file_name.warnings file exists, we read it
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
                    # If we find an occurrence, we add an element to the dictionary
    return sea_warnings


def sea_areas_text(sea, label_groups, seaindex_files, area_attr, area_divisor):
    """ From an input sea (ex: centrarc), and for groups of experiment and observation
    labels (typically customname) sharing the same grid and sea mask file, this function returns
    a line to be printed later on the html atlas page, listing the area of that sea once per
    group (e.g. 'REF1, REF3: 1424; nsidc_cdr_v6: 1376'): the area is not necessarily the same
    from one grid to another, but it is the same for all the datasets of a group. It needs a
    dictionary containing for each variable and each label, its cache file : seaindex_files.
    It also needs the name of the attribute in the cache file in which the sea area is stored :
    'total_area_m2' or 'mean_cell_area_m2' in the outputs from comp_seaiceindex.py at the
    time of writing. area_divisor allows to convert to a printable unit. """
    parts = []
    # We define a list which will hold the parts of the final line to print
    for labels in label_groups:
        # We loop on the groups of experiments and observations sharing a grid
        area = None
        for label in labels:
            # The area is read in the cache of the first dataset of the group which has one
            for variable in ARCTIC_SEAS_VARIABLES:
                # We loop on the Arctic variables just in case some are missing because it is the
                # same land/sea mask for all variables of a given experiment but observations often
                # have only one of the variables
                outfile = seaindex_files.get(variable, dict()).get(label)
                # seaindex_files = { variable : { model_label, cache_file } , variable : ... }
                # We get the cache file corresponding to that label and that variable, None if there
                # isn't any cache. We need to use get to avoid a crash when a variable isn't there.
                if outfile is None or not os.path.exists(outfile):
                    continue # If we don't have any cache file, we go to next variable
                with xr.open_dataset(outfile) as out_ds:
                    if sea in out_ds.data_vars and area_attr in out_ds[sea].attrs:
                    # We look for the area_attr attribute in the cache file we found
                    # It should contain the area of the sea on that grid
                        area = float(out_ds[sea].attrs[area_attr]) / area_divisor
                    # area_divisor is to convert to Thousand km2 or km2
                if area is not None:
                    break
                # If we found the area of the sea for one variable, we skip the next variable
            if area is not None:
                break
            # If we found the area of the sea for one dataset, we skip the next datasets of the group
        parts.append('%s: %s' % (', '.join(labels), 'unknown' if area is None else '%.0f' % area))
        # We add to parts the labels of the group followed by the area if found, unknown if not
    return '; '.join(parts)


def chunked_in_time(da, chunk_time):
    """da read by pieces of chunk_time time steps along its first dimension (time) with dask,
    so that operations on it never load the whole field in RAM; da unchanged if
    chunk_time is None/0, if da has no time dimension (2 dimensions or less, e.g. a single
    (y, x) field) or if dask is not available."""
    if not chunk_time or da.ndim <= 2:
    # da.ndim provides the number of dimensions: with 2 dimensions or less, da is at most a
    # single (y, x) field without any time dimension to split along (and with 0 dimension,
    # da.dims[0] after try would crash)
        return da
    try:
        return da.chunk({da.dims[0]: chunk_time})
    except (ImportError, ValueError) as e:
        print("ArcticSeas: data could not be chunked (is dask installed ?), read at once -> %s" % e)
        return da


# -- Initials of the 12 calendar months, January first
MONTH_INITIALS = 'JFMAMJJASOND'


def season_months(season):
    """Returns a list of numbers from 1 to 12 indicating the position in the year of each month 
    in the argument season, in the order they appear in the argument season, except for
    'ANM' which returns [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] (ANM = annual means for 
    consistency with the other diagnostics of CESMEP). season can be across the
    end of the year: 'DJF' -> [12, 1, 2], 'JJA' -> [6, 7, 8], 'NDJFMA' -> [11, 12, 1, 2, 3, 4],
    'JASONDJFMAMJ' -> July to June. The initials must match one single run of consecutive
    months: 'J', 'M', 'A' alone are ambiguous. Raises ValueError otherwise."""
    season = season.upper() # Turn to upper cases if not already
    if season == 'ANM': # Treated first because it is the only valid season which is not
        # passed through initials of months
        return list(range(1, 13))
    if not 1 <= len(season) <= 12:
        raise ValueError("season '%s' should be 'ANM' or 1 to 12 month initials" % season)
        # The raise ValueError allows to send an error message to the main script which will
        # capture it and print in the log that the season requested by the user was not
        # understood and annual means will be used instead
    # starts holds all possible starting months of 'season' from O for January to 11 for 
    # December (0 to 11 to be python positions in MONTH_INITIALS - different from the
    # months that will be output and range between 1 and 12)
    # Initially, every month is a possible start; then, letter by letter of 'season', we 
    # only keep those whose i-th following month has the i-th letter
    # of season as initial ((start + i) % 12 goes on from December to January).
    # At the end, starts should hold a single value, but it can hold several ('J', 'M' or 'A'
    # alone: ambiguous) or none (initials not following each other in the calendar)
    #
    starts = list(range(12))              
         # initially every month from J to D could be the beginning of 'season'                                        
    for i in range(len(season)):          # i-th letter of season
        kept = []
        for start in starts:              # possible starting months still left
            if MONTH_INITIALS[(start + i) % 12] == season[i]:
            # We find all occurrences of letter 0 of 'season' in MONTH_INITIALS
            # We only keep those positions in starts
            # Then we find all occurrences of letter 1 of 'season' right after
            # an occurrence of letter 0 of 'season' within MONTH_INTIALS and so on
            # Crossing the end of the year is allowed thanks to % 12
                kept.append(start)        
        starts = kept
    if len(starts) != 1:
        # We should end up with only one possibility as a first month for 'season'
        # within MONTH_INITIALS
        raise ValueError("season '%s' is %s" % (season, 'ambiguous' if starts else
                                                 'not a list of consecutive month initials'))
        # Otherwise, we will send an error message to the main script to be printed
        # in the log together with the information that annual means will be used
    return [(starts[0] + i) % 12 + 1 for i in range(len(season))]
    # Once we know for sure which is the number (from 1 to 12) of the first month
    # of 'season', we can list the numbers of all months of the season
    # +1 is because we list the months as 1 to 12 whereas position in python are indexed
    # from 0 (in MONTH_INITIALS)


def without_duplicate_dates(da):
    """Removes the repeated dates of a time series da (time as first dimension). They appear
    when the input files of a simulation overlap (e.g. two sets of files covering the same
    years, concatenated by CliMAF): each year would then hold 24 months instead of 12, and
    seasonal_mean_complete_years would drop it as incomplete. Only the first value of each
    date is kept, and the dates are sorted.
    Returns a tuple: (da without repeated dates, number of values removed, True if the
    removed values were identical to the kept ones)."""
    time_dim = da.dims[0]
    # da has only one dimension : a time dimension
    first = da.drop_duplicates(time_dim, keep='first')
    # first keeps only one value per date, the first found in da
    n_removed = da.sizes[time_dim] - first.sizes[time_dim]
    # we count the number of values removed
    if n_removed == 0:
        return da, 0, True
        # Nothing to remove: da is returned as it is
    last = da.drop_duplicates(time_dim, keep='last')
    # last keeps only one value per date, the last found in da
    identical = bool(((first == last) | (first.isnull() & last.isnull())).all())
    # first == last matches the values date by date and returns an array of True/False if values
    # are the same/different for each date, if first and last both have NaN for a given
    # date, this first test would return False, hence we need a second test to check whether
    # both are NaN with isnull 
    # we obtain an array of True / False if values are the same (including both NaN) / different
    # then .all() returns True if all values of this array are True
    # we change the True / False result from an array to a simple boolean with the bool function
    # identical = True means each single value of last and first is the same
    return first.sortby(time_dim), n_removed, identical


def seasonal_mean_complete_years(da, months):
    """From a monthly time series data array da (time as first dimension), and a
    list of calendar months (between 1 and 12) within the season, this function computes
    seasonal averages and returns a time series in which the time dimension holds one 
    value per year, along a 'year' dimension (integer years). A season with missing months
    (start/end of a series, data gap) is dropped rather than giving a mean skewed towards
    the months present. A season across the end of the year (e.g. DJF) is assigned to the
    year of its last month (DJF 1980 = December 1979 to February 1980)."""
    time_dim = da.dims[0] # time should be the first dimension
    da = da.where(da[time_dim].dt.month.isin(months), drop=True)
    if da.sizes[time_dim] == 0:
        return xr.DataArray(np.array([]), dims='year', coords=dict(year=np.array([], dtype=int)))
        # If there is no occurrence of the months requested by the user in the input da,
        # we return an empty xarray of seasonal means
    # Months coming before January in the season (e.g. December in DJF) should be associated and
    # averaged within the season of the next year
    next_year_months = months[:months.index(1)] if 1 in months else []
    # 'if 1 in months' checks whether January is part of the season, in that case, the season might
    # cross the end of the year and we need to reattribute some months before January to next year
    # [:months.index(1)] takes all the values before the one having value 1 (January)
    date_year = da[time_dim].dt.year
    in_next_year = da[time_dim].dt.month.isin(next_year_months)
    # True/False : if a month is before January in the season, its year will be shifted to next
    # so as to be averaged in the correct season
    season_year = xr.where(in_next_year, date_year + 1, date_year).rename('year')
    # Year of the season of each da value: year + 1 for the months coming before January 
    valid_months = da.notnull().groupby(season_year).sum()
    # season_year is the key to gather data to be averaged together : same season_year, same average 
    # valid_months is the number of months in a season which are not NaN/Inf/None - a time series
    # as long as the seasonal time series computed below and which discard the NaN/Ind/None in its average
    seasonal = da.groupby(season_year).mean()
    return seasonal.where(valid_months == len(months), drop=True)
    # We drop all seasonal means which were computed from less months than the total number of months 
    # in the input argument months


def plotted_series(da, months, period='full'):
    """Returns a tuple (x values, y values) to be plotted from a monthly time series
    da which already holds at this stage the data in the units to be plotted (e.g. /1e12)
    months can be None or a list of numbers between 1 and 12 computed by
    season_months from the string ArcticSeas_season provided in params_ArcticSeas.py,
    If months is None, the function returns the original da monthly values together with
    its original monthly time axis.
    If months is a list of months, this functions computes seasonal averages through
    the seasonal_mean_complete_years function using the list of months provided as argument
    months and returns the list of years those seasons are associated with.
    If period = 'full', we keep all the entire time series (default). If period = 'YYYY-YYYY',
    we only keep the dates falling inside this year interval (empty arrays if there is none)."""
    seasonal = bool(months) # True if months is a list of numbers / False if it is None
    if seasonal:
        da = seasonal_mean_complete_years(da, months)
        # da now holds one value per complete season, along a 'year' dimension of integer years
    time_dim = da.dims[0]
    # da has only one dimension : a time dimension
    if period != 'full':
        first, last = [int(year) for year in period.split('-')]
        # We get the first and last year of the period to be plotted
        years = da[time_dim] if seasonal else da[time_dim].dt.year
        # Year of each value (repeated 12 times per year for monthly data) 
        da = da.where((years >= first) & (years <= last), drop=True)
        # We select the data within this period
    return da[time_dim].values, da.values
    # In case of annual/seasonal means, the x values look like [YYYY, YYYY .... YYYY]
    # In case of monthly time series, they look like [YYYY-MM-DD, YYYY-MM-DD,  .... ]
    # It is a numpy array of datetime64 or cftime for monthly time series, integers for seasonal.


def obs_climatology_curve(da, clim_period, months, sim_times):
    """Computes a climatology of an observed time series da (time as first dimension) over
    clim_period ('full' or 'YYYY-YYYY') and repeats it along the time dimension sim_times 
    of the experiments. For annual or seasonal means, a single mean value is repeated each
    year. For monthly climatologies, the mean of each calendar month is repeated in cycles.
    The input argument months is the list of calendar months of the season (obtained by
    season_months), or None for monthly values (for ArcticSeas_season = 'monthly'). 
    The output is a tupple : 
    (values along sim_times, first year, last year of the climatology), or None if no 
    observation falls within clim_period.
    da is a monthly series: if months is given, its seasonal means are first computed from
    complete seasons only (using seasonal_mean_complete_years)."""
    seasonal = bool(months) # True if months is a list of numbers / False if it is None
    if seasonal:
        da = seasonal_mean_complete_years(da, months)
        # da now holds one value per complete season, along a 'year' dimension of integer years
    time_dim = da.dims[0]
    # time_dim is the name of the time dimension from the input da file (sea ice area/volume)
    years = da[time_dim] if seasonal else da[time_dim].dt.year
    # years holds only the years of each time step (repeated 12 times per year for monthly data)
    if clim_period != 'full':
        first, last = [int(year) for year in clim_period.split('-')]
        # We get the first and last year of the period on which to compute the climatology
        da = da.where((years >= first) & (years <= last), drop=True)
        # We select the data within this period
        years = years[(years >= first) & (years <= last)]
        # years is updated to become the years of each time step within the selected period
    if da.sizes[time_dim] == 0:
        return None
    # The output of this function is None if no data was found inside the clim_period
    # (or, for seasonal means, if no complete season was found)
    if seasonal:
        # If we have seasonal or annual means, we repeat those for as long as the experiment
        values = np.full(len(sim_times), float(da.mean()))
    else:
        monthly_clim = da.groupby(time_dim + '.month').mean().reindex(month=range(1, 13))
        # time_dim.month is the key to gather data by month for groupby to average.
        # If some months are never observed (ex: summer months for CryoSat), 
        # reindex adds those missing months with a NaN
        sim_months = xr.DataArray(sim_times, dims='doesnotmatter').dt.month
        # sim_times is a numpy array of datetime64 values that needs to be converted
        # to an Xarray so that we can use .dt.month to pick the month of each point
        # on the time axis to build sim_months
        values = monthly_clim.sel(month=sim_months).values
        # Knowing the month of each point on the time axis (in sim_months), we can 
        # select the correct monthly value from the monthly climatology monthly_clim
        # so that values is a monthly climatology repeated in cycle
    return values, int(years.min()), int(years.max())


if do_ArcticSeas_timeseries:
    # 
    # ==> -- Open the section and an html file - Function from CLIMAF library
    # -----------------------------------------------------------------------------------------
    index += section("Sea ice area and volume per Arctic sea", level=4)
    #
    # ==> -- Control the size of the thumbnail -> thumbN_size
    # Note: thumbnail_size is None in share/default/default_atlas_settings.py, and can be set
    # in params_ArcticSeas.py (only read for this component); thumbnail_size_global otherwise
    # -----------------------------------------------------------------------------------------
    if thumbnail_size:
        thumbN_size = thumbnail_size
    else:
        thumbN_size = thumbnail_size_global
    #
    # ==> -- If ArcticSeas_season is available in params_ArcticSeas.py:
    # ==> -- 'ANM', month initials such as 'DJF' or 'NDJFMA', or 'monthly' for monthly values
    # ==> -- (season_month_list is then None in that last case). season is the working copy
    # -----------------------------------------------------------------------------------------
    season = str(ArcticSeas_season) if 'ArcticSeas_season' in dir() else 'ANM'
    try:
        season_month_list = None if season.lower() == 'monthly' else season_months(season)
    except ValueError as e:
        print("ArcticSeas: %s, using 'ANM' (annual means) instead" % e)
        season, season_month_list = 'ANM', season_months('ANM')
    season_label = 'monthly' if season_month_list is None else season.upper()
    #
    # ==> -- Southernmost latitude shown on the polar maps of the Arctic seas (ArcticSeas_latcutoff
    # ==> -- in params_ArcticSeas.py; 45 by default, as in check_masks.py)
    # -----------------------------------------------------------------------------------------
    map_latcutoff = float(ArcticSeas_latcutoff) \
        if 'ArcticSeas_latcutoff' in dir() and ArcticSeas_latcutoff is not None else 45.
    #
    # ==> -- Number of time steps read at once when processing the data (None: all at once)
    # -----------------------------------------------------------------------------------------
    chunk_time = ArcticSeas_chunk_time if 'ArcticSeas_chunk_time' in dir() else None
    #
    # Working directory (ArcticSeas_workdir in params_ArcticSeas.py, or ArcticSeas_workfiles
    # under the execution directory by default), created if non-existant already
    # -----------------------------------------------------------------------------------------
    if 'ArcticSeas_workdir' in dir() and ArcticSeas_workdir:
        workdir = ArcticSeas_workdir
    else:
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
                # -- automatically a name from its grid file's name, under ArcticSeas_cache_dir.
                # -- It does not need to exist yet -- ensure_mask() below builds it on demand,
                # -- exactly like any explicitly declared mask file.
                # -- We start from the name of the gridfile without the path and the extension
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
                # could overwrite local variables otherwise, 2. it reads its arguments from the
                # command line and stops with sys.exit() on any error, which would stop the
                # diagnostics otherwise
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
            immediately rename the Arctic one to <mask file name without .nc>_lat<latcutoff>.png,
            so every distinct grid ends up with its own plot name, and a new latcutoff gives a new
            plot. MPLBACKEND=Agg avoids the script's closing plt.show() blocking/failing headless."""
            png = os.path.join(workdir, '%s_lat%g.png' % (os.path.splitext(os.path.basename(maskfile))[0],
                                                          map_latcutoff))
            need = not os.path.exists(png) or os.path.getmtime(png) < os.path.getmtime(maskfile)
            if need and os.path.exists(check_masks_script):
                try:
                    # -- --latcutoff: southernmost latitude shown on the polar map
                    subprocess.run(['python3', check_masks_script,
                                     '--gridfile', gridfile, '--mask', maskfile, '--label', label,
                                     '--latcutoff', str(map_latcutoff)],
                                    cwd=workdir, check=True, env=dict(os.environ, MPLBACKEND='Agg'))
                    generated = os.path.join(workdir, 'check_masks_arctic.png')
                    if os.path.exists(generated):
                        os.replace(generated, png)
                except Exception as e:
                    print("ArcticSeas: check_masks.py failed for %s -> %s" % (label, e))
            return png if os.path.exists(png) else None
            # WARNING : Remove plot for missing seas

        # -- Datasets to process: the simulations of `models` (datasets_setup.py) followed by the
        # -- observations of ArcticSeas_obs (params_ArcticSeas.py)
        # -----------------------------------------------------------------------------------------
        # -- We do not change models so we do not really need that copy but for safety
        Wmodels = copy.deepcopy(models)
        # -- Observational references (params_ArcticSeas.py) go through exactly the same
        # -- processing as the simulations; their treatment is different only when checking for
        # -- missing observations, when computing their climatologies and when plotting
        Wobs = copy.deepcopy(ArcticSeas_obs) if 'ArcticSeas_obs' in dir() and ArcticSeas_obs else []
        # -- How each observation is plotted (see ArcticSeas_obs_display/ArcticSeas_obs_period
        # -- in params_ArcticSeas.py, overridden by 'obs_display'/'obs_period' in each obs
        # -- dictionary of ArcticSeas_obs)
        default_obs_display = ArcticSeas_obs_display if 'ArcticSeas_obs_display' in dir() else 'real'
        default_obs_period = ArcticSeas_obs_period if 'ArcticSeas_obs_period' in dir() else 'sim'
        obs_settings = dict()  # obs label -> (obs_display, obs_period)
        for obs in Wobs:
            obs_display = obs.pop('obs_display', default_obs_display)
            obs_period = str(obs.pop('obs_period', default_obs_period))
            # obs_display and obs_period are stored in the obs_display and obs_period standalone variables 
            # and simultaneously those keys are removed from the obs dictionary to be sent to CLIMAF
            if obs_display not in ('climato', 'real'):
                print("ArcticSeas: obs_display '%s' for %s should be 'climato' or 'real', using 'real'"
                      % (obs_display, model_label_of(obs)))
                obs_display = 'real'
            if obs_period not in ('full', 'sim') and not re.match(r'^\d{4}-\d{4}$', obs_period):
                print("ArcticSeas: obs_period '%s' for %s should be 'full', 'sim' or 'YYYY-YYYY', using 'sim'"
                      % (obs_period, model_label_of(obs)))
                obs_period = 'sim'
            obs_settings[model_label_of(obs)] = (obs_display, obs_period)
        obs_labels = [model_label_of(obs) for obs in Wobs]
        for label in obs_labels:
            if label in [model_label_of(model) for model in Wmodels]:
                print("ArcticSeas: observation %s has the same customname as a simulation; "
                      "their per-simulation settings and cache will be mixed up" % label)
        Wdatasets = Wmodels + Wobs

        # -- Find every distinct grid/mask actually needed by the simulations and
        # -- observations, then build/locate each of them (deduplicated by mask file)
        # -----------------------------------------------------------------------------------------
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
            index += "No Arctic seas mask file could be found or built for any of the " \
                      "simulations or observations (see the job log for details)."
            index += close_line() + close_table()

        # Finally !!! We have everything we need to run the actual diagnostics
        # --------------------------------------------------------------------------------------------
        else:
            # -- Plot of the Arctic seas geometry, one per distinct grid (not one per
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
            # -- set of named seas (netcdf variable names, e.g. 'barentse'), discretized on
            # -- different grids; their 'long_name' attribute (e.g. 'Barents Sea') labels the sections
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
            # seaindex_files is defined to contain the file names holding the sea ice area / volume  
            # organized by simulation and variable
            model_labels = [] # list of model and observations labels (i.e. customnames by default)
            message_failed_comput = []
            # message_failed_comput is a list created to gather error messages as strings of characters
            # which were intended for the logs. Those are captured to be printed on the atlas before
            # the figures in addition to the logs : missing grid, mask, data...
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
                        # model_labels holds all the labels of all models and observations

                    # -- Netcdf variable name to find in output files for this simulation/variable
                    # -----------------------------------------------------------------------------------------
                    real_variable = ArcticSeas_variable_names.get(model_label, dict()).get(variable)
                    if real_variable is None:
                        message = "no ArcticSeas_variable_names entry for %s / %s, skipping" % (model_label, variable)
                        if model_label in obs_labels:
                            # normal for an observation which does not provide this variable
                            # (e.g. no thickness in a concentration product): log only
                            print("ArcticSeas: " + message)
                        else:
                            report_failed_comput(message_failed_comput, message)
                        continue
                        # We skip to next model/obs here if variable not declared in params_ArcticSeas.py

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
                        report_failed_comput(message_failed_comput, "no data found for %s / %s (variable %s), skipping"
                                             % (model_label, variable, real_variable))
                        continue
                        # We skip to next model/obs here if period is not available = no input data found

                    # Build a name for the output files (sea ice volume and sea ice area)
                    cache_label = CACHE_VARIABLE_LABELS[(variable, ArcticSeas_meanORsum)]
                    cache_file = os.path.join(ArcticSeas_cache_dir, sanitize(model_label), cache_label + '.nc')
                    # Finding the last simulation/observation year both in the input and output files
                    simulation_latest_year = latest_year_of_wmodel(wmodel)
                    cache_latest_year = latest_year_in_ncfile(cache_file)
                    up_to_date = (cache_latest_year is not None and simulation_latest_year is not None
                                  and cache_latest_year >= simulation_latest_year)
                    # The output sea ice index file does not need to be updated.

                    if up_to_date:
                        seaindex_files[variable][model_label] = cache_file
                        continue
                        # Here we skip to the next model/obs because we do not need to update the diagnostics

                    # -- Grid and matching sea mask for this simulation or observation 
                    gridfile, maskfile, label = resolve_grid_and_mask(wmodel)

                    if gridfile is None:
                        report_failed_comput(message_failed_comput, "no ArcticSeas_gridfile/ArcticSeas_maskfile entry "
                                             "for %s, skipping" % model_label)
                        continue
                    if not os.path.exists(gridfile):
                        report_failed_comput(message_failed_comput, "grid file not found for %s -> %s"
                                             % (model_label, gridfile))
                        continue
                    if maskfile not in available_masks:
                        report_failed_comput(message_failed_comput, "no sea mask available for %s (mask %s), skipping"
                                             % (model_label, maskfile))
                        continue

                    try:
                        # Gather in dataset dat all informations about model read in datasetup.py 
                        # or in ArcticSeas_obs for the observations together with
                        # the period of the netcdf files computed by get_period_manager
                        dat = ds(**wmodel)
                        # -- No time averaging here (e.g. ccdo yearmean), which would average
                        # -- whatever months are present in each year: comp_seaiceindex.py works
                        # -- on the monthly fields, and seasonal means are computed when plotting,
                        # -- from complete seasons only (see seasonal_mean_complete_years)
                        datafile = cfile(dat)
                        # -- cfile is a CLIMAF function that gathers all the input file into a 
                        # -- single one containing only the variable and period of interest
                        # -- datafile is the complete path to this output file

                        # -- in case sea ice concentration is published as a percentage (0-100) 
                        # -- turn it into a fraction (0-1). correct and write the corrected copy
                        # -- to our own workdir rather than overwriting CliMAF's cache file.
                        # -- The field is read by pieces of ArcticSeas_chunk_time time steps (see
                        # -- chunked_in_time) to limit the RAM used by long model series
                        if variable == 'sic':
                            with xr.open_dataset(datafile) as sic_unit:
                            # datafile is opened only in this block, closed automatically at the end of the block
                            # (it must stay open while the corrected copy is written, since the values
                            # are only read, piece by piece, at that time)
                                sic_field = chunked_in_time(sic_unit[real_variable], chunk_time)
                                is_percent = float(sic_field.max()) > 1.5
                                if is_percent:
                                    fraction_file = os.path.join(
                                        workdir, 'sic_fraction_%s.nc' % sanitize(model_label))
                                    fraction_ds = sic_unit.copy()
                                    fraction_ds[real_variable] = sic_field / 100.
                                    fraction_ds.to_netcdf(fraction_file)
                            if is_percent:
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
                        if chunk_time:
                            cmd += ['--chunk_time', str(chunk_time)]
                        # -- Output captured (and echoed to the log) to keep its 'missing values'
                        # -- warnings next to the cache file, so that they are shown on the atlas
                        # -- page also when the cache is reused later
                        result = subprocess.run(cmd, capture_output=True, text=True)
                        print(result.stdout + result.stderr)
                        # print sends to the log and check_returncode captures the e to
                        # be printed below
                        result.check_returncode()
                        # send each error message occcurence of SEA_WARNING_PATTERN to the 
                        # text file stored next to the cache file
                        write_sea_warnings(cache_file, result.stderr)
                        # cache_file is stored in the seaindex_files dictionary
                        seaindex_files[variable][model_label] = cache_file
                    except Exception as e:
                        report_failed_comput(message_failed_comput, "failed to compute %s (variable %s) for %s -> %s"
                                             % (variable, real_variable, model_label, e))

            # ==> -- Monthly plots (ArcticSeas_season = 'monthly') use the dates of the files as x 
            # ==> -- values: dates in a non-standard calendar (e.g. noleap, 360_day) are cftime objects, 
            # ==> -- which matplotlib only plots with the nc-time-axis package, and cannot mix with 
            # ==> -- standard dates (i.e. in another simulation, or an observation with obs_display = 
            # ==> -- 'real'). Warn in the log and the atlas page if any dataset uses such a calendar 
            # ==> -- (seasonal means use integer years and are not concerned)
            # -----------------------------------------------------------------------------------------
            if season_month_list is None:
                # season_month_list == None means monthly data to be plotted
                non_standard_calendars = dict()  # label -> calendar
                for variable, expe in seaindex_files.items():
                    # variable : sic/sit
                    # expe : {simulation/observation label: its cache file} for this variable
                    for model_label, cache_file in expe.items():
                        if model_label in non_standard_calendars or not os.path.exists(cache_file):
                        # this model/obs is already known to have a non-standard calendar 
                        # or we can not check its calendar because we do not have any cache_file
                            continue
                        with xr.open_dataset(cache_file) as cached:
                            first_var = next(iter(cached.data_vars), None)
                            if first_var is None:
                                continue
                            dates = cached[cached[first_var].dims[0]]
                            # dates contain the time axis of this cache_file
                            if not np.issubdtype(dates.dtype, np.datetime64):
                                # If the type of time-axis is not datetime64, we get in the time axis which 
                                # type it is (360_day, no_leap, ...)
                                non_standard_calendars[model_label] = dates.encoding.get('calendar', 'non-standard')
                                # We store that information to be printed below.
                if non_standard_calendars:
                    # If non_standard_calendars is not empty, we will write a message on the html page because
                    # we can not mix monthly values for that experiment / observation with datetime64 time axes 
                    # of other experiments on the same figure.
                    calendars_text = ', '.join('%s (%s)' % (label, calendar)
                                               for label, calendar in non_standard_calendars.items())
                    # calendars_text = label1 (calendar1) 
                    # or calendars_text = label1 (calendar1), label2 (calendar2), label3 (calendar3) 
                    # with any number of label (calendar)
                    message = "Monthly plots (ArcticSeas_season = 'monthly') with a non-standard calendar: " + \
                              calendars_text + ". Plotting these dates needs the nc-time-axis package, " \
                              "and fails if they are mixed with standard dates (other simulations, " \
                              "observations with obs_display = 'real'): use a season such " \
                              "as 'ANM' instead."
                    print("ArcticSeas: WARNING - " + message)
                    index += open_table()
                    index += start_line('<font color="red">Warning - ' + message + '</font>')
                    index += close_line() + close_table()

            # ==> -- Error messages explaining why sea ice area/volume could not be computed for one
            # ==> -- or more simulations or observations (these are also printed in the log)
            # -----------------------------------------------------------------------------------------
            if message_failed_comput:
                index += open_table()
                index += start_line('<font color="red">Not computed - ' +
                                    '<br>Not computed - '.join(message_failed_comput) + '</font>')
                index += close_line() + close_table()

            # ==> -- Build one table row per sea; for each sea, one plot for the ice area (sic) and
            # ==> -- one for the ice volume (sit), showing all the simulations and observations
            # -----------------------------------------------------------------------------------------
            index += open_table()
            if ArcticSeas_meanORsum == 'sum':
                variable_plot_specs = {
                    'sic': dict(title='Sea ice area',
                                ylabel='Sea ice area (Millions km2)'),
                    'sit': dict(title='Sea ice volume',
                                ylabel='Sea ice volume (Thousand km3)'),
                }
            else:
                variable_plot_specs = {
                    'sic': dict(title='Sea ice concentration',
                                ylabel='Area-weighted mean sea ice concentration (%)'),
                    'sit': dict(title='Sea ice thickness',
                                ylabel='Area-weighted mean sea ice thickness (m)'),
                }
            warnings_per_sea = dict()  # sea -> [ 'model_label (diagnostics): message', ...]
            # Warning messages about missing data are reorganized to be gathered per sea
            # instead of per variable and experiment/observation (initially in a text file
            # next to each cache_file
            for variable, expe in seaindex_files.items():
                # variable : sic/sit
                # expe : {simulation/observation label: its cache file} for this variable
                for model_label, cache_file in expe.items():
                    for warned_sea, message in read_sea_warnings(cache_file).items():
                        # For this set of variable and experiment/observation, if we have a warning
                        # message about the sea warned_sea, we append it to previous messages for
                        # this warned sea (to nothing if there was none before)
                        warnings_per_sea.setdefault(warned_sea, []).append(
                            # if there was no previous message, set default creates an empty entry
                            # for warned_sea in the warning_per_sea dictionary
                            '%s (%s): %s' % (model_label, variable_plot_specs.get(variable, dict(title='unknown diagnostics'))['title'], message))
                            # (variable missing from variable_plot_specs -> 'unknown diagnostics')
                            # The final message to be printed below the sea title will be :
                            # model_label (Sea ice area) : warning message about missing values
                            # model_label (Sea ice volume) : warning message about missing values
                            # etc ... for any number of model_label which has missing values in this sea
                    # -- Repeated dates (overlapping input files) are removed silently before plotting
                    # -- when their values are identical; when they differ, the input files hold two
                    # -- different versions of the series and we add a warning below the sea title 
                    if not os.path.exists(cache_file):
                        # if the cachefile is not there, we go to next model or observation for this variable
                        continue
                    with xr.open_dataset(cache_file) as cached:
                        for sea in cached.data_vars:
                            da, n_removed, identical = without_duplicate_dates(cached[sea].load())
                            # The without_duplicate_dates applied on the cache file checks for this sea
                            # whether there are any duplication of values for the same dates and returns
                            # the number of duplicates if there are any (n_removed) together with 
                            # a True/False flag if both values are the same / different for the same dates
                            if n_removed and not identical:
                                # We found two different values for one or more dates
                                warnings_per_sea.setdefault(sea, []).append(
                                    '%s (%s): %d repeated dates with different values, first value kept'
                                    % (model_label, variable_plot_specs.get(variable, dict(title='unknown diagnostics'))['title'],
                                       n_removed))
            # -- comp_seaiceindex.py stores, as an attribute of each sea in its output file, the
            # -- total area of the sea ('sum') or the mean grid cell area within the sea ('mean')
            # -- on the grid of that simulation/observation
            if ArcticSeas_meanORsum == 'sum':
                area_attr, area_divisor, area_title = 'total_area_m2', 1e9, 'Sea area (Thousand km2)'
            else:
                area_attr, area_divisor, area_title = 'mean_cell_area_m2', 1e6, 'Mean grid cell area (km2)'
            # -- The sea areas depend only on the grid and sea mask: they are listed once per group
            # -- of simulations/observations sharing the same mask file
            # -- available_masks contain : { maskfile : (gridfile, group_of_labels) , ... }, one entry
            # -- per group of simulations/observations sharing the same grid and mask
            area_label_groups = [labels for (gridfile, labels) in available_masks.values()]
            # -- area_label_groups is a list of lists of labels using the same mask (hence the same grid)
            for sea in seas:
                sea_name = sea_display_names.get(sea, sea)
                # If sea_display_names is not defined for sea, sea_name is sea
                # The long_name attributes of seas were gathered earlier in sea_diplay_names
                title = sea_name + '<br><span style="font-size:70%">' + area_title + ': ' + \
                    sea_areas_text(sea, area_label_groups, seaindex_files, area_attr, area_divisor) + '</span>'
                    # sea_areas_text is a function defined earlier which returns a line with 
                    # group of model_labels : area_of_the_sea; group of model_labels area_of_the_sea, etc ...
                if sea in warnings_per_sea:
                    title += '<br><span style="font-size:70%; color:red">Warning - ' + \
                             '<br>Warning - '.join(warnings_per_sea[sea]) + '</span>'
                index += start_line(title)
                for variable, plot_specs in variable_plot_specs.items():
                    # variable : sic/sit
                    # plot_specs : dictionary holding the plot title and ylabel for this variable

                    # -- Gather the data available for this sea/variable first, so we know
                    # -- which simulations (if any) are missing before deciding whether to
                    # -- plot the available ones or to report an error instead
                    # -----------------------------------------------------------------------------------------
                    curves = []  # list of tuples (model_label, time_values, data_values) containing
                                 # the model label and the x and y values to be plotted for each simulation
                                 # which actually has data to be plotted
                    obs_series = dict()  # a dictionary holding for each obs label a DataArray with the
                                 # monthly data for this observational dataset if the sea is in its cache file 
                    present_labels = [] # only the labels of the models or obs which have data for this sea
                    for model_label, outfile in seaindex_files.get(variable, dict()).items():
                        # if variable is not in seaindex_files, empty dictionary so we do not
                        # go through this loop
                        # if variable is in seaindex_files, we loop through the list of (model_label, cachefile) 
                        # for this variable 
                        if not os.path.exists(outfile):
                            continue
                        # if the cachefile is not there, we go to next model or observation for this variable
                        with xr.open_dataset(outfile) as out_ds:
                            if sea not in out_ds.data_vars:
                                continue
                            # if we do not have any data for this sea for this model we go to next model or obs
                            da = out_ds[sea].load()
                        # -- The file out_ds is closed here: da is loaded in RAM, we no longer need out_ds
                        da = without_duplicate_dates(da)[0]
                        # -- We remove duplicated dates/values in da because those repeated dates would trick
                        # -- the seasonal_mean_complete_years function into considering seasons as incomplete because
                        # -- it checks that the number of expected months in a season is indeed there
                        # -- only the first value for each date is kept (a warning was set above, below the
                        # -- sea title, if the repeated values were different)
                        if variable in ('sic', 'sit') and ArcticSeas_meanORsum == 'sum':
                            # -- sea ice area (m^2 -> millions km^2) or volume (m^3 ->
                            # -- thousand km^3), to match the plot's ylabel; purely local
                            # -- to this plot, the cache file on disk is left untouched
                            da = da / 1e12
                        if model_label in obs_labels:
                            # All observational data have been treated as model data and their cachefile are also
                            # stored in seaindex_files so that if model_label is actually an obs_label, da needs 
                            # to be stored in the obs_series dictionary with the obs_label as key and not in the 
                            # curve list which is for model/experiment data
                            obs_series[model_label] = da
                        else:
                            # -- Seasonal means of complete seasons only (ArcticSeas_season in params_ArcticSeas.py)
                            time_values, data_values = plotted_series(da, season_month_list)
                            # We obtain a tuple containing a numpy array of integer years if we computed a seasonal
                            # mean (using only season with no missing data inside the season) as first element
                            # and the corresponding seasonal means as second element
                            # If no seasonal mean was requested (season_month_list = None), the tuple contains
                            # a numpy array of datetime64 or cftime monthly dates and the corresponding monthly values 
                            curves.append((model_label, time_values, data_values))
                            # curves is a list of tuples in which we add the tuple holding this model_label, 
                            # its time_values to be plotted along x, its data_values to be plotted along y
                        present_labels.append(model_label)
                        # only the labels of the model or observations that have data for this sea

                    # -- Observations: either their actual dates within obs_period ('real'), or their
                    # -- climatology over obs_period repeated along the simulation dates ('climato').
                    # -- Without any simulation to follow, observations will be plotted in 'real' mode
                    # -- within obs_period
                    # -----------------------------------------------------------------------------------------
                    obs_curves = []  # a list of tuples : (obs label, legend label, time_values, data_values)
                    sim_times = None
                    if curves:
                        # If we do have model data for this sea, sim_times contains all the dates covered by 
                        # all the simulations of all models/experiments without duplicate
                        sim_times = np.unique(np.concatenate([time_values for (mod_lab, time_values, data) in curves]))
                    for obs_label, da in obs_series.items():
                        # If we do have observations to plot, we loop over those observational datasets
                        obs_display, obs_period = obs_settings[obs_label]
                        # We get back the information stored earlier in obs_settings :
                        # real / climato
                        # full / sim / YYYY-YYYY
                        if obs_period == 'sim':
                            # -- From the first to the last year of the simulations for this diagnostics
                            # -- ('full' if there is no simulation to follow)
                            if sim_times is None:
                                obs_period = 'full'
                            else:
                                sim_years = sim_times if season_month_list \
                                    else xr.DataArray(sim_times, dims='doesnotmatter').dt.year.values
                                # We only create a DataArray to use the dt.year.values and extract the year of each date
                                obs_period = '%d-%d' % (sim_years.min(), sim_years.max())
                                # For option 'sim', we build an obs_period covering all the years of all the simulations
                        if obs_display == 'real' or sim_times is None:
                            # If we do not have any sim_times (no simulation available to be plotted), we plot in mode 'real'
                            # -- Only the points within obs_period (every point if 'full')
                            time_values, data_values = plotted_series(da, season_month_list, obs_period)
                            # We obtain a tuple containing a numpy array of integer years if we computed a seasonal
                            # mean (using only season with no missing data inside the season) as first element
                            # and the corresponding seasonal means as second element
                            # If no seasonal mean was requested (season_month_list = None), the tuple contains
                            # a numpy array of datetime64 or cftime monthly dates and the corresponding monthly values 
                            # If obs_period is 'full', the whole time series is output, if a period is provided
                            # with YYYY-YYYY, only this period is output
                            if len(time_values) == 0:
                                print("ArcticSeas: no %s data within obs_period %s for %s, not plotted"
                                      % (obs_label, obs_period, sea))
                                continue
                            obs_curves.append((obs_label, obs_label, time_values, data_values))
                            # We include in the obs_curves list a new tuple with (obs_customname,
                            # obs_legend, x data to be plotted, y data to be plotted)
                        else:
                            clim = obs_climatology_curve(da, obs_period, season_month_list, sim_times)
                            # obs_climatology_curve computes a climatology of the DataArray da over the
                            # obs_period period and repeats it along the sim_times time axis : a single
                            # value for seasonal mean, a 12-month cycle for monthly means. 
                            # The months within the season are passed through season_month_list. If 
                            # season_month_list is None, we obtain a monthly climatology.
                            # clim is a tuple holding the numpy array of y values to be plotted against 
                            # the sim_times x values, together with the first and last year of the climatology
                            if clim is None:
                                print("ArcticSeas: no %s data within obs_period %s for %s, not plotted"
                                      % (obs_label, obs_period, sea))
                                continue
                            values, first_year, last_year = clim
                            obs_curves.append((obs_label, '%s (clim %d-%d)' % (obs_label, first_year, last_year),
                                               sim_times, values))
                            # We include in the obs_curves list a new tuple with (obs_customname,
                            # obs_legend, x data to be plotted, y data to be plotted)

                    # model_labels holds every simulation listed in datasets_setup.py and every 
                    # observation listed in params_ArcticSeas.py
                    # Any of them missing from present_labels has no data here.
                    missing_labels = [m for m in model_labels
                                      if m not in present_labels and m not in obs_labels]
                    # Observations are left out: a missing observation is just not plotted
                    # if m is in obs_labels, we do not include it in missing_labels even if
                    # there is no data for this observation

                    fig, ax = plt.subplots(figsize=(6, 4))
                    # New figure created for each sea, each variable and associated diagnostics
                    if not curves and not obs_curves:
                        # -- Always an error message here, regardless of ArcticSeas_on_missing_simulations
                        # -- when we don't have any observation nor any model data
                        ax.text(0.5, 0.5, 'No data available for any simulation nor any observation',
                                ha='center', va='center', wrap=True, transform=ax.transAxes)
                    elif missing_labels and ArcticSeas_on_missing_simulations == 'error':
                        # -- As soon as a single simulation is missing, a message is shown here
                        # -- instead of the figure with the other simulations and observations
                        ax.text(0.5, 0.5, 'Data missing for: %s' % ', '.join(missing_labels),
                                ha='center', va='center', wrap=True, transform=ax.transAxes)
                    else:
                        # -- ArcticSeas_on_missing_simulations == 'partial' (or nothing missing):
                        # -- plot the simulations that do have data, silently skipping the rest
                        # -- Simulations first, then observations on top, each with a fixed style
                        # -- (same color for a given simulation on every plot, from its rank among
                        # -- the simulations of datasets_setup.py; same style for a given
                        # -- observation on every plot, from its rank in ArcticSeas_obs)
                        sim_labels = [m for m in model_labels if m not in obs_labels]
                        # sim_labels : all the simulations in the order of datasets_setup.py,
                        # whether or not they have data on this plot
                        for model_label, time_values, data_values in curves:
                            color = SIM_COLORS[sim_labels.index(model_label) % len(SIM_COLORS)]
                            ax.plot(time_values, data_values, lw=1.5, color=color, label=model_label)
                        for obs_label, legend_label, time_values, data_values in obs_curves:
                            color, linestyle = OBS_STYLES[obs_labels.index(obs_label) % len(OBS_STYLES)]
                            ax.plot(time_values, data_values, lw=2, color=color, linestyle=linestyle,
                                    label=legend_label, zorder=3)
                        ax.legend(fontsize=8)

                    ax.set_title("%s - %s (%s)" % (sea_name, plot_specs['title'], season_label))
                    ax.set_xlabel('Time')
                    ax.set_ylabel(plot_specs['ylabel'])
                    fig.tight_layout()

                    png_path = os.path.join(workdir, 'ts_%s_%s.png' % (sanitize(sea), variable))
                    fig.savefig(png_path, dpi=100)
                    plt.close(fig)

                    index += cell(plot_specs['title'], png_path, thumbnail=thumbN_size, hover=hover, **alternative_dir)
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
