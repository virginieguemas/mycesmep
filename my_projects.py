from climaf.api import *
file_patterns = []
file_patterns.append('${root}/${experiment}/${experiment}_*_monthly_${variable}_YYYY-YYYY.nc')
cproject('N4CPL', 'experiment', 'variable','root','model')
dataloc(project = 'N4CPL', organization = 'generic', url = file_patterns)
calias("N4CPL", 'tos','tos',offset=273.15) 
#calias("N4CPL", 'so','so', scale=0.9953067)
calias("N4CPL", 'thetao','thetao',offset=273.15)
calias("N4CPL", 'sic','siconc',scale=100.) #
calias("N4CPL", 'sit','sithic') #
calias("N4CPL", 'hurs','hurs',scale=100.) #


file_patterns = []
file_patterns.append('${root}/${experiment}/${experiment}_*_monthly_${variable}_YYYY-YYYY.nc')
cproject('N3CPL', 'experiment', 'variable','root','model')
dataloc(project = 'N3CPL', organization = 'generic', url = file_patterns)
calias("N3CPL", 'tos','tos',offset=273.15)
calias("N3CPL", 'thetao','thetao',offset=273.15)
#calias("N3CPL", 'so','so', scale=0.9953067)
calias("N3CPL", 'msftyz', 'zomsfatl',scale=1026.e6) 
calias("N3CPL", 'sivol','sivolu') 
calias("N3CPL", 'sit','sivolu') 
calias("N3CPL", 'sic','siconc',scale=100.) #
calias("N3CPL", 'siextentn','extentn',filenameVar='scalari') 
calias("N3CPL", 'siextents','extents',filenameVar='scalari') 
derive("N3CPL", 'hfds','plus','tohfls','hflx_rnf_cea')
calias("N3CPL", 'hurs','hurs',scale=100.) #


# -- Observations: NOAA/NSIDC CDR of passive microwave sea ice concentration, V6,
# -- monthly, Northern Hemisphere, 25 km polar stereographic grid (psn25).
# -- The sensor in the file name changes with time (n07, F08, F11, F13, F17, am2),
# -- but there is never more than one file per month.
# -- Variable: 'cdr_seaice_conc_monthly' (fraction 0-1 once scale_factor is applied).
file_patterns = []
file_patterns.append('/cnrm/ioga/Data/Climate-data-Records/passive_microwave-sea-ice-concentration/'
                     'NOAA/G02202_V6/north/monthly/sic_psn25_YYYYMM_*_v06r00.nc')
cproject('NSIDC_CDR_V6', 'variable')
dataloc(project = 'NSIDC_CDR_V6', organization = 'generic', url = file_patterns)


