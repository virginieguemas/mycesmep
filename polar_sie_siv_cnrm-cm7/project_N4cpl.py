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


