#!/bin/bash
version=dev/climaf_3.0

#export climaf=/home/voldoire/outils/climaf_1.1_dev
export climaf=/cnrm/est/COMMON/climaf/${version}
#export climaf=/cnrm/est/COMMON/climaf/climaf_1.2.13_noslice
export PYTHONPATH=$mafate:$climaf:$PYTHONPATH                                   
export PATH=$climaf/bin:$PATH:/cnrm/est/COMMON/CDFTOOLS_3.0/bin 
export CLIMAF_FIX_NEMO_TIME=yes                                                
export TMPDIR=/sx/d0/Users/voldoire/tmp
    
#ipython notebook  &
