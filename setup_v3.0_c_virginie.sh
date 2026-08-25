#!/bin/bash
#
# Créé à partir du fichier setup_v3.0_c.sh d'Aurore renommé setup_v3.0_c_aurore.sh; ajout de ne définit pas CLIMAF_CACHE, CLIMAF; correction de TMPDIR ; suppression de $mafate pour PYTHONPATH que je n'ai pas défini ; modification du PYTHONPATH pour ne pas utiliser la version dev mais le PYTHONPATH du setenv par défault

pythonversion=3.12
version=dev/climaf_3.0
#version=climaf_3.0

export climaf=/cnrm/est/COMMON/climaf/${version}
#export PYTHONPATH=$climaf:$PYTHONPATH                                   
export PYTHONPATH=/cnrm/est/COMMON/climaf/add_packages/lib/python$pythonversion:/cnrm/ioga/Users/guemas/cesmep/C-ESM-EP/share/cesmep_modules:${climaf}:$PYTHONPATH 
export PATH=$climaf/bin:$PATH:/cnrm/est/COMMON/CDFTOOLS_3.0/bin 
export CLIMAF_FIX_NEMO_TIME=yes                                                
export CLIMAF_CACHE=/cnrm/ioga/Users/guemas/NO_SAVE/CESMEP_clim
export CLIMAF=/cnrm/est/COMMON/climaf/climaf_3.0
export TMPDIR=/sx/d0/Users/guemas/tmp
    
#ipython notebook  &
