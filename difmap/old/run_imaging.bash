#/bin/bash

script=difmap_wrapper.py

# Absolute path of input data
inpath='./data/'

# Absolute path to save outputs
outpath='./difmap_imaging_results/'

# input uvdata
uvfits='hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s'

# imaging parameters
mask_size=95
ALMA_weight=0.5
uvbin=0

# run python wrapper script
python ${script} --inpath ${inpath} --outpath ${outpath} \
       --uvfits ${uvfits} --mask_size ${mask_size} --ALMA_weight ${ALMA_weight} --uvbin ${uvbin} \
