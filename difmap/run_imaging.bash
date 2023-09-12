#/bin/bash

script=difmap_wrapper.py

# Absolute path of input data
# (e.g., inpath='/home/user/current_path/data') -- should exists containing uvfits data
inpath=

# Absolute path to save outputs
# (e.g., outpath='/home/user/current_path/difmap_imaging_results') -- will be generated as defined
outpath=

# Absolute path to the input parameters table and to save the searched parameter table
# (e.g., hdf5in='/home/user/current_path/hdf5s/parms_topset_3599.hdf5') -- should exists
# (e.g., hdf5out='/home/user/current_path/hdf5s/parms_searched_3599.hdf5') -- will be generated as defined
hdf5in=
hdf5out=

# number of processors to use in parallel
nproc=4

# input uvdata
# (e.g., uvfits='hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s')
uvfits='hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s'

#parameters id to image
# param_search_id='all' for imaging all topset parameters
# otherwise, param_search_id='0,1,2,3' for instance
param_search_id='all'


# run python wrapper script
python ${script} --inpath ${inpath} --outpath ${outpath} --hdf5in ${hdf5in} --hdf5out ${hdf5out} --nproc=${nproc} \
       --uvfits ${uvfits} --param_search_id ${param_search_id} \
