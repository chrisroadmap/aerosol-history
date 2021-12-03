# Energy Budget Constraints on the Time History of Aerosol Forcing and Climate Sensitivity
This code demonstrates how we can use constraints of observed surface warming from 1850-2019 and ocean heat update from 1971-2018 to constrain historical aerosol forcing, equilibrium climate sensitivity and transient climate response.

![Comparison of constrained aerosol forcing and CMIP6 models from 1850 to 2019](figures/models_v_constrained.png?raw=true)
Comparison of the energy-budget constrained aerosol forcing best estimate (thick grey line) and 5-95% range (grey shaded band) with CMIP6 model aerosol forcing results (coloured). CMIP6 individual years are coloured points and a 11-year Savitzy-Golay smoothing filter is applied to these model results to produce smoothed model time series estimates (coloured lines). Note this plot is using an 1850 baseline for direct comparison with CMIP6 models, whereas all results in the paper use a 1750 baseline for "pre-industrial" conditions.

## Reference and citation
The paper describing this method is published in Journal of Geophysical Research Atmospheres and available as an open-access paper at https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020JD033622.

Please cite as:

Smith, C. J., Harris, G. R., Palmer, M. D., Bellouin, N., Collins, W., Myhre, G., Schulz, M., Golaz, J.-C., Ringer, M., Storelvmo, T., Forster, P. M., (2021). Energy Budget Constraints on the Time History of Aerosol Forcing and Climate Sensitivity. Journal of Geophysical Research: Atmospheres, 126, e2020JD033622. https://doi.org/10.1029/2020JD033622

## Reproduction
As always, the best way to avoid disaster is to use a virtual environment. My weapon of choice is miniconda3 (https://docs.conda.io/en/latest/miniconda.html). For this code more than ever it is a good idea to use `conda`, as code makes heavy use of `iris` which is difficult(/impossible?) to install through `pip`.

Everything you need should be packed inside the `environment` file. Run this:

    conda env create -f environment.yml

This will create the environment called `aerosol-history`. If developing and you need to add a new package to the environment, run this:

    conda env update -f environment.yml --prune

Once set up, activate the environment:

    conda activate aerosol-history

then run the notebooks in the `notebooks` directory in order. Some notebooks take a while, e.g. the two-layer model calculations for 13 models (notebook 30) and the internal variability calculation (notebook 10).

In the `scripts` directory is the code that runs the APRP decomposition for the 11 climate models used in the paper. However, these scripts won't run because they require the CMIP6 and E3SM data (the file names point to my local paths on the University of Leeds Earth and Environment Linux server, and they are not supplied here due to file sizes, but all are available either from the ESGF portal or https://esgf-node.llnl.gov/projects/e3sm/). The important data that is required to reproduce the results is shipped in the directories `input_data` (if it is a raw input) or `output_data` (if somewhere along the process, the code calculates it).

The only figure in the whole paper and supplement not reproduced by this repository is Figure S3, by Glen Harris. The two-layer results (red curves in fig. S3) can be reproduced by plugging the parameters from each CMIP6 model (`input_data/scmpy2L_calib_n=44_eps=fit_v20200702.txt`) into the two-layer model.
