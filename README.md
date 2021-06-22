# Energy Budget Constraints on the Time History of Aerosol Forcing and Climate Sensitivity
This code demonstrates how we can use constraints of observed surface warming from 1850-2019 and ocean heat update from 1971-2018 to constrain historical aerosol forcing, equilibrium climate sensitivity and transient climate response.

## Reference and citation
The paper describing this method has been accepted in JGRA. The accepted version can be accessed here: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020JD033622

Please cite as:

Smith, C. J., Harris, G. R., Palmer, M. D., Bellouin, N., Collins, W., Myhre, G., et al. (2021). Energy Budget Constraints on the Time History of Aerosol Forcing and Climate Sensitivity. Journal of Geophysical Research: Atmospheres, 126, e2020JD033622. https://doi.org/10.1029/2020JD033622

## Reproduction
As always, the best way to avoid disaster is to use a virtual environment. My weapon of choice is miniconda3 (https://docs.conda.io/en/latest/miniconda.html). For this code more than ever it is a good idea to use `conda`, as code makes heavy use of `iris` which is difficult(/impossible?) to install through `pip`. Therefore, create a new environment, activate it, and use `conda` to grab `iris`:

    conda install -c conda-forge iris

When this is done, grab the repo-specific dependencies:

    pip install -r requirements.txt

Once set up, run the notebooks in the `notebooks` directory in order. Some notebooks take a while, e.g. the two-layer model calculations for 13 models (notebook 02) and the internal variability calculation (notebook 00).

In the `scripts` directory is the code that runs the APRP decomposition for the 11 climate models used in the paper. However, these scripts won't run because they require the CMIP6 and E3SM data (the file names point to my local paths, and they are not supplied here due to file sizes, but all are available either from the ESGF portal or https://esgf-node.llnl.gov/projects/e3sm/). The important data that is required to reproduce the results is shipped in the directories `input_data` (if it is a raw input) or `output_data` (if somewhere along the process, the code calculates it).

The only figure in the whole paper and supplement not reproduced by this repository is Figure S3, by Glen Harris. The two-layer results (red curves in fig. S3) can be reproduced by plugging the parameters from each CMIP6 model (`input_data/scmpy2L_calib_n=44_eps=fit_v20200702.txt`) into the two-layer model.
