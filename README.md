# aerosol-history
Using energy budget constraints to estimate historical aerosol forcing

## Reproduction
As always, the best way to avoid disaster is to use a virtual environment. My weapon of choice is miniconda3 (https://docs.conda.io/en/latest/miniconda.html). For this code more than ever it is a good idea to use `conda`, as code makes heavy use of `iris` which is difficult(/impossible?) to install through `pip`. Therefore, create a new environment, activate it, and use `conda` to grab `iris`:

    conda install -c conda-forge iris

When this is done, grab the repo-specific dependencies:

    pip install -r requirements.txt

In the `scripts` directory is the code that runs the APRP decomposition for the 11 climate models used in the paper. However, these scripts won't run because they require the CMIP6 and E3SM data (the file names point to my local paths, and they are not supplied here due to file sizes, but all are available either from the ESGF portal or https://esgf-node.llnl.gov/projects/e3sm/). The important data that is required is shipped in `output_data`.

The only figure in the whole paper and supplement not reproduced by this repository is Figure S3, by Glen Harris. The two-layer results (red curves in fig. S3) can be reproduced by plugging the parameters from each CMIP6 model (`input_data/scmpy2L_calib_n=44_eps=fit_v20200702.txt`) into the two-layer model.

## Reference
The paper describing this method is under review at the moment. The preprint can be accessed here: https://www.essoar.org/doi/abs/10.1002/essoar.10503977.2
