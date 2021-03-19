# aerosol-history
Using energy budget constraints to estimate historical aerosol forcing

## Reproduction
As always, the best way to avoid disaster is to use a virtual environment. My weapon of choice is miniconda3 (https://docs.conda.io/en/latest/miniconda.html). For this code more than ever it is a good idea to use `conda`, as code makes heavy use of `iris` which is difficult(/impossible?) to install through `pip`. Therefore, create a new environment, activate it, and use `conda` to grab `iris`:

    conda install -c conda-forge iris

When this is done, grab the repo-specific dependencies:

    pip install -r requirements.txt

## Reference
The paper describing this method is under review at the moment. The preprint can be accessed here: https://www.essoar.org/doi/abs/10.1002/essoar.10503977.2
