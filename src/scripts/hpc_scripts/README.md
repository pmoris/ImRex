This directory contains `.pbs` scripts for different experiments. Each script trains models on different datasets and with different architectures.

The different `.pbs` scripts must be run from within their own directory, since they attempt to retrieve the root of the project based on the location from which the script is called.

The output of these scripts will be saved in the `./models` directory. The README in that directory (`./models/README.md`) gives a more in-depth description of the different data/architecture combinations and the naming conventions that were used.

Log files are stored in the same directory as the `.pbs` scripts.

The results of these experiments are available in the associated Zenodo repository: [10.5281/zenodo.3973547](https://doi.org/10.5281/zenodo.3973547).
