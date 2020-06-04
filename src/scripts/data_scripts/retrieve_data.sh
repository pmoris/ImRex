#! /usr/bin/env bash

# retrieve project root (or pass project root from makefile instead?)
# ${parameter:-word} If parameter is unset or null, the expansion of word is substituted. Otherwise, the value of parameter is substituted.
# i.e. if PROJECT_ROOT is not passed through via the Makefile (or when the script is run on its own)
# the oneliner will attempt to move up two directories from the location of this script and set that as the project root.
# A potential problem with this approach is that when PROJECT_ROOT is not a resolved path, VDJDB_RAW_DATA_PATH will not be valid.
# However, normally this script will either be used on its own (and PROJECT_ROOT will be empty initially)
# or it will be invoked by the Makefile, which should pass a fully resolved path.
# To be safe, the PROJECT_ROOT variable is assigned a second time as resolved path.
# Note that "readlink -f/--canonicalize" or "realpath" are not used because it does not exist on all systems.
PROJECT_ROOT=${PROJECT_ROOT:-"$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." >/dev/null 2>&1 && pwd )"}
# resolve path in case the passed variable contains a relative one
PROJECT_ROOT="$(cd ${PROJECT_ROOT} >/dev/null 2>&1 && pwd )"

# define download path based on project root dir
VDJDB_RAW_DATA_PATH="${PROJECT_ROOT}/data/raw/vdjdb/"

# move to project root dir to run all commands
cd ${PROJECT_ROOT}

echo "Project directory: ${PROJECT_ROOT}"
echo "Running from $(pwd)"
echo "Path used for raw vdjdb data directory: $(realpath ${VDJDB_RAW_DATA_PATH})"

echo -e "\nDownloading VDJdb data...\n"

# download VDJdb files
wget --timestamping -P ${VDJDB_RAW_DATA_PATH} https://github.com/antigenomics/vdjdb-db/releases/download/2019-08-08/vdjdb-2019-08-08.zip
# extract files
unzip -u -o "${VDJDB_RAW_DATA_PATH}/vdjdb-2019-08-08.zip" -d "${VDJDB_RAW_DATA_PATH}/vdjdb-2019-08-08"

echo -e "\nDownloaded and extracted the VDJdb GitHub release to ${VDJDB_RAW_DATA_PATH}.\n"

cat << EOF
By default, this project relies on the VDJdb GitHub releases.
If the VDJdb web browser dataset is required for comparison, it must be downloaded manually.
Use the 'Export as' button on https://vdjdb.cdr3.net/search in conjunction with the following filters:
- All species
- Both TRA and TRB chains
- Both MHC types
- No sequence length restriction
- All assay types
- No minimum confidence score
- No spurious CDR3 (non-canonical or unmapped V/J)

The file should be named 'vdjdb-browser.tsv' in order to be compatible with the other scripts in this project.
EOF
