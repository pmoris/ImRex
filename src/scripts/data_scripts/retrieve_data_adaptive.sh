#! /usr/bin/env bash

# retrieve project root (or pass project root from makefile instead?)
# ${parameter:-word} If parameter is unset or null, the expansion of word is substituted. Otherwise, the value of parameter is substituted.
# i.e. if PROJECT_ROOT is not passed through via the Makefile (or when the script is run on its own)
# the oneliner will attempt to move up two directories from the location of this script and set that as the project root.
# A potential problem with this approach is that when PROJECT_ROOT is not a resolved path, RAW_DATA_PATH will not be valid.
# However, normally this script will either be used on its own (and PROJECT_ROOT will be empty initially)
# or it will be invoked by the Makefile, which should pass a fully resolved path.
# To be safe, the PROJECT_ROOT variable is assigned a second time as resolved path.
# Note that "readlink -f/--canonicalize" or "realpath" are not used because it does not exist on all systems.
PROJECT_ROOT=${PROJECT_ROOT:-"$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../.." >/dev/null 2>&1 && pwd )"}
# resolve path in case the passed variable contains a relative one
PROJECT_ROOT="$(cd ${PROJECT_ROOT} >/dev/null 2>&1 && pwd )"

# define download path based on project root dir
RAW_DATA_PATH="${PROJECT_ROOT}/data/raw/"

# move to project root dir to run all commands
cd ${PROJECT_ROOT}

echo "Project directory: ${PROJECT_ROOT}"
echo "Running from $(pwd)"
echo "Path used for raw vdjdb data directory: $(realpath ${RAW_DATA_PATH})"

echo -e "\nPlease download the Adaptive ImmuneCODE dataset (June 25 2020) from the following URL: https://immunerace.adaptivebiotech.com/more-data-and-whats-coming-next/"
echo -e "The download location should read https://enrollimmunerace.adaptivebiotech.com/restricted/ImmuneCODE-Release001.1.zip."
echo -e "The file should be placed in the directory './data/raw/ImmuneCODE-Release001.1.zip'."
echo -e "Once the file is present, press 1 (Yes) to proceed or 2 (No) to cancel this operation."

select yn in "Yes" "No"; do
    case $yn in
        Yes ) echo -e "Unzipping ImmuneCODE zip file..."; break;;
        No ) exit;;
    esac
done

# extract files
unzip -u -o "${RAW_DATA_PATH}ImmuneCODE-Release001.1.zip" -d "${RAW_DATA_PATH}immunecode-adaptive"

echo -e "Extracted data to ${RAW_DATA_PATH}immunecode-adaptive.\nPlease run the processing script next (make preprocess-adaptive)."
