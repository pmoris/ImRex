#################################################################################
# GLOBALS                                                                       #
#################################################################################

# https://stackoverflow.com/questions/322936/common-gnu-makefile-directory-path
# http://andylinuxblog.blogspot.com/2015/06/what-is-colon-equals-sign-in-makefiles.html
# ":=" ensures make ensures the expansion happens immediately  at the start of the makefile,
# instead of when the variable is used in one of the make commands.
PROJECT_ROOT := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
PROJECT_NAME = deepTCR
PYTHON_INTERPRETER = python3

ifeq (,$(shell which conda))
HAS_CONDA=False
else
HAS_CONDA=True
endif

#################################################################################
# SETUP COMMANDS                                                                #
#################################################################################

## Install Python Dependencies
requirements: test_environment
ifeq (True,$(HAS_CONDA))
	@echo ">>> Detected conda, installing requirements inside conda environment."
	conda env create -f environment.yml
	@echo ">>> New conda env created. Activate with:\nsource activate $(PROJECT_NAME)"
else
	@echo ">>> Conda not detected, installing requirements inside virtualenv."
	$(PYTHON_INTERPRETER) -m pip install -U pip setuptools wheel
	virtualenv $(PROJECT_NAME)
	source $(PROJECT_NAME)/bin/activate
	pip install -r requirements.txt
endif

## Test python environment is setup correctly
test_environment:
	$(PYTHON_INTERPRETER) test_environment.py

## Format src code with black
format:
	black src

## Lint using flake8
lint:
	flake8 src

## Run tests
test:
	pytest

## Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

#################################################################################
# ANALYSIS                                                                      #
#################################################################################

## Download dataset
data-vdjdb: requirements
	@echo ">>> Downloading raw data."
	bash src/data_scripts/retrieve_data.sh
	@echo ">>> Creating summary statistics."
	bash src/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-2019-08-08-vdjdb-summary.md data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt
	bash src/data_scripts/vdjdb-content-analyser-slim.sh data/raw/vdjdb/vdjdb-2019-08-08-vdjdb-slim-summary.md data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.slim.txt
	bash src/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-browser-summary.md data/raw/vdjdb/vdjdb-browser.tsv

# Preprocess vdjdb dataset
preprocess-vdjdb:
	# $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-all-species-tra-trb-non-paired.tsv -o data/interim/vdjdb-human.csv --species human
	# $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-all-species-tra-trb-non-paired.tsv -o data/interim/vdjdb-human-trb.csv --tcr-chain TRB --species human
	# $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human.csv --species human
	# $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human-trb.csv --species human --tcr-chain TRB

	# 08-08-2019 release: human TRA+TRB without spurious sequences and without the overabundant KLGGALQAK 10xgenomics epitope
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human.csv --species human --drop-spurious --remove-specific-epitope-reference KLGGALQAK 10xgenomics
	# 08-08-2019 release: human TRB-only without spurious sequences and without the overabundant KLGGALQAK 10xgenomics epitope
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human-trb.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-epitope-reference KLGGALQAK 10xgenomics

	# 08-08-2019 release: human TRB-only without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human-trb-no10x.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics
	# 08-08-2019 release: human TRA+TRB without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human-no10x.csv --species human --tcr-chain all --drop-spurious --remove-specific-reference 10xgenomics

	# 08-08-2019 release: all species TRA+TRB without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-no10x.csv --species all --tcr-chain all --drop-spurious --remove-specific-reference 10xgenomics
	# 08-08-2019 release: all species TRB-only without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-trb-no10x.csv --species all --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics
	# 08-08-2019 release: all species TRB-only without spurious sequences and without the overabundant KLGGALQAK 10xgenomics epitope
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-trb.csv --species all --tcr-chain TRB --drop-spurious --remove-specific-epitope-reference KLGGALQAK 10xgenomics

	# 08-08-2019 release: human TRB-only without spurious sequences and without any 10xgenomics and MHCI only
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human-trb-mhci-no10x.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --mhc MHCI

	# extract 10x data - human TRB-only
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human-trb-10xonly.csv --species human --tcr-chain TRB --drop-spurious --keep-specific-references 10xgenomics
	# extract 10x data - human TRB-only MHCI-only
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human-trb-mhci-10xonly.csv --species human --tcr-chain TRB --drop-spurious --keep-specific-references 10xgenomics --mhc MHCI

	# 08-08-2019 release: human TRB-only without spurious sequences and without any 10xgenomics entries + length filter
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human-trb-no10x-size.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 13
	# 08-08-2019 release: human TRA+TRB without spurious sequences and without any 10xgenomics entries + length filter
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human-no10x-size.csv --species human --tcr-chain all --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 13

# Create interaction-map example
example-figure:
	$(PYTHON_INTERPRETER) src/model_scripts/visualize.py peptide --epitope ELAGIGILTV --cdr3 CASSPGEGLYEQYF --operator absdiff --cmyk True

# Summarise models
model-summary:
	find models/models/ -maxdepth 1 -type d -exec python src/model_scripts/visualize.py metrics {} \;
	# echo "TODO: use find to gather all model directories, then exec {}+ python src/model_scripts/visualize.py metrics"
