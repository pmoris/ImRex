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
	@echo ">>> Conda not detected, installing requirements inside virtualenv. NOTE: cudatoolkit, cudnn and hdf5 need to be installed manually."
	#$(PYTHON_INTERPRETER) -m pip install -U pip setuptools wheel
	virtualenv $(PROJECT_NAME)
	. $(PROJECT_NAME)/bin/activate
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

data-vdjdb-aug-2019: requirements
	@echo ">>> Downloading raw data."
	bash src/data_scripts/retrieve_data.sh
	@echo ">>> Creating summary statistics."
	bash src/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-2019-08-08-vdjdb-summary.md data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt
	bash src/data_scripts/vdjdb-content-analyser-slim.sh data/raw/vdjdb/vdjdb-2019-08-08-vdjdb-slim-summary.md data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.slim.txt
	bash src/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-browser-summary.md data/raw/vdjdb/vdjdb-browser.tsv

preprocess-vdjdb-aug-2019-new:
	mkdir -p data/interim/vdjdb-2019-08-08-new/

	## NO DOWNSAMPLING
	### MHCI
	# 2019-08-08 release: human TRB MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08-new/vdjdb-human-trb-mhci-no10x-size.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# 2019-08-08 release: human TRA MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08-new/vdjdb-human-tra-mhci-no10x-size.csv --species human --tcr-chain TRA --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# 2019-08-08 release: human TRA+B MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08-new/vdjdb-human-tra-trb-mhci-no10x-size.csv --species human --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	## DOWNSAMPLING
	### MHCI
	# 2019-08-08 release: human TRB MHCI without spurious sequences and without any 10xgenomics entries and length restrictions and downsampling
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08-new/vdjdb-human-trb-mhci-no10x-size-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# 2019-08-08 release: human TRA MHCI without spurious sequences and without any 10xgenomics entries and length restrictions and downsampling
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08-new/vdjdb-human-tra-mhci-no10x-size-down.csv --species human --tcr-chain TRA --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# 2019-08-08 release: human TRA+B MHCI without spurious sequences and without any 10xgenomics entries and length restrictions and downsampling
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08-new/vdjdb-human-tra-trb-mhci-no10x-size-down.csv --species human --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

preprocess-vdjdb-aug-2019-pre:
	mkdir -p data/interim/vdjdb-2019-08-08/
	## $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-all-species-tra-trb-non-paired.tsv -o data/interim/vdjdb-human.csv --species human
	## $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-all-species-tra-trb-non-paired.tsv -o data/interim/vdjdb-human-trb.csv --tcr-chain TRB --species human
	## $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human.csv --species human
	## $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-human-trb.csv --species human --tcr-chain TRB

	# 08-08-2019 release: human TRA+TRB without spurious sequences and without the overabundant KLGGALQAK 10xgenomics epitope
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-/vdjdb.txt -o data/interim/vdjdb-human.csv --species human --drop-spurious --remove-specific-epitope-reference KLGGALQAK 10xgenomics
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

## Download dataset
data-vdjdb-jan-2020: requirements
	@echo ">>> Downloading raw data."
	bash src/data_scripts/retrieve_data-2020-01-20.sh
	@echo ">>> Creating summary statistics."
	bash src/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-2020-01-20-vdjdb-summary.md data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt
	bash src/data_scripts/vdjdb-content-analyser-slim.sh data/raw/vdjdb/vdjdb-2020-01-20-vdjdb-slim-summary.md data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.slim.txt
	bash src/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-browser-summary.md data/raw/vdjdb/vdjdb-browser.tsv

## Preprocess vdjdb dataset
preprocess-vdjdb-jan-2020:
	mkdir -p data/interim/vdjdb-2020-01-20/
	# 2020-01-20 release: human TRB-only without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt -o data/interim/vdjdb-2020-01-20/vdjdb-human-trb-no10x.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics
	# 2020-01-20 release: human TRA+TRB without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt -o data/interim/vdjdb-2020-01-20/vdjdb-human-tra-trb-no10x.csv --species human --tcr-chain all --drop-spurious --remove-specific-reference 10xgenomics

	# 2020-01-20 release: human TRB-only without spurious sequences and without any 10xgenomics entries + length filter
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt -o data/interim/vdjdb-2020-01-20/vdjdb-human-trb-no10x-size.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 13
	# 2020-01-20 release: human TRA+TRB without spurious sequences and without any 10xgenomics entries + length filter
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt -o data/interim/vdjdb-2020-01-20/vdjdb-human-tra-trb-no10x-size.csv --species human --tcr-chain all --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 13

	# 2020-01-20 release: all species TRB-only without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt -o data/interim/vdjdb-2020-01-20/vdjdb-allspecies-trb-no10x.csv --species all --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics
	# 2020-01-20 release: all species TRA+TRB without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt -o data/interim/vdjdb-2020-01-20/vdjdb-allspecies-tra-trb-no10x.csv --species all --tcr-chain all --drop-spurious --remove-specific-reference 10xgenomics


	# 2020-01-20 release: human TRB-only without spurious sequences and without any 10xgenomics and MHCI only
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt -o data/interim/vdjdb-2020-01-20/vdjdb-human-trb-mhci-no10x.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --mhc MHCI
	# ?

	# 2020-01-20 release: human TRB-only without spurious sequences and only 10x data
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt -o data/interim/vdjdb-2020-01-20/vdjdb-human-trb-10xonly.csv --species human --tcr-chain TRB --drop-spurious --keep-specific-references 10xgenomics

	# 2020-01-20 release: human TRB-only and MHCI-only without spurious sequences and only 10x data
	# $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt -o data/interim/vdjdb-2020-01-20/vdjdb-human-trb-mhci-10xonly.csv --species human --tcr-chain TRB --drop-spurious --keep-specific-references 10xgenomics --mhc MHCI

## Create interaction-map example
example-figure:
	$(PYTHON_INTERPRETER) src/model_scripts/visualize.py peptide --epitope ELAGIGILTV --cdr3 CASSPGEGLYEQYF --operator absdiff --cmyk True

## Summarise models
model-summary:
	find models/models/ -maxdepth 1 -type d -exec python src/model_scripts/visualize.py metrics {} \;
	# echo "TODO: use find to gather all model directories, then exec {}+ python src/model_scripts/visualize.py metrics"
