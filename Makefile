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

## Calculate metrics across train/test folds and create figures
metrics:
	find models/models -maxdepth 1 -mindepth 1 -type d -exec python src/model_scripts/visualize.py metrics --force True {} \;

## Download data to correct directories
data-vdjdb-aug-2019:
	@echo ">>> Downloading raw data."
	bash src/data_scripts/retrieve_data.sh
	@echo ">>> Creating summary statistics."
	bash src/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-2019-08-08-vdjdb-summary.md data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt
	bash src/data_scripts/vdjdb-content-analyser-slim.sh data/raw/vdjdb/vdjdb-2019-08-08-vdjdb-slim-summary.md data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.slim.txt
	bash src/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-browser-summary.md data/raw/vdjdb/vdjdb-browser.tsv

preprocess-vdjdb-aug-2019:
	mkdir -p data/interim/vdjdb-2019-08-08/

	## NO DOWNSAMPLING
	### MHCI
	# 2019-08-08 release: human TRB MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# 2019-08-08 release: human TRA MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-mhci-no10x-size.csv --species human --tcr-chain TRA --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# 2019-08-08 release: human TRA+B MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhci-no10x-size.csv --species human --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	### MHCII
	# 2019-08-08 release: human TRB MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhcii-no10x-size.csv --species human --tcr-chain TRB --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# 2019-08-08 release: human TRA MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-mhcii-no10x-size.csv --species human --tcr-chain TRA --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# 2019-08-08 release: human TRA+B MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhcii-no10x-size.csv --species human --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	### MHCI+II
	# 2019-08-08 release: human TRB without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-no10x-size.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

		# 2019-08-08 release: human TRB without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-no10x-size-preprint.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 13

	# 2019-08-08 release: human TRA without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-no10x-size.csv --species human --tcr-chain TRA --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# 2019-08-08 release: human TRA+B without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x-size.csv --species human --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	## DOWNSAMPLING
	### MHCI
	# 2019-08-08 release: human TRB MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# 2019-08-08 release: human TRA MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-mhci-no10x-size-down.csv --species human --tcr-chain TRA --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# 2019-08-08 release: human TRA+B MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhci-no10x-size-down.csv --species human --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	### MHCII
	# 2019-08-08 release: human TRB MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhcii-no10x-size-down.csv --species human --tcr-chain TRB --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# 2019-08-08 release: human TRA MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-mhcii-no10x-size-down.csv --species human --tcr-chain TRA --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# 2019-08-08 release: human TRA+B MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhcii-no10x-size-down.csv --species human --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	### MHCI+II
	# 2019-08-08 release: human TRB without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-no10x-size-down.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# 2019-08-08 release: human TRA without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-no10x-size-down.csv --species human --tcr-chain TRA --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# 2019-08-08 release: human TRA+B without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x-size-down.csv --species human --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	## DECOY
	$(PYTHON_INTERPRETER) src/preprocessing/decoy_epitopes.py -i data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size.csv -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-decoy.csv

	$(PYTHON_INTERPRETER) src/preprocessing/decoy_epitopes.py -i data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down-decoy.csv

	$(PYTHON_INTERPRETER) src/preprocessing/decoy_epitopes.py -i data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x-size.csv -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x-size-decoy.csv

	$(PYTHON_INTERPRETER) src/preprocessing/decoy_epitopes.py -i data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x-size-down.csv -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x-size-down-decoy.csv

	## SINGLE EPITOPES
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-size-NLVPMVATV.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --keep-specific-epitopes NLVPMVATV

	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-size-GILGFVFTL.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --keep-specific-epitopes GILGFVFTL

	## REPLACED
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-terminal-replaced.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --terminal_replaced G
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-terminal-replaced-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80 --terminal_replaced G

	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-terminal-only.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --terminal_only
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-terminal-only-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80 --terminal_only

	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-middle-replaced.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --middle_replaced G
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-middle-replaced-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80 --middle_replaced G

	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-middle-only.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --middle_only
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-middle-only-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80 --middle_only

preprocess-vdjdb-aug-2019-old:
	mkdir -p data/interim/vdjdb-2019-08-08-old/
	## $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-all-species-tra-trb-non-paired.tsv -o data/interim/vdjdb-human.csv --species human
	## $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-all-species-tra-trb-non-paired.tsv -o data/interim/vdjdb-human-trb.csv --tcr-chain TRB --species human
	## $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human.csv --species human
	## $(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human-trb.csv --species human --tcr-chain TRB

	# 08-08-2019 release: human TRA+TRB without spurious sequences and without the overabundant KLGGALQAK 10xgenomics epitope
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human.csv --species human --drop-spurious --remove-specific-epitope-reference KLGGALQAK 10xgenomics
	# 08-08-2019 release: human TRB-only without spurious sequences and without the overabundant KLGGALQAK 10xgenomics epitope
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human-trb.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-epitope-reference KLGGALQAK 10xgenomics

	# 08-08-2019 release: human TRB-only without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human-trb-no10x.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics
	# 08-08-2019 release: human TRA+TRB without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human-no10x.csv --species human --tcr-chain all --drop-spurious --remove-specific-reference 10xgenomics

	# 08-08-2019 release: all species TRA+TRB without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-no10x.csv --species all --tcr-chain all --drop-spurious --remove-specific-reference 10xgenomics
	# 08-08-2019 release: all species TRB-only without spurious sequences and without any 10xgenomics entries
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-trb-no10x.csv --species all --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics
	# 08-08-2019 release: all species TRB-only without spurious sequences and without the overabundant KLGGALQAK 10xgenomics epitope
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-trb.csv --species all --tcr-chain TRB --drop-spurious --remove-specific-epitope-reference KLGGALQAK 10xgenomics

	# 08-08-2019 release: human TRB-only without spurious sequences and without any 10xgenomics and MHCI only
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human-trb-mhci-no10x.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --mhc MHCI

	# extract 10x data - human TRB-only
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human-trb-10xonly.csv --species human --tcr-chain TRB --drop-spurious --keep-specific-references 10xgenomics
	# extract 10x data - human TRB-only MHCI-only
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human-trb-mhci-10xonly.csv --species human --tcr-chain TRB --drop-spurious --keep-specific-references 10xgenomics --mhc MHCI

	# 08-08-2019 release: human TRB-only without spurious sequences and without any 10xgenomics entries + length filter
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human-trb-no10x-size.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 13
	# 08-08-2019 release: human TRA+TRB without spurious sequences and without any 10xgenomics entries + length filter
	$(PYTHON_INTERPRETER) src/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08-old/vdjdb.txt -o data/interim/vdjdb-human-no10x-size.csv --species human --tcr-chain all --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 13

## Download dataset
data-vdjdb-jan-2020:
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
