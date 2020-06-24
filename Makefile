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
	pytest src/tests

## Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

#################################################################################
# ANALYSIS                                                                      #
#################################################################################

## Download data to correct directories
data-vdjdb-aug-2019:
	@echo ">>> Downloading raw data."
	bash ./src/scripts/data_scripts/retrieve_data.sh
	@echo ">>> Creating summary statistics."
	bash ./src/scripts/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-2019-08-08-vdjdb-summary.md data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt
	bash ./src/scripts/data_scripts/vdjdb-content-analyser-slim.sh data/raw/vdjdb/vdjdb-2019-08-08-vdjdb-slim-summary.md data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.slim.txt
	bash ./src/scripts/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-browser-summary.md data/raw/vdjdb/vdjdb-browser.tsv

preprocess-vdjdb-aug-2019:
	mkdir -p data/interim/vdjdb-2019-08-08/

	## full untouched dataset for negative data generation through shuffling
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human.csv --species human --drop-spurious

	## NO DOWNSAMPLING
	### MHCI
	# 2019-08-08 release: human TRB MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# 2019-08-08 release: human TRA MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-mhci-no10x-size.csv --species human --tcr-chain TRA --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# 2019-08-08 release: human TRA+TRB MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhci-no10x-size.csv --species human --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# ### MHCII
	# # 2019-08-08 release: human TRB MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhcii-no10x-size.csv --species human --tcr-chain TRB --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# # 2019-08-08 release: human TRA MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-mhcii-no10x-size.csv --species human --tcr-chain TRA --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# # 2019-08-08 release: human TRA+B MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhcii-no10x-size.csv --species human --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# ### MHCI+II
	# # 2019-08-08 release: human TRB without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-no10x-size.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# 	# 2019-08-08 release: human TRB without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-no10x-size-preprint.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 13

	# # 2019-08-08 release: human TRA without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-no10x-size.csv --species human --tcr-chain TRA --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	# # 2019-08-08 release: human TRA+B without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x-size.csv --species human --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11

	## DOWNSAMPLING
	### MHCI
	# 2019-08-08 release: human TRB MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 1000 GILGFVFTL 1000
	# cut -f2 vdjdb-human-trb-mhci-no10x-size.csv -d';' | sort  | uniq -c | sort -n
	# 171 KAFSPEVIPMF
	# 172 RAKFKQLL
	# 191 TPRVTGGGAM
	# 199 VTEHDTLLY
	# 231 LLLGIGILV
	# 316 KRWIILGLNK
	# 404 LLWNGPMAV
	# 883 GLCTLVAML
	# 960 ELAGIGILTV
	# 2856 GILGFVFTL
	# 4387 NLVPMVATV

	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 400 GILGFVFTL 400 ELAGIGILTV 400 GLCTLVAML 400

	# 2019-08-08 release: human TRA MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-mhci-no10x-size-down.csv --species human --tcr-chain TRA --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 300 GILGFVFTL 300
	# cut -f2 vdjdb-human-tra-mhci-no10x-size.csv -d';' | sort  | uniq -c | sort -n
	# 22 KAFSPEVIPMF
	# 25 KLSALGINAV
	# 42 ELAGIGILTV
	# 43 KLVALGINAV
	# 48 NEGVKAAW
	# 69 CINGVCWTV
	# 108 GLCTLVAML
	# 245 LLWNGPMAV
	# 330 LLLGIGILV
	# 2065 NLVPMVATV
	# 2433 GILGFVFTL

	# 2019-08-08 release: human TRA+B MHCI without spurious sequences and without any 10xgenomics entries and length restrictions
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhci-no10x-size-down.csv --species human --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 1000 GILGFVFTL 1000
	# cut -f2 vdjdb-human-tra-trb-mhci-no10x-size.csv -d';' | sort  | uniq -c | sort -n
	# 172 RAKFKQLL
	# 193 KAFSPEVIPMF
	# 197 TPRVTGGGAM
	# 199 CINGVCWTV
	# 199 VTEHDTLLY
	# 322 KRWIILGLNK
	# 561 LLLGIGILV
	# 649 LLWNGPMAV
	# 991 GLCTLVAML
	# 1002 ELAGIGILTV
	# 5289 GILGFVFTL
	# 6452 NLVPMVATV

	# ### MHCII
	# # 2019-08-08 release: human TRB MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhcii-no10x-size-down.csv --species human --tcr-chain TRB --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# # 2019-08-08 release: human TRA MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-mhcii-no10x-size-down.csv --species human --tcr-chain TRA --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# # 2019-08-08 release: human TRA+B MHCII without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhcii-no10x-size-down.csv --species human --mhc MHCII --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# ### MHCI+II
	# # 2019-08-08 release: human TRB without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-no10x-size-down.csv --species human --tcr-chain TRB --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# # 2019-08-08 release: human TRA without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-no10x-size-down.csv --species human --tcr-chain TRA --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	# # 2019-08-08 release: human TRA+B without spurious sequences and without any 10xgenomics entries and length restrictions
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x-size-down.csv --species human --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 0.84 GILGFVFTL 0.80

	## DECOY
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/decoy_epitopes.py -i data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size.csv -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-decoy.csv

	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/decoy_epitopes.py -i data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down-decoy.csv

	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/decoy_epitopes.py -i data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400.csv -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400-decoy.csv

	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/decoy_epitopes.py -i data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhci-no10x-size.csv -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhci-no10x-size-decoy.csv

	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/decoy_epitopes.py -i data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhci-no10x-size-down.csv -o data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-mhci-no10x-size-down-decoy.csv

	# ## SINGLE EPITOPES
	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-size-NLVPMVATV.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --keep-specific-epitopes NLVPMVATV

	# $(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-size-GILGFVFTL.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --keep-specific-epitopes GILGFVFTL

	## REPLACED
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-terminal-replaced.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --terminal_replaced G
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-terminal-replaced-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 1000 GILGFVFTL 1000 --terminal_replaced G

	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-terminal-only.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --terminal_only
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-terminal-only-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 1000 GILGFVFTL 1000 --terminal_only

	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-middle-replaced.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --middle_replaced G
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-middle-replaced-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 1000 GILGFVFTL 1000 --middle_replaced G

	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-middle-only.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --middle_only
	$(PYTHON_INTERPRETER) ./src/scripts/preprocessing/preprocess_vdjdb.py -i data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt -o data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-middle-only-down.csv --species human --tcr-chain TRB --mhc MHCI --drop-spurious --remove-specific-reference 10xgenomics --length-restriction 10 20 8 11 --downsample NLVPMVATV 1000 GILGFVFTL 1000 --middle_only

## Download dataset
data-vdjdb-jan-2020:
	@echo ">>> Downloading raw data."
	bash ./src/scripts/data_scripts/retrieve_data-2020-01-20.sh
	@echo ">>> Creating summary statistics."
	bash ./src/scripts/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-2020-01-20-vdjdb-summary.md data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt
	bash ./src/scripts/data_scripts/vdjdb-content-analyser-slim.sh data/raw/vdjdb/vdjdb-2020-01-20-vdjdb-slim-summary.md data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.slim.txt
	bash ./src/scripts/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-browser-summary.md data/raw/vdjdb/vdjdb-browser.tsv

## Calculate metrics across train/test folds and create figures
metrics:
	find models/models -maxdepth 1 -mindepth 1 -type d -exec python ./src/scripts/evaluate/visualize.py metrics --force True {} \;
	# optionally use --y_lim_loss 2 to reduce y axis to a max of 2 for a better view of loss curves

metrics-compare:
	@echo "Use the following one-liner to compare two directories with trained models."
	@echo 'python ./src/scripts/evaluate/visualize.py compare --force True parent_directory'
	# optionally use --y_lim_loss 2 to reduce y axis to a max of 2 for a better view of loss curves

## Per-epitope evaluation
evaluate_self:
	@echo "Adjust the following one-liners. -name: should contain the name of the models of the same type. --model_type: padded or separated. --features: the features used to construct the padded model, order matters!"
	@echo 'find models/models -maxdepth 1 -mindepth 1 -name "*padded*" -type d -exec python ./src/scripts/evaluate/evaluate_self.py --input {} --model_type padded --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features "atchley1,atchley2,atchley3,atchley4,atchley5" \;'
	@echo 'find models/models -maxdepth 1 -mindepth 1 -name "*nettcr*" -type d -exec python ./src/scripts/evaluate/evaluate_self.py --input {} --model_type separated --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \;'
	@echo $(PROJECT_ROOT)

## Create interaction-map example
example-figure:
	$(PYTHON_INTERPRETER) ./src/model_scripts/visualize.py peptide --epitope ELAGIGILTV --cdr3 CASSPGEGLYEQYF --operator absdiff --cmyk True
