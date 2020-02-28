.PHONY: clean data lint requirements #sync_data_to_s3 sync_data_from_s3

#################################################################################
# GLOBALS                                                                       #
#################################################################################

# https://stackoverflow.com/questions/322936/common-gnu-makefile-directory-path
# http://andylinuxblog.blogspot.com/2015/06/what-is-colon-equals-sign-in-makefiles.html
# ":=" ensures make ensures the expansion happens immediately  at the start of the makefile,
# instead of when the variable is used in one of the make commands.
PROJECT_ROOT := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
# BUCKET = [OPTIONAL] your-bucket-for-syncing-data (do not include 's3://')
# PROFILE = default
PROJECT_NAME = deepTCR
PYTHON_INTERPRETER = python3

ifeq (,$(shell which conda))
HAS_CONDA=False
else
HAS_CONDA=True
endif

#################################################################################
# COMMANDS                                                                      #
#################################################################################

## Install Python Dependencies
requirements: test_environment
	ifeq (True,$(HAS_CONDA))
		@echo ">>> Detected conda, installing requirements inside conda environment."

	# $(PYTHON_INTERPRETER) -m pip install -U pip setuptools wheel
	# $(PYTHON_INTERPRETER) -m pip install -r requirements.txt

## Download dataset
data-vdjdb: requirements
	# download data
	bash src/data_scripts/retrieve_data.sh
	# create summary statistics
	bash src/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-2019-08-08-vdjdb-summary.md data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt
	bash src/data_scripts/vdjdb-content-analyser-slim.sh data/raw/vdjdb/vdjdb-2019-08-08-vdjdb-slim-summary.md data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.slim.txt
	bash src/data_scripts/vdjdb-content-analyser.sh data/raw/vdjdb/vdjdb-browser-summary.md data/raw/vdjdb/vdjdb-browser.tsv

# Preprocess dataset
preprocess-vdjdb: data-vdjdb
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

## Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

## Run tests
test:
	pytest

## Lint using flake8
lint:
	flake8 src

## Upload Data to S3
# sync_data_to_s3:
# ifeq (default,$(PROFILE))
# 	aws s3 sync data/ s3://$(BUCKET)/data/
# else
# 	aws s3 sync data/ s3://$(BUCKET)/data/ --profile $(PROFILE)
# endif

## Download Data from S3
# sync_data_from_s3:
# ifeq (default,$(PROFILE))
# 	aws s3 sync s3://$(BUCKET)/data/ data/
# else
# 	aws s3 sync s3://$(BUCKET)/data/ data/ --profile $(PROFILE)
# endif

## Set up python interpreter environment
create_environment:
ifeq (True,$(HAS_CONDA))
		@echo ">>> Detected conda, creating conda environment."
ifeq (3,$(findstring 3,$(PYTHON_INTERPRETER)))
	conda create --name $(PROJECT_NAME) python=3
else
	conda create --name $(PROJECT_NAME) python=2.7
endif
		@echo ">>> New conda env created. Activate with:\nsource activate $(PROJECT_NAME)"
else
	$(PYTHON_INTERPRETER) -m pip install -q virtualenv virtualenvwrapper
	@echo ">>> Installing virtualenvwrapper if not already installed.\nMake sure the following lines are in shell startup file\n\
	export WORKON_HOME=$$HOME/.virtualenvs\nexport PROJECT_HOME=$$HOME/Devel\nsource /usr/local/bin/virtualenvwrapper.sh\n"
	@bash -c "source `which virtualenvwrapper.sh`;mkvirtualenv $(PROJECT_NAME) --python=$(PYTHON_INTERPRETER)"
	@echo ">>> New virtualenv created. Activate with:\nworkon $(PROJECT_NAME)"
endif

## Test python environment is setup correctly
test_environment:
	$(PYTHON_INTERPRETER) test_environment.py

#################################################################################
# PROJECT RULES                                                                 #
#################################################################################



#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: help
help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
