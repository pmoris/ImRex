# Raw data files

## TCR-epitope sequences

- `./vdjdb`: contains a data dump of VDJdb TCR-epitope pairs ([https://vdjdb.cdr3.net](https://vdjdb.cdr3.net)), see [./VDJdb/README.md](./VDJdb/README.md) for additional information and statistics of this file. This data can be retrieved and unpacked using `Make data-vdjdb-aug-2019`.
- `CDR3_control_sequences.tsv`: 500,000 TRB CDR3 sequences from a study
by Dean et al. ([https://doi.org/10.1186/s13073-015-0238-z](https://doi.org/10.1186/s13073-015-0238-z)).

## Adaptive ImmuneCODE

- `./immunecode-adaptive` contains the SARS-CoV-2 data (June 25 2020) released by Adaptive (source: [https://immunerace.adaptivebiotech.com/more-data-and-whats-coming-next/](https://immunerace.adaptivebiotech.com/more-data-and-whats-coming-next/)). This data can be retrieved and unpacked using `Make data-adaptive`.

## McPAS

- `mcpas.csv` contains human CDR3-epitope sequence pairs, taken from McPAS ([http://friedmanlab.weizmann.ac.il/McPAS-TCR/](http://friedmanlab.weizmann.ac.il/McPAS-TCR/)) and filtered on a number of quality checks.
