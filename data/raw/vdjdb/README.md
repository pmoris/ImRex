# VDJdb dataset retrieval

The files in this directory were retrieved from https://vdjdb.cdr3.net on 12 November 2019. At this time, the last update to VDJdb had happened on 07 August 2019.

## vdjdb-all-species-tra-trb-non-paired.tsv

Contains TRA and TRB chains and their epitope targets for all species in the database (human, monkey and mouse).

Overview of contents:

|Metric|Count|Command|
|---|---|---|
|Total number of records|73224|`tail -n +2 vdjdb-all-species-tra-trb-non-paired.tsv |  wc -l`|
|TRA records|29479|`awk '$2 == "TRA" { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | wc -l`|
|TRB records|43745|`awk '$2 == "TRB" { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | wc -l`|
|Unique TRA sequences|20888|`awk '$2 == "TRA" { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique TRB sequences|33112|`awk '$2 == "TRB" { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique CDR3 sequences|54000|`tail -n +2 vdjdb-all-species-tra-trb-non-paired.tsv | cut -f3 | sort -u | wc -l`|
|Unique epitope sequences|212|`tail -n +2 vdjdb-all-species-tra-trb-non-paired.tsv | cut -f10 | sort -u | wc -l`|
|Unique epitope sequences for TRA records|141|`awk '$2 == "TRA" { print $10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique epitope sequences for TRB records|210|`awk '$2 == "TRB" { print $10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique CDR3-epitope sequence pairs|59072|`tail -n +2 vdjdb-all-species-tra-trb-non-paired.tsv | cut -d $'\t' -f3,10 | sort -u | wc -l`|
|Unique TRA-CDR3-epitope sequence pairs|23578|`awk '$2 == "TRA" { print $3,$10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique TRB-CDR3-epitope sequence pairs|35494|`awk '$2 == "TRB" { print $3,$10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Number of epitope sequences shared between TRA and TRB records|139|`comm -12 <(awk '$2 == "TRA" { print $10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) <(awk '$2 == "TRB" { print $10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) | wc -l`|
|Number of CDR3 sequences shared between TRA and TRB records|0|`comm -12 <(awk '$2 == "TRA" { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) <(awk '$2 == "TRB" { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) | wc -l`|

**Human-only records**
|Metric|Count|Command|
|---|---|---|
|Total number of records|66597|`awk '$6 == "HomoSapiens" { print }' vdjdb-human-tra-trb-non-paired.tsv |  wc -l`|
|TRA records|27295|`awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | wc -l`|
|TRB records|39302|`awk '$2 == "TRB" && $6 == "HomoSapiens"  { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | wc -l`|
|Unique TRA sequences|19637|`awk '$2 == "TRA" && $6 == "HomoSapiens"  { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique TRB sequences|30377|`awk '$2 == "TRB" && $6 == "HomoSapiens"  { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique CDR3 sequences|50014|`awk '$6 == "HomoSapiens" { print }' vdjdb-all-species-tra-trb-non-paired.tsv | cut -f3 | sort -u | wc -l`|
|Unique epitope sequences|177|`awk '$6 == "HomoSapiens" { print }' vdjdb-all-species-tra-trb-non-paired.tsv | cut -f10 | sort -u | wc -l`|
|Unique epitope sequences for TRA records|117|`awk '$2 == "TRA" && $6 == "HomoSapiens" { print $10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique epitope sequences for TRB records|175|`awk '$2 == "TRB" && $6 == "HomoSapiens" { print $10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique CDR3-epitope sequence pairs|54887|`awk '$6 == "HomoSapiens" { print }' vdjdb-all-species-tra-trb-non-paired.tsv | cut -d $'\t' -f3,10 | sort -u | wc -l`|
|Unique TRA-CDR3-epitope sequence pairs|22248|`awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3,$10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique TRB-CDR3-epitope sequence pairs|32639|`awk '$2 == "TRB" && $6 == "HomoSapiens" { print $3,$10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Number of epitope sequences shared between TRA and TRB records|115|`comm -12 <(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) <(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $10 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) | wc -l`|
|Number of CDR3 sequences shared between TRA and TRB records|0|`comm -12 <(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) <(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $3 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) | wc -l`|

These statistics can also be computed using the `vdjdb-content-analyser.sh` bash script.

## vdjdb-human-tra-trb-non-paired.tsv

Contains TRA and TRB chains and their epitope targets for humans.
