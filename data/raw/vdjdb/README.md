# VDJdb dataset retrieval

The files in this directory were retrieved from VDJdb ([https://vdjdb.cdr3.net](https://vdjdb.cdr3.net)] and the associated GitHub release page ([https://github.com/antigenomics/vdjdb-db/releases](https://github.com/antigenomics/vdjdb-db/releases)) on 21 November 2019. At this time, the last update to VDJdb had happened on 07 August 2019 (GitHub release `vdjdb-2019-08-08.zip`).

## All analyses and scripts used in this project start from the normal `vdjdb.txt` file from the 2019-08-08 GitHub release!

See the [Jupyter notebook](../../../notebooks/vdjdb-statistics.ipynb) for a more complete exploration of the differences between all these files.

---

## GitHub release: [MAJOR UPDATE] Aug'19 release

- Source: [https://github.com/antigenomics/vdjdb-db/releases/tag/2019-08-08](https://github.com/antigenomics/vdjdb-db/releases/tag/2019-08-08)
- Expected target folder: `vdjdb-2019-08-08.zip`

**The `vdjdb.txt` file was used for all analyses.** This file contains all epitope-cdr3 pairs as separate entries (i.e. the "unpaired gene export" format when using the web export). The `vdjdb.slim.txt` version contains the same information as the standard file, apart from missing the `cdr3fix`, `meta` data and `method` columns, and having additional `j.start` and `v.end` columns instead (this information is normally contained inside the `cdr3fix` column). The `vdjdb_full.txt` version contains some additional information; among other things it uses a separate column for `cdr3.alpha` and `cdr3.beta` sequences.

<!-- Because of these additional columns, the normal version contains more duplicate entries than the slim version. -->

Note that while these three files all use the same type of header names, the order of the columns differs between the three versions, and not all columns are present in every file.

Lastly, these files contain so called *spurious CDR3* entries, which are filtered out of the VDJdb web results by default.

## VDJdb data browser

A dataset named `vdjdb-browser.tsv` containing CDR3-epitope pairs was retrieved from the VDJdb data browser on 21 November 2019 (Last updated on 07 August, 2019) by using the web interface with the following options:

- All species
- Both TRA and TRB chains
- Both MHC types
- No sequence length restriction
- All assay types
- No minimum confidence score
- No spurious CDR3 (non-canonical or unmapped V/J)

Note that datasets retrieved from the VDJdb data browser contain a different header format than the GitHub release (with additional capitalization and spaces), although the order of the CDR3 and epitope sequences is retained from the normal and full files described above (3 and 10 for CDR3 and epitope respectively).

# Statistics for GitHub release vdjdb.txt

## vdjdb.txt

| Metric                                                         | Count                                                                                                                                       | Command                                                                                                                                        |
|----------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| Total number of records                                        | 75474                                                                                                         | `tail -n +2 vdjdb-2019-08-08/vdjdb.txt  |  wc -l`                                                                                                         |
| TRA records                                                    | 31096                                                                                     | `awk '$2 == "TRA" { print $3 }' vdjdb-2019-08-08/vdjdb.txt | wc -l`                                                                                     |
| TRB records                                                    | 44378                                                                                     | `awk '$2 == "TRB" { print $3 }' vdjdb-2019-08-08/vdjdb.txt | wc -l`                                                                                     |
| Unique TRA sequences                                           | 22248                                                                           | `awk '$2 == "TRA" { print $3 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                           |
| Unique TRB sequences                                           | 33606                                                                           | `awk '$2 == "TRB" { print $3 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                           |
| Unique CDR3 sequences                                          | 55852                                                                                     | `tail -n +2 vdjdb-2019-08-08/vdjdb.txt | cut -f3 | sort -u | wc -l`                                                                                       |
| Unique epitope sequences                                       | 212                                                                                    | `tail -n +2 vdjdb-2019-08-08/vdjdb.txt | cut -f10 | sort -u | wc -l`                                                                                      |
| Unique epitope sequences for TRA records                       | 179                                                                          | `awk '$2 == "TRA" { print $10 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                          |
| Unique epitope sequences for TRB records                       | 219                                                                          | `awk '$2 == "TRB" { print $10 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                          |
| Unique CDR3-epitope sequence pairs                             | 61047                                                                         | `tail -n +2 vdjdb-2019-08-08/vdjdb.txt | cut -d $'\t' -f3,10 | sort -u | wc -l`                                                                           |
| Unique TRA-CDR3-epitope sequence pairs                         | 25044                                                                       | `awk '$2 == "TRA" { print ,vdjdb-2019-08-08-vdjdb-summary.md0 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                        |
| Unique TRB-CDR3-epitope sequence pairs                         | 36000                                                                       | `awk '$2 == "TRB" { print ,vdjdb-2019-08-08-vdjdb-summary.md0 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                        |
| Number of epitope sequences shared between TRA and TRB records | 147   | `comm -12 <(awk '$2 == "TRA" { print $10 }' vdjdb-2019-08-08/vdjdb.txt | sort -u) <(awk '$2 == "TRB" { print $10 }' vdjdb-2019-08-08/vdjdb.txt | sort -u) | wc -l` |
| Number of CDR3 sequences shared between TRA and TRB records    | 2     | `comm -12 <(awk '$2 == "TRA" { print $3 }' vdjdb-2019-08-08/vdjdb.txt | sort -u) <(awk '$2 == "TRB" { print $3 }' vdjdb-2019-08-08/vdjdb.txt | sort -u) | wc -l`   |
| Epitope distribution for the unique CDR3-epitope pairs         |                                 | `tail -n +2 vdjdb-2019-08-08/vdjdb.txt | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20`                                  |

    24639 KLGGALQAK
    6744 NLVPMVATV
    6600 GILGFVFTL
    3255 AVFDRKSDAK
    1654 ELAGIGILTV
    1501 RAKFKQLL
    1190 GLCTLVAML
    1088 IVTDFSVIK
    819 RLRAEAQVK
    685 LLWNGPMAV
    617 LLLGIGILV
    597 TTPESANL
    597 SSLENFRAYV
    586 SSYRRPVGI
    564 FRDYVDRFYKTLRAEQASQE
    548 PKYVKQNTLKLAT
    530 CTPYDINQM
    417 HGIRNASFI
    383 ASNENMETM
    334 KRWIILGLNK

**Human-only records**

| Metric                                                             | Count                                                                                                                                                                          | Command                                                                                                                                                                                                 |
|--------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|   Total number of records                                          | 67766                                                                                                                  | `awk '$6 == "HomoSapiens" { print }' vdjdb-2019-08-08/vdjdb.txt |  wc -l`                                                                                                                                         |
|   TRA records                                                      | 28044                                                                                                 | `awk '$2 == "TRA" && $6 == "HomoSapiens"  { print $3 }' vdjdb-2019-08-08/vdjdb.txt | wc -l`                                                                                                                     |
|   TRB records                                                      | 39722                                                                                                | `awk '$2 == "TRB" && $6 == "HomoSapiens"  { print $3 }' vdjdb-2019-08-08/vdjdb.txt | wc -l`                                                                                                                     |
|   Unique TRA sequences                                             | 20277                                                                                      | `awk '$2 == "TRA" && $6 == "HomoSapiens"  { print $3 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                                                           |
|   Unique TRB sequences                                             | 30720                                                                                      | `awk '$2 == "TRB" && $6 == "HomoSapiens"  { print $3 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                                                           |
|   Unique CDR3 sequences                                            | 50995                                                                                               | `awk '$6 == "HomoSapiens" { print }' vdjdb-2019-08-08/vdjdb.txt | cut -f3 | sort -u | wc -l`                                                                                                                      |
|   Unique epitope sequences                                         | 177                                                                                              | `awk '$6 == "HomoSapiens" { print }' vdjdb-2019-08-08/vdjdb.txt | cut -f10 | sort -u | wc -l`                                                                                                                     |
|   Unique epitope sequences for TRA records                         | 117                                                                                      | `awk '$2 == "TRA" && $6 == "HomoSapiens" { print $10 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                                                           |
|   Unique epitope sequences for TRB records                         | 175                                                                                      | `awk '$2 == "TRB" && $6 == "HomoSapiens" { print $10 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                                                           |
|   Unique CDR3-epitope sequence pairs                               | 55931                                                                                   | `awk '$6 == "HomoSapiens" { print }' vdjdb-2019-08-08/vdjdb.txt | cut -d $'\t' -f3,10 | sort -u | wc -l`                                                                                                          |
|   Unique TRA-CDR3-epitope sequence pairs                           | 22945                                                                                   | `awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3,$10 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                                                       |
|   Unique TRB-CDR3-epitope sequence pairs                           | 32988                                                                                   | `awk '$2 == "TRB" && $6 == "HomoSapiens" { print $3,$10 }' vdjdb-2019-08-08/vdjdb.txt | sort -u | wc -l`                                                                                                       |
|   Number of epitope sequences shared between TRA and TRB records   | 115|`comm -12 <(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $10 }' vdjdb-2019-08-08/vdjdb.txt | sort -u) <(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $10 }' vdjdb-2019-08-08/vdjdb.txt | sort -u) | wc -l`   |
|   Number of CDR3 sequences shared between TRA and TRB records      | 2|`comm -12 <(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3 }' vdjdb-2019-08-08/vdjdb.txt | sort -u) <(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $3 }' vdjdb-2019-08-08/vdjdb.txt | sort -u) | wc -l`     |
|   Epitope distribution for the unique CDR3-epitope pairs           |                             | `awk '$6 == "HomoSapiens" { print }' vdjdb-2019-08-08/vdjdb.txt | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20`                                                                 |

    24367 KLGGALQAK
    6680 NLVPMVATV
    6531 GILGFVFTL
    3222 AVFDRKSDAK
    1642 ELAGIGILTV
    1492 RAKFKQLL
    1165 GLCTLVAML
    1083 IVTDFSVIK
    810 RLRAEAQVK
    676 LLWNGPMAV
    588 LLLGIGILV
    562 FRDYVDRFYKTLRAEQASQE
    532 PKYVKQNTLKLAT
    334 KRWIILGLNK
    240 GLIYNRMGAVTTEV
    227 QARQMVQAMRTIGTHP
    209 KAFSPEVIPMF
    207 TPRVTGGGAM
    203 VTEHDTLLY
    203 CINGVCWTV

## vdjdb.slim

| Metric                                                         | Count                                                                                                                                    | Command                                                                                                                       |
|----------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| Total number of records                                        | 61049                                                                                                        | `tail -n +2 vdjdb-2019-08-08/vdjdb.slim.txt  |  wc -l`                                                                                                 |
| TRA records                                                    | 25051                                                                                    | `awk '$1 == "TRA" { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | wc -l`                                                                             |
| TRB records                                                    | 35998                                                                                    | `awk '$1 == "TRB" { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | wc -l`                                                                             |
| Unique TRA sequences                                           | 22248                                                                          | `awk '$1 == "TRA" { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                   |
| Unique TRB sequences                                           | 33606                                                                          | `awk '$1 == "TRB" { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                   |
| Unique CDR3 sequences                                          | 55852                                                                                    | `tail -n +2 vdjdb-2019-08-08/vdjdb.slim.txt | cut -f2 | sort -u | wc -l`                                                                               |
| Unique epitope sequences                                       | 212                                                                                    | `tail -n +2 vdjdb-2019-08-08/vdjdb.slim.txt | cut -f4 | sort -u | wc -l`                                                                              |
| Unique epitope sequences for TRA records                       | 151                                                                          | `awk '$1 == "TRA" { print $4 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                  |
| Unique epitope sequences for TRB records                       | 212                                                                          | `awk '$1 == "TRB" { print $4 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                  |
| Unique CDR3-epitope sequence pairs                             | 61047                                                                         | `tail -n +2 vdjdb-2019-08-08/vdjdb.slim.txt | cut -d $'\t' -f2,4 | sort -u | wc -l`                                                                   |
| Unique TRA-CDR3-epitope sequence pairs                         | 25051                                                                       | `awk '$1 == "TRA" { print vdjdb-2019-08-08/vdjdb.slim.txt, }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                |
| Unique TRB-CDR3-epitope sequence pairs                         | 35998                                                                       | `awk '$1 == "TRB" { print vdjdb-2019-08-08/vdjdb.slim.txt, }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                |
| Number of epitope sequences shared between TRA and TRB records | 151      | `comm -12 <(awk '$1 == "TRA" { print $4 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u) <(awk '$1 == "TRB" { print $4 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u) | wc -l`  |
| Number of CDR3 sequences shared between TRA and TRB records    | 2      | `comm -12 <(awk '$1 == "TRA" { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u) <(awk '$1 == "TRB" { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u) | wc -l`    |
| Epitope distribution for the unique CDR3-epitope pairs         |                                 | `tail -n +2 vdjdb-2019-08-08/vdjdb.slim.txt | cut -d $'\t' -f2,4 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20`                                  |

  24639 KLGGALQAK
   6744 NLVPMVATV
   6600 GILGFVFTL
   3255 AVFDRKSDAK
   1654 ELAGIGILTV
   1501 RAKFKQLL
   1190 GLCTLVAML
   1088 IVTDFSVIK
    819 RLRAEAQVK
    685 LLWNGPMAV
    617 LLLGIGILV
    597 TTPESANL
    597 SSLENFRAYV
    586 SSYRRPVGI
    564 FRDYVDRFYKTLRAEQASQE
    548 PKYVKQNTLKLAT
    530 CTPYDINQM
    417 HGIRNASFI
    383 ASNENMETM
    334 KRWIILGLNK

**Human-only records**

| Metric                                                             | Count                                                                                                                                                                        | Command                                                                                                                                                                             |
|--------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|   Total number of records                                          | 56592                                                                                                                  | `awk '$3 == "HomoSapiens" { print }' vdjdb-2019-08-08/vdjdb.slim.txt |  wc -l`                                                                                                                              |
|   TRA records                                                      | 23564                                                                                                 | `awk '$1 == "TRA" && $3 == "HomoSapiens"  { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | wc -l`                                                                                                          |
|   TRB records                                                      | 33028                                                                                                | `awk '$1 == "TRB" && $3 == "HomoSapiens"  { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | wc -l`                                                                                                          |
|   Unique TRA sequences                                             | 20841                                                                                      | `awk '$1 == "TRA" && $3 == "HomoSapiens"  { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                                                |
|   Unique TRB sequences                                             | 30757                                                                                      | `awk '$1 == "TRB" && $3 == "HomoSapiens"  { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                                                |
|   Unique CDR3 sequences                                            | 51596                                                                                               | `awk '$3 == "HomoSapiens" { print }' vdjdb-2019-08-08/vdjdb.slim.txt | cut -f2 | sort -u | wc -l`                                                                                                           |
|   Unique epitope sequences                                         | 177                                                                                              | `awk '$3 == "HomoSapiens" { print }' vdjdb-2019-08-08/vdjdb.slim.txt | cut -f4 | sort -u | wc -l`                                                                                                          |
|   Unique epitope sequences for TRA records                         | 125                                                                                      | `awk '$1 == "TRA" && $3 == "HomoSapiens" { print $4 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                                                |
|   Unique epitope sequences for TRB records                         | 177                                                                                      | `awk '$1 == "TRB" && $3 == "HomoSapiens" { print $4 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                                                |
|   Unique CDR3-epitope sequence pairs                               | 56590                                                                                   | `awk '$3 == "HomoSapiens" { print }' vdjdb-2019-08-08/vdjdb.slim.txt | cut -d $'\t' -f2,4 | sort -u | wc -l`                                                                                               |
|   Unique TRA-CDR3-epitope sequence pairs                           | 23564                                                                                   | `awk '$1 == "TRA" && $3 == "HomoSapiens" { print $2,$4 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                                            |
|   Unique TRB-CDR3-epitope sequence pairs                           | 33028                                                                                   | `awk '$1 == "TRB" && $3 == "HomoSapiens" { print $2,$4 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u | wc -l`                                                                                            |
|   Number of epitope sequences shared between TRA and TRB records   | 125|`comm -12 <(awk '$1 == "TRA" && $3 == "HomoSapiens" { print $4 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u) <(awk '$1 == "TRB" && $3 == "HomoSapiens" { print $4 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u) | wc -l` |
|   Number of CDR3 sequences shared between TRA and TRB records      | 2|`comm -12 <(awk '$1 == "TRA" && $3 == "HomoSapiens" { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u) <(awk '$1 == "TRB" && $3 == "HomoSapiens" { print $2 }' vdjdb-2019-08-08/vdjdb.slim.txt | sort -u) | wc -l`   |
|   Epitope distribution for the unique CDR3-epitope pairs           |                             | `awk '$6 == "HomoSapiens" { print }' vdjdb-2019-08-08/vdjdb.slim.txt | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20`                                                                 |

    24639 KLGGALQAK
    6744 NLVPMVATV
    6600 GILGFVFTL
    3255 AVFDRKSDAK
    1654 ELAGIGILTV
    1501 RAKFKQLL
    1190 GLCTLVAML
    1088 IVTDFSVIK
    819 RLRAEAQVK
    685 LLWNGPMAV
    617 LLLGIGILV
    564 FRDYVDRFYKTLRAEQASQE
    548 PKYVKQNTLKLAT
    334 KRWIILGLNK
    244 GLIYNRMGAVTTEV
    229 QARQMVQAMRTIGTHP
    217 TPRVTGGGAM
    209 KAFSPEVIPMF
    206 CINGVCWTV
    203 VTEHDTLLY

**Human-only records**

| Metric                                                             | Count                                                                                                                                                                          | Command                                                                                                                                                                                                 |
|--------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|   Total number of records                                          | 66597                                                                                                                  | `awk '$6 == "HomoSapiens" { print }' vdjdb-browser.tsv |  wc -l`                                                                                                                                         |
|   TRA records                                                      | 27295                                                                                                 | `awk '$2 == "TRA" && $6 == "HomoSapiens"  { print $3 }' vdjdb-browser.tsv | wc -l`                                                                                                                     |
|   TRB records                                                      | 39302                                                                                                | `awk '$2 == "TRB" && $6 == "HomoSapiens"  { print $3 }' vdjdb-browser.tsv | wc -l`                                                                                                                     |
|   Unique TRA sequences                                             | 19637                                                                                      | `awk '$2 == "TRA" && $6 == "HomoSapiens"  { print $3 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                                                           |
|   Unique TRB sequences                                             | 30377                                                                                      | `awk '$2 == "TRB" && $6 == "HomoSapiens"  { print $3 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                                                           |
|   Unique CDR3 sequences                                            | 50014                                                                                               | `awk '$6 == "HomoSapiens" { print }' vdjdb-browser.tsv | cut -f3 | sort -u | wc -l`                                                                                                                      |
|   Unique epitope sequences                                         | 177                                                                                              | `awk '$6 == "HomoSapiens" { print }' vdjdb-browser.tsv | cut -f10 | sort -u | wc -l`                                                                                                                     |
|   Unique epitope sequences for TRA records                         | 117                                                                                      | `awk '$2 == "TRA" && $6 == "HomoSapiens" { print $10 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                                                           |
|   Unique epitope sequences for TRB records                         | 175                                                                                      | `awk '$2 == "TRB" && $6 == "HomoSapiens" { print $10 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                                                           |
|   Unique CDR3-epitope sequence pairs                               | 54887                                                                                   | `awk '$6 == "HomoSapiens" { print }' vdjdb-browser.tsv | cut -d $'\t' -f3,10 | sort -u | wc -l`                                                                                                          |
|   Unique TRA-CDR3-epitope sequence pairs                           | 22248                                                                                   | `awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3,$10 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                                                       |
|   Unique TRB-CDR3-epitope sequence pairs                           | 32639                                                                                   | `awk '$2 == "TRB" && $6 == "HomoSapiens" { print $3,$10 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                                                       |
|   Number of epitope sequences shared between TRA and TRB records   | 115|`comm -12 <(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $10 }' vdjdb-browser.tsv | sort -u) <(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $10 }' vdjdb-browser.tsv | sort -u) | wc -l`   |
|   Number of CDR3 sequences shared between TRA and TRB records      | 0|`comm -12 <(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3 }' vdjdb-browser.tsv | sort -u) <(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $3 }' vdjdb-browser.tsv | sort -u) | wc -l`     |
|   Epitope distribution for the unique CDR3-epitope pairs           |                             | `awk '$6 == "HomoSapiens" { print }' vdjdb-browser.tsv | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20`                                                                 |

  24041 KLGGALQAK
   6514 NLVPMVATV
   6389 GILGFVFTL
   3175 AVFDRKSDAK
   1624 ELAGIGILTV
   1472 RAKFKQLL
   1142 GLCTLVAML
   1067 IVTDFSVIK
    800 RLRAEAQVK
    652 LLWNGPMAV
    577 LLLGIGILV
    523 PKYVKQNTLKLAT
    508 FRDYVDRFYKTLRAEQASQE
    325 KRWIILGLNK
    240 GLIYNRMGAVTTEV
    226 QARQMVQAMRTIGTHP
    201 VTEHDTLLY
    201 CINGVCWTV
    199 TPRVTGGGAM
    193 KAFSPEVIPMF
