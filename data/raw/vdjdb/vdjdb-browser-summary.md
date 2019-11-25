
| Metric                                                         | Count                                                                                                                                       | Command                                                                                                                                        |
|----------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| Total number of records                                        | 73224                                                                                                         | `tail -n +2 vdjdb-browser.tsv  |  wc -l`                                                                                                         |
| TRA records                                                    | 29479                                                                                     | `awk '$2 == "TRA" { print $3 }' vdjdb-browser.tsv | wc -l`                                                                                     |
| TRB records                                                    | 43745                                                                                     | `awk '$2 == "TRB" { print $3 }' vdjdb-browser.tsv | wc -l`                                                                                     |
| Unique TRA sequences                                           | 20888                                                                           | `awk '$2 == "TRA" { print $3 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                           |
| Unique TRB sequences                                           | 33112                                                                           | `awk '$2 == "TRB" { print $3 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                           |
| Unique CDR3 sequences                                          | 54000                                                                                     | `tail -n +2 vdjdb-browser.tsv | cut -f3 | sort -u | wc -l`                                                                                       |
| Unique epitope sequences                                       | 212                                                                                    | `tail -n +2 vdjdb-browser.tsv | cut -f10 | sort -u | wc -l`                                                                                      |
| Unique epitope sequences for TRA records                       | 141                                                                          | `awk '$2 == "TRA" { print $10 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                          |
| Unique epitope sequences for TRB records                       | 210                                                                          | `awk '$2 == "TRB" { print $10 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                          |
| Unique CDR3-epitope sequence pairs                             | 59072                                                                         | `tail -n +2 vdjdb-browser.tsv | cut -d $'\t' -f3,10 | sort -u | wc -l`                                                                           |
| Unique TRA-CDR3-epitope sequence pairs                         | 23578                                                                       | `awk '$2 == "TRA" { print ,vdjdb-browser-summary.md0 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                        |
| Unique TRB-CDR3-epitope sequence pairs                         | 35494                                                                       | `awk '$2 == "TRB" { print ,vdjdb-browser-summary.md0 }' vdjdb-browser.tsv | sort -u | wc -l`                                                                        |
| Number of epitope sequences shared between TRA and TRB records | 139   | `comm -12 <(awk '$2 == "TRA" { print $10 }' vdjdb-browser.tsv | sort -u) <(awk '$2 == "TRB" { print $10 }' vdjdb-browser.tsv | sort -u) | wc -l` |
| Number of CDR3 sequences shared between TRA and TRB records    | 0     | `comm -12 <(awk '$2 == "TRA" { print $3 }' vdjdb-browser.tsv | sort -u) <(awk '$2 == "TRB" { print $3 }' vdjdb-browser.tsv | sort -u) | wc -l`   |
| Epitope distribution for the unique CDR3-epitope pairs         |                                 | `tail -n +2 vdjdb-browser.tsv | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20`                                  |

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
    585 SSLENFRAYV
    577 LLLGIGILV
    562 SSYRRPVGI
    523 PKYVKQNTLKLAT
    511 TTPESANL
    508 FRDYVDRFYKTLRAEQASQE
    500 CTPYDINQM
    398 HGIRNASFI
    345 ASNENMETM
    325 KRWIILGLNK

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

