
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

