
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

