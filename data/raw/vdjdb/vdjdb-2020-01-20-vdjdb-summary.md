
| Metric                                                         | Count                                                                                                                                       | Command                                                                                                                                        |
|----------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| Total number of records                                        | 76164                                                                                                         | `tail -n +2 data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt  |  wc -l`                                                                                                         |
| TRA records                                                    | 31244                                                                                     | `awk '$2 == "TRA" { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | wc -l`                                                                                     |
| TRB records                                                    | 44920                                                                                     | `awk '$2 == "TRB" { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | wc -l`                                                                                     |
| Unique TRA sequences                                           | 22352                                                                           | `awk '$2 == "TRA" { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                           |
| Unique TRB sequences                                           | 33977                                                                           | `awk '$2 == "TRB" { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                           |
| Unique CDR3 sequences                                          | 56327                                                                                     | `tail -n +2 data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | cut -f3 | sort -u | wc -l`                                                                                       |
| Unique epitope sequences                                       | 222                                                                                    | `tail -n +2 data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | cut -f10 | sort -u | wc -l`                                                                                      |
| Unique epitope sequences for TRA records                       | 188                                                                          | `awk '$2 == "TRA" { print $10 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                          |
| Unique epitope sequences for TRB records                       | 227                                                                          | `awk '$2 == "TRB" { print $10 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                          |
| Unique CDR3-epitope sequence pairs                             | 61555                                                                         | `tail -n +2 data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | cut -d $'\t' -f3,10 | sort -u | wc -l`                                                                           |
| Unique TRA-CDR3-epitope sequence pairs                         | 25166                                                                       | `awk '$2 == "TRA" { print ,data/raw/vdjdb/vdjdb-2020-01-20-summary.md0 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                        |
| Unique TRB-CDR3-epitope sequence pairs                         | 36386                                                                       | `awk '$2 == "TRB" { print ,data/raw/vdjdb/vdjdb-2020-01-20-summary.md0 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                        |
| Number of epitope sequences shared between TRA and TRB records | 154   | `comm -12 <(awk '$2 == "TRA" { print $10 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u) <(awk '$2 == "TRB" { print $10 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u) | wc -l` |
| Number of CDR3 sequences shared between TRA and TRB records    | 2     | `comm -12 <(awk '$2 == "TRA" { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u) <(awk '$2 == "TRB" { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u) | wc -l`   |
| Epitope distribution for the unique CDR3-epitope pairs         |                                 | `tail -n +2 data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20`                                  |

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
|   Total number of records                                          | 68067                                                                                                                  | `awk '$6 == "HomoSapiens" { print }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt |  wc -l`                                                                                                                                         |
|   TRA records                                                      | 28183                                                                                                 | `awk '$2 == "TRA" && $6 == "HomoSapiens"  { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | wc -l`                                                                                                                     |
|   TRB records                                                      | 39884                                                                                                | `awk '$2 == "TRB" && $6 == "HomoSapiens"  { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | wc -l`                                                                                                                     |
|   Unique TRA sequences                                             | 20373                                                                                      | `awk '$2 == "TRA" && $6 == "HomoSapiens"  { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                                                           |
|   Unique TRB sequences                                             | 30852                                                                                      | `awk '$2 == "TRB" && $6 == "HomoSapiens"  { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                                                           |
|   Unique CDR3 sequences                                            | 51223                                                                                               | `awk '$6 == "HomoSapiens" { print }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | cut -f3 | sort -u | wc -l`                                                                                                                      |
|   Unique epitope sequences                                         | 186                                                                                              | `awk '$6 == "HomoSapiens" { print }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | cut -f10 | sort -u | wc -l`                                                                                                                     |
|   Unique epitope sequences for TRA records                         | 126                                                                                      | `awk '$2 == "TRA" && $6 == "HomoSapiens" { print $10 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                                                           |
|   Unique epitope sequences for TRB records                         | 184                                                                                      | `awk '$2 == "TRB" && $6 == "HomoSapiens" { print $10 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                                                           |
|   Unique CDR3-epitope sequence pairs                               | 56181                                                                                   | `awk '$6 == "HomoSapiens" { print }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | cut -d $'\t' -f3,10 | sort -u | wc -l`                                                                                                          |
|   Unique TRA-CDR3-epitope sequence pairs                           | 23059                                                                                   | `awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3,$10 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                                                       |
|   Unique TRB-CDR3-epitope sequence pairs                           | 33124                                                                                   | `awk '$2 == "TRB" && $6 == "HomoSapiens" { print $3,$10 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u | wc -l`                                                                                                       |
|   Number of epitope sequences shared between TRA and TRB records   | 124|`comm -12 <(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $10 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u) <(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $10 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u) | wc -l`   |
|   Number of CDR3 sequences shared between TRA and TRB records      | 2|`comm -12 <(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u) <(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $3 }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | sort -u) | wc -l`     |
|   Epitope distribution for the unique CDR3-epitope pairs           |                             | `awk '$6 == "HomoSapiens" { print }' data/raw/vdjdb/vdjdb-2020-01-20/vdjdb.txt | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20`                                                                 |

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
    534 PKYVKQNTLKLAT
    334 KRWIILGLNK
    240 GLIYNRMGAVTTEV
    227 QARQMVQAMRTIGTHP
    209 KAFSPEVIPMF
    207 TPRVTGGGAM
    203 VTEHDTLLY
    203 CINGVCWTV

