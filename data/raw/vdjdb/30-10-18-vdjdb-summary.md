
| Metric                                                         | Count                                                                                                                                    | Command                                                                                                                       |
|----------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| Total number of records                                        | 25420                                                                                                        | `tail -n +2 30-10-18-vdjdb.csv  |  wc -l`                                                                                                 |
| TRA records                                                    | 7487                                                                                    | `awk '$1 == "TRA" { print $2 }' 30-10-18-vdjdb.csv | wc -l`                                                                             |
| TRB records                                                    | 17933                                                                                    | `awk '$1 == "TRB" { print $2 }' 30-10-18-vdjdb.csv | wc -l`                                                                             |
| Unique TRA sequences                                           | 7280                                                                          | `awk '$1 == "TRA" { print $2 }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                   |
| Unique TRB sequences                                           | 17493                                                                          | `awk '$1 == "TRB" { print $2 }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                   |
| Unique CDR3 sequences                                          | 24771                                                                                    | `tail -n +2 30-10-18-vdjdb.csv | cut -f2 | sort -u | wc -l`                                                                               |
| Unique epitope sequences                                       | 182                                                                                    | `tail -n +2 30-10-18-vdjdb.csv | cut -f4 | sort -u | wc -l`                                                                              |
| Unique epitope sequences for TRA records                       | 117                                                                          | `awk '$1 == "TRA" { print $4 }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                  |
| Unique epitope sequences for TRB records                       | 182                                                                          | `awk '$1 == "TRB" { print $4 }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                  |
| Unique CDR3-epitope sequence pairs                             | 25386                                                                         | `tail -n +2 30-10-18-vdjdb.csv | cut -d $'\t' -f2,4 | sort -u | wc -l`                                                                   |
| Unique TRA-CDR3-epitope sequence pairs                         | 7485                                                                       | `awk '$1 == "TRA" { print 30-10-18-vdjdb.csv, }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                |
| Unique TRB-CDR3-epitope sequence pairs                         | 17903                                                                       | `awk '$1 == "TRB" { print 30-10-18-vdjdb.csv, }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                |
| Number of epitope sequences shared between TRA and TRB records | 117      | `comm -12 <(awk '$1 == "TRA" { print $4 }' 30-10-18-vdjdb.csv | sort -u) <(awk '$1 == "TRB" { print $4 }' 30-10-18-vdjdb.csv | sort -u) | wc -l`  |
| Number of CDR3 sequences shared between TRA and TRB records    | 2      | `comm -12 <(awk '$1 == "TRA" { print $2 }' 30-10-18-vdjdb.csv | sort -u) <(awk '$1 == "TRB" { print $2 }' 30-10-18-vdjdb.csv | sort -u) | wc -l`    |
| Epitope distribution for the unique CDR3-epitope pairs         |                                 | `tail -n +2 30-10-18-vdjdb.csv | cut -d $'\t' -f2,4 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20`                                  |

   6689 NLVPMVATV
   4838 GILGFVFTL
   1030 ELAGIGILTV
   1005 GLCTLVAML
    685 LLWNGPMAV
    617 LLLGIGILV
    598 TTPESANL
    597 SSLENFRAYV
    586 SSYRRPVGI
    564 FRDYVDRFYKTLRAEQASQE
    530 CTPYDINQM
    417 HGIRNASFI
    383 ASNENMETM
    363 PKYVKQNTLKLAT
    334 KRWIILGLNK
    264 LSLRNPILV
    237 STPESANL
    215 TPRVTGGGAM
    209 KAFSPEVIPMF
    206 CINGVCWTV

**Human-only records**

| Metric                                                             | Count                                                                                                                                                                        | Command                                                                                                                                                                             |
|--------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|   Total number of records                                          | 20941                                                                                                                  | `awk '$3 == "HomoSapiens" { print }' 30-10-18-vdjdb.csv |  wc -l`                                                                                                                              |
|   TRA records                                                      | 6000                                                                                                 | `awk '$1 == "TRA" && $3 == "HomoSapiens"  { print $2 }' 30-10-18-vdjdb.csv | wc -l`                                                                                                          |
|   TRB records                                                      | 14941                                                                                                | `awk '$1 == "TRB" && $3 == "HomoSapiens"  { print $2 }' 30-10-18-vdjdb.csv | wc -l`                                                                                                          |
|   Unique TRA sequences                                             | 5858                                                                                      | `awk '$1 == "TRA" && $3 == "HomoSapiens"  { print $2 }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                                                |
|   Unique TRB sequences                                             | 14636                                                                                      | `awk '$1 == "TRB" && $3 == "HomoSapiens"  { print $2 }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                                                |
|   Unique CDR3 sequences                                            | 20492                                                                                               | `awk '$3 == "HomoSapiens" { print }' 30-10-18-vdjdb.csv | cut -f2 | sort -u | wc -l`                                                                                                           |
|   Unique epitope sequences                                         | 147                                                                                              | `awk '$3 == "HomoSapiens" { print }' 30-10-18-vdjdb.csv | cut -f4 | sort -u | wc -l`                                                                                                          |
|   Unique epitope sequences for TRA records                         | 91                                                                                      | `awk '$1 == "TRA" && $3 == "HomoSapiens" { print $4 }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                                                |
|   Unique epitope sequences for TRB records                         | 147                                                                                      | `awk '$1 == "TRB" && $3 == "HomoSapiens" { print $4 }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                                                |
|   Unique CDR3-epitope sequence pairs                               | 20928                                                                                   | `awk '$3 == "HomoSapiens" { print }' 30-10-18-vdjdb.csv | cut -d $'\t' -f2,4 | sort -u | wc -l`                                                                                               |
|   Unique TRA-CDR3-epitope sequence pairs                           | 5998                                                                                   | `awk '$1 == "TRA" && $3 == "HomoSapiens" { print $2,$4 }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                                            |
|   Unique TRB-CDR3-epitope sequence pairs                           | 14932                                                                                   | `awk '$1 == "TRB" && $3 == "HomoSapiens" { print $2,$4 }' 30-10-18-vdjdb.csv | sort -u | wc -l`                                                                                            |
|   Number of epitope sequences shared between TRA and TRB records   | 91|`comm -12 <(awk '$1 == "TRA" && $3 == "HomoSapiens" { print $4 }' 30-10-18-vdjdb.csv | sort -u) <(awk '$1 == "TRB" && $3 == "HomoSapiens" { print $4 }' 30-10-18-vdjdb.csv | sort -u) | wc -l` |
|   Number of CDR3 sequences shared between TRA and TRB records      | 2|`comm -12 <(awk '$1 == "TRA" && $3 == "HomoSapiens" { print $2 }' 30-10-18-vdjdb.csv | sort -u) <(awk '$1 == "TRB" && $3 == "HomoSapiens" { print $2 }' 30-10-18-vdjdb.csv | sort -u) | wc -l`   |
|   Epitope distribution for the unique CDR3-epitope pairs           |                             | `awk '$6 == "HomoSapiens" { print }' 30-10-18-vdjdb.csv | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20`                                                                 |

   6689 NLVPMVATV
   4838 GILGFVFTL
   1030 ELAGIGILTV
   1005 GLCTLVAML
    685 LLWNGPMAV
    617 LLLGIGILV
    564 FRDYVDRFYKTLRAEQASQE
    363 PKYVKQNTLKLAT
    334 KRWIILGLNK
    215 TPRVTGGGAM
    209 KAFSPEVIPMF
    206 CINGVCWTV
    201 VTEHDTLLY
    179 RAKFKQLL
    173 GTSGSPIVNR
    166 GTSGSPIINR
    165 ATDALMTGY
    159 LPRRSGAAGA
    158 FLKEKGGL
    153 EIYKRWII

