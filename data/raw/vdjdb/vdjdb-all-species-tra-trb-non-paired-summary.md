
|Metric|Count|Command|
|---|---|---|
|Total number of records| 73224|`tail -n +2 vdjdb-all-species-tra-trb-non-paired.tsv |  wc -l`|
|TRA records|29479|`awk ' == "TRA" { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | wc -l`|
|TRB records|43745|`awk ' == "TRB" { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | wc -l`|
|Unique TRA sequences|20888|`awk ' == "TRA" { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique TRB sequences|33112|`awk ' == "TRB" { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique CDR3 sequences|54000|`tail -n +2 vdjdb-all-species-tra-trb-non-paired.tsv | cut -f3 | sort -u | wc -l`|
|Unique epitope sequences|212|`tail -n +2 vdjdb-all-species-tra-trb-non-paired.tsv | cut -f10 | sort -u | wc -l`|
|Unique epitope sequences for TRA records|141|`awk ' == "TRA" { print 0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique epitope sequences for TRB records|210|`awk ' == "TRB" { print 0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique CDR3-epitope sequence pairs|59072|`tail -n +2 vdjdb-all-species-tra-trb-non-paired.tsv | cut -d $'\t' -f3,10 | sort -u | wc -l`|
|Unique TRA-CDR3-epitope sequence pairs|23578|`awk ' == "TRA" { print ,0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique TRB-CDR3-epitope sequence pairs|35494|`awk ' == "TRB" { print ,0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Number of epitope sequences shared between TRA and TRB records|139|`comm -12 <(awk ' == "TRA" { print 0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) <(awk ' == "TRB" { print 0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) | wc -l`|
|Number of CDR3 sequences shared between TRA and TRB records|0|`comm -12 <(awk ' == "TRA" { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) <(awk ' == "TRB" { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) | wc -l`|

**Human-only records**

|Metric|Count|Command|
|---|---|---|
|Total number of records|66597|`awk ' == "HomoSapiens" { print }' vdjdb-human-tra-trb-non-paired.tsv |  wc -l`|
|TRA records|27295|`awk ' == "TRA" &&  == "HomoSapiens" { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | wc -l`|
|TRB records|39302|`awk ' == "TRB" &&  == "HomoSapiens"  { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | wc -l`|
|Unique TRA sequences|19637|`awk ' == "TRA" &&  == "HomoSapiens"  { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique TRB sequences|30377|`awk ' == "TRB" &&  == "HomoSapiens"  { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique CDR3 sequences|50014|`awk ' == "HomoSapiens" { print }' vdjdb-all-species-tra-trb-non-paired.tsv | cut -f3 | sort -u | wc -l`|
|Unique epitope sequences|177|`awk ' == "HomoSapiens" { print }' vdjdb-all-species-tra-trb-non-paired.tsv | cut -f10 | sort -u | wc -l`|
|Unique epitope sequences for TRA records|117|`awk ' == "TRA" &&  == "HomoSapiens" { print 0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique epitope sequences for TRB records|175|`awk ' == "TRB" &&  == "HomoSapiens" { print 0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique CDR3-epitope sequence pairs|54887|`awk ' == "HomoSapiens" { print }' vdjdb-all-species-tra-trb-non-paired.tsv | cut -d $'\t' -f3,10 | sort -u | wc -l`|
|Unique TRA-CDR3-epitope sequence pairs|22248|`awk ' == "TRA" &&  == "HomoSapiens" { print ,0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Unique TRB-CDR3-epitope sequence pairs|32639|`awk ' == "TRB" &&  == "HomoSapiens" { print ,0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u | wc -l`|
|Number of epitope sequences shared between TRA and TRB records|115|`comm -12 <(awk ' == "TRA" &&  == "HomoSapiens" { print 0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) <(awk ' == "TRB" &&  == "HomoSapiens" { print 0 }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) | wc -l`|
|Number of CDR3 sequences shared between TRA and TRB records|0|`comm -12 <(awk ' == "TRA" &&  == "HomoSapiens" { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) <(awk ' == "TRB" &&  == "HomoSapiens" { print  }' vdjdb-all-species-tra-trb-non-paired.tsv | sort -u) | wc -l`|

