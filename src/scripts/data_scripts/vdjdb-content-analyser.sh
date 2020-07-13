#! /usr/bin/env bash

# retrieve project root (or pass project root from makefile instead?)
# ${parameter:-word} If parameter is unset or null, the expansion of word is substituted. Otherwise, the value of parameter is substituted.
# i.e. if PROJECT_ROOT is not passed through via the Makefile (or when the script is run on its own)
# the oneliner will attempt to move up two directories from the location of this script and set that as the project root.
# A potential problem with this approach is that when PROJECT_ROOT is not a resolved path, VDJDB_RAW_DATA_PATh will not be valid.
# However, normally this script will either be used on its own (and PROJECT_ROOT will be empty initially)
# or it will be invoked by the Makefile, which should pass a fully resolved path.
# To be safe, the PROJECT_ROOT variable is assigned a second time as resolved path.
# Note that "readlink -f/--canonicalize" or "realpath" are not used because it does not exist on all systems.
PROJECT_ROOT=${PROJECT_ROOT:-"$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../.." >/dev/null 2>&1 && pwd )"}

# set default input and output paths
OUTPUT_FILE=${1:-${PROJECT_ROOT}/data/raw/vdjdb/vdjdb-browser-summary.md}
VDJDB_FILE=${2:-${PROJECT_ROOT}/data/raw/vdjdb/vdjdb-browser.tsv}

# check if output file already exists and confirm before overwriting
if [ -f "$OUTPUT_FILE" ]; then
    read -r -p "The file \"$OUTPUT_FILE\" already exists; overwrite it? [y/N] " response
    if [[ ! "$response" =~ ^([yY][eE][sS]|[yY])+$ ]]
    then
        echo "Aborting without writing output..."
        exit 1
    fi
fi

# check if input file exists
if [ ! -f "${VDJDB_FILE}" ]; then
        echo "The input file \"${VDJDB_FILE}\" does not exist!"
        exit 1
fi

cat <<EOF > $OUTPUT_FILE

| Metric                                                         | Count                                                                                                                                       | Command                                                                                                                                        |
|----------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| Total number of records                                        | $(tail -n +2 ${VDJDB_FILE} | wc -l)                                                                                                         | \`tail -n +2 ${VDJDB_FILE}  |  wc -l\`                                                                                                         |
| TRA records                                                    | $(awk '$2 == "TRA" { print $3 }' ${VDJDB_FILE} | wc -l)                                                                                     | \`awk '\$2 == "TRA" { print \$3 }' ${VDJDB_FILE} | wc -l\`                                                                                     |
| TRB records                                                    | $(awk '$2 == "TRB" { print $3 }' ${VDJDB_FILE} | wc -l)                                                                                     | \`awk '\$2 == "TRB" { print \$3 }' ${VDJDB_FILE} | wc -l\`                                                                                     |
| Unique TRA sequences                                           | $(awk '$2 == "TRA" { print $3 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                           | \`awk '\$2 == "TRA" { print \$3 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                           |
| Unique TRB sequences                                           | $(awk '$2 == "TRB" { print $3 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                           | \`awk '\$2 == "TRB" { print \$3 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                           |
| Unique CDR3 sequences                                          | $(tail -n +2 ${VDJDB_FILE} | cut -f3 | sort -u | wc -l)                                                                                     | \`tail -n +2 ${VDJDB_FILE} | cut -f3 | sort -u | wc -l\`                                                                                       |
| Unique epitope sequences                                       | $(tail -n +2 ${VDJDB_FILE} | cut -f10 | sort -u | wc -l)                                                                                    | \`tail -n +2 ${VDJDB_FILE} | cut -f10 | sort -u | wc -l\`                                                                                      |
| Unique epitope sequences for TRA records                       | $(awk '$2 == "TRA" { print $10 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                          | \`awk '\$2 == "TRA" { print \$10 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                          |
| Unique epitope sequences for TRB records                       | $(awk '$2 == "TRB" { print $10 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                          | \`awk '\$2 == "TRB" { print \$10 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                          |
| Unique CDR3-epitope sequence pairs                             | $(tail -n +2 ${VDJDB_FILE} | cut -d $'\t' -f3,10 | sort -u | wc -l)                                                                         | \`tail -n +2 ${VDJDB_FILE} | cut -d $'\t' -f3,10 | sort -u | wc -l\`                                                                           |
| Unique TRA-CDR3-epitope sequence pairs                         | $(awk '$2 == "TRA" { print $3,$10 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                       | \`awk '\$2 == "TRA" { print $3,$10 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                        |
| Unique TRB-CDR3-epitope sequence pairs                         | $(awk '$2 == "TRB" { print $3,$10 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                       | \`awk '\$2 == "TRB" { print $3,$10 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                        |
| Number of epitope sequences shared between TRA and TRB records | $(comm -12 <(awk '$2 == "TRA" { print $10 }' ${VDJDB_FILE} | sort -u) <(awk '$2 == "TRB" { print $10 }' ${VDJDB_FILE} | sort -u) | wc -l)   | \`comm -12 <(awk '\$2 == "TRA" { print \$10 }' ${VDJDB_FILE} | sort -u) <(awk '\$2 == "TRB" { print \$10 }' ${VDJDB_FILE} | sort -u) | wc -l\` |
| Number of CDR3 sequences shared between TRA and TRB records    | $(comm -12 <(awk '$2 == "TRA" { print $3 }' ${VDJDB_FILE} | sort -u) <(awk '$2 == "TRB" { print $3 }' ${VDJDB_FILE} | sort -u) | wc -l)     | \`comm -12 <(awk '\$2 == "TRA" { print \$3 }' ${VDJDB_FILE} | sort -u) <(awk '\$2 == "TRB" { print \$3 }' ${VDJDB_FILE} | sort -u) | wc -l\`   |
| Epitope distribution for the unique CDR3-epitope pairs         |                                 | \`tail -n +2 ${VDJDB_FILE} | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20\`                                  |

$(tail -n +2 ${VDJDB_FILE} | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20)

**Human-only records**

| Metric                                                             | Count                                                                                                                                                                          | Command                                                                                                                                                                                                 |
|--------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|   Total number of records                                          | $(awk '$6 == "HomoSapiens" { print }' ${VDJDB_FILE} |  wc -l)                                                                                                                  | \`awk '\$6 == "HomoSapiens" { print }' ${VDJDB_FILE} |  wc -l\`                                                                                                                                         |
|   TRA records                                                      | $(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3 }' ${VDJDB_FILE} | wc -l)                                                                                                 | \`awk '\$2 == "TRA" && \$6 == "HomoSapiens"  { print \$3 }' ${VDJDB_FILE} | wc -l\`                                                                                                                     |
|   TRB records                                                      | $(awk '$2 == "TRB" && $6 == "HomoSapiens"  { print $3 }' ${VDJDB_FILE} | wc -l)                                                                                                | \`awk '\$2 == "TRB" && \$6 == "HomoSapiens"  { print \$3 }' ${VDJDB_FILE} | wc -l\`                                                                                                                     |
|   Unique TRA sequences                                             | $(awk '$2 == "TRA" && $6 == "HomoSapiens"  { print $3 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                                      | \`awk '\$2 == "TRA" && \$6 == "HomoSapiens"  { print \$3 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                                                           |
|   Unique TRB sequences                                             | $(awk '$2 == "TRB" && $6 == "HomoSapiens"  { print $3 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                                      | \`awk '\$2 == "TRB" && \$6 == "HomoSapiens"  { print \$3 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                                                           |
|   Unique CDR3 sequences                                            | $(awk '$6 == "HomoSapiens" { print }' ${VDJDB_FILE} | cut -f3 | sort -u | wc -l)                                                                                               | \`awk '\$6 == "HomoSapiens" { print }' ${VDJDB_FILE} | cut -f3 | sort -u | wc -l\`                                                                                                                      |
|   Unique epitope sequences                                         | $(awk '$6 == "HomoSapiens" { print }' ${VDJDB_FILE} | cut -f10 | sort -u | wc -l)                                                                                              | \`awk '\$6 == "HomoSapiens" { print }' ${VDJDB_FILE} | cut -f10 | sort -u | wc -l\`                                                                                                                     |
|   Unique epitope sequences for TRA records                         | $(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $10 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                                      | \`awk '\$2 == "TRA" && \$6 == "HomoSapiens" { print \$10 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                                                           |
|   Unique epitope sequences for TRB records                         | $(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $10 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                                      | \`awk '\$2 == "TRB" && \$6 == "HomoSapiens" { print \$10 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                                                           |
|   Unique CDR3-epitope sequence pairs                               | $(awk '$6 == "HomoSapiens" { print }' ${VDJDB_FILE} | cut -d $'\t' -f3,10 | sort -u | wc -l)                                                                                   | \`awk '\$6 == "HomoSapiens" { print }' ${VDJDB_FILE} | cut -d $'\t' -f3,10 | sort -u | wc -l\`                                                                                                          |
|   Unique TRA-CDR3-epitope sequence pairs                           | $(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3,$10 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                                   | \`awk '\$2 == "TRA" && \$6 == "HomoSapiens" { print \$3,\$10 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                                                       |
|   Unique TRB-CDR3-epitope sequence pairs                           | $(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $3,$10 }' ${VDJDB_FILE} | sort -u | wc -l)                                                                                   | \`awk '\$2 == "TRB" && \$6 == "HomoSapiens" { print \$3,\$10 }' ${VDJDB_FILE} | sort -u | wc -l\`                                                                                                       |
|   Number of epitope sequences shared between TRA and TRB records   | $(comm -12 <(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $10 }' ${VDJDB_FILE} | sort -u) <(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $10 }' ${VDJDB_FILE} | sort -u) | wc -l)|\`comm -12 <(awk '\$2 == "TRA" && \$6 == "HomoSapiens" { print \$10 }' ${VDJDB_FILE} | sort -u) <(awk '\$2 == "TRB" && \$6 == "HomoSapiens" { print \$10 }' ${VDJDB_FILE} | sort -u) | wc -l\`   |
|   Number of CDR3 sequences shared between TRA and TRB records      | $(comm -12 <(awk '$2 == "TRA" && $6 == "HomoSapiens" { print $3 }' ${VDJDB_FILE} | sort -u) <(awk '$2 == "TRB" && $6 == "HomoSapiens" { print $3 }' ${VDJDB_FILE} | sort -u)   | wc -l)|\`comm -12 <(awk '\$2 == "TRA" && \$6 == "HomoSapiens" { print \$3 }' ${VDJDB_FILE} | sort -u) <(awk '\$2 == "TRB" && \$6 == "HomoSapiens" { print \$3 }' ${VDJDB_FILE} | sort -u) | wc -l\`     |
|   Epitope distribution for the unique CDR3-epitope pairs           |                             | \`awk '\$6 == "HomoSapiens" { print }' ${VDJDB_FILE} | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20\`                                                                 |

$(awk '$6 == "HomoSapiens" { print }' ${VDJDB_FILE} | cut -d $'\t' -f3,10 | sort -u | cut -f2 | sort | uniq -c | sort -nr | head -20)

EOF
