#!/usr/bin/env bash

db=$1
smiles_file=$2
out_file=$3

search_one() {
    local db="$1"
    local smiles="$2"
    local query_id="$3"

    obabel "$db" \
        -osmi \
        -s "$smiles" \
        -at3 \
        -aa 2>/dev/null |
    awk -v q="$query_id" -v qs="$smiles" '
    BEGIN { OFS="\t" }
    {
        print q, qs, $(NF-1), $1, $NF
    }'
}
export -f search_one
export db

printf "query_id\tquery_SMILES\tcompound_id\tcompound_SMILES\ttanimoto\n" > $out_file

parallel --jobs 32 --bar --joblog $out_file.log --colsep '\t' 'search_one "$db" {1} {2}' :::: $smiles_file >> $out_file