#!/bin/bash

LINES=$(zcat $2 | grep -v "^#" | wc -l | cut -f 1 -d" ")
zcat $2 | \
grep -v "^#" | \
awk -v LINES="${LINES}" '
BEGIN{
    START=1;
    STOP=1;
    CHROM="";
    SPAN=int((LINES - 1)/200 + 1);
    COUNT=0;
}{
    if (NR == FNR) {
        if ($1 == "@SQ") {
            split($2,CHROMNAME,":");
            split($3,CHROMLEN,":");
            CHROM_STOP[CHROMNAME[2]] = CHROMLEN[2];
        }
    } else {
        if ($1 != CHROM) {
            if (CHROM != "") {
                printf "%s\t%i\t%i\n", CHROM, START, STOP+300;
                START=1;
                STOP=CHROM_STOP[CHROM];
                CHROM=$1;
                COUNT=0;
            } else {
                CHROM=$1;
            }
        } else if (COUNT == SPAN) {
            printf "%s\t%i\t%i\n", CHROM, START, STOP+300;
            START=STOP;
            STOP=$2;
            COUNT=1;
        } else { 
            COUNT+=1; 
            STOP=$2;
        }
    }
}END{
    printf "%s\t%i\t%i\n", CHROM, START, CHROM_STOP[CHROM];
}' $1 -