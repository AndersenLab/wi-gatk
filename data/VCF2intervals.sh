#!/bin/bash
# Usage VCF2intervals.sh <genome>.fai <vcf.gz> <breaks>


LINES=$(zcat $2 | grep -v "^#" | wc -l | cut -f 1 -d" ")
zcat $2 | \
grep -v "^#" | \
awk -v LINES="${LINES}" -v BREAKS="$3" '
BEGIN{
    START=1;
    STOP=1;
    CHROM="";
    SPAN=int((LINES - 1)/BREAKS + 1);
    COUNT=0;
    SKIP=TRUE;
}{
    if (NR == FNR) {
        CHROM_STOP[$1] = $2;
    } else {
        if ($1 != CHROM) {
            if (SKIP == FALSE) {
                EFF_STOP = CHROM_STOP[CHROM];
                printf "%s\t%i\t%i\n", CHROM, START, EFF_STOP;
            }
            SKIP=FALSE;
            START=1;
            STOP=1;
            CHROM=$1;
            COUNT=0;
        } else if (COUNT == SPAN) {
            if (STOP + 300 > CHROM_STOP[CHROM]) {
                EFF_STOP = CHROM_STOP[CHROM];
                SKIP = TRUE;
            } else {
                EFF_STOP = STOP + 300;
            }
            printf "%s\t%i\t%i\n", CHROM, START, EFF_STOP;
            START=STOP;
            STOP=$2;
            COUNT=1;
        } else if (SKIP == FALSE) { 
            COUNT+=1; 
            STOP=$2;
        }
    }
}END{
    if (SKIP == FALSE) {
        EFF_STOP = CHROM_STOP[CHROM];
        printf "%s\t%i\t%i\n", CHROM, START, EFF_STOP;
    }
}' $1 -