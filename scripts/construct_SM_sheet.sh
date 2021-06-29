#!/bin/bash

fq_sheet=`mktemp`

#===============================#
# 140905_D00422_0098_AHAJR1ADXX #
#===============================#

seq_folder=140905_D00422_0098_AHAJR1ADXX
prefix=/projects/b1059/data/fastq/WI/dna/processed/$seq_folder
ls -1 $prefix/*1P.fq.gz |\
xargs -n1 basename |\
awk -v prefix=${prefix} -v seq_folder=${seq_folder} '{
    fq1 = $1;
    fq2 = $1;
    gsub("1P.fq.gz", "2P.fq.gz", fq2);
    gsub("_1P.fq.gz", "", ID);
    line = $0;
    gsub("CB4856_CGC_3", "CB4856", $0);
    gsub("CB4857_CGC", "CB4857", $0);
    ID = $1;
    split($0, a, "_");
    SM = a[2];
    gsub("CB4857", "ECA249", SM);
    gsub("N2Baer", "ECA254", SM);
    gsub("-", "", a[3]);
    LB = a[3];
    print SM "\t" ID "\t" LB "\t" prefix "/" fq1 "\t" prefix "/" fq2 "\t" seq_folder;
}' >> ${fq_sheet}

#===================================#
# 151009_D00422_0262_BC7NJ0ANXX-ECA #
#===================================#

seq_folder=151009_D00422_0262_BC7NJ0ANXX-ECA
prefix=/projects/b1059/data/fastq/WI/dna/processed/$seq_folder
ls -1 $prefix/*1P.fq.gz |\
xargs -n1 basename |\
awk -v prefix=${prefix} -v seq_folder=${seq_folder} '{    
    fq1 = $1;
    fq2 = $1;
    gsub("1P.fq.gz", "2P.fq.gz", fq2);
    gsub("_1P.fq.gz", "", ID);
    split($0, a, "_");
    ID = a[1] "_" a[2] "_" a[3];
    gsub("-", "", a[2]);
    LB = a[2];
    SM = a[1];
    gsub("CB4857", "ECA249", SM);
    gsub("N2Baer", "ECA254", SM);
    print SM "\t" ID "\t" LB "\t" prefix "/" fq1 "\t" prefix "/" fq2 "\t" seq_folder;
}' >> ${fq_sheet}

#======================================================#
# 160209_D00422_0276_AC80RDANXX-ECA                    #
# 160512_700819F_0460_BHHHWNBCXX-ECA-14-ln1            #
# 160516_700819F_0461_AHNMLKBCXX-ECA-14-ln2-ECA-15-ln1 #
# 160518_700819F_0462_AHT25JBCXX-ECA-15-ln2            #
# 160518_700819F_0463_BHT25HBCXX-ECA-17-ln1            #  
# 160901_D00422_0374_AC7TWJANXX-ECA                    #
#======================================================#

for seq_may_2016 in 160209_D00422_0276_AC80RDANXX-ECA 160512_700819F_0460_BHHHWNBCXX-ECA-14-ln1 160516_700819F_0461_AHNMLKBCXX-ECA-14-ln2-ECA-15-ln1 160518_700819F_0462_AHT25JBCXX-ECA-15-ln2 160518_700819F_0463_BHT25HBCXX-ECA-17-ln1 160901_D00422_0374_AC7TWJANXX-ECA; do
    out=`mktemp`
    seq_folder=${seq_may_2016}
    prefix=/projects/b1059/data/fastq/WI/dna/processed/$seq_folder
    for i in `ls -1 $prefix/*1P.fq.gz`; do
        bname=`basename ${i}`;
        barcode=`zcat ${i} | grep '@' | cut -f 10 -d ':' | sed 's/+//g' | head -n 20 | uniq -c | sort -k 1,1n | cut -c 9-100 | tail -n 1`
        echo -e "${bname}\t${i}\t${barcode}" >> ${out}
    done;

    cat ${out} |\
    awk -v prefix=${prefix} -v seq_folder=${seq_folder} '{    
        fq1 = $1;
        fq2 = $1;
        LB = $3;
        gsub("1P.fq.gz", "2P.fq.gz", fq2);
        ID = $1;
        gsub("_1P.fq.gz", "", ID);
        split(ID, a, "[-_]")
        SM=a[2];
        print SM "\t" ID "\t" LB "\t" prefix "/" fq1 "\t" prefix "/" fq2 "\t" seq_folder;
    }' >> ${fq_sheet}
done;


#===================================#
# BGI-20161012-ECA21-ECA22          #
#===================================#

out=`mktemp`
seq_folder=BGI-20161012-ECA21-ECA22
prefix=/projects/b1059/data/fastq/WI/dna/processed/$seq_folder
for i in `ls -1 $prefix/*1P.fq.gz`; do
    bname=`basename ${i}`;
    barcode=`zcat ${i} | grep '@' | cut -f 10 -d ':' | sed 's/_//g' | head -n 100 | uniq -c | sort -k 1,1n | cut -c 9-100 | tail -n 1`
    echo -e "${bname}\t${i}\t${barcode}" >> ${out}
done;

cat ${out} |\
awk -v prefix=${prefix} -v seq_folder=${seq_folder} '{
    fq1 = $1;  
    fq2 = $1;
    LB = $3;
    gsub("N", "", LB);
    gsub("1P.fq.gz", "2P.fq.gz", fq2);
    ID = $1;
    gsub("_1P.fq.gz", "", ID);
    split(ID, a, "[-_]")
    SM=a[4];
    gsub("JU1516_MAF", "JU1516", SM);
    print SM "\t" ID "\t" LB "\t" prefix "/" fq1 "\t" prefix "/" fq2 "\t" seq_folder;
}' >> ${fq_sheet}

#===================================#
# BGI-20161012-ECA23                #
#===================================#

out=`mktemp`
seq_folder=BGI-20161012-ECA23
prefix=/projects/b1059/data/fastq/WI/dna/processed/$seq_folder
for i in `ls -1 $prefix/*1P.fq.gz`; do
    bname=`basename ${i}`;
    barcode=`zcat ${i} | grep '@' | cut -f 10 -d ':' | sed 's/_//g' | head -n 100 | uniq -c | sort -k 1,1n | cut -c 9-100 | tail -n 1`
    echo -e "${bname}\t${i}\t${barcode}" >> ${out}
done;

cat ${out} |\
awk -v prefix=${prefix} -v seq_folder=${seq_folder} '{
    fq1 = $1;
    fq2 = $1;
    LB = $3;
    gsub("N", "", LB);
    gsub("1P.fq.gz", "2P.fq.gz", fq2);
    ID = $1;
    gsub("_1P.fq.gz", "", ID);
    split(ID, a, "[-_]")
    SM=a[2];
    print SM "\t" ID "\t" LB "\t" prefix "/" fq1 "\t" prefix "/" fq2 "\t" seq_folder;
}' >> ${fq_sheet}

#===================================#
# original_wi_set                   #
#===================================#

seq_folder=original_wi_set 
prefix=/projects/b1059/data/fastq/WI/dna/processed/${seq_folder}

ls -1 ${prefix}/*1P.fq.gz |\
xargs -n1 basename |\
awk  -F  "-" -v prefix=${prefix} -v seq_folder=${seq_folder} '{
    fq1 = $0;  
    fq2 = $0;
    gsub("1P.fq.gz", "2P.fq.gz", fq2);
    LB = $2;
    ID = $1 "-" $2 "-" $3;
    SM=$3;
    gsub("CB4851", "ECA243", SM);
    gsub("CB4853", "ECA245", SM);
    gsub("CB4857", "ECA249", SM);
    gsub("CB4855", "ECA247", SM);
    gsub("CB4858", "ECA248", SM);
    gsub("PB306", "ECA259", SM);
    print SM "\t" ID "\t" LB "\t" prefix "/" fq1 "\t" prefix "/" fq2 "\t" seq_folder;
}' >> ${fq_sheet}


#===========================================================#
# MAF-JU1249-20170425 - Single Strain from Marie Anne Felix #
#===========================================================#

seq_folder=MAF-JU1249-20170425 
prefix=/projects/b1059/data/fastq/WI/dna/processed/${seq_folder}

echo -e "JU1249\tMAF-JU1259-20170425\tMAF-JU1259-20170425\t${prefix}/JU1249_1P.fq.gz\t${prefix}/JU1249_2P.fq.gz\t${seq_folder}" >> ${fq_sheet}


#===================================#
# Convert Strains to Isotypes; out  #
#===================================#
cat ${fq_sheet} | python strain_to_isotype.py | sort > ../SM_sample_sheet.tsv
