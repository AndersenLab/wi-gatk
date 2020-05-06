# AD_DP Filter
# Author: Daniel E. Cook
# Compile with:
#
# nim c -d:release ad_dp.nim
# 
# This little utility will add a DV/DP filter
# to the VCF.
#
# Previously, we used a DV/DP filter of 0.5,
# DV = # of high quality bases underlying a variant.
# DV has since been renamed AD
# See https://samtools.github.io/bcftools/bcftools.html
#
# 
# FORMAT/DV -  Number of high-quality non-reference bases, (Number=1,Type=Integer)
#
# Changed to:
# 
# FORMAT/AD = Allelic Depths [in GATK]
#
import hts
import strutils

# AD/DP DV/DP annotation script
# This script marks the FT Format field with a ad_dp filter

proc isHomRef(gt: Genotype): bool = 
    for i in gt:
        if i.value() > 0:
            return false
    return true

proc isMissing(gt: Genotype): bool = 
    for i in gt:
        if i.value() < 0:
            return true
    return false

var wtr:VCF
var v:VCF
doAssert(open(v, "-"))
doAssert(open(wtr, "ad_dp.filtered.vcf.gz", mode="w"))
wtr.header = v.header
doAssert(wtr.write_header())

var n_samples = len(v.samples)

for record in v:
    var n_alts = record.ALT.len
    var arr_len = n_samples*(n_alts + 1)
    var ad = new_seq[int32](arr_len)
    var dp = new_seq[int32](arr_len)
    var gts = new_seq[int32](n_samples)
    var ft = new_seq[string](n_samples)
    var gt = record.format.genotypes(gts)
    doAssert record.format.get("AD", ad) == Status.OK
    doAssert record.format.get("DP", dp) == Status.OK
    try:
        doAssert record.format.get("FT", ft) == Status.OK
    except:
        # If FT not present, set values
        for i in 0..<ft.len:
            if gt[i].isMissing():
                ft[i] = "."
            else:
                ft[i] = "PASS"
    for sample_index in 0..<n_samples:
        if gt[sample_index].isHomRef() == false:
            var ad_set: seq[int32]
            for alt_index in 1..n_alts:
                var offset = (sample_index * n_alts) + sample_index + alt_index
                ad_set.add if ad[offset] >= 0: ad[offset] else: 0
            if ad_set.len == 0:
                continue
            if (max(ad_set).float / dp[sample_index].float) < 0.5:
                if ft[sample_index] in ["PASS", "."]:
                    ft[sample_index] = "ad_dp"
                else:
                    if ft[sample_index].split(";").find("ad_dp") == -1:
                        ft[sample_index] = "ad_dp;" & ft[sample_index]
    discard record.format.set("FT", ft)
    doAssert wtr.write_variant(record)

close(v)
close(wtr)

