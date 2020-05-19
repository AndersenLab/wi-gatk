# AD_DP Filter
# Author: Daniel E. Cook
# Compile with:
#
# nim c -d:release ad_dp.nim
# 
# This little utility will polarize variants
#
# FORMAT/DV -  Number of high-quality non-reference bases, (Number=1,Type=Integer)
#
# Changed to:
# 
# FORMAT/AD = Allelic Depths [in GATK]
#
import hts
import strformat
import sequtils
import math

const version = "0.0.1"

var wtr:VCF
var v:VCF
doAssert(open(v, "/dev/stdin"))
doAssert(open(wtr, "/dev/stdout", mode="wb"))
wtr.header = v.header
discard wtr.header.add_format(ID = "HP", Number = "1", Type = "String", Description = fmt"Flag used to mark whether a variant was polarized [AA/BB]; Version {version}")
discard wtr.header.add_format(ID = "HP_VAL", Number = "1", Type = "Float", Description = fmt"-log10(GL-ref/GL-alt); > 2 ALT Polarization; < 2 REF polarization; Version {version}")
doAssert(wtr.write_header())

var n_samples = len(v.samples)

proc between(s: int, a: int, b: int): bool = 
    if a <= s and s <= b:
        return true
    return false


proc rev_phred_to_p(phred: int): float =
     return math.pow(10.0, phred.float / -10.0)

proc is_heterozygous(gt: seq[Allele]): bool = 
    let gt_set = gt.toSeq()
    if (gt_set[0].value() > 0 or
       gt_set[1].value() > 0) and 
       gt_set.len == 2 and 
       gt_set[0].value() != gt_set[1].value():
        return true
    return false

proc gtval(a:int): int32 {.inline.} =
  return (cast[int32](a) + 1) shl 1

var pl_set: seq[int32]
var log_set: seq[float]

for record in v:
    var n_alts = record.ALT.len
    var arr_len = n_samples*(n_alts + 1)
    var pl = new_seq[int32](arr_len)
    var hp = new_seq[string](n_samples)
    var hp_val = new_seq[string](n_samples)
    var gts = new_seq[int32](n_samples)
    var gt = record.format.genotypes(gts)
    doAssert record.format.get("PL", pl) == Status.OK
    # Add test for biallelic
    var idx = 0
    for geno in gt:
        # Only operatate on heterozygous variants
        if geno.is_heterozygous():
            pl_set = pl[ idx * 3 .. (idx * 3) + 2]#.applyIt ( if it < 0: 0 else: it  )
            if pl_set[0].between(-1000, 1000) and pl_set[2].between(-1000, 1000):
                log_set = pl_set.mapIt( rev_phred_to_p(it) )
                var log_score = math.log10(log_set[0] / log_set[2])
                if log_score <= -2.0:
                    gts[(idx*2)] = geno[0].value().gtval()
                    gts[(idx*2) + 1] = geno[0].value().gtval()
                    hp[idx] = fmt"AA"
                elif log_score >= 2.0:
                    gts[(idx*2)] = geno[1].value().gtval()
                    gts[(idx*2) + 1] = geno[1].value().gtval()
                    hp[idx] = fmt"BB"
                else:
                    hp[idx] = "AB"
                hp_val[idx] = fmt"{log_score:<0.3}"
        else:
            # No het polarization fields
            hp[idx] = "."
            hp_val[idx] = "."
        idx += 1
        # Set
        if hp.filterIt( it != "" ).len > 0:
            doAssert record.format.set("HP", hp) == Status.OK
            doAssert record.format.set("HP_VAL", hp_val) == Status.OK
            doAssert record.format.set("GT", gts) == Status.OK
    doAssert wtr.write_variant(record)

close(wtr)
close(v)

