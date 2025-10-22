#!/bin/env python3

import gzip
import sys
import subprocess
import argparse
import re
import vcfpy

indel_re = re.compile(r"(\+|\-)[0-9]+[ACGTNacgtn]+")
carret_re = re.compile(r"\^.")
junk_re = re.compile(r"[.,<>*$]")

def cleanup_pileup(pileup):
    try:
        cleaned =  ''.join(set(c.upper() for c in junk_re.sub("",
                                                   carret_re.sub("",
                                                                 indel_re.sub("", pileup)))))
    except Exception as e:
        print("[parse_vcf.py] :: skipping irregular input:", pileup, file = sys.stderr)
        return None
    # ambiguous sites with multiple contigs mapping to this locus, or reference sites
    if len(cleaned) > 1 or len(cleaned) == 0:
        return None
    return cleaned

def parse_pileup(hap1, hap2):
    hap1_snv = cleanup_pileup(hap1)
    hap2_snv = cleanup_pileup(hap2)
    return(hap1_snv, hap2_snv)

parser = argparse.ArgumentParser(description='Extract SNVs from haplotype resolved genome to genome alignments.')
parser.add_argument('--reference', type = str, nargs = 1, help='Reference genome.')
parser.add_argument('--vcf_out', type = str, nargs = 1, help='Output path for SNV VCF file.')
parser.add_argument('--vcf_template', type = str, nargs = 1, help='Template VCF file.')
parser.add_argument('--sample', type = str, nargs = 1, help='Name of sample being genotyped.', default = "sample")
parser.add_argument('--hap1', type = str, help='First haplotype fasta.', nargs = 1)
parser.add_argument('--hap2', type = str, help = 'Second haplotype fasta.', nargs = 1)
parser.add_argument('--region', type = str, help = 'Subset alignments to region.', nargs = 1)
parser.add_argument('--min_quality', type = int, help='Subset alignments to region.', nargs = 1, default = 20)

args = parser.parse_args()

h1 = subprocess.Popen(
    ['samtools', 'mpileup', '-aa',
     '-q{q}'.format(q=args.min_quality[0]),
     '-r{region}'.format(region = args.region[0]),
     '-f{ref}'.format(ref = args.reference[0]),
     '{bam}'.format(bam = args.hap1[0])],
    stdout = subprocess.PIPE, encoding = 'ascii')

h2 = subprocess.Popen(
    ['samtools', 'mpileup', '-aa',
     '-q{q}'.format(q=args.min_quality[0]),
     '-r{region}'.format(region = args.region[0]),
     '-f{ref}'.format(ref = args.reference[0]),
     '{bam}'.format(bam = args.hap2[0])],
    stdout = subprocess.PIPE, encoding = 'ascii')

i = 0

reader = vcfpy.Reader.from_path(args.vcf_template[0])
reader.header.samples = vcfpy.SamplesInfos([args.sample[0]])
writer = vcfpy.Writer.from_path(args.vcf_out[0], reader.header)

while True:
    hap1_pileup = h1.stdout.readline()
    hap2_pileup = h2.stdout.readline()
    if not hap1_pileup or not hap2_pileup:
        break

    chrom1, pos1, ref1, _, h1_alt, _ = hap1_pileup.split()
    chrom2, pos2, ref2, _, h2_alt, _ = hap2_pileup.split()

    pos1, pos2 = int(pos1), int(pos2)

    assert chrom1 == chrom2 and pos1 == pos2
    a1, a2 = parse_pileup(h1_alt, h2_alt)

    gt = None
    alt = None
    if a1 or a2:
        if a1 == a2:
            gt = "1|1"
            alt = [vcfpy.Substitution("SNV", a1)]
        elif a1 and not a2:
            gt = "1|0"
            alt = [vcfpy.Substitution("SNV", a1)]
        elif not a1  and a2:
            gt = "0|1"
            alt = [vcfpy.Substitution("SNV", a2)]
        elif a1 and a2:
            gt = "1|2"
            alt = [vcfpy.Substitution("SNV", a1), vcfpy.Substitution("SNV", a2)]

        rec = vcfpy.Record(CHROM = chrom1, POS = pos1, ID = ".",
                        REF = ref1, ALT = alt,
                        QUAL = 999, FILTER = ["PASS"], INFO = {}, FORMAT = ["GT"],
                           calls = [vcfpy.Call(sample = args.sample[0], data = vcfpy.OrderedDict(GT = gt))])

        writer.write_record(rec)

    if i % 1e6 == 0:
        print("[parse_snvs.py] :: parsing {chrom}:{start}:{end}".format(chrom = chrom1, start=pos1, end = pos1 + 1e6), file = sys.stderr)
    i = i + 1

writer.close()
