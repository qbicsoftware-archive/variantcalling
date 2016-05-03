import os
import sys
import subprocess
import tempfile
import uuid
import shutil
from datetime import datetime
from os.path import join as pjoin
from os.path import exists as pexists
from xml.etree import ElementTree
import hashlib
import base64
import csv
import glob

configfile: "config.json"
workdir: config['var']

SNAKEDIR = config['src']

try:
    VERSION = subprocess.check_output(
        ['git', 'describe', '--tags', '--always', '--dirty'],
        cwd=SNAKEDIR
    ).decode().strip()
except subprocess.CalledProcessError:
    VERSION = 'unknown'

DATA = config['data']
RESULT = config['result']
LOGS = config['logs']
REF = config['ref']
ETC = config['etc']

def data(path):
    return os.path.join(DATA, path)

def ref(path):
    return os.path.join(REF, path)

def log(path):
    return os.path.join(LOGS, path)

def result(path):
    return os.path.join(RESULT, path)

def etc(path):
    return os.path.join(ETC, path)

if 'params' not in config:
    config['params'] = {}


INPUT_FILES = []
for name in os.listdir(DATA):
    if name.lower().endswith('.sha256sum'):
        continue
    if name.lower().endswith('.bam'):
        if not name.endswith('.bam'):
            print("Extension bam and bam.bai is case sensitive.", file=sys.stderr)
            exit(1)
        #INPUT_FILES.append(os.path.basename(name)[:-len('.bam')])
        INPUT_FILES.append(os.path.basename(name))
    elif name.lower().endswith('bam.bai'):
        if not name.endswith('bam.bai'):
            print("Extension bam and bam.bai is case sensitive.", file=sys.stderr)
            exit(1)
        #INPUT_FILES.append(os.path.basename(name)[:-len('.bam.bai')])
        INPUT_FILES.append(os.path.basename(name))
    elif name.lower().endswith('.bam.gz'):
        if not name.endswith('.bam.gz'):
            print("Extension bam is case sensitive.", file=sys.stderr)
            exit(1)
        #INPUT_FILES.append(os.path.basename(name)[:-len('.bam.gz')])
        INPUT_FILES.append(os.path.basename(name))
    else:
        print("Unknown data file: %s" % name, file=sys.stderr)
        exit(1)

if len(set(INPUT_FILES)) != len(INPUT_FILES):
    print("Some input file names are not unique")
    exit(1)

OUTPUT_FILES = []
OUTPUT_FILES.append("checksums.ok")
OUTPUT_FILES.append(result("variants_merged_sorted.vcf"))

rule all:
    input: OUTPUT_FILES

rule checksums:
    output: "checksums.ok"
    threads: 1
    run:
        out = os.path.abspath(str(output))
        cksums = glob.glob(data("*.sha256sum"))
        if cksums:
            shell("cd %s; "
                  "sha256sum -c *.sha256sum && "
                  "touch %s" % (data('.'), out))
        else:
            shell("touch %s" % out)

#variant calling preparation
rule left_align:
    input: bam=data("aligned.bam"), \
           bai=data("aligned.bam.bai") 
    output: "left_align/{name}.unsorted.bam"
    run:
        shell("BamLeftAlign -in {input.bam} -out {output} -ref %s.fa"
               % ref(config['params']['fasta']))


rule re_align_intervals:
    input: bam="left_align/{name}.unsorted.bam", \
           bai="left_align/{name}.unsorted.bam.bai"
    output: "left_align/{name}.intervals"
    run:
        shell(("java -jar %s/GenomeAnalysisTK.jar "  # TODO Lite
               "-T RealignerTargetCreator "
               "-I {input.bam} "
               "-R %s.fa "
               "-o {output}")
              % (os.environ['GATK_JARS'],
                 ref(config['params']['fasta'])))
        #"--intervals {config['intervals']}")


rule re_align:
    input: bam="left_align/{name}.unsorted.bam", \
           bai="left_align/{name}.unsorted.bam.bai", \
           intervals="left_align/{name}.intervals"
    output: "realigned/{name}.unsorted.bam"
    run:
        shell(("java -jar %s/GenomeAnalysisTK.jar "  # TODO Lite
               "-T IndelRealigner "
               "-I {input.bam} "
               "-R %s.fa "
               "--targetIntervals {input.intervals} "
               "-o {output}")
              % (os.environ['GATK_JARS'],
                 ref(config['params']['fasta'])))


rule clip_overlap:
    input: bam="realigned/{name}.unsorted.bam", \
           bai="realigned/{name}.unsorted.bam.bai"
    output: "clipped/{name}.unsorted.bam"
    run:
        shell("BamClipOverlap -in {input.bam} -out {output}")

##index new bam file
rule merge_sam_index:
    input:  "clipped/{name}.unsorted.bam"
    output: "clipped/{name}.unsorted.bam.bai"
    run:
        shell("samtools index {input}")

## create basic variant calls ,here with freebayes
rule freebayes:
    input: bam="clipped/{name}.unsorted.bam" \
           bai="clipped/{name}.unsorted.bam.bai" 
    output: "vcf_base/merged_base.vcf"
    run:
        fasta = ref(config['params']['fasta'])
        shell("freebayes --min-base-quality 20 --min-alternate-qsum 90 -f %s.fa -b {input.bam} > {output}" % fasta)

## filter variants according to variant quality>5 , alternate observations>=2
rule filter_variants:
    input: "vcf_base/merged_base.vcf"
    output: "vcf_filtered/merged_filtered.vcf"
    run:
        shell("vcffilter -f 'QUAL > 5 & AO > 2' {input} > {output}")

## split complex variants to primitives
## this step has to be performed before vcfbreak multi - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail
rule split_complex_variants:
    input: "vcf_filtered/merged_filtered.vcf"
    output: "vcf_prim/merged_prim.vcf"
    run:
        shell("vcfallelicprimitives -kg {input} > {output}")

rule vcffixup:
    input: "vcf_prim/merged_prim.vcf"
    output: "vcf_prim/merged_prim_fixed.vcf"
    run:
        shell("vcffixup {input} > {output}")

## split multi-allelic variants: If multiple alleles are specified in a single record, break the record into multiple lines, preserving allele-specific INFO fields.
rule split_multiallelic_variants:
    input: "vcf_prim/merged_prim_fixed.vcf"
    output: "vcf_split/merged_split.vcf"
    run:
        shell("vcfbreakmulti {input} > {output}")

## align INDELs to the left. Note that complex indels and multi-allelic deletions are not shifted!
rule re_align_INDELS:
    input: "vcf_split/merged_split.vcf"
    output: "vcf_aligned/merged_aligned.vcf"
    run:
        shell("bcftools norm -f %s.fa {input} -o {output}"
                % ref(config['params']['fasta']))

## sort variants by genomic position. 
rule sort_variants:
    input: "vcf_aligned/merged_aligned.vcf"
    output: result("variants_merged_sorted.vcf")
    run:
        shell("cat {input} | vcf-sort > {output}")

