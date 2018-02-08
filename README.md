# variantcalling
Repository for variant calling pipeline using FreeBayes

Note that you have to incorporate the path to the reference genome manually:

Add
```
"fasta": "/lustre_cfc/qbic/reference_genomes/Homo_sapiens/DNA/2016.01.21.UCSC/hg19/Sequence/BWAIndex/hg19/hg19"
```
in
```
params.json
```
and feed it as input for qproject:

```
qproject create -t . -w github:qbicsoftware/variantcalling --group qbicgrp --params params.json
```
