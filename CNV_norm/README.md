---
title: "Normalization of cancer Hi-C dataLOIC - CAIC"
date: "05/06/2018"
output: html_document
---

This document describes how to run the LOIC / CAIC normalizations method implemented in the [iced](https://pypi.org/project/iced/#files) python package


# Python wrapper

The ice_cnv.py script allows to run the LOIC/CAIC normalization on Hi-C data generated by the HiC-Pro pipeline

```python
python ice_cnv.py -h
usage: ice_cnv.py [-h] [--remove-all-zeros-loci]
                  [--filter_low_counts_perc FILTER_LOW_COUNTS_PERC]
                  [--classic] [--bed-file BED_FILE] [--seg-file SEG_FILE]
                  [--max_iter MAX_ITER] [--remove-cnv] [--outfile OUTFILE]
                  [--type TYPE] [--dense]
                  filename

positional arguments:
  filename

optional arguments:
  -h, --help            show this help message and exit
  --remove-all-zeros-loci
  --filter_low_counts_perc FILTER_LOW_COUNTS_PERC
  --classic
  --bed-file BED_FILE, -b BED_FILE
  --seg-file SEG_FILE, -s SEG_FILE
  --max_iter MAX_ITER, -m MAX_ITER
                        Maximum number of iterations
  --remove-cnv, -r
  --outfile OUTFILE, -o OUTFILE
  --type TYPE, -t TYPE  Can be 'independant', 'joint', or 'chrjoint'
  --dense
```

# Python API

### Importing packages
+ **datasets** contains dataset example to practice 
+ **normalizations** contains normalization functions
+ **io** contains methods for loading and managing sparse data

```python
from iced import datasets
from iced import normalization
from iced import io
```


### Data loading
3 data files are needed to process CAIC normalization 

+ **counts** : count matrix can be in sparse format or matrix format
+ **lengths** : array with each chromosome length
+ **cnv** : smooth profil line array (for each bin) generated with `segment_hic_data.sh` (*file_seg.bed* $4^{th}$ column)
 
```python
counts, lengths, cnv = datasets.load_sample_cancer()

counts = io.load_counts(mat, lengths=None)
```

### ICE Normalization
To compute ICE, provide the **counts** matrix
```python
ice_normed = normalization.ICE_normalization(counts, counts_profile=None)
```

### LOIC Normalization
To compute LOIC Normalization, function requires to provide the profil line **cvn** and the matrix **counts** 
```python
loic_normed = normalization.ICE_normalization(counts, counts_profile=cnv)
```

### CAIC Normalization
CAIC Normalization requires to compute LOIC Normalization before. 
Then, call *estimate_block_biases* with **counts** matrix, chr **lengths** file and **cnv** profile line.
Final step is to divide, **loic_normed** matrix by the **block bias** matrix estimated
```python
loic_normed = normalization.ICE_normalization(counts, counts_profile=cnv)
block_biases = normalization.estimate_block_biases(counts, lengths, cnv)
caic_normed = loic_normed / block_biases
```

