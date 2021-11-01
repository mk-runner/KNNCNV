# KNNCNV

KNNCNV: A k-nearest neighbor based method for detection of copy number variations using NGS data

------

# Usage

## 1. Input

- bam_path (str, for the six real blood samples in Section 3.2.1, the *.bam file can be obtained from the [1000 Genomes Project](http://www.1000genomes.org). The three cancer samples in Section 3.2.2 can be downloaded from the [European Genome-Phenome Archive](https://ega-archive.org/). Additionally, these *.bam files can be download from the [Baidu Netdisk](https://pan.baidu.com/s/1Ja4XH2wZupeAcwc9qhZn8A), and its extraction code is `29to`).
- fa_path (str, the reference genome that can be in the fa format or the fasta format).
- gt_path (str, optional, the confirmed CNVs file of the *.bam file. if this file is provided, some performance metrics including precision and sensitivity can be calculated. Note that the file can be obtained from the  [database of genomic variants](http://dgv.tcag.ca/dgv/app/home)).

## 2. Output of the entire process

The experimental result on the NA12878 is NA12878.txt that is a txt text delimited by commas and represents the detailed descriptions of CNVs predicted by KNNCNV. More specifically, the first line of the text is the column descriptions including chr, start, end, variant type, and RD, where chr denotes the chromosome ID, and start and end represent the start and end positions of declared CNVs, respectively. The variant type contains deletion and duplication, and the RD is the read depth of declared CNVs.

## 3. API

- [ ] ### API for the entire process

- open the file `knncnv.py` and modify the variables `bam_path`, `fa_path` inside; 

```python
if __name__ == '__main__':
    # Local path of the *.bam file
    bam_path = r"/real_data/NA12878.bam"
    
    # Local path of the *.fasta file or the *.fa file
    fa_path = r"./data/chr21.fa"
    
    # Local path of the ground truth (i.e., confirmed CNVs) for the *.bam file.
    gt_path = r"./data/NA12878.gt"
    
    # parameter setting of the preprocessing
    bin_size = 1000  # the bin size ('1000' by default)

    knncnv(bam_path, fa_path, gt_path=gt_path, iter_num=20, bin_size=bin_size)
```

- run the `knncnv.py`;

- output of the entire process is NA12878.txt.

- [ ] ### API only for the VBGMM

```python
from ihybcnv import vbgmm

labels = vbgmm(scores)
# 0 stands for inliers and 1 for outliers(CNVs).
```

------

# Required Dependencies

Python 3.8            

- biopython     1.78
- numpy         1.18.5
- pandas        1.0.5
- pysam         0.16.0.1
- pyod          0.8.4
- rpy2          3.4.2
- scikit-learn  0.23.1
- scipy         1.5.0

R 3.4.4

- DNAcopy
