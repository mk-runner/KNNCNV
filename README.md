# KNNCNV
KNNCNV: A k-nearest neighbor based method for detection of copy number variations using NGS data </br>

# Usage
1. API for the entire process
   
(1) open the file `knncnv.py`, and modify the variables `bam_path`, `fa_path` inside; </br>

Note that the *.fa file is the reference genome that can be in the fa format or the fasta format. The ground truth file is the confirmed CNVs of the *.bam file, which can be obtained from the [database of genomic variants](http://dgv.tcag.ca/dgv/app/home) and help us calculate some performance metrics, such as precision, sensitivity, and F1-score. Additionally, for the six real blood samples in Section 3.2.1, the *.bam file can be obtained from the [1000 Genomes Project](http://www.1000genomes.org). The three cancer samples in Section 3.2.2 can be downloaded from the [European Genome-Phenome Archive](https://ega-archive.org/).
   
```python
if __name__ == '__main__':
    # Local path of the *.bam file
    bam_path = r"/real_data/NA19238.chrom21.SLX.maq.SRP000032.2009_07.bam"
    
    # Local path of the *.fasta file or the *.fa file
    fa_path = r"./data/chr21.fa"
    
    # Local path of the ground truth (i.e., confirmed CNVs) for the *.bam file.
    gt_path = r"./data/NA19238.gt"
    
    # parameter setting of the preprocessing
    bin_size = 1000  # the bin size ('1000' by default)

    knncnv(bam_path, fa_path, gt_path=gt_path, iter_num=20, bin_size=bin_size)
```
   (2) run the `knncnv.py`;
   
   (3) check the `knncnv` function in the file `knncnv.py` for more information;
   
   (4) output of the entire process (e.g., output of NA12878.txt). </br>
The output format is a txt text delimited by commas, which represents the detained descriptions of CNVs predicted by KNNCNV. More specifically, the first line of the text is the column descriptions including chr, start, end, variant type, and RD, where chr denotes the chromosome ID, and start and end represent the start and end positions of predicted CNVs, respectively. The variant type contains deletion and duplication, and the RD is the read depth of predicted CNVs.

2. API only for the VBGMM

```python
from ihybcnv import vbgmm

labels = vbgmm(scores)
# 0 stands for inliers and 1 for outliers(CNVs).
```

# Real Datasets
- the real blood samples: clink this link：https://pan.baidu.com/s/1Ja4XH2wZupeAcwc9qhZn8A extraction code：29to or [1000 Genomes Project](https://www.internationalgenome.org/)
- The genome-wide samples can be obtained from the [European Genome-Phenome Archive](https://ega-archive.org/)

# Required Dependencies
1. Python 3.8            
    - biopython     1.78
    - numpy         1.18.5
    - pandas        1.0.5
    - pysam         0.16.0.1
    - pyod          0.8.4
    - rpy2          3.4.2
    - scikit-learn  0.23.1
    - scipy         1.5.0
2. R 3.4.4
    - DNAcopy
