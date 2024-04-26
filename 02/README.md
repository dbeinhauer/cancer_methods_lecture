# Task 2: Cancer Genomics Data Analysis

First, if not indexed then index the reference genome (computationally demanding task):

```bash
bwa index data/hg19.fa
```

Then, execute the analysis pipeline for both tumor and wildtype samples (for further details see [pipeline.sh](pipeline.sh)):

```bash
./pipeline.sh data/hg19.fa data/tu.r1.fq data/tu.r2.fq "tu" tumor_analysis
./pipeline.sh data/hg19.fa data/wt.r1.fq data/wt.r2.fq "wt" wildtype_analysis
```

Finally, the read depth plot analysis can be seen in [read_depth.ipynb](read_depth.ipynb).