# Batch Effects Correction With Unknown Subtypes for scRNA-seq Data (BUSseq)

## Contents

- [Overview](#overview)
- [Implementation](#implementation）
- [Remarks](#remarks)

# Overview

Single-cell RNA-sequencing (scRNA-seq) technologies enable the measurement of the transcriptome of individual cells, which provides unprecedented opportunities to discover cell types and understand cellular heterogeneity. Despite their widespread applications, single-cell RNA-sequencing (scRNA-seq) experiments are still plagued by batch effects and dropout events.

One of the major tasks of scRNA-seq experiments is to identify cell types for a population of cells. Therefore, the cell type of each individual cell is always unknown and is the target of inference. However, most existing methods for batch effects correction, such as Combat and the surrogate variable analysis (SVA), are designed for bulk experiments and require knowledge of the subtype information, which corresponds to cell type information for scRNA-seq data, of each sample a priori.
  
Here, our proposed `BUSseq` method fits an interpretable Bayesian hierarchical model---the Batch Effects Correction with Unknown Subtypes for scRNA seq Data (BUSseq)---to correct batch effects in the presence of unknown cell types. BUSseq is able to simultaneously correct batch effects, clusters cell types, and takes care of the count data nature, the overdispersion, the dropout events, and the cell-specific sequencing depth of scRNA-seq data. After correcting the batch effects with BUSseq, the corrected value can be used for downstream analysis as if all cells were sequenced in a single batch. BUSseq can integrate read count matrices obtained from different scRNA-seq platforms and allow cell types to be measured in some but not all of the batches as long as the experimental design fulfills the conditions.

# Implementation

1. Installing RngStreams: Please download [RngStreams](http://statmath.wu.ac.at/software/RngStreams/) and install it by
```
wget http://statmath.wu.ac.at/software/RngStreams/rngstreams-1.0.1.tar.gz
tar zxvf rngstreams-1.0.1.tar.gz
cd rngstreams-1.0.1
./configure --prefix=<prefix_path> #prefix_path is your home directory/any directory you want to install at
make
make install
```
2. Compiling BUSseq: Please download all source code in *BUSseq-1.0* into your working directory and compile by
```
make
```
We use [OMPRNG](https://homepage.divms.uiowa.edu/~mbognar/omprng/) to generate random number in parallel.

3. Prepare your data files: Here we use the hematopoietic study as an example. There are N = 4,649 cells and G = 2,470 genes profiled in 2 batches. The first batch contains 2729
cells, while the second batch consists of 1920 cells. Then, the G by N read count matrix should be put in the `count_data_hemat_v1.txt` file, and `dim_hemat_v1.txt` contains the dimension information as the following:
```
4649
2470
2
2729
1920
```
For generate cases, the count data matrix and the dimension files should be named as `count_data_[your study name]_v[version].txt` and `dim_[your study name]_v[version].txt`.

4. Open the terminal, change the current diretory to your working directory and input the following two commands for MCMC sampling and posterior inference, respectivly.
```
cd /your/working/directory/
/your/BUSseq/source/code/directory/BUSseq -d./ -r/your/count/data/directory/ -p [your study name] -v [version] -K [the number of cell types] -i [the number of iterations] -o [the number of iterations for each output] -s [seed] -c [the number of cores for parallel]
/your/BUSseq/source/code/directory/BUSseq_inference -d./ -r/your/count/data/directory/ -p [your study name] -v [version] -K [the number of cell types] -i [the number of iterations] -b [the number of burn-ins] -c [the number of cores for parallel]
```
where 
   - *-d* for the directory to save the posterior sampling and inference result;
   - *-r* for the directory of read count data and dimension information created in the step 3;
   - *-p* for the project name, which shoud be consistent with the name of read count data and dimension file;
   - *-v* for the version number;
   - *-K* for the total number of cell types;
   - *-i* for the number of iterations in the MCMC sampling;
   - *-o* for the number of iterations printed into the hard disk to control memory usage;
   - *-s* for the seed of MCMC sampling to let the results be reproducible;
   - *-c* for the number of cores for parallel computing;
   - *-b* for the number of burn-ins in the posterior inference.

# Results
There are two folders created to store the results of BUSseq. `MCMC_sampling_K[the number of cell types]` stores the posterior sampling of MCMC algorithm, and `Inference_K[the number of cell types]` saves the estimation of all parameters. 

# Remarks
Please find the more details of how to running *BUSseq* on the real studies in the [BUSseq-1.0_implementation](https://github.com/songfd2018/BUSseq-1.0_implementation) folder.








