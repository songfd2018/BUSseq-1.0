# Batch Effects Correction With Unknown Subtypes for scRNA-seq Data (BUSseq)

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Implementation](#implementation)
- [Remarks](#remarks)

# Overview

Single-cell RNA-sequencing (scRNA-seq) technologies enable the measurement of the transcriptome of individual cells, which provides unprecedented opportunities to discover cell types and understand cellular heterogeneity. Despite their widespread applications, single-cell RNA-sequencing (scRNA-seq) experiments are still plagued by batch effects and dropout events.

One of the major tasks of scRNA-seq experiments is to identify cell types for a population of cells. Therefore, the cell type of each individual cell is always unknown and is the target of inference. However, most existing methods for batch effects correction, such as Combat and the surrogate variable analysis (SVA), are designed for bulk experiments and require knowledge of the subtype information, which corresponds to cell type information for scRNA-seq data, of each sample a priori.
  
Here, our proposed `BUSseq` method fits an interpretable Bayesian hierarchical model---the Batch Effects Correction with Unknown Subtypes for scRNA seq Data (BUSseq)---to correct batch effects in the presence of unknown cell types. BUSseq is able to simultaneously correct batch effects, clusters cell types, and takes care of the count data nature, the overdispersion, the dropout events, and the cell-specific sequencing depth of scRNA-seq data. After correcting the batch effects with BUSseq, the corrected value can be used for downstream analysis as if all cells were sequenced in a single batch. BUSseq can integrate read count matrices obtained from different scRNA-seq platforms and allow cell types to be measured in some but not all of the batches as long as the experimental design fulfills the conditions.

# System Requirements

## OS Requirements

The C++ source codes support *Linux* operating systems. It has been tested on the following systems:

Linux: Ubuntu 18.04 and CentOS 7

## Hardware Requirements

The C++ source codes can work on a standard personal computer (PC). The runtimes reported below were generated on 4 Intel i7-6700HQ CPU @ 2.60GHz Cores.

# Installation

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


# Implementation

Here we introduce the workflow of BUSseq on a small demo. We prepare the input files and draw figures in *R*, so please find the Rscript in the [demo](demo) folder.
 
1. Prepare the input files: we generate a synthetic dataset in *R*, please run
```
cd /your/working/directory/demo/
R --vanilla --slave < simulate_data.R
```
In this demo, there are N = 450 cells and G = 1,000 genes profiled in 3 batches, each of which has 150 cells. Then, the G by N read count matrix should be put in the `count_data_demo_v1.txt` file, and `dim_demo_v1.txt` contains the dimension information as the following:
```
450
1000
3
150
150
150
```
As the same time, `metadata_demo_v1.txt` stores the metadata of the demo study, including the batch and cell type information of each cell. For general cases, the count data matrix, dimension and metadata files should be named as `count_data_[your study name]_v[version].txt`, `dim_[your study name]_v[version].txt` and `metadata_[your study name]_v[version].txt`, respectively.

In this demo, I also draw the t-SNE plot of the synthetic count data colored by batch labels and cell type labels.
![raw_colored_by_batch](/demo/Image/tsne_simulation_uncorrected_by_batch.jpeg =250x) 
![raw_colored_by_celltype](/demo/Image/tsne_simulation_uncorrected_by_celltype.jpeg =250x) 
![legend](/demo/Image/legend.pdf) 

2. Run MCMC sampling: 
```
R --vanilla --slave < run_MCMC.R
```
When we do not the number of true cell type number *K*, we run MCMC sampling for different *K*s and select the optimal *K* to obtain the minimum BIC. Thus, we generate the scatter plot of the number of cell type *K* versus BIC values in the [Image](demo/Image) folder. 

For general cases, change the working directory and run the following two commands in the terminal to conduct MCMC sampling and posterior inference directly:
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
   - *-b* for the number of burn-in iterations in the posterior inference.

After running MCMC algorithm, there are two folders created to store the results. `MCMC_sampling_K[the number of cell types]` stores the posterior sampling of MCMC algorithm, and `Inference_K[the number of cell types]` saves the posterior inference of all parameters. 

3. Correct batch effects by running:
```
R --vanilla --slave < correct_batch_effects.R
```
In the posterior inference, we correct the raw read count data as a version of corrected count data stored in the `x_corrected.txt` file by the quantile matching approach. At the same time, we draw the t-SNE plot of the corrected count data colored by batch labels and cell type labels in the [Image](demo/Image) folder. As a result, we can find that these cells are clustered by their cell types after corrected.


# Remarks
Please find the more details of how to running *BUSseq* on a simulation dataset and two case studies in the [BUSseq-1.0_implementation](https://github.com/songfd2018/BUSseq-1.0_implementation) directory.








