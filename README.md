# Benchmarking of computational error-correction methods for next-generation sequencing data


[![Preprint Available](https://img.shields.io/badge/Preprint-online-green.svg)](https://doi.org/10.1101/642843) [![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)


This project contains the links to the datasets and the code that was used for our study : ["Benchmarking of computational error-correction methods for next-generation sequencing data"](https://doi.org/10.1101/642843)

**Table of contents**

* [How to cite this study](#how-to-cite-this-study)
* [Reproducinf results](#reproducing-results)
  * [Tools](#tools)
  * [Data](#data)
  * [Notebooks](#notebooks-and-figures)
* [License](#license)
* [Contact](#contact)


# How to cite this study

> Mitchell, Keith, et al. "Benchmarking of computational error-correction methods for next-generation sequencing data" bioRxiv, doi: https://doi.org/10.1101/642843


# Reproducing results

## Tools

We have evaluated 10 error correction tools: BFC, Bless, Coral, Fiona, Lighter, Musket, Pollux, Reckoner, Racer and SGA. Details about the tools and instructions for running can be found in our ["paper"](https://doi.org/10.1101/642843).

We have prepared ["wrappers"](https://github.com/Mangul-Lab-USC/benchmarking_error_correction/tree/master/scripts/wrappers) in order to run each of the respective tools as well as create standardized log files.


We have also prepared ["scripts"](https://github.com/Mangul-Lab-USC/benchmarking_error_correction/tree/master/scripts/evaluation) to perform the evaluation of the error correction methods.


## Data
[![Hosted on Figshare-DOI: 10.6084/m9.figshare.11776413](https://img.shields.io/badge/Hosted%20on%20Figshare-DOI:%2010.6084/m9.figshare.7738901-blue.svg)](https://dx.doi.org/10.6084/m9.figshare.11776413)

We have evaluated 5 datasets composed of raw reads and their respective true reads (gold standard).

* **D1 dataset**: D1 was produced by computational simulations using a customized version of the tool WgSim. We generated simulated data mimicking the WGS human data using a customized version of the tool WgSim. Read coverage varied between 1 and 32. The WgSim fork is available at https://github.com/mandricigor/wgsim.



* **D2 dataset**: Raw reads corresponding to 8 samples (SRR1543964, SRR1543965, SRR1543966, SRR1543967, SRR1543968, SRR1543969, SRR1543970, and SRR1543971) were downloaded from https://www.ncbi.nlm.nih.gov/. The error-free (true) reads for the D2 dataset were generated using a UMI-based high-fidelity sequencing protocol, also known as safe-SeqS.



* **D3 dataset**: We generated simulated data mimicking the TCR-Seq data using the T cell receptor alpha chain (TCRA). Samples have read lengths of 100bp and read coverage varied between 1 and 32.



* **D4 dataset**: D4 corresponds to HIV population sequencing of an infected patient. The error-free (true) reads for the D4 dataset were generated using a UMI-based high-fidelity sequencing protocol.



* **D5 dataset**: We prepared the viral dataset D5 using real sequencing data from NCBI with the accession number SRR961514. Each read was assigned to the reference with which it has a minimum number of mismatches. The original error rate in the dataset was 1.44%. We modified these reads as follows: first, we corrected the corresponding portion of errors with a corresponding reference nucleotides to obtain different levels of errors in the datasets (1.44%, 0.33%, 0.1%, 0.033%, 0.01% , 0.0033%, 0.001%, 0.00033%, 0.0001%); We also created datasets with mixtures of two haplotypes with the original 1.44% error rate but with different levels of diversity between haplotypes (Hamming distance=5.94%, 0.29%, 0.02%). We applied a haplotype-based error correction protocol to eliminate sequencing errors from the D5 dataset.


## Notebooks and Figures

We have prepared Jupyter Notebooks that utilize the raw data described above to reproduce the results and figures presented in our [manuscript](https://www.biorxiv.org/content/early/2018/10/25/452532).

* [D1 dataset (E.coli) Jupyter Notebook](https://github.com/Mangul-Lab-USC/benchmarking_error_correction/blob/master/notebooks/D1_WGS_E.coli.ipynb)
* [D1 dataset (Human) Jupyter Notebook](https://github.com/Mangul-Lab-USC/benchmarking_error_correction/blob/master/notebooks/D1_WGS_human.ipynb)
* [D2 dataset Jupyter Notebook](https://github.com/Mangul-Lab-USC/benchmarking_error_correction/blob/master/notebooks/D2_TCRA_real.ipynb)
* [D3 dataset Jupyter Notebook](https://github.com/Mangul-Lab-USC/benchmarking_error_correction/blob/master/notebooks/D3_TCRA_simulated.ipynb)
* [D4 dataset Jupyter Notebook](https://github.com/Mangul-Lab-USC/benchmarking_error_correction/blob/master/notebooks/D4_HIV.ipynb)
* [D5 dataset Jupyter Notebook](https://github.com/Mangul-Lab-USC/benchmarking_error_correction/blob/master/notebooks/D5_HIV.ipynb)


# License

This repository is under MIT license. For more information, please read our [LICENSE.md](./LICENSE.md) file.


# Contact

Please do not hesitate to contact us (mangul@usc.edu) if you have any comments, suggestions, or clarification requests regarding the study or if you would like to contribute to this resource.
