# Simulation of Transcription-Coupled DNA Supercoiling in Bacteria

The aim of this package is to simulate the Transcription-Coupled DNA Supercoiling (TCDS) and its impact in gene expression in Bacteria, we used python programming language with the help of a Python-based scientific ecosystem for scientific computing. We developed this stochastic model based on [Meyer et al.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003785) article. In order to calculate the variation of supercoiling density and on experiments done in Vitro ([Chong et al.](http://www.cell.com/cell/abstract/S0092-8674(14)00739-9)) and in Vivo ([Moulin et al.](https://www.ncbi.nlm.nih.gov/pubmed/15659173)) to verify and refine the results. We improved the model later on by including the effect of the protein barriers such as Nucleoid Associated Proteins (NAPs) and Topoisomerases.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

The following library are required in order for the script to run successfully:

* [Numpy](http://www.numpy.org/) - Fundamental package for scientific computing with Python.
* [Pandas](https://pandas.pydata.org/) - A library providing high-performance, easy-to-use data structures and data analysis tools.
* [matplotlib](https://matplotlib.org/) - A plotting library for the Python programming language and its numerical mathematics extension NumPy.
* [dnaplotlib](https://github.com/VoigtLab/dnaplotlib) - A library that enables highly customizable visualization of individual genetic constructs and libraries of design variants.

### Installation

The installation is simple:

1. Download the package & open the terminal in the script main directory.
2. Install the package using `pip`

```
pip install .
```

## Inputs, Outputs & Parameters files :

### Input files :

#### The General Feature Format file (GFF):
Each line of the  [GFF](https://www.ensembl.org/info/website/upload/gff.html) file contains 9 columns of data : *Seqname*, *Source*, *Feature*, *Start*, *End*, *Score*, *Strand*, *Frame* and *Attribute*.

```
##gff-version 3
#!gff-spec-version 1.20
#!processor NCBI annotwriter
##sequence-region chrom1genes 1 14000
chrom5genes RefSeq  region  1 14000 . + . ID=id0;Name=chrom1genes
chrom1genes RefSeq  gene  1201  13200 . + . ID=gene0;Name=gene0
```

#### The Transcription Start Site file (TSS):
The TSS file contains 4 information :
* *TUindex* : The transcription unit index.
* *TUorient* : The TU orientation ("+" for forward or "-" for reverse)
* *TSS_pos* : The TSS position in the chromosome.
* *TSS_strength* : The promoter strength.

```
TUindex TUorient  TSS_pos TSS_strength
0 - 4250  0.001
1 + 5150  0.015
2 + 16150 0.02
3 - 22250 0.02
```

#### The Transcription Termination Site file (TSS):
The TTS file also provides the *TUindex*, *TUorient*, *TTS_pos* and :
* *TTS_proba_off* : The probability that the transcription will end at this TTS.

```
TUindex TUorient  TTS_pos TTS_proba_off
0 - 3150  1.
1 + 6250  .6
2 + 17250 .4
3 - 21150 1.
```

### Output files :

The structure of the output files :

```
output
├── all_res
│   ├── save_RNAPs_info.npz
│   ├── save_sigma_info.npz
│   └── save_tr_info.npz
├── resume_sim
│   ├── resume_sim_Barr.npz
│   ├── resume_sim_RNAPs.npz
│   └── resume_sim_tr.npz
├── save_tr_nbr.csv
└── save_tr_times.csv
```

  * *all_res* : The *all_res* directory contains the files in which the information are saved **at each step during the whole simulation**.
      * *save_RNAPs_info.npz* : We store in this file RNA Polymerase IDs which are bound to  the transcripts and their positions along the genome.
      * *save_sigma_info.npz* : This file contains information related to the supercoiling such as the size of each region (regions are delimited by the protein barriers) and the supercoiling density ($\sigma$) in each one of them and $\sigma$ on the whole genome too.
      * *save_tr_info.npz* :  This file contains information about transcripts such as the number of transcripts and the initiation rate ($k_{on}$)
  * *resume_sim* : The *resume_sim* directory contains the files in which the last information are saved **directly before the end of the simulation**, those information can be used to resume the simulation later on.
      * *resume_sim_Barr.npz* : This file contains information related to barriers such as the barriers positions and types, the regions size and the supercoiling density in each one of them and finally the remaining steps for the ARN Polymerase in order to reach the end of the transcript (this information will help us if we want to resume the simulation).
      * *resume_sim_RNAPs.npz* : This file contains the RNA polymerases IDs which are bound to the transcripts, their positions along the genome and the unhooked RNA polymerases IDs (those are free RNA Pol, they can be bond and start transcribing).
      * *resume_sim_tr.npz* :  In this file we save the number of transcripts and the initiation rate ($k_{on}$)
  * *save_tr_nbr.csv* : A human readable csv file in which the number of transcripts made after each transcription are stored.
  * *save_tr_times.csv* : This file contains the time of creation of each transcript (the moment when the ARN Polymarase finishes the transcription process).
 
 NOTE : `.npz` files can be opened only by using [`numpy.load`](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.load.html) function.

### Parameters file :

The user can edit the following parameters :

1. The path to the *GFF*, *TSS*, *TTS* and the *Protein barrier* file.
2. The simulation parameters such as : 
  * *SIGMA_0* : The basal supercoiling density, which can be either specified explecitly or calculated based on Topoisomerases concentration.
    * *DELTA_X* : The spatial discretization in bp.
    * *RNAPS_NB* : The number of RNA Polymerases.
    * *ITERATIONS_NB* : Simulation time (in seconds).
    * *OUTPUT_STEP* : Time interval at which a graphical and/or text output is given (in seconds).
    * *GYRASE_CONC* : Gyrase concentration (in micromolar)
    * *TOPO_CONC* : Topoisomerase I concentration (micromolar)

Example of the parameters file: 
```
[INPUTS]
# GFF file Path
GFF = input/chrom5genes.gff
# TSS file Path
TSS = input/TSS.dat
# TSS file Path
TTS = input/TTS.dat
# Protein barrier file Path
BARR_FIX = input/prot.dat

[GLOBAL]
# The supercoiling generated in each step
# (equivalent to J_0 in the previous model)
genSC_RNAP = 0.1
# the sequence-dependent threshold of
# supercoiling-induced destabilization
sigma_t = -0.04
# Constants (from Meyer article)
epsilon = 0.01
m = 2.20

[SIMULATION]
SIGMA_0 = 0
# Spatial discretization in bp
DELTA_X = 60
# Time step corresponding to 25~30nt/s
DELTA_T = 2
# Number of RNAP
RNAPS_NB = 5
# Number of iterations (s)
ITERATIONS_NB = 6000
# Time interval at which a graphical and/or text output is given
OUTPUT_STEP = 1
# Gyrase concentration (micromolar)
GYRASE_CONC = 0.1
# Topo I concentration (micromolar)
TOPO_CONC = 0.041
## Constant defining the efficiency of topoisomerases
# Gyrase efficiency
GYRASE_EFFICIENCY = 0.09
# TopoI efficiency
TOPO_EFFICIENCY = 0.31
# TopoI constant 0.0005
TOPO_CTE = 0.0005
# Gyrase constant
GYRASE_CTE = 0.001
k_GYRASE = 50
x0_GYRASE = 0.016
k_TOPO = 80
x0_TOPO = -0.04
```

## Code Example

You can test the example by placing the `example` directory (provided with the package) in your *current working directory* and type the following command :

```
# Execute the script by providing the parameter file and the input directory (the output directory path is optional)
# python start_simulation.py path/to/the/params.ini output
python start_simulation.py params.ini 
```

NOTE : The input directory should contain the *GFF*, *TSS*, *TTS*  and the *Protein barrier* file.

If the simulation ended successfully, you'll get the output files as described [above](#output-files).

You can use the script that read the npz files and show the simulation by using the `film_from_npz.py` script and specifying the parameters and the output files path from which the information will be read :
```
python film_from_npz.py params.ini output
```  

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* So many thanks goes to *MEYER Sam* for his guidance and suggestion during the development.
* Many thanks to the scientists that did experiments on which we were based to refine and calibrate the model.
* And finally to everyone who contributed directly or indirectly to the realization of this  project.

## Things to improve

* Improve the documentation
* Refine the model
* Make different config files for each type of the parameters (put the global unchanged parameters separately from the simulation parameters).
* ...

#### Contact
* [sam.meyer@insa-lyon.fr](sam.meyer@insa-lyon.fr)
* [bilal.elhoudaigui@gmail.com](bilal.elhoudaigui@gmail.com).
