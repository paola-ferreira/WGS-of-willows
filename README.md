# Whole-genome sequencing data of Willows
This is a step-by-step tutorial how it was analyzed the WGS data of willows from Illumina sequencing. Please notice that all analyses must be run in a cluster due to the amount of data generated.
This tutorial is supposed to be friendly-user but I will be delight to answer any questions.

### #Before you start
We used our Whole-genome sequencing data in order to propose a evolutionary hypothesis for the enigmatic genus Salix. To do this, we used the alignments from Wagner et al. 2020 as targets in our dataset. The wagner's dataset comprised 23,393 loci for 133 taxa.

1) Download the DNA matrices and import in Geneious
2) Select all alignments and create a consensus for every aligment by clicking on 
``` Tools -> Generate Sequencing Consensus ```
4) Export all the consensus sequences (fasta format) to a folder in your computer
5) Open the folder in a terminal and join all the consensus sequences into a single file 
``` cat * > reference.fasta ```

The file generate here is will be used as the target.

#### 1) Importing data
Go to the folder comprising all the raw data reads generated from the last step and upload all of them to the metacentrum cluster:

``` 
scp * paolaferreira@nympha.zcu.cz:/storage/plzen1/home/paolaferreira/1.raw_data
```

You should also upload the fasta file reference in a new folder:

``` 
scp reference.fasta paolaferreira@nympha.zcu.cz:/storage/plzen1/home/paolaferreira/0.reference
```

#### 2) Check how many reads we have in every file. 
Note our willow data is paired end, which means we sequenced sequence both ends of a fragment and generate high-quality, alignable sequence data. Therefore for every sample we have SampleA_R1.fq.gz SampleA_R2.fq.gz and they must have the same amount of reads. (Read more about Illumina Paired-End Sequencing at https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html)

``` for i in *.fq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done```

The results should look like this:
``` 
ACU2_1.fq.gz
48960515
ACU2_2.fq.gz
48960515
```

Lastly, save all the information of the number of reads because this is usually reported in the supplementary materials of any NGS publication.

#### 3) Merge raw reads data
In order to increase the depth coverage of your samples, we have sequenced S. glaucosericea (GSR7), S. lapponum (LAP1), S. mielichhoferi (MIE5), S. myrsinifolia (MYS5) twice. Before proceeding with the further analyses, let's merge the sequences. For example, in S. myrsinifolia you will see four reads (two in every sequencing)
``` 
MYS5_1.fq.gz
46385628
MYS5_2_1.fq.gz
39397922
MYS5_2_2.fq.gz
39397922
MYS5_2.fq.gz
46385628
``` 
Merge the files MYS5_1.fq.gz and MYS5_2_1.fq.gz which represents both R1 sequencing:

``` cat  MYS5_1.fq.gz MYS5_2_1.fq.gz > MYS5_R1.fq.gz ``` 

Do the same with the R2

``` cat  MYS5_2.fq.gz MYS5_2_2.fq.gz > MYS5_R2.fq.gz ``` 

Check if the merging was performed correctly by checking the number of reads. 

``` for i in MYS5*.fq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done```

The results should be the sum of both reads:

``` 
MYS5_R1.fq.gz
85783550
MYS5_R2.fq.gz
85783550
``` 

### #Installation and Setup

For this project we need to install the pipeline (SECAPR; Andermann et al. 2018) and an additional program for trimming the alignments (Gblocks; Castresana, 2002) 

#### 1. Installing SECAPR pipeline
In our willow's project we used the pipeline entitled Sequence Capture Processor (hereafter SECAPR; Andermann et al. 2018). Please notice that the installation and setup here is just a copy of the original SECAPR repository on github. The complete documentation and installation files can be found at: http://antonellilab.github.io/seqcap_processor/ 

#### SECAPR data analyses overview
SECAPR is a pipeline written in phyton by Tobias Andermann. It comprises a number of steps tha can drive you from the raw sequencing data to multiple sequencing alignment. It also has a number of data visualization codes where you can decided if you need to re-adjust your parameters or go forward.

![image](https://user-images.githubusercontent.com/88035938/130038396-5e2cceed-edce-4e03-a400-4735c33c8795.png)

SECAPR analytical workflow. The flowchart shows the basic SECAPR functions, which are separated into two separate steps (colored boxes). Blue box (1. reference library from raw data): in this step the raw reads are cleaned and assembled into contigs (de novo assembly); Orange box (2. reference based assembly with custom reference library): the contigs from the previous step are used for reference-based assembly, enabling allele phasing and additional quality control options, e.g. concerning read-coverage. Black boxes show SECAPR commands and white boxes represent the input and output data of the respective function. Boxes marked in grey represent multiple sequence alignments (MSAs) generated with SECAPR, which can be used for phylogenetic inference. **Extracted from SECAPR Github.


#### 1.1. Install conda

Download the Python2.7 version of Miniconda and install it by executing the downloaded sh-file (see commands below).

Download conda (MacOS 64bit):
``` 
wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh 
``` 

Download conda (Linux 64bit):
``` 
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
``` 

Download conda (Linux 32bit):
``` 
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86.sh
``` 

Install conda:
``` 
sh Miniconda2-latest-*.sh
``` 

Activate conda (on metacentrum cluster):
``` 
module load conda-modules-py37
``` 

Add Bioconda channels (containing bioinformatics software):
``` 
conda config --add channels defaults; conda config --add channels conda-forge; conda config --add channels bioconda; conda config --add channels https://conda.anaconda.org/faircloth-lab 
``` 

#### 1.2. Install the SECAPR environment

Conda automatically downloads and installs all necessary software dependencies. We strongly recommend to install SECAPR and all it's dependencies in a separate virtual environment, in order to not interfer with potentially already installed verisons of the software dependencies.

Install SECAPR in virtual environment (here named secapr_env):
``` 
wget https://raw.githubusercontent.com/AntonelliLab/seqcap_processor/master/recipe/install_secapr_env.sh
sh install_secapr_env.sh
``` 


#### 1.3. Activate the environment

To activate the newly created environment, type:

Activate environment:
``` 
source activate secapr_env
``` 

Check the installed version:
``` 
secapr_env --version
``` 

When the environment is activated, all the necessary software dependencies will be available in the standarad path, e.g. when you type samtools the samtools version required for SECAPR will be executed. After you are done using secapr, you can deactivate the environment to switch back to your standard environment with this command:

``` 
source deactivate
``` 

#### 1.4. Check active environment

Check if you are connected to the correct environment (there should eb a star in front of secapr_env in the output of this command):

Check active environment:
``` 
conda info --envs 
``` 

#### 1.5. Install SECAPR development version
For the purpose of this project, we are using the SECAPR version 2.0.2. Therefore we need to install the version by following the instructions below:

Connect to your secapr environment 

``` 
source activate secapr_env
```

Remove the current secapr installation 

``` 
pip uninstall secapr
```

Download the version 2.0.2 (Source code.zip) from github 

``` 
wget https://github.com/AntonelliLab/seqcap_processor/releases/tag/v2.0.2
```

Unzip the downloaded file 

``` 
unzip seqcap_processor-2.0.2.zip
``` 

Move the unzipped directory to a safe location on your computer (in our case to the metacentrum cluster), i.e. not on your Desktop or Download folder, since this will be the path where secapr will be executed from in the future

Enter the unzipped secapr directory 

``` 
cd seqcap_processor-2.0.2
``` 

Install secapr from the folder 

``` 
python -m pip install -e .
``` 

Check the version installed 

``` 
secapr --version
``` 

#### 2. Install Gblocks
Gblocks is a software that eliminates poorly aligned positions and divergent regions of an alignment of DNA or protein sequences (Castresana, 2000; Talavera & Castresana, 2007).

For this project we are using Gblocks version 0.91. You can download this version through conda.

First activate conda:

``` 
conda install -c bioconda gblocks
``` 

Try to see if the Gblock is working by typing 

``` 
Gblocks
``` 

Alternatively you can download the software through its own website at http://molevol.cmima.csic.es/castresana/Gblocks.html


### #Data analyses




#### 1) Rename and unzip your files
Since we are going to use the SECAPR pipeline, we need to have our files in a way that the pipeline will be able to recognize it.
First unzip our .fastqc files using the gunzip program:

``` 
gunzip *
``` 

Second rename your files: A simple sample ID is enough followed by _R1 for the forward reads and _R2 for the backward reads


#### 2) Quality Check your raw reads
In our case, the sequencing company was responsible for removing the adapters and do a quality control of the reads. However, most of the time you will receive the raw data. To convince yourself that the data is in a good quality, let's double check our reads:

``` 
secapr quality_check --input /storage/plzen1/home/paolaferreira/1.raw_data --output /storage/plzen1/home/paolaferreira/2.Checking_quality_raw_data
``` 

- SECAPR will produce two plots using the R-script to show summary statistics for each individual test. The test names carry 3-letter acronyms, and the corresponding full test-name can be found by opening one of the html files. The first plot shows how many occurrences of each test-result (fail,pass,warn) were found for each test among all samples (per-test basis). The second plot shows for each sample (y-axis) which test had which result (per-sample basis).
- If you are satisfyed with the quality of your reads, you should go for the next step. 

#### 3) Assembling contigs
After checking the quality (and/or cleaning and trimming the reads) we are now ready to use the fastq-reads for de-novo contig assembly. In this step the overlap between fastq reads is being used to build long, uninterrupted sequences, with no a priori knowledge of the correct sequence or order of those fragments. We will use the contigs to hopefully find the target regions that were selected for during sequence capture.
For our willows project we used the development version of SECAPR which includes the ABYSS assembler. This is a much faster program compared to ABySS and Trinity (mostly used for transcriptome analysis) and also allows you to test several k-mer sizes at once. But first create a folder for each sample follow by ____clean___  and transfer both fastaq to there. This is necessary since we skipped the cleaning step, otherwise SECAPR would create a folder like that to you.
We ran our assemble using: 

``` 
secapr assemble_reads --input /storage/plzen1/home/paolaferreira/1.raw_data  --output /storage/plzen1/home/paolaferreira/3.assembling --assembler spades --kmer 21,33,55,77,99,127 --contig_length 50 --cores 16
``` 

#### 4) Mapping
Mapping is the process to align contigs (or reads) obtained by through high-throughput genome sequencing to a reference (genome, genes, cds, etc). To do this we will be using the consensus sequence from the multiple sequence alignment uploaded in the __before you start section.  

![image](https://user-images.githubusercontent.com/88035938/130071327-a8ebf5eb-8239-4309-b9cd-6277f53caf98.png)

Figure 1: Illustration of the mapping process. The input consists of a set of reads and a reference genome. In the middle, it gives the results of mapping: the locations of the reads on the reference genome. The first read is aligned at position 100 and the alignment has two mismatches. The second read is aligned at position 114. It is a local alignment with clippings on the left and right. The third read is aligned at position 123. It consists of a 2-base insertion and a 1-base deletion. Illustration and Figure legend extracted from https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/mapping/tutorial.html

We mapped our contigs using the following commands:

``` 
secapr find_target_contigs --contigs /storage/plzen1/home/paolaferreira/3.assembling  --output /storage/plzen1/home/paolaferreira/4.mapping --reference /storage/plzen1/home/paolaferreira/0.reference 
``` 

#### 4) Alignment
Now we have our selected contigs which were mapped against the RADseq references and we need to perform a Multiple Sequence Alignments (MSAs) using mafft or muscle. This function builds a separate alignment for each locus with matching contigs for ≥3 samples.

``` 
secapr align_sequences --sequences /storage/plzen1/home/paolaferreira/4.mapping/combined_joined_unphased_fastas.fasta --output /storage/plzen1/home/paolaferreira/4.mapping/5.Alignment --aligner mafft --output-format fasta --no-trim --cores 12
``` 

#### 5) Trimming
Once we get our multiple sequence alignments usually we need to remove poor aligned and divergent regions that may not be homologous or satured by multiple substitutions.

Go to the folder where your aligments are and create a loop in bash using default parameters:

``` 
for i in *fasta; do Gblocks ${i} -t=y -p=y; done
``` 

This step is going to create two outputs:
your_file_alignment.fasta-gb - This file is our trimmed multiple sequence alignment 
your_file_alignment.fasta-gb.htm - This file comprises the trimming statistics for each alignment.

Let's create a new folder and move the cleaned alignments:

``` 
paolaferreira@skirit:/storage/plzen1/home/paolaferreira$ mkdir 6.cleaned_alignments
mv *gb ../6.cleaned_alignments
``` 

#### 7) Mapping
Usually an analysis of most NGS datasets comprises the previous steps that we have gone through (cleaning and trimming of the reads, de-novo assembly of contigs, mapping, multiple sequence alignments and trimming). Now we will use those mutiple sequence alignments in order to generate new reference libraries and assemble the clean reads:

``` 
secapr reference_assembly --reads /storage/plzen1/home/paolaferreira/1.raw_data --reference_type alignment-consensus --reference /storage/plzen1/home/paolaferreira/6.cleaned_alignments --output /storage/plzen1/home/paolaferreira/8.reference_based_mapping --cores 16 --min_coverage 4
``` 








#### References:
Andermann, T., Cano, Á., Zizka, A., Bacon, C., & Antonelli, A. (2018). SECAPR-a bioinformatics pipeline for the rapid and user-friendly processing of targeted enriched Illumina sequences, from raw reads to alignments. PeerJ, 6, e5175. https://doi.org/10.7717/peerj.5175

Castresana, J. (2000). Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Molecular Biology and Evolution 17, 540-552. https://doi.org/10.1093/oxfordjournals.molbev.a026334 

Talavera, G., and Castresana, J. (2007). Improvement of phylogenies after removing divergent and ambiguously aligned blocks from protein sequence alignments. Systematic Biology 56, 564-577. https://doi.org/10.1080/10635150701472164

Wagner, N. D., He, L., & Hörandl, E. (2020). Phylogenomic Relationships and Evolution of Polyploid Salix Species Revealed by RAD Sequencing Data. Frontiers in plant science, 11, 1077. https://doi.org/10.3389/fpls.2020.01077 

