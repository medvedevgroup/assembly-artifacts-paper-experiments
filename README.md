# Assembly Artifacts Analyses


## Link to code and scripts

 Tool    | Version / Git commit | Link / Source |
| ------ | -------------------- | ---- |
| Modified-directed-bcalm (for directed unitig generation)  |   [e4d67f](https://github.com/medvedevgroup/directed-unitigs/tree/e4d67f87e5a6f774f30427535a5eb1012f8b64bb)         |   [GitHub link](https://github.com/medvedevgroup/directed-unitigs)  |
| Custom python simulator (simu_rand.py) |        -            |[GitHub link](https://github.com/medvedevgroup/assembly-artifacts-paper-experiments)     |
|All other code and scripts| - | [GitHub link](https://github.com/medvedevgroup/assembly-artifacts-paper-experiments) |  

## Other softwares used 




 Tool   | Purpose | Version / Git commit | Link / Source |
| ------ | -------------------- | ---- | ---- |
| bcalm2 |   Bidirected unitig generation |   [b8cde9](https://github.com/GATB/bcalm/tree/b8cde9cd4c5c17303b792e4017ce0728c6688666)                |   [GitHub link](https://github.com/GATB/bcalm)   |
| Jellyfish 2       |   Canonical and non-canonical k-mer counting| [029416](https://github.com/gmarcais/Jellyfish/tree/029416926c5c9c80ae3f62d94497af8166012855)         |   [GitHub link](https://github.com/gmarcais/Jellyfish/)   |
| ABySS  | Assembler | 2.0.2                | Conda: `conda install -c bioconda abyss`         |
| SPAdes | Assembler | 3.15.3                    |  [GitHub link](https://github.com/ablab/spades/releases/tag/v3.15.3)             |
| MEGAHIT       |   Assembler |  1.2.9                  |  [Binary link](https://github.com/voutcn/megahit/releases/tag/v1.2.9)             |
| GATB-Minia  | Assembler | [831ba4](https://github.com/GATB/gatb-minia-pipeline/tree/831ba4e8812aaed9f4ae1f00f9d8bf77439ba13f)                 | [GitHub link](https://github.com/GATB/gatb-minia-pipeline)          |
| Bowtie 2 | Alignment | [c04cfb](https://github.com/BenLangmead/bowtie2/tree/c04cfb8f027b1bc98a8800044dfb71a360a73a96) | [Github link](https://github.com/BenLangmead/bowtie2)|
| ART       |   Read simulator | 2.6.0         |   [Link](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)   |
| CAMISIM       |   Metagenomic simulation | [001a33](https://github.com/CAMI-challenge/CAMISIM/tree/001a3397a6a7c1ba4e0424988e9cf63067a34752)         |   [GitHub link](https://github.com/CAMI-challenge/CAMISIM/)   |



# Instructions to reproduce all results in paper

Reproducing most of the results in the paper is pretty straightforward. The readers should be able to reproduce the results following the descriptions provided in the paper with their own scripts. Still, here we provide the scripts we used along with instructions. The scripts and instructions are included to serve as guidelines only, and by no means are complete pipelines to reproduce the results. Some scripts might need minor changes, i.e. you may need to put correct path of certain files. 


## Reproducing results in Table 1 

The presence of unsafe and mis-assembled unitigs in human chromosome 1 is evaluated using simulated error-free reads. For the scripts, check `scripts/exp1` folder.

### Dataset

| Dataset | Download | 
| -------- | -------- | 
| T2T v1 chr 1    | [link](https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.2)     | 


### Method

From reference fasta file, remove all N and uppercase and make a file called "ref.sset" containing a single string of alphabet A, C, G, T only.
`cat chr1.fa | sed '/^>/d' > ref.sset`


The file `pipe_invalid_d_read.sh` simulates read at a certain coverage C (simple random model), then outputs sequenced segments, reads, unitigs (mis-assembled, unsafe and all others).    
`pipe_invalid_d_read.sh [folder containing ref.sset] [k-mer-size] [coverage]`


If the reference file `ref.sset` is in the folder `t2t_chr1`, then for coverage 20x, k=71, generate the unitigs using following command:
`pipe_invalid_d_read.sh t2t_chr1/ 71 20`

This outputs the following 4 important files.
The file "unitigs.sset" contains all unitigs, one unitig per line.
The file "segfile.sset" contains sequenced segments, one segment per line.
The file "invalids_true" contains mis-assembled unitigs, one unitig per line.
The file "invalids" contains unsafe unitigs, one unitig per line. 


To obtain Figure 3 cases, 

* First get a file named "gtouch" containing all unitigs and touching sequenced segments or reads.

Syntax: 
`bash gtouch_maker.sh <k-mer-size> <reference-string> <output-directory> <unitigs:one-unitig-per-line> <reads/segments:one-read-per-line>`

Example:
`bash gtouch_maker.sh 71 t2t_chr1/ref.sset out/ t2t_chr1/unitigs.sset t2t_chr1/segfile.sset`

An example to show the contents of "gtouch" file (generated with k=3): 

```bash
AAATTTA TTTTTTTTAAATTTA
AAATTTA GGGGAAATTTA
GGG GGGGAAAAAAA
CCAAA CCA
```

Reverse complements are not considered. All k-mers in unitigs and all segment/read are substrings of reference.
 
Contains 2 columns separated by " " (space).   
Column 1: unitig (directed).   
Column 2: a substring (i.e. segment/simulated read) of reference genome that contains at least one k-mer of the unitig (column 1).   
 
Grouped by column 1: That means if you want all the substrings touching unitig AAAAAAA, then start with first occurrence of AAAAAAA in column 1, and then scan subsequent lines untill you get a different unitig (i.e. GGG in the example). That gives you list of all substrings touching at least a k-mer from the unitig.


* Now, execute the command: 
`python pretty_printer.py 71 gtouch invalids_true ref.sset`

The output will indicate how many Figure 3 cases (labelled as "number of misjoins") are there. 

## Reproducing Table 3 Results

Compile `kmer_operations.cpp` : `g++ kmer_operations.cpp -o kmer_operations`
Then, execute:
`kmer_operations <k> <cov> <invalid_spades> <invalid_bcalm> <con/uni>`
<k> denotes k-mer size, <cov> denotes coverage, <invalid_spades> denotes path to misassembled unitigs/contigs generated by spades (or, other assembler), <invalid_bcalm> denotes path to misassembled unitigs predicted by our theorem. <con/uni> denotes label of output, use it to distinguish between contig and unitig. 
The output file includes number of k-mers shared between each pair of misassembled unitigs (or, contigs) generated by assembler and misassembled unitigs predicted by our theorem. Filter out the rows with numbers less than threshold using awk to generate the numbers in table 3 and 4.


## Reproducing Table 5 Results

Use bcalm2 to generate bidirected unitigs, and modified-directed-bcalm to generated double directed unitigs. Then, generate the table following the definitions in the paper. 

## Reproducing Table 6 Results

#### Dataset Generation

We use chr4 from grch38.

| Dataset | Download | 
| -------- | -------- | 
| grch38   | [link](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/)     | 


To generate synthetic dataset from chr4: 
`python artificial_palin.py ref.fa`


This generates ~7 million bp long contig from chr4 where around 700K locations were picked randomly and palindromes were introduced. 

#### 

Generate reads of length 100 with 1x coverage, and then build unitigs. For k=31,   
`art_sp_invalid.sh t2t_chr1 31 1 100`

## CAMI simulated dataset generation and reproducing Tables 2 and 4

Dataset: Download the 30 reference genomes by using the accessions listed in supplementary of the paper. 

1. Follow the instructions to download and install [CAMISIM](https://github.com/CAMI-challenge/CAMISIM). 

2. CAMISIM requires a configuration file (the one we used is located in cami_config directory). Now, execute the simulation
`python metagenomesimulation.py spades_config.ini`. This results in paired-end reads from 30 reference genomes listed in CAMI low complexity dataset in the supplementary of metaSPAdes paper. 

3. Concatenate all paired-end reads and convert to a fasta file named `read.fa`.
```sh
cat *.fq > read.fq
fq_to_fa.sh read.fq > read.fa
```
You can get number of total reads using this command: `echo "$(cat read.fa | wc -l) / 2" | bc`

4. Now, generate double directed unitigs using the following command. For k=55,
`dd_unitigs.sh read.fa 55`. This generates unitigs in a file named `unid55.fa`.

4. Concatenate all references to `ref.fa`, then build bowtie2 index from it
`bowtie2-build -p -f ref.fa refindex`

5. Align unitig to reference. We use options to align exactly (no errors or mismatch tolerated). The ones that are unaligned are `misassembled`
    `bowtie2 -f -x refindex  -U unid55.fa -S unid55_ref.sam --un misassembled55.fa --score-min 'C,0,-1' -p 8`

6. Align reads to references and save the alignment sam file (`read_ref.sam`) Extract only the relevant columns from sam file using: `cat read_to_ref.sam | grep -v @ | cut -f3,4,10 > read_to_ref.awkresult; cp read_to_ref.awkresult read_to_ref.samcut; cat read_to_ref.samcut | awk '{print $0 " " 150}' > read_to_ref.awkresult`. Here 150 is the read length. Now, generate `segments` from the reads for k=55 using:
`python segment.py read_to_ref.awkresult ref.fa 55 segs.sset`. This outputs a file called `segs.sset` which contains one segment per line.

7. Convert the file containing one segment per line to a fasta file. Then build a bowtie2 index from it.
`sset_to_fa.sh segs.sset > segs.fa`
`bowtie2-build -p -f segs.fa segindex`

8. Align unitig to segments. We use options to align exactly (no errors or mismatch tolerated). The ones that are unaligned are `unsafe`
    `bowtie2 -f -x segindex  -U unid55.fa -S unid55_seg.sam --un unsafe55.fa --score-min 'C,0,-1' -p 8`

Now, use instructions to reproduce Table 1 and 3, and apply it to cami dataset to generate data for  Table 2 and 4 respectively.







