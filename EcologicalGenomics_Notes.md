# Ecological Genomics Notes



##  P. Alexander Burnham



# January 23, 2017: 

## Info Update: Mellisa

* Advantages in Seq tech. 

* Ranage of applications:

     - WGS (whole genome sequencing)
     - RNAseq (seqencing RNA conveted to cDNA)
     - Targeted caputure seq. (string of probes mixed with sample, pulls immune related genes from organisim and seq just those)
     - Chipseq (chromatin imunoprecipitation seq, recognizezes and antibody and puls out all DNA bound to that portein)

**Why one or the other?**

Genetic variation


      - phenotypes

​      
number of samples 


      - population
      - individual 
      - comparitive studies 
      - model or note


Demographic history

Adaptive geneitc varaitin

gene expression car.

length of reads

number of reads

distribution

Reads:
​      
​      
      - short = 50bp
      - long 100 bp, 150 bp, 300 bp (miseq) 
      - 10,000-60,000bp = SMRT 


Single vs. paired end

*    General library Prep. Workflow  
     ​    
     - extraction (DNA, RNA -> to cDNA)
     - fragment sample
     - ligate adaptors (indvidual barcodes)
     - add seq. adaptors

     *Reduced Rep*  

     - RNA -> coding
     - GBS/RAD-seq
     - near restriciton sites

*    Sequencing-by-synthesis (SBS)   

     - bridge amp
     - cluster gen.
     - labeled dNTP (ATCG)

*    Other Technologies   

*    Learning Activity   


**Human Genome Project (2001-2003)**


      - ABI = Sanger 
      - 15 years  
      - 1 geneome (one person)    
      - $3 billion    


Uses PCR and sequences broken by faulty base pairs to work backwards

**2014 X-Ten releases**


      - Hiseq by Illumina (look up video of how it works)
      - 1 day  
      - 45 whole genomes    
      - $1000 bucks each  


Sheet of glass with 8 lanes with flow cells...(look up)



**Take home messages**

* Likely using Illumina seq. (just library tech. changes) usually SBS

* adaptors are markers (barcodes) used to identify samples during sequencing 

     - first thing that's sample and gives ID barcode = seq adaptor 
     - alago is a sequence of DNA that is attatched to plate binds to sample

* Model vs. non Model:

     - short reads (assembly to create long sequnce based on short reads that shift)
     - denavo assembly -> computer program (added variability with mixed sample)
     - 15% error for SMRT can be reduced to less than 1 with repeated passes

     - illumina is much smaller but more accurate (0.05% error)

     - combine the two to have a higher degree of confidnece 


## Paper Discussion: 

### Three advances in biology:

* Modern Synthesis (evolution (Darwin) and population genetics (Mendel))

* Watson and Crick (molecular biology DNA)

* Omics Era (genomics, preteomics etc.)


### What do we think about this?

* Phylogenies -> reducing error bars or reshaping question?   

* Do you throw out experimetnal design and scientific method for large scale shot gun blast sequencing? - > can still be used to do hypothesis driven science, but also something to be said for sending a teleoscope into deep space.   

* Most journals and funding sources require all seq data be made public - leeds to a storage space issure -> genbank ran out of space and needed emergency funding from congress.   

* data are not reviewed very well some gen bank sequences are kinda rough   

* reference genomes are a sample size of 1 so it might not be representative

On Wednesday two papers to talk about:

1) where do genomics data have limitations and why we use them

2) discussion leaders (sign up via blackboard)

3) next week talk about some of library preps (four of them) four differnt update people





![Illumina Sequencing](https://upload.wikimedia.org/wikipedia/commons/thumb/6/65/Cluster_Generation.png/440px-Cluster_Generation.png)



# January 25, 2017: 

## Info Update: Steve

### Outline:

* What are QTN

* Quantitative Gentics
    - theory of adaptive traits
        - V~a~
        - h^2^
* Methods
    - linkage mapping
    - GWAS
    - selection scans

### Notes on topics in outline:

* QTN = "quantitative trait necleotides"

    - Flowering times (quantitative traits)
    - Flower color (more descrete - mendalian major effects traits)
    - THermal tolerance 
    - venom potancy 
    - defense compounds
    - toxin tolerance
    - drought tolerance
    - altitiude tolerance

**Looking at haldane and fishers work on Quant genetics:**   
    - addative effect at each locus and addative effect average our Additive genetic variance or V~a~   
    - the additive portion of the addedative variance or phenotypic variance equals heritiability $$(V_a/V_p = h^2)$$   
    - the differnce between pheotypic variance between loci is addative effect and average out to average affect which is usually a linear functiion   
    - most populations are at the value of trait where fitness is max (local addpatition)
    - mutation with small effect on trait (fisher thought most important mutation 50% of being beneficial)  
    - large effect mutation -> tend not to be selected -> to extreme and usually cause problems
    - best traits climb to the fitness peak through small incremental steps

## Three main methods:

 ### 1) QTL mapping (linkage mapping)
    - assumed that two parents are homozygous (for continuos trats) 
    - F~1~ is heterozygous but no recombinance
    - mate the F~1~ generation and get one generation of recombination for these traits
    - **The Idea of QTL** use markers linked to chromosomes, crossing (over many generations)


```
linkage
|           QTN
|             o
|         o       o
|       o           o
| 
|   o                 o
|o                        o
--------|-----------|-----------
0                             5
          Chromsome

```

### 2) GWAS:

$$y_{trait} = mu_{intercept} + beta_{effectsize} * SNP_i + covariates_{population size}$$

collect data run regressions (millions)

**Manhatten plot
```
log_p
|                       .
|                       .
|                       ..
|                     ....
|.....  ......  ..........
|  ...............  ........
--------------------------------
1    2      3      4     5     6

```
### 3) Selection SCANS (based ons selection sweeps):

```

x= SNPs
selction favors mutation next to O

pop1                                  pop1
________x_________                  _xO__________________

___x______x______x                  _xO__________________
                          time
_xO__________________     ------>   _xO__________________  

______________x____                 _xO__________________


pop2
________x_________ 

___x______x______
                      
_xO__________________ 

______________x____


- lose variation as everything favors the mutation
- diversity is higher in pop 2
- FST is populations struction high FST will be for pop2 low for pop 1 after selction

```
- we study the few large effect size mutations for QTL and QTN that we study because they are easy to find   
- majority are small but hard to find

## PAPER DISCUSSION:

forward vs reverse approach

apply scientific method and don't assume that evrything is representing the sample of arge mutations we have mostly been studying

1) large effect allele mutation (mendelian)

2) basic theory in terms of gen. theoary neither requires or predicts the importance of large mutations

3) making assumption that the little parts we find act the same way the big parts act the sam way and thats probably not true

They will argue that there will be knowledge to found using this information

##DISCUSION OF GROUP PROJECT
### Notes on Seastar Wasting Disease:

**BACKGROUND ON DISEASE**   

* High mortaility (between 70-100%) 
* many species impacted (not species specific) 
* loss of legs withing hours or days, turns inside out!
* turn into goo
* Alaska to Baha
* First report in 2012, less severe currently (reported in the 1970s and 80s in isolation)
* study (Hewson et al., 2014) linked Densovirus implicated but not very confidnet 
* potentially something in microbiom can lower host immunity and allow for opportunistic infection by pathogen
* [Seastar Disease Outbreak Distribution](http://data.piscoweb.org/marine1/seastardisease.html)

**INFO ON WHAT DATA WE HAVE**   

* focusing on *Pisaster ercraceous*
* RNA extractions and polly a tail and sequenced mRNA on 3 illumina hiseq lanes
* also pulled our ribosmal 16s structure
* none of samples pos for densavirus (qPCR) positive control
* DNA virus
* is it cyclical
    - 1) never been so extreme before
    - 2) do go through booms and busts though
* All arrived healthy
* put into aquaria with artifial seawater (kept at 12^o^C)
* From same site
* were no infected; some developed SWD on there own 
* first and only time course study of this kind


## Hypothesis ideas:

### My idea: (Envirnomental Componant)

* Intertidal VS Subtidal
    - genetic differences (local adapt.) related to suseptabilty    
    - gene expression (immune related genes higher in SUB group)   

### Other ideas:

* Resistance genes (stayed healthy or recovered)
* w/in Intertidal -> gen. differences between three groups
* both groups, Microbiome -> H vs S (differences (yes or no))
* Microbiome through time 
* h^2^ of microbiome
    - day as replicate   
    - relatedness VS microbiome trait
* immune related genes (Sam) 
    - what kind of pathogen based which genes are expressed   
    - which are the immune relatedness genes (conserved in other species)
    - potential problems -> viruses that are opportunistic and enter at tail end of infection and immune repsonse genes still expressed (confounded?)


# January 30, 2017: 

## Info Update: Group

**2 ideas i'm interested in**

* Intertidal VS Subtidal
    - genetic differences (local adapt.) related to suseptabilty    
    - gene expression (immune related genes higher in SUB group) 
    - same population but handeling was different (intertidal was more stressful)
    - gene expression differences
    - microbiome differences
* Immune-related gene expression
    - reverse pathology
    - looking for specific classes of genes
    - a-priori tests for resitance genes (compaire indv. that got sick to healthy)
    - looking at the S -> H transition 


**To Do:** Team Sherlock 

Chose reverse gene exression 

One page proposal for ideas for our project (1 per group) for next monday. Use handout as outline. Word and email. For wednesday, bioblitz (library prep types)



# February 1, 2017

## Outline Today:

* Announcements:
  * send linkts to github notebooks to Andrew
  * signups for discussion and info updates
  * project proposals due by email next monday
  * we start transcritomics next week 
* Info update blitz on library preps
* UNIX turorial 



### Whole Genome Sequencing 

* Applications

  * High Power and resolution
  * population genetics
  * local addpatatoins
  * plastic responses of the environment
  * inbreeding

* Methods

  * Somtimes you need a reference genome
    * No: De Novo Assembly (gene expression, Adaptations)
    * Yes: select important variation (DNA protein, epigenetic mod)
  * Need money (servers, large storage etc)
  * How to work with a server (command line)
  * Python or Perl

* Limitations

  * Polymorphic genes + core genes (highly conserved)
  * panalogs  
  * Rapidly evolving genes  (poor resolution)
  * large genes families (poor resolution)

* Usually using one indivudal (could also pool samples)

  * it is impossible to sequence everything (highly repetative regions and centro/telomeres etc. )

  * somtimes DNA is bound up by cromotin so enzymes can't get to it

  * genome will be a working hypothysis not a set entity 

    ​

* Seq. Platforms:

  * short reads 
    * illumina seq. 150bp 
    * SOLID 50bp
  * longer reads 
    * Pacific bioscience 5 kb
    * Ion torrent (500bp)
    * illumine moleculo (10kb)

* knowledge on organism

  * genome size is important (k-ner approach short unique elements of DNA seq. of length K)

* Wet lab procedures

  * high quality DNA (tissue)
  * no energetic tissue (muscle) contains too much mitchondrial DNA
  * no gut and skin: other DNA that isn't from organism
  * Quantity: 1mg -> 6ug (short)

* Library Prep.

  ->         ->     <-

  ​      <-    ->    <-      single end

  **->	       ->**

  **<-          ->**             paired end   (bold is paired)

  <—————————> mate pair

  **->	       ->**

  **<-          ->**             if gaps

  contig        contig

  <—— ……..——>

  Scaffold

* statistics: N50 = 50% of the assembly comes from contigs (# of base pairs in middle contig)

* Annotate your sequence (use a related genome) how related to simalar organism's genome

  * could be automated 
  * or manually
  * publish the genome 
  * find them on NCBI



### RNA SEQ.

* advantage

  * differential genome expression

    * between populations
    * treatments in experiment 

  * allele specific expression

    * envirnomental response
    * adaptation

  * functional relavent subset of them genome

  * RNA Seq. VS Micro-array:

    | RNA SEQ                              | Micro-array                              |
    | :----------------------------------- | :--------------------------------------- |
    | wide range of expression values      | x                                        |
    | x                                    | saturation of analog-type flourscent signal |
    | gives information on splicing events | x                                        |

    ​

* Limitation

  * can tell you about differential gene expression but not protein expression so difficult to say something about the function

* work flow

  * set-up
  * wet lab
  * seq. stratagy
  * bio-info
  * statistical measures





**Set Up**

* purpose:
  * coding or regualtory - non coding
  * reference genome?
  * alternative splicing?
  * technology
  * population or specific treatments?
* stat: biological replication
* choice of tissue
  * circaidan rythme 
  * embrilogical expression in animals 
  * pooling tissue with small tissue

**Wet Lab**

* RNA Extration
  * RNase free envirnment
  * DNase treament
  * Get rid of Ribosmal RNA and mRNA 
  * enrich mRNA with poly-A tails
* cDNA
  * reverse transcribe RNA to cDNA
* library
  * single end and paired end

**Seq. Strategy**

* platform
  * pyrosequencing by Roche
  * Ion torrent
  * GA/hiseq by illumina 
* Error profiles
* seq coverage:
  * greater than 1000000 bp
* programming:
  * unix
  * python 
  * R!!!

### Amplicon SEQ.

* Methods
  * library prep  ->  Extract -> PCR1 -> clean up product (gel) -> PCR2 (add barcodes + adaptors) -> SEQ
  * sequencing 
  * Data analysis 
    * learning activity 
  * Applications 
* Data Analysis 
  * trim adaptor 
  * allign sequences
  * compare genes with the same gene to look for mutations etc. 

### STEVE  GBS rad-seq



WGS					RNAseq			**GBS rad-seq.**		     amplicon

<————————————————————————————————————>   

everything			gene space					    				one gene      (comleteness of data)

​					expressed											

______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

single indv.			fewer indv.									many indv.     (sampling tradeoffs)															



GBS = genotyping by sequencing (genotype by sequencing at each loci)

RAD-seq = restriction assisted DNA sequencing (restriction enzymes - cuts double strand DNA - sticky ends add barcodes and addaptors like in all other library preps)

 

* lots of individuals

* lots of SNPs across the genome

* dont care about sepcific genes

* don't need complete genome

  ​


## Comand line



Server log in:

ssh pburnham@pbio381.uvm.edu

my personal net id password:



"top" shows how many cpus are being used and who is on the server

"q" = quit

everyone has a home directory ~/mainfolder 

"pwd" finds working directory 

"ll" = list in long format

make new folder called "mydata" (mkdir mydata)

cd "directory name" is change directory

shared space goes to most basel directory in the entire machine /data/

cp is copy

how to peak into big files without opening them "head" and -n "number" allows you to specify number of rows

tail is the last 10

">" new file.txt" makes new file

grep generalized regular expression (search and replace) grep "what im looking for" file.txt

mv is move

"anything after * is included" called a wild card

unix help manual "man" then command 

rm = remove 

-a shows hidden files

vim is built in text editor vim then file name

i is insert mode and allows for us to edit 

alias sets comman equal to another command (-i asks for conformation)

:w save and :q is quit



Used 38_6-21_H_0_R2.fq.gz as my file for class I missed 



```
38_6-21_H_0_R1.fq.gz
38_6-21_H_0_R2.fq.gz
```



## February 13, 2017

Lauren info update:

Transcriptanomics RNAseq experimental design 

OUTLINE:

* background
* issues
* R excersize 
* general rules of thumb 

**Background:**

* enables differential expression examination 
  * disease resistane
  * mating behavior
  * addaptive significane 
* Molecular mech. -> phenotypic/behavioral plasticity (migration patterns) 
* Some limitations:
  * reference genome quality
  * gene annotation availablity 
  * expensive for sample library prep.

**Issues:**

* underutilization of biological replicates (usually differnt indivuals instead of same sample which is a tech replicate)
  * reqiuire independant replications
  * doesn't include pooled samples
  * only 23/158 studies used greater than 3 biological replications 
  * most studies derive broad biological conclusions with little to no bio rep.
* studies usally prioritize seq. depth over replication 
  * coverage is...
  * depth is the amount of reads per sample
* wide dynamic range of RNA seq data makes it noise (unexplained variation in the data set)
  * poisson counting error
  * non poisson technical variance (storage, processing, differences in lanes etc.)
  * biological variance (natural variance due to difference in genetic or envirnmental variance)
* power prob of rejecting a FALSE hypothesis (want 80% power)

**Rules of Thumb:**

1) More bio reps instead of increasing depths

2) sequence depth > 10 reads pre transcript (sufficient)  

3) ~ 10-20 mil mapped reads per sample is suffient 

4) 3 biological replicates per condition 

5) conduct a pilot experiment 

- what is the best and most powerful experiment that I can afford
- what is the smallest fold change I can detect 



# February 15, 2017

## Sam Info Update:

**SNPs and population Genomics.**

SNPs (single nucleatide polymorphism - different nucs at same locus)

(Tissue -> sequence -> clean/trim -> assembly ->)(SNP detection/validation)

===> practical applications



1) Tissue: 

* breadth of tissue
* different stages (exon skipping)

2) Pool and sequence libraries 

* ~30-100 paired end long reads

3) Process raw sequence data 

* for errors
* important for SNP detection

4) digital normalization 

* remove high coverage reads and thus errors contained therein
* reduces sampling variation
* not quantatative info anymore

5) Assemble clean paired end reads

6) prune reads to reduce: 

* DNA contamination
* non coding RNA 
* gene fragments 

7) Assembly Evaluation 

* reference genome 
* cogs (conserved genes in other eukaryotic species)

**8) SNP Detection and Evaluation**

* Find your SNPs
* Software - constant patterns of sequence variations
* probability of sequence error - eliminate SNPs of low frequency 
* Artifacts caused by insertions and deletions (InDels) -> filter out SNP near indels
* use quality score

**Practical Application**

* differines in pop structure
* natural selection acting on a particular loci?

**Using outliers:**

* for a given locus, what's the level of differentiation compared to differances acrros genome (FST)
* FST -> look up
* Plotted FST in a a histogram is ~N()
  * outliers in this plot can be used to look at differtional selction acting on the populaton of loci

**Not Outlier approach:**

* takes high FST values and tests them for other features associated with seection 
  * how those features asssociated with inv. fitness
  * enrichment for certain functional roles



## Coding in Unix for today - working with my sam file I made last class:

#### Let’s check out our .sam files! Try `head` and `tail`.

```
tail -n 100 YOURFILENAME.sam > tail.sam
vim tail.sam

:set nowrap
```

#### Let’s see how many of our reads map uniquely.

Why is it important to consider whether a read maps uniquely (i.e., to one place in the transcriptome) for gene expression studies?

```
$ grep -c XT:A:U YOURFILENAME.sam 
1177827

$ grep -c X0:i:1 YOURFILENAME.sam
1182952
```

You can check a number of other elements, total number of reads, search for the various flags…

## Extract read counts from the .sam file from each sample

We will use a custom python script (by my friend Dan Barshis and published with the Simple Fool’s Guide to Population Genomics) called **countxpression.py**. This script will take any number of input *.sam files and, for each .sam file, extract the number of reads that map to each gene (i.e. the “counts”). It will also generate a summary output of useful information including proportion of quality read alignments. The script requires 4 input variables: mapqualitythreshold, lengththreshold, outputstatsfilename, anynumberofinputfiles.

```
cd /data/scripts
cp countxpression_pe.py ~/scripts      #or copy to your directory with the .sam file

python countxpression_pe.py 20 35 countstatssummary.txt YOURFILENAME.sam
```

This python script will generate two files: a .txt file you named (3rd argument you passed the script) and a counts .txt file that includes the number of uniquely mapped reads to each gene in our transcriptome.

Below are what the files should look like:

```
$ head NC_AD4_M3_bwaaln_counts.txt
ContigName  UniqueTotReads  MultiTotReads   totalreadsgoodmapped
OTAU000001-RA   11  207 218
OTAU000002-RA   982 49  1031
OTAU000003-RA   867 0   867
OTAU000004-RA   338 0   338
OTAU000005-RA   154 0   154
OTAU000006-RA   26  0   26
OTAU000007-RA   17  0   17
OTAU000008-RA   1017    55  1072
OTAU000009-RA   1984    0   1984
```

Once we have all the read counts extracted from each .sam file and in one directory, we can stitch them together with some bash scripting that I wrote.

```
# This loop takes the second column of data and renames the file to a shorter version of itself
for filename in *counts.txt; do
    myShort=`echo $filename | cut -c1-10` 
    echo "$myShort" > $myShort"_uniqmaps.txt"    
    cut -f 2 "$filename" > $myShort"_uniqmaps.txt"  
done 
# makes many individual files, but they don't have the header inserted

# This loop uses the tail command to get rid of the the first line
for filename in *_uniqmaps.txt; do
    tail -n +2 -- "$filename" > $filename"_uniqmapsNH.txt"  
done 

# This loop inserts the shortened version of the filename as the first line using the echo (print) and cat functions
for filename in *_uniqmapsNH.txt; do (myShort=`echo $filename | cut -c1-10`;echo "$myShort"; cat $filename) > tmp; mv tmp $filename; done

# This combines all the single column datafiles into one!
paste *_uniqmapsNH.txt > allcountsdata.txt

# Add row/gene names to table
cut -f 1 WA_PP1_M3__bwaaln_counts.txt | paste - allcountsdata.txt > allcountsdataRN.txt

vim allcountsdata.txt

# to view with tabs aligned.
:set nowrap  

# clean up files, get rid of intermediate files
rm *uniqmaps*

```

## Processes going on in the background

1. Making new assembly with more sequence data
2. Testing assembly quality several ways:

- stats (N50, number of contigs, etc.)
- proportion with high quality blastp match to (a) uniprot, (b) Patiria miniata, and (c) Strongylocentrotus purpuratus
- Proportion single-copy core eukaryotic genes represented.

1. Cleaning the second half of the samples/reads that just finished downloading 2/14!
2. Map all reads to new assembly
3. Extract read counts using custom python script.

Hopefully all of this will be done within a week by our next class session - next Wednesday

------

Clean up file names. Run script for the paired and unpaired files:

```
for file in *.fq.gz_*_clean_paired.fq ; do mv $file ${file//.fq.gz_*_clean_paired.fq/.cl.pd.fq} ; done

for file in *.fq.gz_*_clean_unpaired.fq ; do mv $file ${file//.fq.gz_*_clean_unpaired.fq/.cl.un.fq} ; done
```

```
cat 08_5-08_H_0_R1.cl.pd.fq 08_5-11_S_1_R1.cl.pd.fq 08_5-14_S_1_R1.cl.pd.fq 08_5-17_S_2_R1.cl.pd.fq 08_5-20_S_3_R1.cl.pd.fq 10_5-08_H_0_R1.cl.pd.fq 10_5-11_H_0_R1.cl.pd.fq 10_5-14_H_0_R1.cl.pd.fq 10_5-17_H_0_R1.cl.pd.fq 10_5-20_S_2_R1.cl.pd.fq 35_6-12_H_0_R1.cl.pd.fq 35_6-15_H_0_R1.cl.pd.fq 35_6-18_H_0_R1.cl.pd.fq 35_6-21_H_0_R1.cl.pd.fq 36_6-12_S_1_R1.cl.pd.fq 36_6-15_S_2_R1.cl.pd.fq 36_6-18_S_3_R1.cl.pd.fq > 08-11-35-36_R1.cl.pd.fq

cat 08_5-08_H_0_R2.cl.pd.fq 08_5-11_S_1_R2.cl.pd.fq 08_5-14_S_1_R2.cl.pd.fq 08_5-17_S_2_R2.cl.pd.fq 08_5-20_S_3_R2.cl.pd.fq 10_5-08_H_0_R2.cl.pd.fq 10_5-11_H_0_R2.cl.pd.fq 10_5-14_H_0_R2.cl.pd.fq 10_5-17_H_0_R2.cl.pd.fq 10_5-20_S_2_R2.cl.pd.fq 35_6-12_H_0_R2.cl.pd.fq 35_6-15_H_0_R2.cl.pd.fq 35_6-18_H_0_R2.cl.pd.fq 35_6-21_H_0_R2.cl.pd.fq 36_6-12_S_1_R2.cl.pd.fq 36_6-15_S_2_R2.cl.pd.fq 36_6-18_S_3_R2.cl.pd.fq > 08-11-35-36_R2.cl.pd.fq
```

Install latest Trinity and bowtie2, based on directions [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing%20Trinity):

```
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.4.0.zip
unzip Trinity-v2.4.0.zip
rm Trinity-v2.4.0.zip
cd trinityrnaseq-Trinity-v2.4.0
make
make plugins
```

```
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.0/bowtie2-2.3.0-linux-x86_64.zip
unzip bowtie2-2.3.0-linux-x86_64.zip
rm bowtie2-2.3.0-linux-x86_64.zip
ln -s /data/popgen/bowtie2-2.3.0/bowtie2 /usr/local/bin
ln -s /data/popgen/bowtie2-2.3.0/bowtie2-build /usr/local/bin
```

Start screen and run trinity.

```
screen -r 30308.pts-1.pbio381

/data/popgen/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --max_memory 50G \
    --left /data/project_data/fastq/cleanreads/08-11-35-36_R1.cl.pd.fq \
    --right /data/project_data/fastq/cleanreads/08-11-35-36_R2.cl.pd.fq \
    --CPU 20 > run.log 2>&1 &
```

### Let’s install DESeq2 in R studio and look at a script and example data file.

[Here’s](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) the package website with installation instructions, manual, tutorials, etc.

Love MI, Huber W and Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” *Genome Biology*, **15**, pp. 550. [doi: 10.1186/s13059-014-0550-8](http://doi.org/10.1186/s13059-014-0550-8).

```
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```
# February 22, 2017

## R Script from today's work:

```

```

**Terminal Work**

Moved files from server to my desktop to use in my WD in R

```
scp pburnham@pbio381.uvm.edu:/data/project_data/DGE/* .
```





