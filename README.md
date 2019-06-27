# IslandPath-DIMOB

IslandPath-DIMOB is a standalone software to predict genomic islands in bacterial and archaeal genomes based on the presence of dinucleotide biases and mobility genes.

Genomic islands (GIs) are clusters of genes in prokaryotic genomes of probable horizontal origin. 
GIs are disproportionately associated with microbial adaptations of medical or environmental interest.

This version here is a modification of the original version (https://github.com/brinkmanlab/islandpath) that allows IslandPath-DIMOB to process assembly draft genomes in multiple contigs and not using a reference genome.

## Install

Although the original version has multiple instalation options: docker, github releases... for this particular version we would only rely on github clone. 

Clone the latest code from github:

    git clone https://github.com/JFsanchezherrero/islandpath
    
**_Dependencies_**

No additional dependencies were added to this new implementation. IslandPath-DIMOB remains with the same original dependencies. 

1. Though IslandPath-DIMOB should work with any OS, it has only been tested with linux. 

2. Perl version 5.18.x or higher  
The latest version of Perl can be obtained from http://www.cpan.org

3. The following Perl libraries are also required:
    - Data::Dumper
    - Log:Log4perl
    - Config::Simple
    - Moose
    - MooseX::Singleton
    - Bio::Perl

4. A working installation of HMMER3  
HMMER can be obtained from http://hmmer.org/  
"hmmscan" must be within your executable path.


## Run

IslandPath-DIMOB v1.0.1_b takes as input an annotated complete/draft genome as a genbank (.gbk) or an embl (.embl) file.

	perl Dimob.pl <genome.gbk> <output_name> [cutoff_dinuc_bias] [min_length]

	Default values:
		cutoff_dinuc_bias = 8
		min_length = 8000

	Example:
		perl Dimob.pl example/NC_003210.gbk NC_003210_GIs
		perl Dimob.pl example/NC_003210.gbk NC_003210_GIs 6 10000
		perl Dimob.pl example/NC_000913.embl NC_000913_GIs 6 10000
        
## Citation

[Bertelli and Brinkman, 2018](https://doi.org/10.1093/bioinformatics/bty095)  
[Hsiao et al., 2005](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0010062)


## Questions? Comments? Bugs?

Email islandpick-mail@sfu.ca (contact person: Claire Bertelli) for the original version.

Email jsanchez@igtp.cat (contact person: Jose F. Sanchez-Herrero) or via github issue for the modification version.


## Copyright and License

IslandPath-DIMOB is distributed under the GNU General Public License. See also the LICENSE file included with this package.

Give credit to the original version of this software available at https://github.com/brinkmanlab/islandpath


## Versions - New features

### 27/06/2019 - IslandPath-DIMOB v1.0.1_b
Several new implementations were performed in order to add multicontig functionality.

We additionally fix some bugs and warning messages: 
- Fix smartmacth experimental warning message
- Fix dinuc bias loop iteration bug https://github.com/brinkmanlab/islandpath/issues/8#issue-459228372

We increase the input and output options
- Use GI min_length as a variable
- Use cutoff_genes_dinuc as a variable
- Output csv information for dinucleotide bias.
- Provide additional information such as annotation of the genes within each GI

### 23/12/2016 - IslandPath-DIMOB v1.0.0  
Increased recall and precision in the prediction of genomic islands based on the presence of dinucleotide bias and mobility genes. Standardization of input file types, and automatic generation of the other file types required by IslandPath-DIMOB.  
Input: gbk or embl file  
Publication: [Bertelli and Brinkman, 2018](https://doi.org/10.1093/bioinformatics/bty095)  

### 2008 - IslandPath-DIMOB
Improvement and assessment of IslandPath-DIMOB predictions by Morgan Langille  
Input files: ffn, faa, ptt  
Publication: [Langille et al., 2008](http://www.biomedcentral.com/1471-2105/9/329)

### 2005 - IslandPath-DIMOB 
Second version developed by Will Hsiao  
Further studies used dinucleotide sequence composition bias and the presence of mobility genes to develop a data set of GIs (IslandPath DIMOB) for multiple organisms and revealed that these genomic regions contain higher proportions of novel genes.  
Input files: ffn, faa, ptt  
Publication: [Hsiao et al., 2005](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0010062)

### 2003 - IslandPath-DINUC
IslandPath-DINUC developed by Will Hsiao  
IslandPath was originally designed to aid to the identification of prokaryotic genomics islands (GIs), by visualizing several common characteristics of GIs such as abnormal sequence composition or the presence of genes that functionally related to mobile elements (termed mobility genes).  
Publication: [Hsiao et al., 2003](http://bioinformatics.oxfordjournals.org/content/19/3/418.short)  


## Authors

IslandPath-DIMOB was written and updated by several members of the [Brinkman Laboratory](http://www.brinkman.mbb.sfu.ca/) at Simon Fraser University, Burnaby, BC, Canada

2015 - present:     Claire Bertelli    claire.bertelli@sfu.ca  
2009 - 2015:    Matthew Laird    lairdm@sfu.ca  
2007 - 2009: Morgan Langille    morgan.g.i.langille@dal.ca  
2003 - 2007: Will Hsiao william.hsiao@bccdc.ca

