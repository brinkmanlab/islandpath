# IslandPath-DIMOB

IslandPath-DIMOB is a standalone software to predict genomic islands in bacterial and archaeal genomes based on the presence of dinucleotide biases and mobility genes.

Genomic islands (GIs) are clusters of genes in prokaryotic genomes of probable horizontal origin. 
GIs are disproportionately associated with microbial adaptations of medical or environmental interest.

The latest IslandPath-DIMOB version is integrated in [IslandViewer 4](http://www.pathogenomics.sfu.ca/islandviewer/browse/), the leading integrated web interface for genomic island prediction.

## Install

A pre-built Docker image is available in the [brinkmanlab docker hub](https://hub.docker.com/r/brinkmanlab/islandpath/). Using this pre-installed version of IslandPath-DIMOB ensures the software runs according to expectations.

Users wishing to install locally IslandPath-DIMOB can download a [release](http://github.com/brinkmanlab/islandpath/releases/) and install required perl libraries and HMMER listed below:

Alternatively, you can also clone the latest code from github:

    git clone https://github.com/brinkmanlab/islandpath
    
Please note that IslandPath-DIMOB predictions should only take a couple of minutes per bacterial genome. It was recently reported that IslandPath-DIMOB was extremely slow on Mac OS X with a conda installation of perl libraries. While we investigate the reason for this issue, we recommend using the Docker image.
    
  
**_Dependencies_**

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

IslandPath-DIMOB v1.0.0 takes as input an annotated complete genome as a genbank (.gbk) or an embl (.embl) file.

    # gbk file
    ./Dimob.pl example/NC_003210.gbk NC_003210_GIs.txt
    # embl file
    ./Dimob.pl example/NC_000913.embl NC_000913_GIs.txt


## Citation

[Bertelli and Brinkman, 2018](https://doi.org/10.1093/bioinformatics/bty095)  
[Hsiao et al., 2005](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0010062)


## Questions? Comments? Bugs?

Email islandpick-mail@sfu.ca (contact person: Claire Bertelli) and we'd be happy to help.

If you find a bug, please report it to islandpick-mail@sfu.ca along with the
following information:

* version of Perl (output of 'perl -V' is best)
* version of IslandPath-DIMOB
* operating system type and version
* exact text of error message or description of problem

If we don't have access to a system similar to yours, you may be asked to insert some debugging lines and report back on the results. The more help and information you can provide, the better!


## Copyright and License

IslandPath-DIMOB is distributed under the GNU General Public License. See also the LICENSE file included with this package.


## Versions - New features

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

