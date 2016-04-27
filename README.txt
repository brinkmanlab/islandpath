IslandPath-DIOMB 0.1 README
---------------------

OVERVIEW

IslandPath-DIMOB is a computational genomic island prediction method. 
It looks for genomic regions that contain dinucleotide bias and contain at least one mobility gene such as transposases and integrases. 

Genomic islands (GIs) are clusters of genes in prokaryotic genomes of probable horizontal origin. 
GIs are disproportionately associated with microbial adaptations of medical or environmental interest.

IslandPath (see Hsiao et al., 2003) was originally designed to aid to the identification of prokaryotic genomics islands (GIs), by visualizing several common characteristics of GIs such as abnormal sequence composition or the presence of genes that functionally related to mobile elements (termed mobility genes). 

Further studies (see Hsiao et al., 2005) used dinucleotide sequence composition bias and the presence of mobility genes to develop a data set of GIs (IslandPath DIMOB) for multiple organisms and revealed that these genomic regions contain higher proportions of novel genes. 

For more information, see http://www.pathogenomics.sfu.ca/islandviewer and  http://www.pathogenomics.sfu.ca/islandpath 

-----------------------------------------------------------------------------
PREREQUISITES

1. Though IslandPath-DIMOB should work with any OS, it has only
been tested with linux. 

2. Perl version 5.005_03 or higher (5.6.x or higher recommended)
The latest version of Perl can be obtained from http://www.cpan.org

3. The Bioperl library version 1.4 or higher
While IslandPath-DIMOB should work with Bioperl 1.2, 1.4 or above is highly
recommended.  Bioperl can be obtained from www.bioperl.org

4. A working installation of HMMER2 or HMMER3 
HMMER can be obtained from http://hmmer.janelia.org/
"hmmpfam" or "hmmscan" must be within your executable path.

-----------------------------------------------------------------------------
INSTALLATION

1)Extract the archive file
tar xvzf IslandPath_DIMOB_0.1.tar.gz

2)The program is run using the dimob.pl script along with the 3 input files
./dimob.pl your_file.faa your_file.ffn your_file.ptt

If using HMMER2 instead of HMMER3 then you must use the --hmmer2 option.
The output is send to the screen with the start and end coordinates for each predicted GI seperated by a tab (\t)

You can try running the program using the genome files in the example directory
./dimob.pl example/NC_000913.faa example/NC_000913.ffn example/NC_000913.ptt > islands.txt


-----------------------------------------------------------------------------
QUESTIONS? PROBLEMS? COMMENTS?

Email islandpick-mail@sfu.ca (contact person: Morgan Langille) and we'd be happy to help.

-----------------------------------------------------------------------------
BUG REPORTING

If you find a bug, please report it to  islandpick-mail@sfu.ca along with the
following information:

    * version of Perl (output of 'perl -V' is best)
    * version of IslandPath-DIMOB
    * operating system type and version
    * exact text of error message or description of problem

If we don't have access to a system similar to yours, you may be asked to
insert some debugging lines and report back on the results. The more help
and information you can provide, the better!

-----------------------------------------------------------------------------
COPYRIGHT AND LICENCE

IslandPick is distributed under the GNU General Public License. See also
the LICENSE file included with this package.

-----------------------------------------------------------------------------
AUTHOR INFORMATION

IslandPath-DIMOB was written by Will Hsiao and updated by Morgan Langille of the Brinkman Laboratory 
at Simon Fraser University, Burnaby, BC, Canada
http://www.pathogenomics.sfu.ca/brinkman




