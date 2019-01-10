README for file ChlamydomonasTranscriptNameConversionBetweenReleases.Mch12b.txt

4/10/2014

File contents
=============

This file contains transcript identifiers from Chlamydomonas annotation versions
organized in rows of equivalent transcripts allowing batch conversion of transcript
identifiers between versions.

This version of the conversion file is Mch 12b

The file starts with a header line like so
#Correspondence table from v5.5

The next line of the file contains column headings, starting with a comment character
('#'). Columns are space-padded to 25 characters.
These are the column headings, in order, together with an explanation of what version they correspond to
5.5       JGI v5.5 in Phytozome v10 
3.1       JGI v3.1 (published in genome paper Merchant et al., 2007)
Genbank   Genbank submission of genome and annotations from Merchant et al. (2007)
4         JGI v4 annotations 
4.3       JGI v4.3 (based on Augustus u10.2 annotations)
u5        Augustus u5 annotations
u9        Augustus u9 annotations
5.3.1     JGI v5.3.1 in Phytozome v9.1


Conversion algorithms used to generate the file
================================================

Gene mapping from previous releases was done by Simon Prochnik at the Joint Genome Institute using PASA (Haas, 2003) to map annotations on genome v3 and v4 and a locus mapping tool for the others (This tool is JGI software that uses BLAT (Kent, 2002) and BLAST (Altschul, 1990) to map equivalent gene regions using information from adjacent loci to distinguish paralogous regions). Multiple transcripts in different annotation versions based on the v5 assembly are matched up with each other's intron-exon coordinates using the bedtools jaccard tool. The mappings available on this site will be updated as new versions of the genome and/or annotations are released.

Below is a summary of the different annotation versions that are available for interconversion.

From JGI v4.3 (Augustus u10.2) onwards, stable locus identifiers were generated. These are of the form Cre01.g123450.t2.3 where Cre01 means Chlamydomonas reinhardtii chromosome 1; 123450 is a gene serial number, starting at the beginning of chromosome 1; 't2' means this is the second transcript at this locus and '.3' indicates it is the third version of the transcript. The locus identifier (e.g. Cre01.g006450') or the transcript identifier (e.g. 'Cre01.g006450.t1') may be used as well as the full identifier including the version number (e.g. 'Cre01.g006450.t1.1') can be used in this search tool. For more details on the annotations, see Blaby et al (2014) Trends in Plant Science (submitted). 

Annotation Version         Version Details	                  
JGI v3.1                   protein IDs of JGI v3.1 annotations published in 2007 genome paper
Genbank	                   protein IDs of JGI v3.1 annotation submitted to GenBank based on genome paper
JGI v4                     protein IDs of JGI annotations on v4 genome assembly
Augustus u5                Augustus update 5 (u5) annotations produced on v4 assembly
Augustus u9                Augustus update 9 (u9) produced on v4 assembly
JGI v4.3 (Phytozome 8)     Augustus update 10.2 (u10.2) produced on v4 assembly. Released as JGI v4.3 in Phytozome 8
JGI v5.3.1 (Phytozome 9/9.1)  Augustus update 11.6 (u11.6)-based annotations on v5 assembly, released as JGI v5.3.1 in Phytozome 9/9.1
JGI v5.5 (Phytozome 10)	   Augustus update 11.6 (u11.6)-based annotations on v5 assembly, released as JGI v5.5 in Phytozome 10

Batch conversion tool
======================

There is a batch ID conversion tool on the UCLA Algal Functional Annotation Tool site at
http://pathways.mcdb.ucla.edu/algal/id_conversion.html
