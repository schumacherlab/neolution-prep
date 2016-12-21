### Script for preparing input for neolution pipeline

1. Make project folder and navigate to it in the Terminal
2. Do 'git clone http://gitlab.nki.nl/l.fanchi/neolution-prep.git .' (the dot matters)
3. In the project folder, create directory '1a\_variants' and copy VCF files to it
4. Optional: Create directory '1b\_rnaseq\_data' in project folder
	* create sub-directory 'bam' and copy BAM and BAI files to it
	* create sub-directory 'processed' and copy expression level data to it (e.g. cufflinks output)
5. Create a TSV file with sample information. See below for format & required info

|  patient_id  |     dna\_data\_prefix     |      rna\_data\_prefix       |  hla\_a\_1  |  hla\_a\_2  |  hla\_b\_1  |  hla\_b\_2  |  hla\_c\_1  |  hla\_c\_2  |
|:------------:|:-------------------------:|:----------------------------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|
|  ID #1       | 4152\_1\_CF8585\_GATAGACA |  4153\_1\_CF8597\_TTAGGCA\_  |  A03:01   |  A01:01   |  B08:01   |  B16:01   |    NA     |    NA     |
|  ID #2       | 4152\_2\_CF8714\_GCCACATA |  4153\_2\_CF8716\_ACTTGAA\_  |  A02:01   |  A09:01   |  B36:03   |  B52:01   |    NA     |    NA     |

* The '_prefix' columns should contain the (unique) beginnings of the dna and rna input filenames  
* Empty cells or NAs can be used to exclude alleles 