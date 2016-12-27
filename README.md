### Script for preparing input for neolution pipeline

1. Make project directory and navigate to it in the Terminal
2. Do 'git clone http://gitlab.nki.nl/l.fanchi/neolution-prep.git .' (the dot matters) **NOTE: directory must be EMPTY**
3. In the project directory, create sub-directory '1a\_variants' 
	* in it, create sub-directory 'vcf' and copy VCF file(s) to it  
4. Optional: Create directory '1b\_rnaseq\_data' in project folder
	* create sub-directory 'processed' and copy expression level data to it (e.g. cufflinks output)
	* create sub-directory 'bam' and copy BAM and BAI file(s) to it (necessary in case you want to determine expression of mutant allele)
5. Edit the sample_info.tsv file. **Make sure to leave tab separation between data!**  
Example data:

| patient_id | dna\_data\_prefix | rna\_data\_prefix | hla\_a\_1 | hla\_a\_2 | hla\_b\_1 | hla\_b\_2 | hla\_c\_1 | hla\_c\_2 |
|:-----:|:-------------------------:|:------------------------:|:-----:|:-----:|:-----:|:-----:|:----:|:----:|
| ID #1 | 4152\_1\_CF8585\_GATAGACA | 4153\_1\_CF8597\_TTAGGCA | A0301 | A0101 | B0801 | B1601 |  NA  |  NA  |
| ID #2 | 4152\_2\_CF8714\_GCCACATA | 4153\_2\_CF8716\_ACTTGAA | A0201 | A0901 | B3603 | B5201 |  NA  |  NA  |

* The '_prefix' columns should contain the (unique) beginnings of the dna and rna input filenames
* Don't use special characters in HLA allele names (no asterix '*' or colon ':')  
* Empty cells or NAs can be used to exclude alleles 