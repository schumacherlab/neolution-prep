## Script prepares files for neolution pipeline
# make sure to edit runConfig before running this script

# rootDirectory = path.expand(getwd())
rootDirectory = path.expand('~/projects/mypetproject')
setwd(rootDirectory)

source('./helperFunctions.R')
source('./runConfig.R')

# parse and clean VCF
parseAndExtractFieldsFromVcf()

# find evidence for variants in RNAseq data (skipped if no RNAseq data is found)
findRnaReadLevelEvidenceForVariants(quant_mode = 'salmon')

# run variant context generation
performVarcontextGeneration()

# prepare Neolution input files (also allows optional merging with RNAseq data; skipped if no RNAseq data found)
prepareNeolutionInput(rna_path = '1b_rnaseq_data/processed_salmon',
											rna_file_suffix = 'salmon-quant-by-ensg\\.tsv',
											expression_unit = 'tpm')

# generate SnpEff output
runSnpEff()
