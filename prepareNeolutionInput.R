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
findRnaReadLevelEvidenceForVariants()

# run variant context generation
performVarcontextGeneration()

# prepare Neolution input files (also allows optional merging with RNAseq data; skipped if no RNAseq data found)
prepareNeolutionInput()

# generate SnpEff output
runSnpEff()
