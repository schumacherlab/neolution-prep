## Script prepares files for neolution pipeline
# make sure to edit runConfig before running this script

# rootDirectory = path.expand(getwd())
rootDirectory = path.expand('~/projects/mypetproject')
setwd(rootDirectory)

source('./helperFunctions.R')
source('./runConfig.R')

# extract fields from VCF files and run variant context generation
performVarcontextGeneration()

# prepare Neolution input files (also allows optional merging with RNAseq data)
prepareNeolutionInput()

# find evidence for variants in RNAseq data (optional)
findRnaReadLevelEvidenceForVariants()

# generate SnpEff output
runSnpEff()
