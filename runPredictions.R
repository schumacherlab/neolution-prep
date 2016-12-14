## Script runs neolution pipeline
rootDirectory = getwd()

setwd(rootDirectory)
source('./helperFunctions.R')
source('./runConfig.R')

registerDoMC(runOptions$neolution$numberOfWorkers)

# inventorize input files
contextLists = list.files(path = file.path(rootDirectory, "3_neolution"),
													pattern = "varcontext_rna.tsv$",
													full.names = TRUE)

# load HLA typing info
sampleInfo = fread(input = 'sample_info.tsv', na.strings = c('', 'NA', 'N.A.'))

sampleInfo$filepath = sapply(sampleInfo$dna_data_prefix,
														 function(x) {
														 	contextLists[grep(pattern = x,
														 										x = contextLists)]
														 }, USE.NAMES = FALSE)

sampleHlaTypes = as.data.table(sampleInfo %>%
																	 	gather(hla_allele, hla_type, c(4:9)))
sampleHlaTypes = sampleHlaTypes[naturalorder(sampleHlaTypes$dna_data_prefix)]
sampleHlaTypes = unique(sampleHlaTypes[!is.na(hla_type)],
														by = c("filepath", "hla_type"))
sampleHlaTypes[, hla_type := gsub(pattern = 'HLA-|\\*|\\:', replacement = '', x = sampleHlaTypes$hla_type)]

# optional: exclude C alleles
# sampleHlaTypes = sampleHlaTypes[hla_allele != 'hla_c_1' & hla_allele != 'hla_c_2', ]

if (any(nchar(sampleHlaTypes$hla_type) < 5)) {
	stop('Please check HLA type input: HLA types with irregular name(s) found.\n', paste(sampleHlaTypes[nchar(sampleHlaTypes$hla_type) < 5]$hla_type, collapse = ', '))
}

# start predictions
setwd("/home/NFS/users/l.fanchi/dev_environments/neolution-live/")

x = foreach(i = 1:nrow(sampleHlaTypes)) %dopar% {
	invisible(sapply(seq(1, length(runOptions$neolution$xmer)),
									 function(y) {
									 	system(command = paste("Rscript performPredictions.R",
									 												 "-f", sampleHlaTypes$filepath[i],
									 												 "-m", sampleHlaTypes$hla_type[i],
									 												 # "-a", sampleHlaTypes$affinity_cutoff[i],
									 												 "-r", runOptions$neolution$rankCutoff,
									 												 "-p", runOptions$neolution$processingCutoff,
									 												 "-l", runOptions$neolution$xmer[y],
									 												 "--selfsim"))
									 }))
}
