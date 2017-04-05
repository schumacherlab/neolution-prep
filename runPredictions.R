## Script runs neolution pipeline
project_directory = getwd()

setwd(project_directory)
source('./helperFunctions.R')
source('./runConfig.R')

user_options = parseCommandlineArguments()
if ('threads' %in% names(user_options)) {runOptions$neolution$numberOfWorkers = as.numeric(user_options$threads)}

registerDoMC(runOptions$neolution$numberOfWorkers)

# inventorize input files
context_files = list.files(path = file.path(project_directory, '3_neolution'),
													 pattern = 'varcontext\\.tsv$',
													 full.names = TRUE)

# load HLA typing info
sample_info = fread(input = 'sample_info.tsv', na.strings = c('', 'NA', 'N.A.'))

sample_info[, filepath := sapply(dna_data_prefix, function(x) {context_files[grep(pattern = x,
																																									x = context_files)]}, USE.NAMES = FALSE)]

samples_by_hla = as.data.table(sample_info %>% gather(hla_allele, hla_type, c(4:9)) %>% dplyr::filter(!is.na(hla_type)))
samples_by_hla = unique(x = samples_by_hla[naturalorder(samples_by_hla$dna_data_prefix)],
												by = c('filepath', 'hla_type'))
samples_by_hla[, hla_type := gsub(pattern = 'HLA-|\\*|\\:', replacement = '', x = hla_type)]

# optional: exclude C alleles
# samples_by_hla[hla_allele != 'hla_c_1' & hla_allele != 'hla_c_2']

if (any(nchar(samples_by_hla$hla_type) < 5)) {
	stop('Please check HLA type input: HLA types with irregular name(s) found.\n', paste(samples_by_hla[nchar(hla_type) < 5, hla_type], collapse = ', '))
}

# switch to run on mulitple servers
switch(EXPR = system(command = 'hostname', intern = TRUE),
			 # 'steroid' = {subset = seq(round(nrow(samples_by_hla)/2) + 1, nrow(samples_by_hla))},
			 'void' = {subset = seq(1, round(nrow(samples_by_hla)/2))},
			 'paranoid' = {subset = seq(round(nrow(samples_by_hla)/2) + 1, nrow(samples_by_hla))})

samples_by_hla = samples_by_hla[subset]

# start predictions
setwd('/home/NFS/users/l.fanchi/dev_environments/neolution-live/')

invisible(foreach(i = 1:nrow(samples_by_hla)) %dopar% {
	invisible(sapply(seq(1, length(runOptions$neolution$xmer)),
									 function(y) {
									 	system(command = paste('Rscript performPredictions.R',
									 												 '-f', samples_by_hla$filepath[i],
									 												 '-m', samples_by_hla$hla_type[i],
									 												 # '-a', samples_by_hla$affinity_cutoff[i],
									 												 # '-r', runOptions$neolution$rankCutoff,
									 												 # '-p', runOptions$neolution$processingCutoff,
									 												 '-d', runOptions$neolution$model_cutoff,
									 												 '-l', runOptions$neolution$xmer[y],
									 												 if (runOptions$neolution$selfsim_filter_mode == 'simple') {
									 												 	'--selfsim'
									 												 } else if (runOptions$neolution$selfsim_filter_mode == 'extended') {
									 												 	'--extselfsim'
									 												 },
									 												 if (runOptions$neolution$selflist) {
									 												 	'--selflist'
									 												 })
									 	)
									 })
	)
})
