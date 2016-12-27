# load required packages
if (!require("pacman")) install.packages("pacman")

required_packages = c('data.table', 'gtools', 'utils', 'optparse', 'RMySQL', 'compiler', 'naturalsort', 'parallel', 'doMC', 'stringr', 'tidyr', 'foreach')

library(pacman)
pacman::p_load(char = required_packages)

# helper functions
performVarcontextGeneration = function(vcf_path = file.path(rootDirectory, '1a_variants', 'vcf'), vcf_fields = c('ID', 'CHROM', 'POS', 'REF', 'ALT')) {
	registerDoMC(runOptions$varcontext$numberOfWorkers)

	# extract relevant fields from VCFs
	message('Step 1: Extracting "', paste(vcf_fields, collapse = ' '), '" fields from VCF')

	dir.create(file.path(rootDirectory, '1a_variants', 'extr_fields'), showWarnings = F)

	extractFieldsFromVCF(vcf_path = vcf_path,
											 vcf_fields = vcf_fields)

	# make list of input files
	variantLists = list.files(path = vcf_path,
														pattern = '\\.tsv$',
														recursive = FALSE,
														full.names = TRUE)

	# remove 'chr' prefix, if present
	invisible(sapply(variantLists,
									 function(x) {
									 	system(command = paste("perl -p -i -e 's/chr//g'", x))
									 }))

	# generate variant contexts
	message('Step 2: Generating context for variants')
	generateVarcontext(input_list = variantLists)

	# move files
	system(command = paste('cd', vcf_path, ';', 'mv -f *.tsv ../extr_fields'))
}

extractFieldsFromVCF = function(vcf_path, vcf_fields = c('ID', 'CHROM', 'POS', 'REF', 'ALT')) {
	extractFields = paste("export SHELL=/bin/bash;",
												"find", vcf_path, "-name \'*.vcf\' -print0",
												"| sed s/\\\\.vcf//g",
												"| xargs -0 -I % echo java -jar -Xmx2G", file.path(runOptions$snpeff$path, "SnpSift.jar"), "extractFields -s \\\",\\\" -e \\\".\\\" %.vcf", paste(vcf_fields, collapse = ' '),
												"\\> %.tsv",
												"| cat | parallel -j8")
	system(command = extractFields)
}

generateVarcontext = function(input_list) {
	if (length(input_list) < 1) {
		stop('No input lists found...')
	}
	dir.create(file.path(rootDirectory, '2_varcontext'),
						 showWarnings = FALSE)
	dir.create(file.path(rootDirectory, '2_varcontext', 'varcontext_logs'),
						 showWarnings = FALSE)

	setwd(runOptions$varcontext$varcontextDirectory)

	invisible(foreach(i = 1:length(input_list)) %dopar% {
		filename = sub(pattern = '\\.[^.]+$', replacement = '', x = basename(input_list[i]))

		runStart = format(Sys.time(), '%Y%m%d')

		# write run info to log
		write(x = paste0(Sys.time(),' - Varcontext run start\n\n',
										 'Branch:\t\t\t', system('git symbolic-ref --short -q HEAD', intern = TRUE), '\n',
										 'Commit hash:\t\t', system('git rev-parse HEAD', intern = TRUE), '\n\n',
										 'Input file:\t\t', input_list[i], '\n',
										 'Ensembl API:\t\t', runOptions$varcontext$ensemblApi, '\n\n',
										 'Field separator:\t', runOptions$varcontext$fieldSeparator, '\n',
										 'Canonical transcripts:\t', runOptions$varcontext$canonicalOnly, '\n',
										 'Peptide context:\t', runOptions$varcontext$peptideContext,'\n',
										 'NMD status:\t', runOptions$varcontext$nmdStatus, '\n'),
					file = file.path(rootDirectory, '2_varcontext', 'varcontext_logs',
													 paste(filename,
													 			'runInfo.txt',
													 			sep = "_")),
					append = FALSE)

		system(command = paste0('export ENSEMBLAPI="', runOptions$varcontext$ensemblApi, '";',
														'export PERL5LIB="$PERL5LIB:', runOptions$varcontext$perlLibs,'";',
														'perl ', file.path(runOptions$varcontext$varcontextDirectory, 'varcontext/create_context.pl'), ' ',
														'--separator=', runOptions$varcontext$fieldSeparator, ' ',
														ifelse(runOptions$varcontext$canonicalOnly, '--canonical ', ''),
														ifelse(runOptions$varcontext$peptideContext, '--peptide ', ''),
														ifelse(runOptions$varcontext$nmdStatus, '--nmd ', ''),
														'"', input_list[i], '"',
														' 1> "', file.path(rootDirectory, '2_varcontext', paste(filename, 'varcontext.tsv"', sep = '_')),
														' 2> "', file.path(rootDirectory, '2_varcontext', 'varcontext_logs', paste(filename, 'warnings.log"', sep = '_'))),
					 intern = FALSE)
	})

	setwd(rootDirectory)
}

prepareNeolutionInput = function(varcontext_path = file.path(rootDirectory, '2_varcontext'),
																 rna_path = file.path(rootDirectory, '1b_rnaseq_data', 'processed'),
																 sample_info_path = file.path(rootDirectory, 'sample_info.tsv'),
																 rna_file_suffix = 'genes\\.fpkm_tracking',
																 expression_unit = 'FPKM') {
	varcontext_data = lapply(list.files(path = varcontext_path,
																			pattern = 'varcontext\\.tsv',
																			full.names = TRUE),
													 fread,
													 colClasses = list(character = c('chromosome', 'nmd_remark')))
	varcontext_data = setNames(object = varcontext_data,
														 nm = list.files(path = varcontext_path,
														 								pattern = 'varcontext\\.tsv'))

	dir.create(file.path(rootDirectory, '3_neolution'),
						 showWarnings = FALSE)

	if (file.exists(sample_info_path) & file.exists(rna_path)) {
		message('Step 3a: Merging RNA expression data')

		rnaseq_files = list.files(path = rna_path,
															pattern = rna_file_suffix,
															full.names = TRUE)

		if (length(rnaseq_files) < 1) {
			stop('RNAseq directory does not contain files with suffix: ', rna_file_suffix)
		}

		rnaseq_data = lapply(rnaseq_files,
												 function(x) {
												 	data = fread(x,
												 							 colClasses = list(character = c(expression_unit)))
												 	data[[expression_unit]] = as.numeric(data[[expression_unit]])
												 	return(data)
												 })

		rnaseq_data = setNames(object = rnaseq_data,
													 nm = list.files(path = rna_path,
													 								pattern = rna_file_suffix))

		sample_info = fread(sample_info_path,
												na.strings = c('', 'NA', 'N.A.'))

		sample_combinations = data.table(variants = sapply(sample_info$dna_data_prefix, function(x) grep(pattern = x,
																																																		 x = list.files(path = varcontext_path,
																																																		 							 pattern = 'varcontext\\.tsv',
																																																		 							 full.names = FALSE),
																																																		 value = T),
																											 USE.NAMES = FALSE),
																		 rna_expression_data = sapply(sample_info$rna_data_prefix, function(x) grep(pattern = x,
																		 																																					 x = list.files(path = rna_path,
																		 																																					 							 pattern = rna_file_suffix,
																		 																																					 							 full.names = FALSE),
																		 																																					 value = T),
																		 														 USE.NAMES = FALSE)
		)

		# check if both varcontext and rna_expression data are present before merge
		sample_combinations = sample_combinations[sapply(1:nrow(sample_combinations),
																										 function(x) file.exists(file.path(varcontext_path, sample_combinations$variants[x]))
																										 &
																										 	file.exists(file.path(rna_path, sample_combinations$rna_expression_data[x])))]

		if (nrow(sample_combinations) < 1) {
			stop('No samples left after file checking, please check sample combinations')
		}

		prediction_input = lapply(seq(1, nrow(sample_combinations)),
															function(x) {
																mergeByEnsemblId(variant_table = varcontext_data[[sample_combinations[x, variants]]],
																								 expression_table = rnaseq_data[[sample_combinations[x, rna_expression_data]]])
															})
		prediction_input = setNames(object = prediction_input, nm = sub("[.][^.]*$", "", sample_combinations$variants))

		rna_coverage_summary = data.table(sample = names(prediction_input),
																			percent_ensg_coverage = sapply(prediction_input, function(x) round(x = length(which(!is.na(x[[expression_unit]]))) / nrow(x) * 100,
																																																				 digits = 1)),
																			percent_no_expression = sapply(prediction_input, function(x) round(x = length(which(x[[expression_unit]] == 0)) / nrow(x) * 100,
																																																				 digits = 1))
		)

		write.table(x = rna_coverage_summary,
								file = file.path(rootDirectory, '3_neolution', 'rna_expression_coverage_info.tsv'),
								sep = '\t',
								quote = FALSE,
								append = FALSE)

		message('Step 3b: Generating Neolution pipeline input')
	} else {
		message('Step 3: Generating Neolution pipeline input (no RNAseq data)')

		prediction_input = varcontext_data
	}

	prediction_input = lapply(prediction_input,
														function(x) {
															if (all(is.na(x[, rna_expression]))) return(x[, !names(x) == 'rna_expression', with = F])
														})

	invisible(mapply(FUN =
									 	function(x, y) {
									 		write.table(x = x,
									 								file = file.path(rootDirectory, '3_neolution', paste0(names(prediction_input)[y], '_rna.tsv')),
									 								sep = '\t',
									 								row.names = F)
									 	},
									 prediction_input,
									 seq(1, length(prediction_input))))
}

mergeByEnsemblId = function(variant_table, expression_table, expression_unit = 'FPKM') {
	if (is.null(expression_table)){
		return(variant_table[, rna_expression := NA])
	} else if ('gene_id' %in% names(variant_table) && 'gene_id' %in% names(expression_table)) {
		table_merged = merge(x = variant_table,
												 y = expression_table[!duplicated(gene_id), c('gene_id', expression_unit), with = FALSE],
												 by = 'gene_id',
												 all.x = TRUE)
		return(table_merged)
	} else {
		warning('gene_id column missing in either variant table or expression table')
	}
}

findRnaReadLevelEvidenceForVariants = function(neolution_input_path = file.path(rootDirectory, '3_neolution'),
																							 rna_path = file.path(rootDirectory, '1b_rnaseq_data/bam'),
																							 sample_info_path = file.path(rootDirectory, 'sample_info.tsv')) {
	if (!file.exists(rna_path)) {
		stop('Please put RNAseq BAM files in "./1b_rnaseq_data/bam" subdirectory')
	}

	# parse neolution input data
	input_data = lapply(list.files(path = neolution_input_path,
																 pattern = '\\.tsv',
																 full.names = TRUE),
											fread,
											colClasses = list(character = c('chromosome', 'nmd_remark')))
	input_data = setNames(object = input_data,
												nm = list.files(path = neolution_input_path,
																				pattern = '\\.tsv'))

	# load sample info
	if (file.exists(sample_info_path)) {
		sample_info = fread(sample_info_path, sep = '\t', header = T, na.strings = c('','NA', 'N.A.'))
	} else {
		stop('Sample info file missing, please provide path in argument to findRnaReadLevelEvidenceForVariants')
	}

	# make list of all unique variants found in varcontext (mapping to an exon)
	snv_positions = lapply(input_data,
												 function(x) {
												 	setorder(x = unique(x[grepl(pattern = '^[0-9]{1,2}$|^[XY]$',
												 															x = x$chromosome) &
												 													!grepl(pattern = '[rg]s[0-9]+',
												 																 x = x$variant_id) &
												 													nchar(x$ref_allele) == 1 &
												 													nchar(x$alt_allele) == 1,
												 												.(chromosome, start_position), ]),
												 					 chromosome,
												 					 start_position)
												 })

	# generate position list files to feed to samtools mpileup
	dir.create(path = file.path(rootDirectory, '1a_variants', 'poslist'), showWarnings = F)

	invisible(sapply(seq(1, length(snv_positions)),
									 function(x) write.table(x = snv_positions[[x]],
									 												file = file.path(rootDirectory, '1a_variants', 'poslist', paste0(sub("\\.[^.]*$", "", names(snv_positions)[x]), '_poslist.tsv')),
									 												quote = FALSE,
									 												col.names = FALSE,
									 												row.names = FALSE)))

	# make overview of sample input files
	sample_combinations = data.table(variants = sapply(sample_info$dna_data_prefix, function(x) grep(pattern = x,
																																																x = list.files(path = file.path(rootDirectory, '1a_variants', 'poslist'),
																																																							 pattern = '_poslist\\.tsv',
																																																							 full.names = TRUE),
																																																value = T),
																										 USE.NAMES = FALSE),
																	 rna_bam_file = sapply(sample_info$rna_data_prefix, function(x) grep(pattern = x,
																	 																																				 x = list.files(path = rna_path,
																	 																																				 							 pattern = '\\.bam$',
																	 																																				 							 full.names = TRUE),
																	 																																				 value = T),
																	 												USE.NAMES = FALSE)
	)

	# # perform pileups on whole bams
	# invisible(sapply(seq(1, nrow(sample_combinations)),
	# 								 function(x) {
	# 								 	if (!file.exists(file.path(rootDirectory, '1b_rnaseq_data', 'pileups', paste0(sub('[.][^.]*$', '', basename(sample_combinations$rna_bam_file[x])), '_mpil.tsv')))) {
	# 								 		performSamtoolsPileup(bam_file = sample_combinations$rna_bam_file[x])
	# 								 	}})
	# )
	#
	# # determine bam read depth summary
	# pileup_data = lapply(list.files(path = file.path(rootDirectory, '1b_rnaseq_data', 'pileups'),
	# 																pattern = '_mpil.tsv',
	# 																full.names = TRUE),
	# 										 fread,
	# 										 colClasses = c(V1 = 'character'),
	# 										 col.names = c('chromosome', 'start_position', 'ref_base', 'number_of_reads', 'rna_read_bases', 'base_quality'))
	#
	# pileup_data = lapply(pileup_data,
	# 										 function(x) {
	# 										 	x$rna_read_bases = gsub(pattern = '[^atgcATGC]',
	# 										 													replacement = '',
	# 										 													x = x$rna_read_bases)
	# 										 	x$rna_read_bases = ifelse(test = is.na(x$rna_read_bases) | x$rna_read_bases == '',
	# 										 														yes = NA,
	# 										 														no = x$rna_read_bases)
	#
	# 										 	x$rna_total_read_count = ifelse(test = is.na(x$rna_read_bases),
	# 										 																	yes = NA,
	# 										 																	no = nchar(x$rna_read_bases))
	# 										 })
	#
	# read_depth_summary = lapply(pileup_data,
	# 														function(x) summary(x$rna_total_read_count))
	#
	# read_depth_summary = setNames(object = read_depth_summary, nm = list.files(path = file.path(rootDirectory, '1b_rnaseq_data', 'pileups'),
	# 																																					 pattern = '_mpil.tsv',
	# 																																					 full.names = FALSE))

	# perform pileups on variant locations
	invisible(sapply(seq(1, nrow(sample_combinations)),
									 function(x) {
									 	if (!file.exists(file.path(rootDirectory, '1b_rnaseq_data', 'pileups', paste0(sub('[.][^.]*$', '', basename(sample_combinations$rna_bam_file[x])), '_mpil_loc.tsv')))) {
									 		performSamtoolsPileup(locations_file = sample_combinations$variants[x], bam_file = sample_combinations$rna_bam_file[x])
									 	}})
						)

	pileup_loc_data = lapply(list.files(path = file.path(rootDirectory, '1b_rnaseq_data', 'pileups'),
																			pattern = '_mpil_loc.tsv',
																			full.names = TRUE),
													 fread,
													 colClasses = c(V1 = 'character'),
													 col.names = c('chromosome', 'start_position', 'ref_base', 'number_of_reads', 'rna_read_bases', 'base_quality'))

	pileup_loc_data = setNames(object = pileup_loc_data, nm = list.files(path = file.path(rootDirectory, '1b_rnaseq_data', 'pileups'),
																																			 pattern = '_mpil_loc.tsv',
																																			 full.names = FALSE))

	input_pileup_merge = lapply(seq(1, length(input_data)),
															function(index_input_data) {
																index_prefix_dna = which(sapply(sample_info$dna_data_prefix, function(y) grepl(pattern = y, x = names(input_data)[index_input_data], fixed = T), USE.NAMES = F))
																if (index_prefix_dna > 0) {
																	index_pileup_data = which(grepl(pattern = sample_info$rna_data_prefix[index_prefix_dna], x = names(pileup_loc_data), fixed = T))
																}
																if (index_pileup_data > 0) {
																	merged_data = merge(x = input_data[[index_input_data]],
																											y = pileup_loc_data[[index_pileup_data]][, .(chromosome, start_position, rna_read_bases)],
																											by = c('chromosome', 'start_position'),
																											all.x = TRUE)
																	setorder(merged_data, variant_id, chromosome, start_position)
																}})

	input_pileup_merge = lapply(input_pileup_merge,
															function(x) {
																x$rna_read_bases = gsub(pattern = '[^atgcATGC]',
																												replacement = '',
																												x = x$rna_read_bases)
																x$rna_read_bases = ifelse(test = is.na(x$rna_read_bases) | x$rna_read_bases == '',
																													yes = NA,
																													no = x$rna_read_bases)

																x$rna_ref_read_count = str_count(string = toupper(x$rna_read_bases),
																																 pattern = toupper(x$ref_allele))

																x$rna_alt_read_count = str_count(string = toupper(x$rna_read_bases),
																																 pattern = toupper(x$alt_allele))

																x$rna_total_read_count = ifelse(test = is.na(x$rna_read_bases),
																																yes = NA,
																																no = nchar(x$rna_read_bases))

																x$rna_vaf = str_count(string = toupper(x$rna_read_bases),
																											pattern = toupper(x$alt_allele)) /
																	nchar(x$rna_read_bases)

																x$rna_read_bases = NULL

																return(x)
															})

	input_pileup_merge = lapply(input_pileup_merge,
															function(x) {
																x$rna_alt_expression = sapply(seq(1, nrow(x)),
																															function(y){
																																if (is.na(x$rna_ref_read_count[y]) | is.na(x$rna_alt_read_count[y]) | is.na(x$rna_total_read_count[y])) return(NA)
																																if (x$rna_ref_read_count[y] >= 5
																																		| x$rna_alt_read_count[y] >= 5
																																	# x$rna_total_read_count[y] >= 5 &
																																	# x$rna_total_read_count[y] >= summary(unique(x = x,
																																	#																							by = c('chromosome', 'start_position'))$rna_total_read_count,
																																	#																			 na.rm = T)[["1st Qu."]]
																																)
																																{
																																	if (x$rna_alt_read_count[y] < 1) {
																																		return(FALSE)
																																	} else {
																																		return(TRUE)
																																	}
																																} else {
																																	return(NA)
																																}
																															},
																															USE.NAMES = F)
																return(x)
															})

	input_pileup_merge = lapply(input_pileup_merge,
															function(x) {
																order = c('variant_id', 'chromosome', 'start_position', 'end_position', 'variant_strand', 'ref_allele' , 'alt_allele',
																					'dna_ref_read_count', 'dna_alt_read_count', 'dna_vaf', 'rna_ref_read_count', 'rna_alt_read_count', 'rna_total_read_count', 'rna_vaf', 'rna_alt_expression',
																					'gene_id', 'transcript_id', 'transcript_strand')
																setcolorder(x = x,
																						neworder = c(order, names(x)[-match(x = order, table = names(x))]))
																return(x)
															})

	invisible(mapply(FUN =
									 	function(x, y) {
									 		write.table(x = x,
									 								file = file.path(rootDirectory, '3_neolution', names(input_data)[y]),
									 								sep = '\t',
									 								row.names = F)
									 	},
									 input_pileup_merge,
									 seq(1, length(input_pileup_merge))))
}

performSamtoolsPileup = function(locations_file = NULL, bam_file) {
	dir.create(file.path(rootDirectory, '1b_rnaseq_data', 'pileups'), showWarnings = FALSE)

	system(command = paste('samtools mpileup',
												 ifelse(test = is.null(locations_file),
												 			 yes = '',
												 			 no = paste('-l', locations_file)),
												 bam_file, '>',
												 file.path(rootDirectory, '1b_rnaseq_data', 'pileups', paste0(sub('\\.[^.]*$', '', basename(bam_file)),
												 																														 ifelse(test = is.null(locations_file),
												 																														 			 yes = '_mpil.tsv',
												 																														 			 no = '_mpil_loc.tsv')))),
				 intern = FALSE,
				 wait = TRUE)

	Sys.sleep(time = 1)
}


runSnpEff = function(vcf_path = file.path(rootDirectory, '1a_variants', 'vcf')) {
	#registerDoMC(2)

	message('Step 4: Running snpEff')
	dir.create(path = file.path(rootDirectory, '4_snpEff'),
						 showWarnings = FALSE)

	# make list of input files
	variantLists = list.files(path = vcf_path,
														pattern = '\\.vcf$',
														recursive = FALSE,
														full.names = TRUE)

	invisible(foreach(i = 1:length(variantLists)) %do% {
		command = paste('java -Xmx4g -jar', file.path(runOptions$snpeff$path, 'SnpSift.jar'), 'filter " (ID "\'!\'"~ \'[gr]s*\') "', variantLists[i],
										'| java -Xmx4g -jar', file.path(runOptions$snpeff$path, 'snpEff.jar'), '-nodownload -canon -stats',
										file.path(rootDirectory,'4_snpEff', paste0(sub('\\.[^.]*$','',basename(variantLists[i])), '-snpEff_summary.html')),
										runOptions$snpeff$build, '>', file.path(rootDirectory,'4_snpEff', paste0(sub('\\.[^.]*$','',basename(variantLists[i])), '-snpEff.vcf')))
		system(command = command,
					 wait = TRUE)
	})
}
