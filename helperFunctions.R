# load required packages
if (!require("pacman")) install.packages("pacman")

required_packages = c('data.table', 'gtools', 'utils', 'optparse', 'RMySQL', 'compiler', 'pbapply', 'naturalsort', 'parallel', 'doMC', 'stringr', 'tidyr', 'foreach', 'pander')

library(pacman)
pacman::p_load(char = required_packages)

# drop NA from vector
dropNa = function(vector) { vector[!is.na(vector)] }

# helper functions for input file generation
parseAndExtractFieldsFromVcf = function(vcf_path = file.path(rootDirectory, '1a_variants', 'vcf'), extract_fields = NULL, write = TRUE) {
  # extract relevant info from VCF
  message('Step 1: Parsing & extracting fields from VCF')

  # make list of VCF files
  vcf_files = list.files(path = vcf_path,
                         pattern = '\\.vcf$',
                         recursive = FALSE,
                         full.names = TRUE)

  vcf_files = setNames(object = vcf_files,
                       nm = sub("[.][^.]*$", "",
                                list.files(path = vcf_path,
                                           pattern = '\\.vcf$')))

  if (length(vcf_files) < 1) stop('No VCF files found in ', vcf_path)

  vcf_data = pblapply(vcf_files,
  										function(file_path) parseVcf(vcf_path = file_path,
  																								 sample_tag = 'TUMOR',
  																								 extract_fields = extract_fields))

  vcf_data = lapply(vcf_data,
                    function(vcf) {
                      data = vcf[grepl('^[gr]s\\d+$', variant_id) # always include SNPs
                      					 | sapply(1:nrow(vcf), function(index) nchar(vcf$ref_allele[index]) != nchar(vcf$alt_allele[index])) # include all indels, since sample order is inconsistent (can't be sure we're looking at NORMAL or TUMOR sample)
                      					 | genotype != '0/0', # exclude tumor-specific variants which are ref
                                 .(variant_id, chromosome, start_position, ref_allele, alt_allele, dna_ref_read_count, dna_alt_read_count, dna_total_read_count, dna_vaf)]
                      return(data)
                    })

  if (write) {
  	dir.create(file.path(rootDirectory, '1a_variants', 'parsed'), showWarnings = F)

  	invisible(mapply(function(vcf, index) {
  		write.table(x = vcf,
  								file = file.path(rootDirectory, '1a_variants', 'parsed', paste0(names(vcf_data)[index], '.tsv')),
  								sep = '\t',
  								append = FALSE,
  								quote = FALSE,
  								row.names = FALSE)
  	},
  	vcf_data,
  	seq(1, length(vcf_data))))
  } else {
  	return(vcf_data)
  }

}

parseVcf = function(vcf_path, sample_tag, extract_fields = NULL) {
	require(data.table)
	require(stringr)

	if (!file.exists(vcf_path)) stop('VCF not found')

	# read data
	vcf_lines = readLines(vcf_path)
	# find table header line & read table
	vcf_dt = fread(input = vcf_path,
								 skip = max(grep(pattern = '^#CHROM', x = vcf_lines)) - 1,
								 sep = '\t',
								 colClasses = list(character = c('#CHROM')))
	setnames(x = vcf_dt,
					 old = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'),
					 new = c('chromosome', 'start_position', 'variant_id', 'ref_allele', 'alt_allele', 'quality', 'filter', 'info', 'format'))

	vcf_parsed = vcf_dt[, .(chromosome, start_position, variant_id, ref_allele, alt_allele)]

	# replace 'chr' prefix in chromosome, if present
	vcf_parsed[, chromosome := gsub('^chr', '', chromosome)]

	# extract various tags from sample
	vcf_parsed[, genotype := extractDataFromVcfField(vcf_table = vcf_dt,
																									 sample_tag = sample_tag,
																									 format_tag = 'GT')]

	vcf_parsed[, genotype_quality := as.numeric(extractDataFromVcfField(vcf_table = vcf_dt,
																																			sample_tag = sample_tag,
																																			format_tag = 'GQ'))]

	vcf_parsed[, variant_allele_quality := as.numeric(extractDataFromVcfField(vcf_table = vcf_dt,
																																						sample_tag = sample_tag,
																																						format_tag = 'VAQ'))]

	# vcf_parsed[, average_base_quality := extractDataFromVcfField(vcf_table = vcf_dt,
	# 																														 sample_tag = sample_tag,
	# 																														 format_tag = 'BQ')]

	vcf_parsed[, average_mapping_quality := as.numeric(extractDataFromVcfField(vcf_table = vcf_dt,
																																						 sample_tag = sample_tag,
																																						 format_tag = 'MQ'))]

	# extract sample base read counts
	vcf_parsed = cbind(vcf_parsed, extractVariantReadCountsFromVcf(vcf_table = vcf_dt,
																																 sample_tag = sample_tag,
																																 count_tag = 'BCOUNT'))

	# pick alt allele with highest read count (in case two alt's are given)
	vcf_parsed[, alt_allele := sapply(seq(1, nrow(vcf_parsed)),
																		function(row_index) {
																			if (any(is.na(vcf_parsed[, .(A,C,G,T)][row_index])) | nchar(vcf_parsed$alt_allele[row_index]) == 1) {
																				return(vcf_parsed$alt_allele[row_index])
																			} else {
																				alts = unlist(str_split(vcf_parsed$alt_allele[row_index], ','))
																				counts = sapply(alts,
																												function(alt) {
																													count = vcf_parsed[[alt]][row_index]
																												})
																				base = names(which.max(counts))
																				return(base)
																			}
																		})]

	# extract ref and alt read counts
	dna_read_counts = extractVariantSupportingReadCountsFromVcf(vcf_table = vcf_dt,
																															sample_tag = sample_tag)

	vcf_parsed[, dna_ref_read_count := dna_read_counts[[1]]]
	vcf_parsed[, dna_alt_read_count := dna_read_counts[[2]]]

	# vcf_parsed[, dna_ref_read_count := sapply(seq(1, nrow(vcf_parsed)),
	# 																										 function(row_index) {
	# 																										 	if (any(is.na(vcf_parsed[, .(A,C,G,T)][row_index])) | nchar(vcf_parsed$ref_allele[row_index]) != 1) {
	# 																										 		return(NA)
	# 																										 	} else {
	# 																										 		return(vcf_parsed[[vcf_parsed$ref_allele[row_index]]][row_index])
	# 																										 	}
	# 																										 })]
	# vcf_parsed[, dna_alt_read_count := sapply(seq(1, nrow(vcf_parsed)),
	# 																										 function(row_index) {
	# 																										 	if (any(is.na(vcf_parsed[, .(A,C,G,T)][row_index])) | nchar(vcf_parsed$alt_allele[row_index]) != 1) {
	# 																										 		return(NA)
	# 																										 	} else {
	# 																										 		return(vcf_parsed[[vcf_parsed$alt_allele[row_index]]][row_index])
	# 																										 	}
	# 																										 })]

	# extract total read counts
	vcf_parsed[, dna_total_read_count := as.numeric(extractDataFromVcfField(vcf_table = vcf_dt,
																																					sample_tag = sample_tag,
																																					format_tag = 'DP'))]

	# compute dna_vaf
	vcf_parsed[, dna_vaf := sapply(seq(1, nrow(vcf_parsed)),
																 function(row_index) {
																 	if (is.numeric(vcf_parsed$dna_alt_read_count[row_index]) & is.numeric(vcf_parsed$dna_total_read_count[row_index])) {
																 		return(vcf_parsed$dna_alt_read_count[row_index] / vcf_parsed$dna_total_read_count[row_index])
																 	} else {
																 		return(NA)
																 	}
																 })]

	return(vcf_parsed)
}

extractDataFromVcfField = function(vcf_table, sample_tag, format_tag) {
	# if sample_tag not found in columns, return NA
	if (!any(grepl(pattern = sample_tag, x = names(vcf_table)))) {
		return(NA)
	}

	tag_positions = unlist(sapply(str_split(vcf_table$format, pattern = ':'),
																function(x) {
																	match(x = format_tag, table = x, nomatch = NA)
																}))

	vcf_tumor_split = str_split(vcf_table[[sample_tag]], pattern = ':')

	tag_data = sapply(seq(1, length(tag_positions)),
										function(x) {
											if (grepl(pattern = '^[gr]s\\d+$', x = vcf_table$variant_id[x]) & is.na(tag_positions[[x]])) {
												snp_split = unlist(str_split(string = vcf_table$info[x], pattern = ';'))
												tag_data = grep(pattern = paste0('^', format_tag, '='), x = snp_split, value = TRUE)

												if (length(tag_data) < 1) return(NA)

												return(sub(pattern = paste0('^', format_tag, '='), replacement = '', x = tag_data))
											}

											vcf_tumor_split[[x]][tag_positions[x]]
										})

	return(tag_data)
}

extractVariantReadCountsFromVcf = function(vcf_table, sample_tag, count_tag) {
	# if sample_tag not found in columns, return NA
	if (!any(grepl(pattern = sample_tag, x = names(vcf_table)))) {
		return(NA)
	}

	count_data = extractDataFromVcfField(vcf_table = vcf_table,
																			 sample_tag = sample_tag,
																			 format_tag = count_tag)

	count_data = lapply(str_split(string = count_data, pattern = ','), as.numeric)

	count_data = lapply(count_data,
											function(counts) {
												if (any(!is.na(counts))) {
													data = setNames(object = counts, nm = c('A', 'C', 'G', 'T'))
													data = as.data.table(as.list(data))
												} else {
													data = as.data.table(as.list(rep(x = NA, 4)))
													data = setNames(object = data, nm = c('A', 'C', 'G', 'T'))
												}
												return(data)
											})

	count_data = rbindlist(count_data, use.names = TRUE)
}

extractVariantSupportingReadCountsFromVcf = function(vcf_table, sample_tag) {

	##### sample order is not consistent for somaticIndelCaller (seems alphabetically ordered by filename, instead of NORMAL-TUMOR)
	##### end result is we don't know which column is TUMOR, therefore don't process indels at this moment in time

	# # SC field: "counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads"
	# count_data_indels = extractDataFromVcfField(vcf_table = vcf_table,
	# 																						sample_tag = sample_tag,
	# 																						format_tag = 'SC')
	#
	# # reverse 'counts_data_indels' so that order is identical to 'counts_data_snvs'
	# count_data_indels = lapply(str_split(string = count_data_indels, pattern = ','),
	# 													 function(counts) {
	# 													 	as.numeric(counts[c(3, 4, 1, 2)])
	# 													 })

	# if sample_tag not found in columns, return NA
	if (!any(grepl(pattern = sample_tag, x = names(vcf_table)))) {
		return(list(NA, NA))
	}

	# DP4 field: "# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases"
	count_data_snvs = extractDataFromVcfField(vcf_table = vcf_table,
																						sample_tag = sample_tag,
																						format_tag = 'DP4')

	# pre-process 'count_data_snvs' for merge
	count_data_snvs = lapply(str_split(string = count_data_snvs, pattern = ','),
													 as.numeric)

	# # merge count data from indels and snvs
	# count_data = lapply(seq(1, length(count_data_snvs)),
	# 										function(index) {
	# 											if (any(is.na(count_data_snvs[[index]]))) return(count_data_indels[[index]]) else return(count_data_snvs[[index]])
	# 										})

	count_data = count_data_snvs

	ref_supporting_read_counts = sapply(count_data, function(counts) sum(counts[c(1,2)]))
	alt_supporting_read_counts = sapply(count_data, function(counts) sum(counts[c(3,4)]))

	return(list(ref_supporting_read_counts, alt_supporting_read_counts))
}

performVarcontextGeneration = function(variant_path = file.path(rootDirectory, '1a_variants', 'parsed'), filter_rna_alt_expression = TRUE, vcf_fields = c('ID', 'CHROM', 'POS', 'REF', 'ALT')) {
	registerDoMC(runOptions$varcontext$numberOfWorkers)

	dir.create(file.path(rootDirectory, '2_varcontext'),
						 showWarnings = FALSE)
	dir.create(file.path(rootDirectory, '2_varcontext', 'input_lists'),
						 showWarnings = FALSE)
	dir.create(file.path(rootDirectory, '2_varcontext', 'varcontext_logs'),
						 showWarnings = FALSE)

	# make list of input files
	variant_lists = list.files(path = variant_path,
														pattern = '\\.tsv$',
														recursive = FALSE,
														full.names = TRUE)
	variant_data = lapply(variant_lists,
											 fread,
											 colClasses = list(character = c('chromosome')))
	variant_data = setNames(object = variant_data,
												 nm = list.files(path = variant_path,
												 								pattern = '\\.tsv$'))

	# remove variants without rna_alt_expression
	if (filter_rna_alt_expression) {
		variant_data = lapply(variant_data,
													function(variants) {
														return(variants[rna_alt_expression == TRUE | is.na(rna_alt_expression)])
													})
	}

	invisible(mapply(function(x, y) {
		write.table(x = x,
								file = file.path(rootDirectory, '2_varcontext', 'input_lists', names(variant_data)[y]),
								sep = '\t',
								row.names = F)
	},
	variant_data,
	seq(1, length(variant_data))))

	# generate variant contexts
	message('Step 2: Generating context for variants')
	message('Using gene build: ', runOptions$varcontext$ensemblApi)

	generateVarcontext(input_list = list.files(path = file.path(rootDirectory, '2_varcontext', 'input_lists'),
																						 pattern = '\\.tsv',
																						 full.names = TRUE))

	# move files
	# system(command = paste('cd', variant_path, ';', 'mv -f *.tsv ../extr_fields'))
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
																																																		 x = names(varcontext_data),
																																																		 value = T),
																											 USE.NAMES = FALSE),
																		 rna_expression_data = sapply(sample_info$rna_data_prefix, function(x) grep(pattern = x,
																		 																																					 x = names(rnaseq_data),
																		 																																					 value = T),
																		 														 USE.NAMES = FALSE)
		)

		if (nrow(sample_combinations) < 1) stop('Input files not found, check directory structure and/or sample_info.tsv')

		prediction_input = lapply(seq(1, nrow(sample_combinations)),
															function(x) {
																if (length(sample_combinations[x, variants]) > 0 & length(sample_combinations[x, rna_expression_data]) > 0) {
																	message('Merging RNA expression data for: ', sample_combinations[x, variants], ' & ', sample_combinations[x, rna_expression_data])
																	mergeByEnsemblId(variant_table = varcontext_data[[sample_combinations[x, variants]]],
																									 expression_table = rnaseq_data[[sample_combinations[x, rna_expression_data]]])
																} else if (length(sample_combinations[x, variants]) > 0) {
																	message('Adding empty rna_expression column for ', sample_combinations[x, variants])
																	mergeByEnsemblId(variant_table = varcontext_data[[sample_combinations[x, variants]]],
																									 expression_table = NULL)
																} else {
																	warning('No input data present for ', sample_combinations[x, variants], ' & ', sample_combinations[x, rna_expression_data])
																}
															})
		prediction_input = setNames(object = prediction_input, nm = sub("[.][^.]*$", "", sample_combinations$variants))

		rna_coverage_summary = data.table(sample = names(prediction_input),
																			percent_ensg_coverage = sapply(prediction_input, function(x) round(x = length(which(!is.na(x[[expression_unit]]))) / nrow(x) * 100,
																																																				 digits = 1))
		)

		rna_expression_summary = data.table(sample = names(prediction_input),
																				percent_no_expression = sapply(rnaseq_data, function(x) round(x = length(which(x[[expression_unit]] == 0)) / nrow(x) * 100,
																																																			digits = 1))
		)

		write.table(x = rna_expression_summary,
								file = file.path(rootDirectory, '1b_rnaseq_data', 'rna_expression_info.tsv'),
								sep = '\t',
								row.names = FALSE,
								quote = FALSE,
								append = FALSE)

		write.table(x = rna_coverage_summary,
								file = file.path(rootDirectory, '3_neolution', 'rna_coverage_info.tsv'),
								sep = '\t',
								row.names = FALSE,
								quote = FALSE,
								append = FALSE)

		message('Step 3b: Generating Neolution pipeline input')
	} else {
		message('Step 3: Generating Neolution pipeline input (no RNAseq data)')

		prediction_input = varcontext_data
	}

	invisible(mapply(FUN =
									 	function(x, y) {
									 		write.table(x = x,
									 								file = file.path(rootDirectory, '3_neolution', paste0(names(prediction_input)[y], '.tsv')),
									 								sep = '\t',
									 								row.names = F)
									 	},
									 prediction_input,
									 seq(1, length(prediction_input))))
}

mergeByEnsemblId = function(variant_table, expression_table, expression_unit = 'FPKM') {
	if (is.null(expression_table)) {
		variant_table[[expression_unit]] = NA
		return(variant_table)
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

findRnaReadLevelEvidenceForVariants = function(vcf_input_path = file.path(rootDirectory, '1a_variants', 'parsed'),
																							 rna_path = file.path(rootDirectory, '1b_rnaseq_data', 'bam'),
																							 fasta_genome_ref = runOptions$samtools$fastaGenomeRef,
																							 sample_info_path = file.path(rootDirectory, 'sample_info.tsv')) {

	registerDoMC(runOptions$samtools$numberOfWorkers)

	if (!dir.exists(rna_path)) {
		message('RNAseq data directory does not exist at ', rna_path)
		user_choice = menu(choices = c('yes', 'no'), title = 'Do you want to continue without RNAseq data?')

		if (user_choice == 2) {
			q()
		} else if (user_choice == 1) {
			# insert code to add empty rna_* columns to tsv files
		} else {
			stop('Please enter a valid option (1 or 2)')
		}
	}

	# parse input data
	input_data = lapply(list.files(path = vcf_input_path,
																 pattern = '\\.tsv',
																 full.names = TRUE),
											fread,
											colClasses = list(character = c('chromosome')))
	input_data = setNames(object = input_data,
												nm = list.files(path = vcf_input_path,
																				pattern = '\\.tsv'))

	# load sample info
	if (file.exists(sample_info_path)) {
		sample_info = fread(sample_info_path, sep = '\t', header = T, na.strings = c('','NA', 'N.A.'))
	} else if (sample_info$patient_id[1] == 'place_holder') {
		stop('Please fill in sample_info.tsv')
	} else {
		stop('Sample info file missing at "', sample_info_path, '", please provide correct path in argument to findRnaReadLevelEvidenceForVariants')
	}

	# make list of all unique variants found in parsed VCF
	snv_positions = lapply(input_data,
												 function(x) {
												 	setorder(x = unique(x[grepl(pattern = '^[0-9]{1,2}$|^[XY]$',
												 															x = x$chromosome) &
												 													!grepl(pattern = '^[rg]s\\d+$',
												 																 x = x$variant_id) &
												 													nchar(x$ref_allele) == 1 &
												 													nchar(x$alt_allele) == 1,
												 												.(chromosome, start_position), ]),
												 					 chromosome,
												 					 start_position)
												 })

	# generate position list files to feed to samtools mpileup
	dir.create(path = file.path(rootDirectory, '1a_variants', 'poslist'), showWarnings = F)

	invisible(sapply(1:length(snv_positions),
									 function(x) write.table(x = snv_positions[[x]],
									 												file = file.path(rootDirectory, '1a_variants', 'poslist', paste0(sub("\\.[^.]*$", "", names(snv_positions)[x]), '_poslist.tsv')),
									 												quote = FALSE,
									 												col.names = FALSE,
									 												row.names = FALSE)))

	# make overview of sample input files
	sample_combinations = data.table(locations_file = sapply(sample_info$dna_data_prefix, function(x) grep(pattern = x,
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

	if (nrow(sample_combinations) < 1) {
		stop('No valid sample combinations found, please check dna/rna_prefixes and rna bamfile filenames')
	}

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
	# invisible(sapply(seq(1, nrow(sample_combinations)),
	# 								 function(x) {
	# 								 	if (!file.exists(file.path(rootDirectory, '1b_rnaseq_data', 'pileups', paste0(sub('[.][^.]*$', '', basename(sample_combinations$rna_bam_file[x])), '_mpil_loc.tsv')))) {
	# 								 		performSamtoolsPileup(locations_file = sample_combinations$variants[x], bam_file = sample_combinations$rna_bam_file[x])
	# 								 	}})
	# 					)

	x = foreach(i = 1:nrow(sample_combinations)) %dopar% {
		if (!file.exists(file.path(rootDirectory, '1b_rnaseq_data', 'pileups', paste0(sub('\\.[^.]*$', '', basename(sample_combinations$rna_bam_file[i])),
																																									ifelse(test = is.null(sample_combinations$locations_file[i]),
																																												 yes = '_mpil.tsv',
																																												 no = '_mpil_loc.tsv'))))) {
			if (is.null(fasta_genome_ref) | !file.exists(fasta_genome_ref)) {
				stop('Please provide valid fasta DNA reference (set location in runConfig.R)')
			} else {
				message('Performing pileup with genomic reference: ', fasta_genome_ref, '\n')
			}

			performSamtoolsPileup(bam_file = sample_combinations$rna_bam_file[i],
														locations_file = sample_combinations$locations_file[i],
														fasta_reference = fasta_genome_ref)
		}
	}

	pileup_loc_data = lapply(list.files(path = file.path(rootDirectory, '1b_rnaseq_data', 'pileups'),
																			pattern = '_mpil(_loc){0,1}\\.tsv',
																			full.names = TRUE),
													 fread,
													 colClasses = c(V1 = 'character'),
													 col.names = c('chromosome', 'start_position', 'ref_base', 'number_of_reads', 'rna_read_bases', 'base_quality'))

	pileup_loc_data = setNames(object = pileup_loc_data, nm = list.files(path = file.path(rootDirectory, '1b_rnaseq_data', 'pileups'),
																																			 pattern = '_mpil(_loc){0,1}\\.tsv',
																																			 full.names = FALSE))

	input_pileup_merge = lapply(1:length(input_data),
															function(index_input_data) {
																index_prefix_dna = which(sapply(sample_info$dna_data_prefix, function(y) grepl(pattern = y, x = names(input_data)[index_input_data], fixed = T), USE.NAMES = F))
																if (index_prefix_dna > 0) {
																	index_pileup_data = which(grepl(pattern = sample_info$rna_data_prefix[index_prefix_dna], x = names(pileup_loc_data), fixed = T))
																}
																if (index_pileup_data > 0) {
																	merged_data = merge(x = input_data[[index_input_data]],
																											y = pileup_loc_data[[index_pileup_data]][, .(chromosome, start_position, ref_base, rna_read_bases)],
																											by = c('chromosome', 'start_position'),
																											all.x = TRUE)
																	setorder(merged_data, variant_id, chromosome, start_position)
																}})

	input_pileup_merge = lapply(input_pileup_merge,
															function(dt) {
																dt$rna_read_bases = gsub(pattern = '[^atgcATGC,\\.]',
																												 replacement = '',
																												 x = toupper(dt$rna_read_bases))
																dt$rna_read_bases = ifelse(test = is.na(dt$rna_read_bases) | dt$rna_read_bases == '',
																													 yes = NA,
																													 no = dt$rna_read_bases)

																dt$rna_ref_read_count = str_count(string = dt$rna_read_bases,
																																	pattern = '[\\.,]')

																dt$rna_alt_read_count = mapply(FUN = function(ref_allele, alt_allele, ref_base, read_bases)
																{
																	if (identical(ref_allele, ref_base)) str_count(string = read_bases,
																																								 pattern = alt_allele)
																	else if (identical(chartr(old = 'ATGC',
																														new = 'TACG',
																														x = ref_allele), ref_base)) str_count(string = read_bases,
																																																	pattern = chartr(old = 'ATGC',
																																																									 new = 'TACG',
																																																									 x = alt_allele))
																	else NA
																},
																dt$ref_allele,
																dt$alt_allele,
																dt$ref_base,
																dt$rna_read_bases,
																SIMPLIFY = T)

																dt$rna_total_read_count = ifelse(test = is.na(dt$rna_read_bases),
																																 yes = NA,
																																 no = nchar(dt$rna_read_bases))

																dt$rna_read_bases = NULL
																dt$ref_base = NULL

																return(dt)
															})

	input_pileup_merge = lapply(input_pileup_merge,
															function(dt) {
																dt$rna_alt_expression = sapply(1:nrow(dt),
																															 function(row_index){
																															 	if (is.na(dt$rna_ref_read_count[row_index]) | is.na(dt$rna_alt_read_count[row_index]) | is.na(dt$rna_total_read_count[row_index])) return(NA)
																															 	if (dt$rna_ref_read_count[row_index] >= 5
																															 			| dt$rna_alt_read_count[row_index] >= 5
																															 	) {
																															 		if (dt$rna_alt_read_count[row_index] > 0) {
																															 			return(TRUE)
																															 		} else if (dt$rna_alt_read_count[row_index] < 1) {
																															 			return(FALSE)
																															 		} else {
																															 			return(NA)
																															 		}
																															 	} else {
																															 		return(NA)
																															 	}
																															 },
																															 USE.NAMES = F)

																dt$rna_vaf = sapply(1:nrow(dt),
																										function(row_index){
																											if (is.na(dt$rna_alt_read_count[row_index]) | is.na(dt$rna_total_read_count[row_index])) return(NA)
																											if (dt$rna_ref_read_count[row_index] >= 5
																													| dt$rna_alt_read_count[row_index] >= 5
																											) {
																												return(dt$rna_alt_read_count[row_index] / dt$rna_total_read_count[row_index])
																											} else {
																												return(NA)
																											}
																										},
																										USE.NAMES = F)
																return(dt)
															})

	input_pileup_merge = lapply(input_pileup_merge,
															function(dt) {
																order = c('variant_id', 'chromosome', 'start_position', 'ref_allele' , 'alt_allele',
																					'dna_ref_read_count', 'dna_alt_read_count', 'dna_total_read_count', 'dna_vaf', 'rna_ref_read_count', 'rna_alt_read_count', 'rna_total_read_count', 'rna_vaf', 'rna_alt_expression')
																setcolorder(x = dt,
																						neworder = order)
																return(dt)
															})

	invisible(mapply(FUN =
									 	function(x, y) {
									 		write.table(x = x,
									 								file = file.path(rootDirectory, '1a_variants', 'parsed', names(input_data)[y]),
									 								sep = '\t',
									 								row.names = F)
									 	},
									 input_pileup_merge,
									 seq(1, length(input_pileup_merge))))
}

performSamtoolsPileup = function(bam_file, locations_file = NULL, fasta_reference = NULL) {
	dir.create(file.path(rootDirectory, '1b_rnaseq_data', 'pileups'), showWarnings = FALSE)

	system(command = paste('samtools',
												 'mpileup',
												 ifelse(test = is.null(locations_file),
												 			 yes = '',
												 			 no = paste('-l', locations_file)),
												 ifelse(test = is.null(fasta_reference),
												 			 yes = '',
												 			 no = paste('-f', fasta_reference)),
												 bam_file, '>',
												 file.path(rootDirectory, '1b_rnaseq_data', 'pileups', paste0(sub('\\.[^.]*$', '', basename(bam_file)),
												 																														 ifelse(test = is.null(locations_file),
												 																														 			 yes = '_mpil.tsv',
												 																														 			 no = '_mpil_loc.tsv')))),
				 intern = FALSE,
				 wait = TRUE)

	Sys.sleep(time = 1)
}


runSnpEff = function(vcf_path = file.path(rootDirectory, '1a_variants', 'vcf'), filter_snps = TRUE) {
	#registerDoMC(2)

	message('Step 4: Running snpEff')
	message('Using genome & gene builds: ', runOptions$snpeff$build)
	dir.create(path = file.path(rootDirectory, '4_snpEff'),
						 showWarnings = FALSE)

	# make list of input files
	variantLists = list.files(path = vcf_path,
														pattern = '\\.vcf$',
														recursive = FALSE,
														full.names = TRUE)

	invisible(foreach(i = 1:length(variantLists)) %do% {
		command_snpsift = paste('java -Xmx4g -jar', file.path(runOptions$snpeff$path, 'SnpSift.jar'),
														'filter " (ID "\'!\'"~ \'[gr]s*\') "',
														variantLists[i])
		command_snpeff = paste('java -Xmx4g -jar', file.path(runOptions$snpeff$path, 'snpEff.jar'), '-nodownload -canon -stats',
													 file.path(rootDirectory,'4_snpEff', paste0(sub('\\.[^.]*$','',basename(variantLists[i])), '-snpEff_summary.html')),
													 runOptions$snpeff$build, '>', file.path(rootDirectory,'4_snpEff', paste0(sub(pattern = '\\.[^.]*$',
													 																																						 replacement = '',
													 																																						 x =  basename(variantLists[i])),
													 																																				 ifelse(test = filter_snps,
													 																																				 			 yes = '',
													 																																				 			 no = '-with_snps'),
													 																																				 '-snpEff.vcf')))
		command = ifelse(test = filter_snps,
										 yes = paste(command_snpsift, command_snpeff, sep = '|'),
										 no = paste(paste('cat', variantLists[i]), command_snpeff, sep = '|'))

		system(command = command,
					 wait = TRUE)
	})
}

# helper functions for analysis document
"%nin%" = Negate("%in%")

parseEpitopePredictions = function(path, pattern = '_epitopes\\.csv') {
	require(stringr)
	require(data.table)

	# get predictions paths
	files = list.files(path = path,
										 pattern = pattern,
										 full.names = TRUE)

	# extract filenames
	short_names = sub(pattern = '.+[/]', replacement = '', x = files)
	short_names = sub(pattern = '_epitopes\\.csv', replacement = '', x = short_names)

	# get data in & sort by tumor peptide affinity
	predictions = lapply(seq(1, length(files)),
											 function(i) {
											 	data = fread(files[i])
											 	data[, sample_prefix := sub(pattern = '_mg.+|_L[0-9]{3}.+', replacement = '', x = short_names[i])]

											 	return(data)
											 })

	sapply(predictions,
				 function(x) setkeyv(x = x, grep(pattern = 'tumor_.+affinity', x = names(x), value = TRUE)))

	# set table names
	predictions = setNames(object = predictions,
												 nm = short_names)

	return(predictions)
}

applyCutoffs = function(predictions, rank = 3, processing = 0.5, expression = 0, selfsim = TRUE) {
	# take subsets
	subsets = lapply(seq(1, length(predictions), 1),
									 function(x) {
									 	data_subset = subset(x = predictions[[x]],
									 											 subset = predictions[[x]][[grep(pattern = 'tumor.+rank',
									 											 																x = colnames(predictions[[x]]),
									 											 																value = T)]] <= rank)
									 	data_subset = data_subset[tumor_processing_score >= processing & (rna_expression > expression | is.na(rna_expression)) & different_from_self == selfsim]
									 	# data_subset = data_subset[rna_expression > expression]
									 	# data_subset = data_subset[different_from_self == selfsim]
									 	return(data_subset)
									 })

	# set table names
	subsets = setNames(object = subsets,
										 nm = names(predictions))

	return(subsets)
}

prepareEpitopeLists = function(list_of_predictions, split_by= c('9mer', '10mer', '11mer')) {
	# split all predictions by xmer
	split_by_xmer = lapply(split_by,
												 function(x) {
												 	join_by_xmer = rbindlist(list_of_predictions[grepl(pattern = x,
												 																										 x = names(list_of_predictions))])
												 	subset_by_xmer = subset(x = join_by_xmer,
												 													select = c('sample_prefix', 'hla_allele', 'tumor_peptide'))
												 	subset_by_xmer[, hla_allele := gsub(pattern = "*",
												 																			replacement = "\\*",
												 																			x = subset_by_xmer$hla_allele,
												 																			fixed = TRUE)]
												 	return(subset_by_xmer)
												 })

	# add info required by peptide synth facility (add xmer, n-1 peptides and ni column)
	split_by_xmer_info = lapply(split_by_xmer,
															function(x) {
																xmer = nchar(x$tumor_peptide)
																tumor_peptide_trunc = substr(x = x$tumor_peptide,
																														 start = 1,
																														 stop = nchar(x$tumor_peptide) - 1)
																tumor_peptide_resin = substr(x = x$tumor_peptide,
																														 start = nchar(x$tumor_peptide),
																														 stop = nchar(x$tumor_peptide))
																table = as.data.table(cbind(x, tumor_peptide_trunc, tumor_peptide_resin, xmer))
																setkey(x = table, tumor_peptide_resin, hla_allele)
																return(table)
															})
	return(split_by_xmer_info)
}

