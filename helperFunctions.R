# install/load required packages
if (!require('pacman')) install.packages('pacman')
library(pacman)

required_packages = c('BSgenome', 'compiler', 'cowplot', 'data.table', 'deconstructSigs', 'doMC', 'dplyr', 'dtplyr', 'foreach', 'ggplot2', 'gtools', 'koRpus',
                      'naturalsort', 'optparse', 'pander', 'parallel', 'pbapply', 'rtracklayer', 'stringr', 'tidyr', 'utils')

if (!p_isinstalled('BSgenome')) {
  source('https://bioconductor.org/biocLite.R')
  biocLite('BSgenome', suppressUpdates = T)
}
if (!p_isinstalled('rtracklayer')) {
  source('https://bioconductor.org/biocLite.R')
  biocLite('rtracklayer', suppressUpdates = T)
}

pacman::p_load(char = required_packages)


## Convenience functions ---------------------------------------------------
"%nin%" = Negate("%in%")

# drop NA from vector
dropNa = function(vector) { vector[!is.na(vector)] }

# wrapper for executing commands (or not)
commandWrapper = function(command, nice = 19, wait = TRUE, execute) {
  if (is.numeric(nice)) {command = paste('nice -n', nice, command)}

  if (!execute) {command = paste('nohup', command, '2> nohup.out &\n')}

  if (execute) {
    system(command = command,
           intern = FALSE,
           wait = wait)
  } else {
    message(command)
  }
}

# parse commandline arguments, return list of arguments and values
parseCommandlineArguments = function() {
  cmdln_args = commandArgs(TRUE)
  if (length(cmdln_args) < 1) {return()}

  if (any(!grepl('^-{1,2}', cmdln_args))) {message('Commandline arguments must begin with single ("-") or double ("--") dashes\nStopping execution'); q()}
  if (any(!grepl('[A-Za-z]+=[A-Za-z0-9.,]+', cmdln_args))) {message('Commandline arguments and values must be separated by an equal sign ("=")\nStopping execution'); q()}

  cmdln_args = gsub('^-{1,2}', '', cmdln_args)
  cmdln_args = strsplit(cmdln_args, '=', TRUE)

  parsed_args = lapply(cmdln_args,"[[",2 )
  names(parsed_args) = lapply(cmdln_args, "[[", 1)

  return(parsed_args)
}


## VCF parsing ---------------------------------------------------
# helper functions for input file generation
parseAndExtractFieldsFromVcf = function(vcf_path = file.path(rootDirectory, '1a_variants', 'vcf'), vcf_regex = '\\.vcf$', normal_tag = 'NORMAL', tumor_tag = 'TUMOR', extract_fields = NULL, write = TRUE) {
  # extract relevant info from VCF
  message('Step 1: Parsing & extracting fields from VCF')

  # make list of VCF files
  vcf_files = list.files(path = vcf_path,
                         pattern = vcf_regex,
                         recursive = FALSE,
                         full.names = TRUE)

  vcf_files = setNames(object = vcf_files,
                       nm = sub(regexPatterns$file_extension, '',
                                list.files(path = vcf_path,
                                           pattern = vcf_regex,
                                           recursive = FALSE)))

  if (length(vcf_files) < 1) stop('No VCF files found in ', vcf_path)

  progress_bar = startpb(min = 0, max = length(vcf_files))
  on.exit(closepb(progress_bar))

  vcf_data = mapply(function(file, n, t, i) {
    data = parseVcf(vcf_path = file,
                    n_tag = n,
                    t_tag = t,
                    extract_fields = extract_fields)
    setpb(progress_bar, i)
    return(data)
  },
  vcf_files,
  normal_tag,
  tumor_tag,
  1:length(vcf_files),
  SIMPLIFY = F)

  vcf_data = lapply(vcf_data,
                    function(vcf) {
                      # subset data
                      data = vcf[grepl(regexPatterns$snp_identifier, variant_id) # always include SNPs
                                 | sapply(1:nrow(vcf), function(index) nchar(vcf$ref_allele[index]) != nchar(vcf$alt_allele[index])) # legacy (somSniper pipeline) support: include all indels, since sample order (NORMAL/TUMOR) is inconsistent for indel calls
                                 | genotype != '0/0' | is.na(genotype), # safety check, shouldn't happen: exclude tumor-specific variants which are ref
                                 .(variant_id, chromosome, start_position, ref_allele, alt_allele, dna_ref_read_count, dna_alt_read_count, dna_total_read_count, dna_vaf)]

                      # sort data
                      data_sorted = rbindlist(list(data %>%
                                                     filter(!grepl(regexPatterns$gs_identifier, variant_id)) %>%
                                                     .[naturalorder(.$chromosome)],
                                                   data %>%
                                                     filter(grepl(regexPatterns$gs_identifier, variant_id)) %>%
                                                     .[naturalorder(.$chromosome)])
                                              )

                      return(data_sorted)
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

parseVcf = function(vcf_path, n_tag, t_tag, extract_fields = NULL) {
  require(data.table)
  require(stringr)

  if (!file.exists(vcf_path)) stop('VCF not found')

  # read data
  vcf_lines = readLines(vcf_path)
  # find table header line & read table
  vcf_dt = fread(input = vcf_path,
                 skip = max(grep(pattern = '#CHROM', x = vcf_lines)) - 1,
                 sep = '\t',
                 colClasses = list(character = c('#CHROM')),
                 fill = TRUE)
  setnames(x = vcf_dt,
           old = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'),
           new = c('chromosome', 'start_position', 'variant_id', 'ref_allele', 'alt_allele', 'quality', 'filter', 'info', 'format'))

  if (n_tag %nin% names(vcf_dt) | t_tag %nin% names(vcf_dt)) {
    stop('"',n_tag, '" and/or "', t_tag, '" not found in VCF header, please doublecheck normal & tumor sample tags')
  }

  if (all(vcf_dt[grepl(regexPatterns$gs_identifier, variant_id), n_tag, with = FALSE] == '')) {
    vcf_dt[grepl(regexPatterns$gs_identifier, variant_id), (n_tag) := vcf_dt[grepl(regexPatterns$gs_identifier, variant_id), t_tag, with = FALSE]]
    vcf_dt[grepl(regexPatterns$gs_identifier, variant_id), (t_tag) := NA]
  }

  # check if IndelGenotyper is used: if so, set genotype to NA, as sample_tag order is not consistent and can therefore not be extracted with confidence
  indel_genotyper = 'IndelGenotyperV2' %in% sort(gsub(pattern = '##source=', '', grep(pattern = '##source', x = vcf_lines, value = TRUE)))

  # enumerate tags
  tags = list(filter_tags = sort(gsub(pattern = 'ID=', '', x = str_extract(string = grep(pattern = '##FILTER', x = vcf_lines, value = TRUE),
                                                                           pattern = 'ID=[^,]+'))),
              format_tags = sort(gsub(pattern = 'ID=', '', x = str_extract(string = grep(pattern = '##FORMAT', x = vcf_lines, value = TRUE),
                                                                           pattern = 'ID=[^,]+'))),
              info_tags = sort(gsub(pattern = 'ID=', '', x = str_extract(string = grep(pattern = '##INFO', x = vcf_lines, value = TRUE),
                                                                         pattern = 'ID=[^,]+'))))

  # take subset of full vcf and start parsing/adding data
  vcf_parsed = vcf_dt[, .(chromosome, start_position, variant_id, ref_allele, alt_allele)]

  # replace 'chr' prefix in chromosome, if present
  vcf_parsed[, chromosome := gsub('^chr', '', chromosome)]

  # extract genotype tag info for variants
  if ('GT' %in% tags$format_tags) {
    vcf_parsed[!grepl(regexPatterns$gs_identifier, variant_id), genotype := extractSingleDataFromVcfField(vcf_table = vcf_dt[!grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                          sample_tag = t_tag,
                                                                                                          format_tag = 'GT')]
    vcf_parsed[grepl(regexPatterns$gs_identifier, variant_id), genotype := extractSingleDataFromVcfField(vcf_table = vcf_dt[grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                         sample_tag = n_tag,
                                                                                                         format_tag = 'GT')]
  }

  if (indel_genotyper) {vcf_parsed[which(grepl(pattern = 'SOMATIC', x = vcf_dt$info)), genotype := NA]}

  # extract genotype/variant allele/mapping quality info for variants
  if ('GQ' %in% tags$format_tags) {
    vcf_parsed[!grepl(regexPatterns$gs_identifier, variant_id), genotype_quality := as.numeric(extractSingleDataFromVcfField(vcf_table = vcf_dt[!grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                             sample_tag = t_tag,
                                                                                                                             format_tag = 'GQ'))]
    vcf_parsed[grepl(regexPatterns$gs_identifier, variant_id), genotype_quality := as.numeric(extractSingleDataFromVcfField(vcf_table = vcf_dt[grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                            sample_tag = n_tag,
                                                                                                                            format_tag = 'GQ'))]
  }
  if ('VAQ' %in% tags$format_tags) {
    vcf_parsed[!grepl(regexPatterns$gs_identifier, variant_id), variant_allele_quality := as.numeric(extractSingleDataFromVcfField(vcf_table = vcf_dt[!grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                   sample_tag = t_tag,
                                                                                                                                   format_tag = 'VAQ'))]
    vcf_parsed[grepl(regexPatterns$gs_identifier, variant_id), variant_allele_quality := as.numeric(extractSingleDataFromVcfField(vcf_table = vcf_dt[grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                  sample_tag = n_tag,
                                                                                                                                  format_tag = 'VAQ'))]
  }
  if ('MQ' %in% tags$format_tags) {
    vcf_parsed[!grepl(regexPatterns$gs_identifier, variant_id), average_mapping_quality := as.numeric(extractSingleDataFromVcfField(vcf_table = vcf_dt[!grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                    sample_tag = t_tag,
                                                                                                                                    format_tag = 'MQ'))]
    vcf_parsed[grepl(regexPatterns$gs_identifier, variant_id), average_mapping_quality := as.numeric(extractSingleDataFromVcfField(vcf_table = vcf_dt[grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                   sample_tag = n_tag,
                                                                                                                                   format_tag = 'MQ'))]
  }

  # extract base quality scores
  if ('QSS' %in% tags$format_tags) {
    vcf_parsed[!grepl(regexPatterns$gs_identifier, variant_id), c('ref_base_quality_score', 'alt_base_quality_score') := extractSplitDataFromVcf(vcf_table = vcf_dt[!grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                                 sample_tag = t_tag,
                                                                                                                                                 format_tag = 'QSS')]
    vcf_parsed[grepl(regexPatterns$gs_identifier, variant_id), c('ref_base_quality_score', 'alt_base_quality_score') := extractSplitDataFromVcf(vcf_table = vcf_dt[grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                                sample_tag = n_tag,
                                                                                                                                                format_tag = 'QSS')]
  }

  # extract sample base read counts
  if ('BCOUNT' %in% tags$format_tags) {
    vcf_parsed[!grepl(regexPatterns$gs_identifier, variant_id), c('A', 'C', 'G', 'T') := extractVariantReadCountsFromVcf(vcf_table = vcf_dt[!grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                         sample_tag = t_tag,
                                                                                                                         count_tag = 'BCOUNT')]
    vcf_parsed[grepl(regexPatterns$gs_identifier, variant_id), c('A', 'C', 'G', 'T') := extractVariantReadCountsFromVcf(vcf_table = vcf_dt[grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                        sample_tag = n_tag,
                                                                                                                        count_tag = 'BCOUNT')]
  }

  # pick alt allele with highest read count (in case two alt's are given)
  if ('BCOUNT' %in% tags$format_tags) {
    vcf_parsed[, alt_allele := unlist(mclapply(seq(1, nrow(vcf_parsed)),
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
                                               }, mc.cores = 20))]
  }

  # extract ref and alt read counts
  if ('DP4' %in% tags$format_tags) {
    vcf_parsed[!grepl(regexPatterns$gs_identifier, variant_id), c('dna_ref_read_count', 'dna_alt_read_count') := extractSplitDataFromVcf(vcf_table = vcf_dt[!grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                         sample_tag = t_tag,
                                                                                                                                         format_tag = 'DP4')]
    vcf_parsed[grepl(regexPatterns$gs_identifier, variant_id), c('dna_ref_read_count', 'dna_alt_read_count') := extractSplitDataFromVcf(vcf_table = vcf_dt[grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                        sample_tag = n_tag,
                                                                                                                                        format_tag = 'DP4')]
  } else if ('AD' %in% tags$format_tags) {
    vcf_parsed[!grepl(regexPatterns$gs_identifier, variant_id), c('dna_ref_read_count', 'dna_alt_read_count') := extractSplitDataFromVcf(vcf_table = vcf_dt[!grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                         sample_tag = t_tag,
                                                                                                                                         format_tag = 'AD')]
    vcf_parsed[grepl(regexPatterns$gs_identifier, variant_id), c('dna_ref_read_count', 'dna_alt_read_count') := extractSplitDataFromVcf(vcf_table = vcf_dt[grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                        sample_tag = n_tag,
                                                                                                                                        format_tag = 'AD')]
  }

  # extract total read counts
  if ('DP' %in% tags$format_tags) {
    vcf_parsed[!grepl(regexPatterns$gs_identifier, variant_id), dna_total_read_count := as.numeric(extractSingleDataFromVcfField(vcf_table = vcf_dt[!grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                 sample_tag = t_tag,
                                                                                                                                 format_tag = 'DP'))]
    vcf_parsed[grepl(regexPatterns$gs_identifier, variant_id), dna_total_read_count := as.numeric(extractSingleDataFromVcfField(vcf_table = vcf_dt[grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                                sample_tag = n_tag,
                                                                                                                                format_tag = 'DP'))]
  }

  # extract/compute dna_vaf
  if ('AF' %in% tags$format_tags) {
    vcf_parsed[!grepl(regexPatterns$gs_identifier, variant_id), dna_vaf := as.numeric(extractSingleDataFromVcfField(vcf_table = vcf_dt[!grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                    sample_tag = t_tag,
                                                                                                                    format_tag = 'AF',
                                                                                                                    split_by = ','))]
    vcf_parsed[grepl(regexPatterns$gs_identifier, variant_id), dna_vaf := as.numeric(extractSingleDataFromVcfField(vcf_table = vcf_dt[grepl(regexPatterns$gs_identifier, variant_id)],
                                                                                                                   sample_tag = n_tag,
                                                                                                                   format_tag = 'AF',
                                                                                                                   split_by = ','))]
  } else {
    vcf_parsed[, dna_vaf := unlist(mclapply(seq(1, nrow(vcf_parsed)),
                                            function(row_index) {
                                              if (is.numeric(vcf_parsed$dna_alt_read_count[row_index]) & is.numeric(vcf_parsed$dna_total_read_count[row_index])) {
                                                return(vcf_parsed$dna_alt_read_count[row_index] / vcf_parsed$dna_total_read_count[row_index])
                                              } else {
                                                return(NA)
                                              }
                                            }, mc.cores = 20))]
  }

  # # if multiple variant_ids are present:
  # # (1) if cosmic_id is present, replace variant_id with '.' (i.e. tumor-specific variant)
  # # (2) if no cosmic_id is present, rs_id take preference over gs_id
  # # (3) if no rs_id is present, return first gs_id (should only be one anyway)
  # vcf_parsed$variant_id = unlist(lapply(str_split(string = vcf_parsed$variant_id, pattern = ';'), # for some reason data.table's [, := ] notation didn't work here..
  # 																			function(vctr) {
  # 																				if (length(vctr) >= 1) {
  # 																					# if COSMIC_id match found, return '.', as variant is tumor variant
  # 																					if (any(grepl(pattern = regexPatterns$cosmic_identifier, x = vctr))) {match = '.'}
  # 																					# if no COSMIC_id, match rs_id and return first hit
  # 																					else if (any(grepl(pattern = regexPatterns$rs_identifier, x = vctr))) {match = grep(pattern = regexPatterns$rs_identifier, x = vctr, value = T)[1]}
  # 																					# if no rs_id match, return first value
  # 																					else {match = vctr[1]}
  #
  # 																					return(match)
  # 																				} else {
  # 																					return(vctr)
  # 																				}
  # 																			}))

  return(vcf_parsed)
}

extractSingleDataFromVcfField = function(vcf_table, sample_tag, format_tag, split_by = NULL) {
  # if sample_tag not found in columns, return NA
  if (!any(grepl(pattern = sample_tag, x = names(vcf_table)))) {
    return(NA)
  }

  tag_positions = unlist(mclapply(str_split(vcf_table$format, pattern = ':'),
                                  function(x) {
                                    match(x = format_tag, table = x, nomatch = NA)
                                  }, mc.cores = 20))

  vcf_tumor_split = str_split(vcf_table[[sample_tag]], pattern = ':')

  tag_data = unlist(mclapply(seq(1, length(tag_positions)),
                             function(x) {
                               # if we cannot find tag position & variant is a gs SNP, look in INFO field
                               if (is.na(tag_positions[[x]]) & grepl(pattern = regexPatterns$gs_identifier, x = vcf_table$variant_id[x])) {
                                 snp_split = unlist(str_split(string = vcf_table$info[x], pattern = ';'))
                                 tag_data = grep(pattern = paste0('^', format_tag, '='), x = snp_split, value = TRUE)

                                 if (length(tag_data) < 1) return(NA)

                                 tag_data = sub(pattern = paste0('^', format_tag, '='), replacement = '', x = tag_data)

                                 if (!is.null(split_by)) {
                                   tag_data = str_split(string = tag_data, pattern = split_by)
                                   tag_data = unlist(mclapply(tag_data, function(data) data[1], mc.cores = 20))
                                 }

                                 return(tag_data)
                               }

                               vcf_tumor_split[[x]][tag_positions[x]]

                             }, mc.cores = 20))

  return(tag_data)
}

extractVariantReadCountsFromVcf = function(vcf_table, sample_tag, count_tag) {
  count_data = extractSingleDataFromVcfField(vcf_table = vcf_table,
                                             sample_tag = sample_tag,
                                             format_tag = count_tag)

  count_data = mclapply(str_split(string = count_data, pattern = ','), as.numeric, mc.cores = 20)

  count_data = mclapply(count_data,
                        function(counts) {
                          if (any(!is.na(counts))) {
                            data = setNames(object = counts, nm = c('A', 'C', 'G', 'T'))
                            data = as.data.table(as.list(data))
                          } else {
                            data = as.data.table(as.list(rep(x = NA_integer_, 4)))
                            data = setNames(object = data, nm = c('A', 'C', 'G', 'T'))
                          }
                          return(data)
                        }, mc.cores = 20)

  count_data = lapply(rbindlist(count_data), as.numeric)
  return(count_data)
}

extractSplitDataFromVcf = function(vcf_table, sample_tag, format_tag, split_by = ',') {
  # if sample_tag not found in columns, return NA
  if (!any(grepl(pattern = sample_tag, x = names(vcf_table)))) {
    return(list(NA))
  }

  tag_data = extractSingleDataFromVcfField(vcf_table = vcf_table,
                                           sample_tag = sample_tag,
                                           format_tag = format_tag)

  tag_data = mclapply(str_split(string = tag_data, pattern = split_by),
                      as.numeric, mc.cores = 20)

  if (format_tag == 'DP4') { # for allele count data
    split_one = unlist(mclapply(tag_data, function(data) sum(data[c(1,2)]), mc.cores = 20)) # sum forward & reverse strand counts for ref supporting reads
    split_two = unlist(mclapply(tag_data, function(data) sum(data[c(3,4)]), mc.cores = 20)) # sum forward & reverse strand counts for alt supporting reads
  } else {
    split_one = unlist(mclapply(tag_data, function(data) data[1], mc.cores = 20))
    split_two = unlist(mclapply(tag_data, function(data) data[2], mc.cores = 20))
  }

  return(list(split_one, split_two))
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


# Variant allele expression -----------------------------------------------
findRnaReadLevelEvidenceForVariants = function(variant_input_path = file.path(rootDirectory, '1a_variants', 'parsed'),
                                               variant_regex = '\\.tsv$',
                                               rna_path = file.path(rootDirectory, '1b_rnaseq_data', 'bam'),
                                               quant_mode = 'salmon',
                                               pileup_mode = 'samtools',
                                               fasta_genome_ref = runOptions$samtools$fastaGenomeRef,
                                               sample_info_path = file.path(rootDirectory, 'sample_info.tsv')) {

  registerDoMC(runOptions$samtools$numberOfWorkers)

  if (!dir.exists(rna_path)) {
    message('RNAseq data directory does not exist at ', rna_path)
    user_choice = menu(choices = c('yes', 'no'), title = 'Do you want to continue without RNAseq data?')

    if (user_choice == 2) {
      stop('Please copy RNAseq data to ', rna_path)
    } else if (user_choice == 1) {
      # insert code to add empty rna_* columns to tsv files
    } else {
      stop('Please enter a valid option (1 or 2)')
    }
  }

  # parse input data
  input_data = lapply(list.files(path = variant_input_path,
                                 pattern = variant_regex,
                                 full.names = TRUE),
                      fread,
                      colClasses = list(character = c('chromosome')))
  input_data = setNames(object = input_data,
                        nm = list.files(path = variant_input_path,
                                        pattern = variant_regex))

  # load sample info
  if (file.exists(sample_info_path)) {
    sample_info = fread(sample_info_path, sep = '\t', header = T, na.strings = c('','NA', 'N.A.'))
  } else {
    stop('Sample info file missing at "', sample_info_path, '", please provide correct path in argument to findRnaReadLevelEvidenceForVariants')
  }

  if (sample_info$patient_id[1] == 'place_holder') {
    stop('Please fill in sample_info.tsv')
  }

  # make list of all unique variants found in parsed VCF
  snv_positions = lapply(input_data,
                         function(x) {
                           setorder(x = unique(x[grepl(pattern = '^[0-9]{1,2}$|^[XY]$',
                                                       x = x$chromosome) &
                                                   !grepl(pattern = regexPatterns$gs_identifier, # rs_id's are now regarded as tumor-specific variants
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
                                           file = file.path(rootDirectory, '1a_variants', 'poslist', paste0(sub(regexPatterns$file_extension, '', names(snv_positions)[x]), '_poslist.tsv')),
                                           quote = FALSE,
                                           col.names = FALSE,
                                           row.names = FALSE)))

  # make overview of sample input files
  regex_pattern = if (quant_mode == 'cufflinks') {'accepted_hits\\.bam$'} else if (quant_mode == 'salmon') {'Aligned\\.sortedByCoord\\.out\\.bam$'} else {stop('Please define valid quant_mode (cufflinks/salmon)')}

  sample_combinations = data.table(locations_file = sapply(sample_info$dna_data_prefix, function(x) grep(pattern = x,
                                                                                                         x = list.files(path = file.path(rootDirectory, '1a_variants', 'poslist'),
                                                                                                                        pattern = paste0(gsub(regexPatterns$file_extension, '', variant_regex), '_poslist\\.tsv'),
                                                                                                                        full.names = TRUE),
                                                                                                         value = T),
                                                           USE.NAMES = FALSE),
                                   rna_bam_file = sapply(sample_info$rna_data_prefix, function(x) grep(pattern = x,
                                                                                                       x = list.files(path = rna_path,
                                                                                                                      pattern = regex_pattern,
                                                                                                                      recursive = TRUE,
                                                                                                                      full.names = TRUE),
                                                                                                       value = T),
                                                         USE.NAMES = FALSE)
  )

  if (nrow(sample_combinations) < 1) {
    stop('No valid sample combinations found, please check dna/rna_prefixes and rna bamfile filenames')
  }

  if (is.null(fasta_genome_ref) | !file.exists(fasta_genome_ref)) {
    stop('Please provide valid fasta DNA reference (set location in runConfig.R)')
  } else {
    message('Performing pileup with genomic reference: ', fasta_genome_ref, '\n')
  }

  invisible(foreach(i = 1:nrow(sample_combinations)) %dopar% {
    if (!file.exists(file.path(rootDirectory, '1b_rnaseq_data', 'pileups',
                               paste0(sub(regexPatterns$file_extension, '', basename(sample_combinations$rna_bam_file[i])),
                                      if (is.null(sample_combinations$locations_file[i])) {'_mpil.tsv'} else {'_mpil_loc.tsv'})))) {

      if (pileup_mode == 'samtools') {
        performSamtoolsPileup(bam_file = sample_combinations$rna_bam_file[i],
                              locations_file = sample_combinations$locations_file[i],
                              fasta_reference = fasta_genome_ref)
      } else if (pileup_mode == 'sambamba') {
        performSambambaPileup(bam_file = sample_combinations$rna_bam_file[i],
                              locations_file = sample_combinations$locations_file[i],
                              fasta_reference = fasta_genome_ref)
      } else {
        stop('Please specify valid argument for "pileup_mode" (either "samtools" or "sambamba")')
      }
    }
  })

  pileup_loc_data = lapply(list.files(path = file.path(rootDirectory, '1b_rnaseq_data', 'pileups'),
                                      pattern = '_mpil_loc\\.tsv',
                                      full.names = TRUE),
                           fread,
                           colClasses = c(V1 = 'character'),
                           col.names = c('chromosome', 'start_position', 'ref_base', 'number_of_reads', 'rna_read_bases', 'base_quality'))

  pileup_loc_data = setNames(object = pileup_loc_data, nm = list.files(path = file.path(rootDirectory, '1b_rnaseq_data', 'pileups'),
                                                                       pattern = '_mpil_loc\\.tsv',
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

                                  # order data
                                  setorder(merged_data, chromosome, start_position)
                                  merged_data = rbindlist(list(merged_data %>%
                                                                 filter(!grepl(regexPatterns$gs_identifier, variant_id)) %>%
                                                                 .[naturalorder(chromosome)],
                                                               merged_data %>%
                                                                 filter(grepl(regexPatterns$gs_identifier, variant_id)) %>%
                                                                 .[naturalorder(chromosome)])
                                  )
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
                                                                 if (sum(dt$rna_ref_read_count[row_index], dt$rna_alt_read_count[row_index]) >= 7 |
                                                                     (dt$rna_ref_read_count[row_index] >= 5 | dt$rna_alt_read_count[row_index] >= 5)) {
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

performSamtoolsPileup = function(bam_file, locations_file = NULL, fasta_reference = NULL, optional_args = c('-E'), execute = TRUE) {
  dir.create(file.path(rootDirectory, '1b_rnaseq_data', 'pileups'), showWarnings = FALSE)

  command = paste(runOptions$samtools$samtoolsPath,
                  'mpileup',
                  if (!is.null(locations_file)) {paste('-l', locations_file)},
                  if (!is.null(fasta_reference)) {paste('-f', fasta_reference)},
                  optional_args,
                  bam_file,
                  '-o', file.path(rootDirectory, '1b_rnaseq_data', 'pileups', paste0(sub(regexPatterns$file_extension, '', basename(bam_file)),
                                                                                     if (is.null(locations_file)) {'_mpil.tsv'} else {'_mpil_loc.tsv'})))

  commandWrapper(command = command, nice = NULL, execute = execute)

  Sys.sleep(time = 1)
}

performSambambaPileup = function(bam_file, locations_file = NULL, fasta_reference = NULL, execute = TRUE) {
  dir.create(file.path(rootDirectory, '1b_rnaseq_data', 'pileups'), showWarnings = FALSE)

  command = paste(runOptions$samtools$sambambaPath,
                  'mpileup',
                  '-t', 6,
                  '-o', file.path(rootDirectory, '1b_rnaseq_data', 'pileups', paste0(sub(regexPatterns$file_extension, '', basename(bam_file)),
                                                                                     if (is.null(locations_file)) {'_mpil.tsv'} else {'_mpil_loc.tsv'})),
                  bam_file,
                  '--samtools',
                  if (!is.null(locations_file)) {paste('-l', locations_file)},
                  if (!is.null(fasta_reference)) {paste('-f', fasta_reference)})

  commandWrapper(command = command, nice  = NULL, execute = execute)

  Sys.sleep(time = 1)
}


# Varcontext generation ---------------------------------------------------
performVarcontextGeneration = function(variant_path = file.path(rootDirectory, '1a_variants', 'parsed'),
                                       variant_regex = '\\.tsv$',
                                       filter_rna_alt_expression = TRUE,
                                       vcf_fields = c('ID', 'CHROM', 'POS', 'REF', 'ALT'),
                                       execute = TRUE) {
  registerDoMC(runOptions$varcontext$numberOfWorkers)

  dir.create(file.path(rootDirectory, '2_varcontext'),
             showWarnings = FALSE)
  dir.create(file.path(rootDirectory, '2_varcontext', 'input_lists'),
             showWarnings = FALSE)
  dir.create(file.path(rootDirectory, '2_varcontext', 'varcontext_logs'),
             showWarnings = FALSE)

  # make list of input files
  variant_lists = list.files(path = variant_path,
                             pattern = variant_regex,
                             recursive = FALSE,
                             full.names = TRUE)
  variant_data = lapply(variant_lists,
                        fread,
                        colClasses = list(character = c('chromosome')))
  variant_data = setNames(object = variant_data,
                          nm = list.files(path = variant_path,
                                          pattern = variant_regex))

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
                                             pattern = variant_regex,
                                             full.names = TRUE),
                     execute = execute)
}

generateVarcontext = function(input_list, execute = TRUE) {
  if (length(input_list) < 1) {
    stop('No input lists found...')
  }

  setwd(runOptions$varcontext$varcontextDirectory)

  commandList = sapply(1:length(input_list),
                       function(i) {
                         filename = sub(pattern = regexPatterns$file_extension, replacement = '', x = basename(input_list[i]))

                         command = paste0('export ENSEMBLAPI="', runOptions$varcontext$ensemblApi, '";',
                                          'export PERL5LIB="$PERL5LIB:', runOptions$varcontext$perlLibs,'";',
                                          'perl ', file.path(runOptions$varcontext$varcontextDirectory, 'varcontext/create_context.pl'), ' ',
                                          '--separator=', runOptions$varcontext$fieldSeparator, ' ',
                                          ifelse(runOptions$varcontext$canonicalOnly, '--canonical ', ''),
                                          ifelse(runOptions$varcontext$peptideContext, '--peptide ', ''),
                                          ifelse(runOptions$varcontext$nmdStatus, '--nmd ', ''),
                                          '"', input_list[i], '"',
                                          ' 1> "', file.path(rootDirectory, '2_varcontext', paste(filename, 'varcontext.tsv"', sep = '_')),
                                          ' 2> "', file.path(rootDirectory, '2_varcontext', 'varcontext_logs', paste(filename, 'warnings.log"', sep = '_')))
                       })

  if (execute) {
    invisible(foreach(i = 1:length(input_list)) %dopar% {
      filename = sub(pattern = regexPatterns$file_extension, replacement = '', x = basename(input_list[i]))

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

      command = commandList[i]

      commandWrapper(command = command, nice = NULL, wait = FALSE, execute = execute)
    })
  } else {
    sapply(commandList,
           function(msg) {
             message(msg, '\n')
           })
  }

  setwd(rootDirectory)
}


# Neolution input generation ----------------------------------------------
prepareNeolutionInput = function(varcontext_path = file.path(rootDirectory, '2_varcontext'),
                                 rna_path = file.path(rootDirectory, '1b_rnaseq_data', 'processed_salmon'),
                                 sample_info_path = file.path(rootDirectory, 'sample_info.tsv'),
                                 rna_file_suffix = 'salmon-quant-by-ensg\\.tsv',
                                 expression_unit = 'tpm') {

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
                                                   expression_table = rnaseq_data[[sample_combinations[x, rna_expression_data]]],
                                                   expression_unit = expression_unit)
                                } else if (length(sample_combinations[x, variants]) > 0) {
                                  message('Adding empty rna_expression column for ', sample_combinations[x, variants])
                                  mergeByEnsemblId(variant_table = varcontext_data[[sample_combinations[x, variants]]],
                                                   expression_table = NULL,
                                                   expression_unit = expression_unit)
                                } else {
                                  warning('No input data present for ', sample_combinations[x, variants], ' & ', sample_combinations[x, rna_expression_data])
                                }
                              })
    prediction_input = setNames(object = prediction_input, nm = sub(regexPatterns$file_extension, '', sample_combinations$variants))

    rna_coverage_summary = data.table(sample = names(prediction_input),
                                      percent_ensg_coverage = sapply(prediction_input, function(x) round(x = length(which(!is.na(x[[expression_unit]]))) / nrow(x) * 100,
                                                                                                         digits = 1))
    )

    rna_expression_summary = data.table(sample = names(rnaseq_data),
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


# SnpEff report generation ------------------------------------------------
runSnpEff = function(vcf_path = file.path(rootDirectory, '1a_variants', 'vcf'), vcf_regex = '\\.vcf$', filter_snps = TRUE, canon_only = TRUE, execute = TRUE) {
  #registerDoMC(2)

  message('Step 4: Running snpEff')
  message('Using genome & gene builds: ', runOptions$snpeff$build)
  dir.create(path = file.path(rootDirectory, '4_snpEff'),
             showWarnings = FALSE)

  # make list of input files
  variantLists = list.files(path = vcf_path,
                            pattern = vcf_regex,
                            recursive = FALSE,
                            full.names = TRUE)

  invisible(foreach(i = 1:length(variantLists)) %do% {
    command_snpsift = paste('java -Xmx4g -jar', file.path(runOptions$snpeff$path, 'SnpSift.jar'),
                            'filter " (ID "\'!\'"~ \'gs[0-9]+\') "',
                            variantLists[i])

    command_snpeff = paste('java -Xmx4g -jar', file.path(runOptions$snpeff$path, 'snpEff.jar'),
                           '-nodownload',
                           if (canon_only) {'-canon'},
                           '-stats',
                           file.path(rootDirectory,'4_snpEff', paste0(sub(pattern = regexPatterns$file_extension,
                                                                          replacement = '',
                                                                          x = basename(variantLists[i])),
                                                                      if (!filter_snps) {'-with_snps'},
                                                                      '-snpEff_summary.html')),
                           runOptions$snpeff$build, '>', file.path(rootDirectory,'4_snpEff', paste0(sub(pattern = regexPatterns$file_extension,
                                                                                                        replacement = '',
                                                                                                        x =  basename(variantLists[i])),
                                                                                                    if (!filter_snps) {'-with_snps'},
                                                                                                    '-snpEff.vcf')))

    command = if (filter_snps) {paste(command_snpsift, command_snpeff, sep = '|')} else {paste(paste('cat', variantLists[i]), command_snpeff, sep = '|')}


    commandWrapper(command = command, execute = execute)

  })
}


# Mutation analysis -------------------------------------------------------

mutationalSignatureAnalysis = function(table, genome_build = 'GRCh38') {
  required_colnames = c('patient_id', 'chromosome', 'start_position', 'end_position', 'ref_allele', 'alt_allele')
  if (any(required_colnames %nin% names(table))) {
    stop('Please make sure all required cols are present in input table. \n', paste(required_colnames, collapse = ', '))
  }

  if (genome_build == 'GRCh38') {
    if (!require('BSgenome.Hsapiens.NCBI.GRCh38')) {
      library(BiocInstaller)
      biocLite('BSgenome.Hsapiens.NCBI.GRCh38', suppressUpdates = T)
    }

    genome = getBSgenome('BSgenome.Hsapiens.NCBI.GRCh38')
  } else if (genome_build == 'hg19' | genome_build == 'GRCh37') {
    if (!require('BSgenome.Hsapiens.UCSC.hg19')) {
      library(BiocInstaller)
      biocLite('BSgenome.Hsapiens.UCSC.hg19', suppressUpdates = T)
    }

    genome = getBSgenome('BSgenome.Hsapiens.UCSC.hg19')
  }

  # exclude germline SNPs from analysis
  table_subset = table[(nchar(ref_allele) == 1 & nchar(alt_allele) == 1)] # make sure we take only SNVs

  sigs_input = mut.to.sigs.input(mut.ref = table_subset,
                                 sample.id = 'patient_id',
                                 chr = 'chromosome',
                                 pos = 'start_position',
                                 ref = 'ref_allele',
                                 alt = 'alt_allele',
                                 bsg = genome)

  output_sigs = setNames(object = lapply(row.names(sigs_input),
                                         function(name) whichSignatures(tumor.ref = sigs_input,
                                                                        signatures.ref = signatures.nature2013,
                                                                        sample.id = name,
                                                                        contexts.needed = TRUE)),
                         nm = row.names(sigs_input))

  return(output_sigs)
}

# override deconstructSigs plotting function
plotSignaturesSmallLabels = function(sigs.output, sub = "")
{
  op <- graphics::par()
  tumor <- sigs.output[["tumor"]]
  product <- sigs.output[["product"]]
  diff <- sigs.output[["diff"]]
  weights <- sigs.output[["weights"]]
  y_limit <- 1.2 * max(tumor, product)
  tumor_plotting <- formatContexts(tumor)
  product_plotting <- formatContexts(product)
  error_summed <- round(sqrt(sum(diff * diff)), digits = 3)
  diff_plotting <- formatContexts(diff)
  name <- unique(tumor_plotting$sample.id)
  subtype <- sub
  tmp <- which(weights != 0)
  c <- paste(colnames(weights)[tmp[1]], " : ", round(weights[tmp[1]],
                                                     3), sep = "")
  if (length(tmp) > 1) {
    for (i in tmp[2:length(tmp)]) {
      c <- paste(c, " & ", colnames(weights)[i], " : ",
                 round(weights[i], 3), sep = "")
    }
  }
  if (subtype == "") {
    top.title <- name
  }
  if (subtype != "") {
    top.title <- paste(name, " -- ", subtype, sep = "")
  }
  graphics::par(mfrow = c(3, 1), xpd = FALSE, mar = c(5, 4,
                                                      4, 2), oma = c(0, 0, 0, 4))
  grDevices::palette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2"))
  graphics::barplot(tumor_plotting$fraction, names.arg = tumor_plotting$full_context,
                    cex.names = 0.7, las = 2, col = NA, ylim = c(0, y_limit),
                    border = NA, xaxt = "n", ann = FALSE, yaxt = "n", space = 0.3)
  graphics::box()
  x = graphics::par("usr")
  graphics::abline(h = seq(from = 0, to = y_limit, by = 0.01),
                   col = "#d3d3d350", lty = 1)
  graphics::abline(v = seq(from = x[1], to = x[2], by = 1),
                   col = "#d3d3d350", lty = 1)
  graphics::barplot(tumor_plotting$fraction, names.arg = tumor_plotting$full_context,
                    cex.names = 0.7, las = 2, col = factor(tumor_plotting$mutation),
                    ylim = c(0, y_limit), border = NA, space = 0.3, main = top.title,
                    ylab = "fraction", add = TRUE)
  graphics::barplot(product_plotting$fraction, names.arg = product_plotting$full_context,
                    cex.names = 0.7, las = 2, col = NA, ylim = c(0, y_limit),
                    border = NA, xaxt = "n", ann = FALSE, yaxt = "n", space = 0.3)
  graphics::box()
  x = graphics::par("usr")
  graphics::abline(h = seq(from = 0, to = y_limit, by = 0.01),
                   col = "#d3d3d350", lty = 1)
  graphics::abline(v = seq(from = x[1], to = x[2], by = 1),
                   col = "#d3d3d350", lty = 1)
  graphics::barplot(product_plotting$fraction, names.arg = product_plotting$full_context,
                    cex.main = 1, cex.names = 0.7, las = 2, col = factor(product_plotting$mutation),
                    ylim = c(0, y_limit), border = NA, space = 0.3, main = c,
                    ylab = "fraction", add = TRUE)
  graphics::barplot(diff_plotting$fraction, names.arg = diff_plotting$full_context,
                    cex.names = 0.7, las = 2, col = NA, ylim = c(-y_limit,
                                                                 y_limit), border = NA, xaxt = "n", ann = FALSE,
                    yaxt = "n", space = 0.3)
  graphics::box()
  x = graphics::par("usr")
  graphics::abline(h = seq(from = -y_limit, to = y_limit,
                           by = 0.02), col = "#d3d3d350", lty = 1)
  graphics::abline(v = seq(from = x[1], to = x[2], by = 1),
                   col = "#d3d3d350", lty = 1)
  graphics::barplot(diff_plotting$fraction, names.arg = diff_plotting$full_context,
                    cex.names = 0.7, las = 2, col = factor(diff_plotting$mutation),
                    ylim = c(-y_limit, y_limit), border = "black", space = 0.3,
                    main = paste("error = ", error_summed, sep = ""), ylab = "fraction",
                    add = TRUE)
  graphics::par(fig = c(0, 1, 0, 1), oma = c(1, 1, 1, 1),
                mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n",
                 yaxt = "n")
  graphics::legend("right", legend = unique(tumor_plotting$mutation),
                   col = c("#999999", "#E69F00", "#56B4E9", "#009E73",
                           "#F0E442", "#0072B2"), bty = "n", ncol = 1, inset = c(-0,
                                                                                 0), pch = 15, xpd = TRUE, pt.cex = 2.5)
  on.exit(suppressWarnings(graphics::par(op)))
}


# Peptide order -----------------------------------------------------------
parseEpitopePredictions = function(path, sample_table = sample_info, pattern = '_epitopes\\.csv') {
  require(stringr)
  require(data.table)

  # get predictions paths
  files = list.files(path = path,
                     pattern = pattern,
                     full.names = TRUE)

  # extract filenames
  short_names = sub(pattern = '.+[/]', replacement = '', x = files)
  short_names = sub(pattern = pattern, replacement = '', x = short_names)

  # get data in & sort by tumor peptide affinity
  predictions = lapply(seq(1, length(files)),
                       function(i) {
                         data = fread(files[i])
                         data[, sample_prefix := sub(pattern = regexPatterns$seqdata_prefix,
                                                     replacement = '',
                                                     x = short_names[i],
                                                     perl = T)]
                         data[, patient_id := sample_table[dna_data_prefix == data[, unique(sample_prefix)], patient_id]]

                         return(data)
                       })

  sapply(predictions,
         function(x) setkeyv(x = x, grep(pattern = 'tumor_.+affinity', x = names(x), value = TRUE)))

  # set table names
  predictions = setNames(object = predictions,
                         nm = short_names)

  return(predictions)
}

applyCutoffs = function(predictions, model = NULL, rank = NULL, affinity = NULL, processing = NULL, expression = NULL, selfsim = NULL, invert = FALSE) {
  # do checks for necessary cutoff values here
  if (is.numeric(model) & is.numeric(expression) & is.logical(selfsim)) {
    message('Applying cutoffs -> model >= ', model, ' | expression > ', expression, ' | selfsim: ', selfsim)
  } else if (is.numeric(rank) & is.numeric(processing) & is.numeric(expression) & is.logical(selfsim)) {
    message('Applying cutoffs -> rank <= ', rank, ' | processing >= ', processing, ' | expression > ', expression, ' | selfsim: ', selfsim)
  } else if (is.numeric(affinity) & is.numeric(processing) & is.numeric(expression) & is.logical(selfsim)) {
    message('Applying cutoffs -> affinity <= ', rank, ' | processing >= ', processing, ' | expression > ', expression, ' | selfsim: ', selfsim)
  } else {
    stop('Please provide required cutoff values as arguments')
  }

  # take subsets
  if (!invert) {
    subsets = lapply(seq(1, length(predictions), 1),
                     function(x) {
                       if (is.numeric(model)) {
                         data_subset = predictions[[x]][model_prediction >= model & (rna_expression > expression | is.na(rna_expression))]
                         if (selfsim) {data_subset = data_subset[different_from_self == TRUE]}

                       } else if (is.numeric(rank)) {
                         data_subset = subset(x = predictions[[x]],
                                              subset = predictions[[x]][[grep(pattern = 'tumor.+rank',
                                                                              x = names(predictions[[x]]),
                                                                              value = T)]] <= rank)
                         data_subset = data_subset[tumor_processing_score >= processing & (rna_expression > expression | is.na(rna_expression))]
                         if (selfsim) {data_subset = data_subset[different_from_self == TRUE]}

                       } else if (is.numeric(affinity)) {
                         data_subset = subset(x = predictions[[x]],
                                              subset = predictions[[x]][[grep(pattern = 'tumor.+affinity',
                                                                              x = names(predictions[[x]]),
                                                                              value = T)]] <= affinity)
                         data_subset = data_subset[tumor_processing_score >= processing & (rna_expression > expression | is.na(rna_expression))]
                         if (selfsim) {data_subset = data_subset[different_from_self == TRUE]}

                       } else {
                         stop('No model, rank or affinity cutoff provided')
                       }

                       return(data_subset)
                     })
  } else {
    message('Inverting cutoffs; returning peptides not making cutoffs')
    subsets = lapply(seq(1, length(predictions), 1),
                     function(x) {
                       if (is.numeric(model)) {
                         data_subset = predictions[[x]][model_prediction < model & rna_expression <= expression]
                         if (selfsim) {data_subset = data_subset[different_from_self == FALSE | is.na(different_from_self)]}

                       } else if (is.numeric(rank)) {
                         data_subset = subset(x = predictions[[x]],
                                              subset = predictions[[x]][[grep(pattern = 'tumor.+rank',
                                                                              x = names(predictions[[x]]),
                                                                              value = T)]] > rank)
                         data_subset = data_subset[tumor_processing_score < processing & rna_expression <= expression]
                         if (selfsim) {data_subset = data_subset[different_from_self == FALSE | is.na(different_from_self)]}

                       } else if (is.numeric(affinity)) {
                         data_subset = subset(x = predictions[[x]],
                                              subset = predictions[[x]][[grep(pattern = 'tumor.+affinity',
                                                                              x = names(predictions[[x]]),
                                                                              value = T)]] > affinity)
                         data_subset = data_subset[tumor_processing_score < processing & rna_expression <= expression]
                         if (selfsim) {data_subset = data_subset[different_from_self == FALSE | is.na(different_from_self)]}

                       } else {
                         stop('No model, rank or affinity cutoff provided')
                       }

                       return(data_subset)
                     })
  }

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
                                                   select = c('patient_id', 'sample_prefix', 'hla_allele', 'tumor_peptide'))
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

  split_by_xmer_info = setNames(object = split_by_xmer_info,
                                nm = split_by)

  return(split_by_xmer_info)
}

