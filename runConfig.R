## set runtime configuration

# regexes
regexPatterns = list(file_extension = '\\.[^.]+$', # match file extension (everything after last dot, inclusive)
                     snp_identifier = '[gr]s\\d+', # for matching SNPs
                     gs_identifier = 'gs\\d+', # for matching germline snps annotated by our pipeline, keep boundless (no '^' or '$')
                     rs_identifier = 'rs\\d+', # for matching snps found in dbSNP, keep boundless (no '^' or '$')
                     cosmic_identifier = 'COSM\\d+', # for matching variants found in COSMIC coding muts database, keep boundless
                     seqdata_prefix = '(^.+?[ATCG]{7,8})(.+)', # for isolating GCF prefix
                     allele_exclusion = 'C[0-9]{4}') # for excluding particular alleles from analysis

# run options
commonPaths = list(libs_path = '/DATA/users/l.fanchi/libs',
                   resources_path = '/DATA/resources')

runOptions = list(general = list(gtf_annotation = file.path(commonPaths$resources_path,
                                                            'ensembl_88/gtf/Homo_sapiens.GRCh38.88.gtf')),

                  # set samtools mpileup options
                  samtools = list(samtoolsPath = file.path(commonPaths$libs_path,
                                                           'samtools-1.5/bin/samtools'),
                                  sambambaPath = file.path(commonPaths$libs_path,
                                                           'sambamba-0.6.6-Linux_x86_64/sambamba'),
                                  numberOfWorkers = 10,
                                  fastaGenomeRef = file.path(commonPaths$resources_path,
                                                             'ensembl_88/fasta_dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa')),

                  # set varcontext options
                  varcontext = list(number_of_workers = 10,

                                    field_separator = '"\t"',
                                    print_overlapping_bases = TRUE,
                                    mutant_allele_expression_filter = TRUE,

                                    ensembl_build = 88,
                                    assembly_build = 38,
                                    canonical_only = FALSE,
                                    cdna_context = FALSE,
                                    cdna_context_size = 30,
                                    peptide_context = FALSE,
                                    protein_context = TRUE,
                                    nmd_status = TRUE,

                                    varcontext_directory = file.path(commonPaths$resources_path, 'pipeline/varcontext'),
                                    ensembl_api = file.path(commonPaths$libs_path, 'ensembl_88/'),
                                    perl_libs = paste(c(file.path(commonPaths$libs_path,
                                                                 c('perl5/lib/perl5', 'bioperl-live'))),
                                                     collapse = ':')),

                  # set neolution options
                  neolution = list(path = file.path(commonPaths$resources_path, 'pipeline/neolution-live'),
                                   numberOfWorkers = 3, # number of sample threads for predictions
                                                        # (each thread spawns 'machine_cores/6' children during self-sim checking)

                                   model_cutoff = 0.02, # use 0.01 for TIL screens
                                                        # (more inclusive; we picked up low magnitude TIL hits at low prob scores)
                                                        # use 0.02 for PBMC screens
                                                        # (more stringent; we're unlikely to pick up low magnitude responses with low prob scores)
                                   rank_cutoff = NA,
                                   processing_cutoff = NA,
                                   expression_cutoff = 0,
                                   random_forest_model = TRUE,
                                   xmer = c(9:11),
                                   selfsim_filter = TRUE,
                                   selfsim_filter_mode = 'simple',
                                   selflist = FALSE,
                                   panversion = 4),

                  # set snpEff locations
                  snpeff = list(path = file.path(commonPaths$libs_path, 'snpEff'),
                                build = 'GRCh38.86'))
