## set runtime configuration

# regexes
regexPatterns = list(file_extension = '\\.[^.]+$', # match file extension (everything after last dot, inclusive)
                     snp_identifier = '[gr]s\\d+', # for matching SNPs
                     gs_identifier = 'gs\\d+', # for matching germline snps annotated by our pipeline, keep boundless (no '^' or '$')
                     rs_identifier = 'rs\\d+', # for matching snps found in dbSNP, keep boundless (no '^' or '$')
                     cosmic_identifier = 'COSM\\d+', # for matching variants found in COSMIC coding muts database, keep boundless
                     seqdata_prefix = '(^.+[ATCG]{7,8})(.+)', # for isolating GCF prefix
                     allele_exclusion = 'C[0-9]{4}') # for excluding particular alleles from analysis

# run options
userPaths = list(libs_path = '/DATA/users/l.fanchi/libs',
                 resources_path = '/DATA/users/l.fanchi/resources',
                 home_path = '/home/l.fanchi')

runOptions = list(general = list(gtf_annotation = file.path(userPaths$resources_path,
                                                            'ensembl_90/gtf/Homo_sapiens.GRCh38.90.gtf')),

                  # set samtools mpileup options
                  samtools = list(samtoolsPath = file.path(userPaths$libs_path,
                                                           'samtools-1.5/bin/samtools'),
                                  sambambaPath = file.path(userPaths$libs_path,
                                                           'sambamba-0.6.6-Linux_x86_64/sambamba'),
                                  numberOfWorkers = 10,
                                  fastaGenomeRef = file.path(userPaths$resources_path,
                                                             'ensembl_88/fasta_dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa')),

                  # set varcontext options
                  varcontext = list(number_of_workers = 10,

                                    field_separator = '"\t"',
                                    trim_overlapping_bases = FALSE,

                                    ensembl_build = 90,
                                    assembly_build = 38,
                                    canonical_only = FALSE,
                                    cdna_context = FALSE,
                                    cdna_context_size = 54,
                                    peptide_context = FALSE,
                                    protein_context = TRUE,
                                    nmd_status = TRUE,

                                    varcontext_directory = file.path(userPaths$home_path, 'stable_environments/varcontext'),
                                    ensembl_api = file.path(userPaths$libs_path, 'ensembl_90/'),
                                    perl_libs = paste(c(file.path(userPaths$libs_path,
                                                                 c('perl5/lib/perl5', 'bioperl-live'))),
                                                     collapse = ':')),

                  # set neolution options
                  neolution = list(path = file.path(userPaths$home_path, 'stable_environments/neolution-live'),
                                   numberOfWorkers = 3, # number of sample threads for predictions
                                                        # (each thread spawns 'machine_cores/6' children during self-sim checking)

                                   rank_cutoff = NA,
                                   processing_cutoff = NA,
                                   model_cutoff = 0.02, # use 0.01 for TIL screens
                                                        # (more inclusive; we picked up low magnitude TIL hits at low prob scores)
                                                        # use 0.02 for PBMC screens
                                                        # (more stringent; unlikely to pick up low magnitude responses with low prob scores)
                                   expression_cutoff = 0,
                                   random_forest_model = TRUE,
                                   xmer = c(9:11),
                                   selfsim_filter = FALSE,
                                   selfsim_filter_mode = 'simple',
                                   selflist = FALSE),

                  # set snpEff locations
                  snpeff = list(path = file.path(userPaths$libs_path, 'snpEff'),
                                build = 'GRCh38.86'))



