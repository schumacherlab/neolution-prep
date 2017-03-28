## set runtime configuration

# regexes
regexPatterns = list(file_extension = '\\.[^.]+$', # match file extension (everything after last dot, inclusive)
										 snp_identifier = '[gr]s\\d+', # for matching SNPs
										 gs_identifier = 'gs\\d+', # for matching snps not found in dbSNP, keep boundless (no '^' or '$')
										 rs_identifier = 'rs\\d+', # for matching snps found in dbSNP, keep boundless (no '^' or '$')
										 seqdata_prefix = '_mg.+|_S\\d_L.+', # for isolating GCF prefix
										 allele_exclusion = 'C[0-9]{4}') # for excluding particular alleles from analysis

# run options
runOptions = list(general = list(),

									# set samtools mpileup options
									samtools = list(samtoolsPath = '/home/NFS/users/l.fanchi/libs/samtools-1.3/bin/samtools',
																	sambambaPath = '/home/NFS/users/l.fanchi/libs/sambamba-0.6.6-Linux_x86_64/sambamba',
																	numberOfWorkers = 8,
																	fastaGenomeRef = '/home/NFS/users/l.fanchi/resources/ensembl_87/fasta_dna/Homo_sapiens.GRCh38.87.dna.primary_assembly.fa'),

									# set varcontext options
									varcontext = list(numberOfWorkers = 10,

																		fieldSeparator = '"\t"',
																		canonicalOnly = FALSE,
																		peptideContext = FALSE,
																		nmdStatus = TRUE,

																		varcontextDirectory = '/home/NFS/users/l.fanchi/dev_environments/varcontext',
																		ensemblApi = '/home/NFS/users/l.fanchi/libs/ensembl_81/',
																		perlLibs = paste('/home/NFS/users/l.fanchi/perl5/lib/perl5','/home/NFS/users/l.fanchi/libs/bioperl-live',
																										 '/home/NFS/users/l.fanchi/libs/ensembl/modules','/home/NFS/users/l.fanchi/libs/ensembl-variation/modules', sep = ':')),

									# set neolution options
									neolution = list(numberOfWorkers = 3, # number of threads for predictions (each thread spawns 16 children during self-sim checking)

																	 rank_cutoff = NA,
																	 processing_cutoff = NA,
																	 model_cutoff = 0.02, # use 0.02 for TIL screens (find largest responses), use 0.01 for PBMC screens (more inclusive)
																	 expression_cutoff = 0,
																	 random_forest_model = TRUE,
																	 xmer = c(9:11),
																	 selfsim_filter = FALSE,
																	 selfsim_filter_mode = 'simple',
																	 selflist = FALSE),

									# set snpEff locations
									snpeff = list(path = '/home/NFS/users/l.fanchi/libs/snpEff',
																build = 'GRCh38.86'))



