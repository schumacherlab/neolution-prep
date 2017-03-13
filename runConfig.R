## set runtime configuration

# regex for matching SNPs
regex_snps = '[gr]s\\d+'

# regex for isolating GCF prefix
regex_prefix = '_mg.+|_S\\d_L.+'

# regex for excluding particular alleles from analysis
exclusion_pattern = 'C[0-9]{4}'

# run options
runOptions = list(general = list(),

									# set samtools mpileup options
									samtools = list(numberOfWorkers = 8,
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
																	 model_cutoff = 0.02,
																	 expression_cutoff = 0,
																	 random_forest_model = TRUE,
																	 xmer = c(9:11),
																	 selfsim_filter = FALSE,
																	 selfsim_filter_mode = 'simple',
																	 selflist = FALSE),

									# set snpEff locations
									snpeff = list(path = '/home/NFS/users/l.fanchi/libs/snpEff',
																build = 'GRCh38.86'))



