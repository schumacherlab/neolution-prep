## set runtime configuration
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

																	 rankCutoff = 3,
																	 processingCutoff = 0.5,
																	 expressionCutoff = 0,
																	 modelCutoff = 0.02,
																	 xmer = c(9:11)),

									# set snpEff locations
									snpeff = list(path = '/home/NFS/users/l.fanchi/libs/snpEff',
																build = 'GRCh38.86'))



