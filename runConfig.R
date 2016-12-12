# set root directory
# rootDirectory = '/home/NFS/users/l.fanchi/projects/opacin'

runOptions = list(general = c(),
									varcontext = c(),
									neolution = c(),
									snpeff = c())

# set general options
runOptions$general$numberOfWorkers = 3

# set varcontext options
runOptions$varcontext$fieldSeparator = '"\t"';
runOptions$varcontext$canonicalOnly = FALSE;
runOptions$varcontext$peptideContext = FALSE;
runOptions$varcontext$nmdStatus = TRUE;

runOptions$varcontext$varcontextDirectory = '/home/NFS/users/l.fanchi/dev_environments/varcontext'
runOptions$varcontext$ensemblApi = '/home/NFS/users/l.fanchi/libs/ensembl_81/'
runOptions$varcontext$perlLibs = paste('/home/NFS/users/l.fanchi/perl5/lib/perl5','/home/NFS/users/l.fanchi/libs/bioperl-live',
																			 '/home/NFS/users/l.fanchi/libs/ensembl/modules','/home/NFS/users/l.fanchi/libs/ensembl-variation/modules', sep = ':')

# set neolution options
runOptions$neolution$rankCutoff = 3
runOptions$neolution$processingCutoff = 0.5
runOptions$neolution$xmer = c(9:11)

# set snpEff locations
runOptions$snpeff$path = '/home/NFS/users/l.fanchi/libs/snpEff'
runOptions$snpeff$build = 'GRCh38.81'
