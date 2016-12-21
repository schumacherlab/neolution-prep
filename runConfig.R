## set runtime configuration
runOptions = list(general = c(),
                  varcontext = c(),
                  neolution = c(),
                  snpeff = c())

# set varcontext options
runOptions$varcontext$numberOfWorkers = 10

runOptions$varcontext$fieldSeparator = '"\t"';
runOptions$varcontext$canonicalOnly = FALSE;
runOptions$varcontext$peptideContext = FALSE;
runOptions$varcontext$nmdStatus = TRUE;

runOptions$varcontext$varcontextDirectory = '/home/NFS/users/l.fanchi/dev_environments/varcontext'
runOptions$varcontext$ensemblApi = '/home/NFS/users/l.fanchi/libs/ensembl_81/'
runOptions$varcontext$perlLibs = paste('/home/NFS/users/l.fanchi/perl5/lib/perl5','/home/NFS/users/l.fanchi/libs/bioperl-live',
																			 '/home/NFS/users/l.fanchi/libs/ensembl/modules','/home/NFS/users/l.fanchi/libs/ensembl-variation/modules', sep = ':')

# set neolution options
runOptions$neolution$numberOfWorkers = 3 # number of threads for predictions (each thread spawns 16 children during self-sim checking)

runOptions$neolution$rankCutoff = 3
runOptions$neolution$processingCutoff = 0.5
runOptions$neolution$xmer = c(9:11)

# set snpEff locations
runOptions$snpeff$path = '/home/NFS/users/l.fanchi/libs/snpEff'
runOptions$snpeff$build = 'GRCh38.81'
