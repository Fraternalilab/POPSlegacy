#! /usr/bin/R
#===============================================================================
# run POPS ovar all XML files of PDB databank
# (C) 2018 Jens Kleinjung
#===============================================================================

#_______________________________________________________________________________
## initialise resources
## synchronise XML database
#system("rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/XML/ ./XML");

## link POPS and database
#system("ln -s ~/POPS/src/pops");
#system("ln -s ~/database/XML/");

#_______________________________________________________________________________
## configure runs
## get input names of all subdirectories and 'xml.gz' files
dirnames = list.dirs(path = "./XML", full.names = FALSE);
filenames = sapply(dirnames, function(x) {
			list.files(path = paste("./XML/", x, sep = ""),
			full.names = FALSE, pattern = 'xml\\.gz$');
})

## create output structure
## each input file will have its own output directory to accommodate multiple
##   output files from POPS
#sapply(names(filenames), function(x) {
#  if (! identical(x, "")) { 
#    dir.create(paste("./JSON", x, sep = "/"));
#    sapply(1:length(filenames[[x]]), function(y) {
#      dir.create(paste("./JSON", x, filenames[[x]][y], sep = "/"));
#    });
#  }
#})

#_______________________________________________________________________________
## run all
sapply(names(filenames)[2:10], function(x) {
	if (! identical(x, "")) { 
		sapply(1:length(filenames[[x]]), function(y) {
			infile = paste("./XML", x, filenames[[x]][y], sep = "/");
			## shell command for POPSing current input file
			command = paste("./pops --pdbml", infile, "--zipped --jsonOut || exit 1"); 
			#print(command);
			try(system(command));
		});
	}
})

#===============================================================================


