#! /usr/bin/R
#===============================================================================
# run POPS ovar all XML files of PDB databank
# (C) 2018 Jens Kleinjung
#===============================================================================

#_______________________________________________________________________________
#' sadata: An S4 class input/output control.
#' @slot dirnames : array of directory names
#' @slot filenames : list of input files, list elements are the directories
#' @slot filenames : data frame of input files, first column contains directories
#' @slot outpath : list of output paths
#' @slot command : command lines for shell
ioctrl <- setClass(
  "ioctrl",
  
  slots = c(
    dirnames = "character",
    filenames = "list",
    filenames_df = "data.frame",
    outpath = "list",
    command = "list"
  )
)

ioctrl_o = ioctrl();

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
ioctrl_o@dirnames = list.dirs(path = "./XML", full.names = FALSE);
ioctrl_o@filenames = sapply(dirnames, function(x) {
			list.files(path = paste("./XML/", x, sep = ""),
			          full.names = FALSE, pattern = 'xml\\.gz$');
})

#_______________________________________________________________________________
## function to coerce list of directories and filenames to data frame
coerce_filenames = function(ioctrl_o) {
  dir_file_l = sapply(names(ioctrl_o@filenames), function(x) {
    if (! identical(x, "")) {
      sapply(1:length(ioctrl_o@filenames[[x]]), function(y) {
        dir_file = c(x, ioctrl_o@filenames[[x]][y]);
        return(dir_file);
      });
    }
  })
  dir_file_df = unlist(dir_file_l);
  dim(dir_file_df) = c(2, (length(dir_file_df) / 2));
  rownames(dir_file_df) = c("dir", "file");
  return(dir_file_df);
}

## create POPS command
make_command = function(x, y) {
  
}

dirfile = coerce_filenames(ioctrl_o);

make_command = function(z) {
  infile = paste("./XML", z$x, z$y, sep = "/");
  outdir = paste("./JSON", z$x, z$y, sep = "/");
  ## shell command for POPSing current input file
  command = paste("./pops --pdbml", infile, "--outDirName", outdir, "--zipped --jsonOut --silent || exit 1"); 
  return(command);
}

tt = make_command(traverse_dir(dirnames, filenames));

#_______________________________________________________________________________
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
run_results_l = sapply(names(filenames), function(x) {
	if (! identical(x, "")) { 
		sapply(1:length(filenames[[x]]), function(y) {
			infile = paste("./XML", x, filenames[[x]][y], sep = "/");
			outdir = paste("./JSON", x, filenames[[x]][y], sep = "/");
			## shell command for POPSing current input file
			command = paste("./pops --pdbml", infile, "--outDirName", outdir, "--zipped --jsonOut --silent || exit 1"); 
			#print(command);
			try(system(command));
		});
	}
})

#_______________________________________________________________________________
## validate all

#_______________________________________________________________________________
## upload all

#===============================================================================


