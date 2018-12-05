#! /usr/bin/R
#===============================================================================
# run POPS ovar all XML files of PDB databank
# (C) 2018 Jens Kleinjung
#===============================================================================

#_______________________________________________________________________________
#' sadata: An S4 class input/output control.
#' @slot dirnames : array of directory names
#' @slot filenames : list of input files, list elements are the directories
#' @slot filenames : matrix of directory and file names in rows
#' @slot outpath : list of output paths
#' @slot command : command lines for shell
ioctrl <- setClass(
  "ioctrl",
  
  slots = c(
    dirnames = "character",
    filenames = "list",
    filenames_df = "matrix",
    outpath = "character",
    command = "character"
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
## function to coerce list of directories and filenames to matrix
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

## populate matrix of directory and file names
ioctrl_o@filenames_df = coerce_filenames(ioctrl_o);

#_______________________________________________________________________________
## create output directory structure
## each input file will have its own output directory to accommodate multiple
##   output files from POPS
ioctrl_o@outpath = apply(ioctrl_o@filenames_df, 2, function(x) {
  outpath = paste("./JSON", x[1], x[2], sep = "/");
  dir.create(paste(outpath, showWarnings = FALSE, recursive = TRUE));
  return(outpath);
})

#_______________________________________________________________________________
## create POPS commands
ioctrl_o@command

make_command = function(z) {
  infile = paste("./XML", z$x, z$y, sep = "/");
  outdir = paste("./JSON", z$x, z$y, sep = "/");
  ## shell command for POPSing current input file
  command = paste("./pops --pdbml", infile, "--outDirName", outdir, "--zipped --jsonOut --silent || exit 1"); 
  return(command);
}

tt = make_command(traverse_dir(dirnames, filenames));

#_______________________________________________________________________________
## run all command lines
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


