#!/usr/bin/perl

#===============================================================================
# pops.pl : Retrieve complex PDB file from identifier and POPS it;
#			should be re-written (Jens 2009)
#			POPS server Version 1.0.6 (using POPSc Version 1.5.3)
# (C) 2002-2012 Franca Fraternali and Jens Kleinjung
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#===============================================================================

#______________________________________________________________________________
# main function
#_______________________________________________________________________________
# read the command line argument: sequence/structure number
# default settings
$query  = "";
# MODE OPTIONS
$mulmod = 0; # multi-model calculation
$resmod = 0; # residue-only calculation 
$radmod = 0; # probe radius
# OUTPUT OPTIONS
$comout = 0; # molecule composition
$typout = 0; # atom types
$topout = 0; # molecule topology
$atoout = 0; # POPS area per atom
$resout = 0; # POPS area per residue
$chaout = 0; # POPS area per chain
$neiout = 0; # neighbour list
$parout = 0; # parameter list

$inputFlags = "";
$local  = 1; #
$error  = 0; #

$myhome = "/srv/www/htdocs/mathbio/tool/POPS";
$pwddir = $ENV{'PWD'};
@pwdpath = split('/',$pwddir);
$jobId = $pwdpath[-1];
$myjobdir = "../job/" . "$jobId";

#_______________________________________________________________________________
# initialise status
system("chmod 311 /tmp/job/$jobId");
system("echo '-3' > status");

#_______________________________________________________________________________
# parameters from website
&read_web_parameters;

# load PDB file
# from PDB resource
if (($pdbid ne "") && ($local == 1)) {
    $pdbid =~ tr/A-Z/a-z/; # translate to lower case
    $query = $pdbid; # define query name
    &pdb_download(); # load PDB file from PDB resource 
}

# from user as uploaded file
if ($local == 0) {
    $query = $pdbup; # define query name
    &pdb_upload(); # PDB file uploaded by the user
}

# calculate POPS
$inputFlags .= " --pdb $query";
&pops_it();

# write result page
&print_result;

# change permission of job directory and parameter files
system("chmod 444 parameters pops.out result.html stderr stdout");

# flag up success
system("echo '0' > status");

#_______________________________________________________________________________
# Subroutines
#_______________________________________________________________________________

#_______________________________________________________________________________
# read WEB parameters
sub read_web_parameters
{
    # file containing parameters from WEB site
    open(PARAM, "parameters");
    
    while(<PARAM>) {
		# PHP ID
		if ($_ =~ /mbjob.uniqid/) {
			chomp $_;
			($tag, $php_id) = split(/\=/, $_);
			print("PHPSESSID = $php_id \n");
		}
		# pasted PDB identifier
		if ($_ =~ /^PDBID/) {
			chomp $_;
			($tag, $pdbid) = split(/\=/, $_);
			print ("PDBID = $pdbid\n");
			if ($pdbid ne "") { $local = 1; }
		}
		# pasted PDB upload file name, check only if no PDBID specified
		if ($pdbid eq "") {
			if ($_ =~ /^file.PDBUP.name/) {
				chomp $_;
				($tag, $pdbup) = split(/\=/, $_);
				print ("PDBUP = $pdbup\n");
				if ($pdbup ne "") { $local = 0; }
			}
		}
		# multi model switch
		#if ($_ =~ /^MULMOD/) {
		#	chomp $_;
		#	($tag, $mulmod) = split(/\=/, $_);
		#	if ($mulmod) { $inputFlags .= " --multiModel"; }
		#	print ("MULMOD = $mulmod\n");
		#}
		# coarse grained calculation
		if ($_ =~ /^RESMOD/) {
			chomp $_;
			($tag, $resmod) = split(/\=/, $_);
			if ($resmod) { $inputFlags .= " --coarse"; }
			print ("RESMOD = $resmod\n");
		}
		# probe radius 
		if ($_ =~ /^RADMOD/) {
			chomp $_;
			($tag, $radmod) = split(/\=/, $_);
			if ($radmod) { $inputFlags .= " --rProbe $radmod"; }
			print ("RADMOD = $radmod\n");
		}
		# molecule composition
		if ($_ =~ /^COMOUT/) {
			chomp $_;
			($tag, $comout) = split(/\=/, $_);
			if ($comout) { $inputFlags .= " --compositionOut"; }
			print ("COMOUT = $comout\n");
		}
		# atom types
		if ($_ =~ /^TYPOUT/) {
			chomp $_;
			($tag, $typout) = split(/\=/, $_);
			if ($typout) { $inputFlags .= " --typeOut"; }
			print ("TYPOUT = $typout\n");
		}
		# molecule topology
		if ($_ =~ /^TOPOUT/) {
			chomp $_;
			($tag, $topout) = split(/\=/, $_);
			if ($topout) { $inputFlags .= " --topologyOut"; }
			print ("TOPOUT = $topout\n");
		}
		# POPS area per atom
		if ($_ =~ /^ATOOUT/) {
			chomp $_;
			($tag, $atoout) = split(/\=/, $_);
			if ($atoout) { $inputFlags .= " --atomOut"; }
			print ("ATOOUT = $atoout\n");
		}
		# POPS area per residue
		if ($_ =~ /^RESOUT/) {
			chomp $_;
			($tag, $resout) = split(/\=/, $_);
			if ($resout) { $inputFlags .= " --residueOut"; }
			print ("RESOUT = $resout\n");
		}
		# POPS area per chain
		if ($_ =~ /^CHAOUT/) {
			chomp $_;
			($tag, $chaout) = split(/\=/, $_);
			if ($chaout) { $inputFlags .= " --chainOut"; }
			print ("CHAOUT = $chaout\n");
		}
		# neighbour list
		if ($_ =~ /^NEIOUT/) {
			chomp $_;
			($tag, $neiout) = split(/\=/, $_);
			if ($neiout) { $inputFlags .= " --neighbourOut"; }
			print ("NEIOUT = $neiout\n");
		}
		# parameter list
		if ($_ =~ /^PAROUT/) {
			chomp $_;
			($tag, $parout) = split(/\=/, $_);
			if ($parout) { $inputFlags .= " --parameterOut"; }
			print ("PAROUT = $parout\n");
		}
    }
}
    
#_______________________________________________________________________________
# download PDB file from PDB resource
sub pdb_download
{
    # link to pdb structure and unzip
    $pdbquery = "$query" . ".pdb.gz";
    $unzipquery = "$query" . ".pdb";
    print("Retrieving $pdbquery\n");
    $fileUrl = "http://www.rcsb.org/pdb/files/" . $pdbquery;
    $saveAs = "./" . $pdbquery;
    print("wget -v -O $saveAs $fileUrl\n");
    $error = system("export http_proxy=http://proxy.nimr.mrc.ac.uk:8000; wget -v -O $saveAs $fileUrl");
	print("Download error status: $error\n");
    if ($error > 0) {
		print STDERR ("Warning: Failed to download $query!\n");
		system("echo '1' > status");
		exit(1);
    }
    system("gunzip $pdbquery");
    system("mv $unzipquery $query");
}

#_______________________________________________________________________________
# uploaded PDB file
sub pdb_upload
{
    if (-e "$pdbup") { print("Found $pdbup\n"); }
    else {
		print STDERR ("Warning: cannot find $pdbup\n");
		system("echo '1' > status");
		exit(1);
	}
}

#_______________________________________________________________________________
# POPS calculation
sub pops_it
{
    # execute POPS
    print("$myhome/pops $inputFlags\n");
    $error = system("($myhome/pops $inputFlags > stdout) >& stderr");
    
    if ($error > 0) {
		print STDERR ("Warning: Failed calculation on $query!\n");
		system("echo '2' > status");
		exit(1);
    }
}

#______________________________________________________________________------
# write WEB results page
sub print_result
{
    open(POPOUT, ">result.html") || die "Cannot create 'result.html'\n";

    print POPOUT "<font face=\"Arial,Helvetica\" color=\"#000066\" size=+1>\n",
				"<b>Job Details</b></font><br>\n",
				"POPS Server Version: <b>1.0.6</b>\n", 
				"POPSc Version: <b>1.5.3</b>\n",
				"Job ID: <b>$php_id</b>\n";

	open(DATE,"date |") || print POPOUT "CGI-ERROR - date?\n";
	$line = <DATE>;
	chomp($line);
	close(DATE);

	print POPOUT "Date: <b>$line</b>",
				"<p><font color='green'> INPUT OPTIONS </font>\n",
				"&nbsp; PDB $query\n</p>", 
				"<p><font color='green'> MODE OPTIONS </font>\n",
				"&nbsp; coarse grained calculation $resmod\n",
				"&nbsp; probe radius $radmod\n</p>",
				"<p><font color='green'> OUTPUT OPTIONS </font>\n",
				"&nbsp; molecule composition $comout\n",
				"&nbsp; atom types $typout\n",
				"&nbsp; molecule topology $topout\n",
				"&nbsp; POPS area per atom $atoout\n",
				"&nbsp; POPS area per residue $resout\n",
				"&nbsp; POPS area per chain $chaout\n",
				"&nbsp; neighbour list $neiout\n",
				"&nbsp; parameter list $parout</p>\n",
				"<hr>";

    print POPOUT "<font face=\"Arial,Helvetica\" color=\"#000066\" size=+1>\n",
				"<b>Result Details</b></font><br>\n",
				"<p><font color='green'> OUTPUT FILES </font>\n",
				"&nbsp; <a href=$myjobdir/pops.out>POPS* output</a>\n",
				"&nbsp; <a href=$myjobdir/stdout>stdout (infos)</a>\n",
				"&nbsp; <a href=$myjobdir/stderr>stderr (warnings, errors)</a>\n";

    # print the main result line
	# read the result file and transform to HTML format
	print POPOUT "</p><br>",
				 "<p><font color='green'>OUTPUT TEXT</font>\n",
				 "<font face=\"Courier\">\n";
	open(DAT, "pops.out") || die "Cannot read pops.out\n";
	while (<DAT>) {
		chomp;
		s/\t/        /g;
		s/ /&nbsp;/g;
		print POPOUT "$_\n";
	}
	close(DAT);
	print POPOUT "</font></p><hr>\n"; 

	close(POPOUT);
}

