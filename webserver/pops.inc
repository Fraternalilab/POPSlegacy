<?php
/*============================================================================*/
/* POPS website layout and functions                                          */
/* (C) 2007-2010 Jens Kleinjung, Franca Fraternali, Alessandro Pandini        */
/*============================================================================*/

$scriptdir = "/srv/www/htdocs/wiki/php";
$wikibaseurl = "mathbio.nimr.mrc.ac.uk/wiki";

/* variables */
$pops_input = <<<EOF
<br>
<form method="post" enctype="multipart/form-data">
<p>
	<font color='green'> INPUT OPTIONS<br></font>
	&nbsp; Enter a PDB identifier:
	  <input name="PDBID" value="" size="4" type="text"><br>
	&nbsp; <i><b>OR</b></i> Upload a PDB file:
	  <input name="PDBUP" value="" accept="*" size="50" type="file"><br>
	<br>
	<font color='green'> MODE OPTIONS<br></font>
	&nbsp; coarse-grained calculation <input name="RESMOD" value="1" type="checkbox"><br>
	&nbsp; probe radius               <input name="RADMOD" value="1.4" size="4" type="text"><br>
	<br>
	<font color='green'> OUTPUT OPTIONS<br></font>
	&nbsp; molecule composition       <input name="COMOUT" value="1" type="checkbox"><br>
	&nbsp; atom types                 <input name="TYPOUT" value="1" type="checkbox"><br>
	&nbsp; molecule topology          <input name="TOPOUT" value="1" type="checkbox"><br>
	&nbsp; POPS* area per chain        <input name="CHAOUT" value="1" type="checkbox"><br>
	&nbsp; POPS* area per residue      <input name="RESOUT" value="1" type="checkbox"><br>
	&nbsp; POPS* area per atom         <input name="ATOOUT" value="1" type="checkbox"><br>
	&nbsp; neighbour list             <input name="NEIOUT" value="1" type="checkbox"><br>
	<br>
</p>
<p>
	<input name="action" value="Submit" type="submit">
	<input name="mbjob[method]" value="pops" type="hidden">
	<a href="#Options">Help on Options</a>
</p>
</form>
EOF;

/*___________________________________________________________________________*/
/* overview */
$pops_desc = <<<EOF
<br>
<h2><a name="Overview">Overview</a></h2>
<p>
	POPS* is a fast algorithm to calculate solvent accessible surface areas (SASAs)
	of proteins and nucleic acids at atomic (default) and 
	residue (coarse-grained) level. Atomic and residue area parameters 
	have been optimised versus an accurate all-atom method (Naccess).
	Residue areas (coarse-grained) are simulated with a single sphere centered 
	at the C-alpha atom for amino acids and at the P atom for nucleotides.
</p>
<p>
	An analytical formula is used for the calculation of solvent accessibilities. 
	The formula is simple, easily derivable and fast to compute, therefore, 
	practical for use in molecular dynamics (MD) simulations as an approximation 
	to the first solvation shell.
	The default output contains:
	list of hydrophobic, hydrophilic and total contributions to the accessible 
	surface area of the entire molecule (all atoms), first model only.
</p>
<p>
	After submission a holding page provides access to the result page.
	Otherwise results (also of previous runs) can be accessed through the
	<a href="http://$wikibaseurl/Online_Tools#Retrieve">Retrieve function</a>.
</p>
<p>
	The POPS* program code can be downloaded 
	<a href="http://mathbio.nimr.mrc.ac.uk/wiki/Software#POPS.2A">here</a>.<br>
	Users publishing results obtained with the POPS* program should acknowledge 
	its use by citing <a href="#Fraternali_2002">Fraternali and Cavallo (2002)</a>
	and for the POPS* server <a href="#Cavallo_2003">Cavallo et al. (2003)</a>.
</p>
EOF;

$pops_parameters = <<<EOF
<br>
<h2><a name="Parameters">Parameters</a></h2>
The <b>optimised parameters</b> are listed here:<br>
<ul>
	<li>
		<a href="http://mathbio.nimr.mrc.ac.uk/tool/POPS/sasa_data_atom">atom parameters</a><br>
	</li>
</ul>
<ul>
	<li>
		<a href="http://mathbio.nimr.mrc.ac.uk/tool/POPS/sasa_data_coarse">residue parameters (coarse-grained)</a><br>
	</li>
</ul>
EOF;

/*___________________________________________________________________________*/
/* help on options */
$pops_help_options = <<<EOF
<br>
<h2><a name="Options">Help on Options</a></h2>
<dl>
	<dt><b>PDB structure file</b></dt>
	<dd>
		Entering a PDB identifier (lower case or upper case without the extension
	   .pdb) in the first text window will automatically perform the POPS*
	   calculation on the corresponding structure. Uploading of structure files 
	   by the user can be performed by specifying the file path and name in the 
	   second text window. Please ensure that atom and residue names conform to
	   the <a href="http://www.wwpdb.org/docs.html">PDB format</a>. 
	</dd>
	<dt><b>coarse-grained</b></dt>
	<dd>
		Calculation of accessible surface area on residue level.
	</dd>
	<dt><b>probe radius</b></dt>
	<dd>
	   The radius of the surface probe (solvent molecule) in Angstrom.
	</dd>
	<dt><b>multiple models (option currently disabled)</b></dt>
	<dd>
		By default the POPS* calculation will be performed on the first model only, 
		if a structure file containing multiple models is chosen/uploaded. Activation 
		of the 'Multiple models' button triggers POPS* calculation of all models.
	</dd>
	<dt><b>molecule composition</b></dt>
	<dd>
		Number of chains, standard residues, total residues and atoms.
	</dd>
	<dt><b>atom types</b></dt>
	<dd>
		List of atom types and their POPS* parameters.
	</dd>
	<dt><b>molecule topology</b></dt>
	<dd>
		Number of bonds, angles, torsions and non-bonded interactions.
	</dd>
	<dt><b>POPS* area per atom</b></dt>
	<dd>
		List of atom areas [A<sup>2</sup>], number of overlaps and atom grouping 
		(1: positive, 2: negative, 3: polar, 4: aromatic, 5: aliphatic).
	</dd>
	<dt><b>POPS* area per residue</b></dt>
	<dd>
		List of residue areas [A<sup>2</sup>] (hydrophilic, hydrophobic, total) and number of overlaps.
	</dd>
	<dt><b>POPS* area per chain</b></dt>
	<dd>
		List of chain areas [A<sup>2</sup>] (hydrophilic, hydrophobic, total).
	</dd>
	<dt><b>neighbour list</b></dt>
	<dd>
		List of neighbours per atom, provided in a separate output file.
	</dd>
</dl>
EOF;

$pops_help_parameters = <<<EOF
<br>
<h2><a name="Background Information">Background Information</a></h2>
The area equation is defined by a product &#928; of terms that estimate
the reduction of SASA of atom i by the overlap with its neighbours j (Hasel at al., 1988):<br><br>
&#928;<sub>i=1</sub><sup>N</sup> [ 1 - (p<sub>i</sub> p<sub>ij</sub> b<sub>ij</sub>(r<sub>ij</sub>) / S<sub>i</sub>)&nbsp;].<br><br>
i is the atom for which the POPS* area is computed, j is any of N neighbour atoms.<br>
p<sub>i</sub> is an atom type specific SASA parameter.<br>
p<sub>ij</sub> is a sphere overlap parameter depending on the degree of bonding between i and j (1-2, 1-3, 1-4, 1-(>4)).<br>
b<sub>ij</sub> is a geometric construct based on the radii and distance (r<sub>ij</sub>)
	of i and j.<br>
S<sub>i</sub> is the SASA of the free atom i (no neighbours). 
<p>
	The atom specific parameters (radii, SASAs) are listed in the 
	<a href="#Parameters">Parameters</a> files for the atoms of all standard PDB residues,
	followed by the connectivity parameters (p<sub>ij</sub>, b<sub>ij</sub>) and the
	solvent radius for water.

</p>
EOF;

$pops_references = <<<EOF
<br>
<h2><a name="References">References</a></h2>
<ol>
	<li>
		<a name="Fraternali_1996"></a>Fraternali, F. and van Gunsteren, W.F.<br>
		An efficient mean solvation force model for use in molecular dynamics 
		simulations of proteins in aqueous solution.<br>
		<i>Journal of Molecular Biology</i> <b>256</b> (1996) 939-948. &nbsp;
		<a href="http://www.ncbi.nlm.nih.gov/pubmed/8601844">[PubMed Abstract]</a>
	</li>
	<p></p>
	<li>
		<a name="Fraternali_2002"></a>Fraternali, F. and Cavallo, L.<br>
		Parameter optimized surfaces (POPS*): analysis of key interactions and 
		conformational changes in the ribosome.<br>
		<i>Nucleic Acids Research</i> <b>30</b> (2002) 2950-2960. &nbsp;
		<a href="http://nar.oupjournals.org/cgi/reprint/30/13/2950">[PDF]</a>
	</li>
	<p></p>
	<li>
		<a name="Cavallo_2003"></a>Cavallo, L., Kleinjung, J. and Fraternali, F.<br>
		POPS: A fast algorithm for solvent accessible surface areas at atomic 
		and residue level.<br>
		<i>Nucleic Acids Research</i> <b>31</b> (2003) 3364-3366. &nbsp; 
		<a href="http://nar.oxfordjournals.org/cgi/reprint/31/13/3364">[PDF]</a>
	</li>
	<p></p>
	<li>
		<a name="Kleinjung_2005"></a>Kleinjung, J. and Fraternali, F.<br>
		POPSCOMP: an automated interaction analysis of biomolecular complexes.<br>
		<i>Nucleic Acids Research</i> <b>33</b> (2005) W342-W346. &nbsp;
		<a href="http://nar.oxfordjournals.org/cgi/reprint/33/suppl_2/W342">[PDF]</a>
	</li>
</ol>
<br>
EOF;

$pops_footer = <<<EOF
<hr>
The POPS server was developed by Franca Fraternali, the POPS logo was designed 
by Domenico Fraternali. &nbsp;
<a href="mailto:ffranca@kcl.ac.uk,jkleinj@nimr.mrc.ac.uk">Contact Authors</a>
EOF;

$pops_warning = <<<EOF
<p>
	<font color='red'>
	PLEASE NOTE THAT THE POPS WEB SERVER IS BEING UPDATED CURRENTLY AND MIGHT 
	NOT WORK PROPERLY UNTIL THIS NOTE DISAPPEARS!
	</font>
</p>
EOF;

$pops_version = <<<EOF
<p>
	Server Version 1.0.6; POPS* Version 1.5.3
</p>
EOF;

$pops_message = <<<EOF
<p>
	<font color='red'>
	</font>
</p>
EOF;


/* functions */
function parse_results($jobid, $wwwdir)
{
	// parse POPS result page
	$myFile = $wwwdir . "/job/" . $jobid . "/result.html";
	$output = file($myFile);

	// append POPS output to result string
	$resultstring .= join("<br>\n",$output);

	return $resultstring;
}

function check_input($_FILES, $_POST)
{
	return false;
}

?>
