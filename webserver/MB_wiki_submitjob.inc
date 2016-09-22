<?php
/*============================================================================*/
/*                                                                            */
/* MB_wiki_submitjob                                                          */
/*                                                                            */
/* (C) 2007 Jens Kleinjung and Alessandro Pandini                             */
/*                                                                            */
/*============================================================================*/

$phpdir = "/srv/www/htdocs/wiki/php";
$mathbiodir = "http://mathbio.crick.ac.uk";

include("$phpdir/MB_submitjob.inc");

/*___________________________________________________________________________*/
/* variables scope */
global $wgOut;
global $action;

/*___________________________________________________________________________*/
/* page layout variables */
$results_name = <<<EOD
<i><font color="#4a96ad">Results </font></i><font color="#940c63">Page</font>
EOD;

$results_page_title = <<<EOD
<div align="center">
    <h1>$results_name</h1>
</div>
EOD;

$results_desc = <<<EOF
<p>The server has generated a temporary <i><b>$job</b></i> subdirectory with your submitted data and results.<br>
These can be retrieved by providing the following unique job ID. Record this job ID for later access to the page.<br></p>
EOF;

/*___________________________________________________________________________*/
/* Action Submit */
if ($action == "Submit"){

	$returned_value = false;

	/* check input consistency */
	$method = $_POST['mbjob']['method'];
	include($phpdir . "/" . $method . ".inc");
	$returned_value = check_input($_FILES, $_POST);
	
	/* submit job request and capture results only if input is correct */
	if (!$returned_value)
	{
		$returned_value = SubmitJob();
	}

	$results_input = <<<EOD
	<form method="POST" action="Retrieve Results">
   	<h4>        
      	<br>Enter job ID:
      	<br><input name="JOBID" type="text" value="$uniqId"/>
      	<br><input type="submit" name="action" value="Show Results"/>
      </h4>
   </form>
EOD;
		
	if (!$returned_value)
	{
		$outstring = "$results_page_title<p>$results_desc</p><p>Job ID: $uniqId</p><hr><p>$results_input</p>";
	}else{
		$outstring = "$results_page_title<p>Error: $returned_value</p>Please resubmit your job.";
	}
}

/*___________________________________________________________________________*/
/* Action Submit Example */
if ($action == "Submit Example"){

	$jobDir = $wikidir . "/job";

	//unique job identifier
	$uniqId=md5(uniqId(time()));

	// create working dir - if required
	if(!is_writeable($jobDir))
	{
		$err="Problem with job directory - " . $jobDir;
		return $err;
	}

	$workingDir= $jobDir . "/" . $uniqId;
	umask(000);
	if(!is_dir($workingDir) && (mkdir($workingDir,0777))==FALSE)
	{
		$err="Failed to create working dir $workingDir";
		return $err;
	}	  

	// set the method from the form
	if (isset($_POST['mbjob']['method']))
	{ 
		$job = $_POST['mbjob']['method'];	
	}

	// execute example script
	chdir($workingDir);
	exec("echo '1' > status");
	exec("perl $phpdir/script/$job.example.pl >& output &");

	$results_input = <<<EOD
	<form method="POST" action="Retrieve Results">
   	<h4>        
      	<br>Enter job ID:
      	<br><input name="JOBID" type="text" value="$uniqId"/>
      	<br><input type="submit" name="action" value="Show Results"/>
      </h4>
   </form>
EOD;

	$outstring = "$results_page_title<p>$results_desc</p><p>Job ID: $uniqId</p><hr><p>$results_input</p>";
}

/*___________________________________________________________________________*/
/* Action Retrieve Results */
if ($action == "Retrieve Results"){
	$results_input = <<<EOD
	<form method="POST" action="Retrieve Results">
   	<h4>        
      	<br>Enter job ID:
      	<br><input name="JOBID" type="text" value=""/>
      	<br><input type="submit" name="action" value="Show Results"/>
      </h4>
   </form>
EOD;
	$outstring = "$results_page_title<p><p>$results_input</p>";
}

/* parse errors */
function parse_errors($jobid, $wikidir)
{
    // parse POPS result page
    $myFile = $wikidir . "/job/" . $jobid . "/stderr";
    $erroroutput = file($myFile);

    // append error output to result string
    $resultstring .= join("<br>\n",$erroroutput);

    return $resultstring;
}

/*___________________________________________________________________________*/
/* Action Show Results */
if ($action == "Show Results"){
	
	$jobid = $_POST['JOBID'];	

	$fh=fopen($wikidir . "/job/" . $jobid . "/parameters", 'r');
	while($line = fgets($fh)){
	   list($tag, $value) = split("=", $line);
	   if ($tag == "mbjob.method")
	   {
	   	$method = rtrim($value);
	   }		 
	}	   
	fclose($fh);

	include($phpdir . "/" . $method . ".inc");
	$resultstring = parse_results($jobid, $wikidir);
	$status_file = $wikidir . "/job/" . $jobid . "/status";


	/* passthru("grep 2 $status_file >& /dev/null", $greperr); */
	$grepval = `cat $status_file`;
	if ($grepval == 0){
		$progress_string = "Job successfully completed!";
	}
	if ($grepval < 0) {
		$progress_string = "Job computation running: reload the page later again";
	}
	if ($grepval > 0){
		$error_output = $mathbiodir . "/job/" . $jobid . "/stderr";
		$progress_string = "Job failed! Please read the error report:";
		$resultstring = parse_errors($jobid, $wikidir); 
	}

	$outstring = "$results_page_title<p><center><h2>$progress_string</h2></center></p><p>$resultstring</p>";
}

/*___________________________________________________________________________*/
/* Webpublish the output */
$wgOut->addHTML($outstring);
?>

