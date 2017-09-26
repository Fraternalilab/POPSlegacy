<?php
//=============================================================================
// MB_submitjob.inc : script to process web server jobs 
// Mathematical Biology Division, NIMR, London
//=============================================================================

global $uniqId;
global $workingDir;
global $method;
global $script;
global $type;
global $jobDir;
global $job;
global $wikidir;
global $phpdir;

$debug = 1;
$number = 0;
$message = "";
$wikidir = "/srv/www/htdocs/wiki";
$phpdir = "/srv/www/htdocs/wiki/php";

//_____________________________________________________________________________
function setError($message, $number = 0)
{
	$message = $message;
	$number = $number;
	if ($debug) 
		echo $message;
	return $number;
}

//_____________________________________________________________________________
function writeParameterArray($fh, $parameters, $prefix = "")
{
	// recursive function to write out array of parameters
	//	in the form key1.key2=value
  
	foreach($parameters As $key => $value)
	{
		if (is_array($value)) 
			writeParameterArray($fh,&$value,$prefix . $key . ".");
		else
			fputs($fh,$prefix . $key . "=" . $value . "\n");
	}
}

//_____________________________________________________________________________
function stageIn()
{
	// 1. create working directory
	// 2. write a file 'parameters' containing the contents of $_POST
	// 3. upload any files
  
	global $uniqId;
	global $workingDir;
	global $jobDir;
	global $job;
	global $wikidir;

	$jobDir = $wikidir . "/job";

	//unique job identifier
	if(isset($_POST['mbjob']['uniqid']))
	{ 
		$uniqId = $_POST['mbjob']['uniqid'];
	} else {
		$uniqId=md5(uniqId(time()));
	}

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

	// write parameters
	chdir($workingDir);

	if (isset($_POST['mbjob']['append']) && $_POST['mbjob']['append'] == "true")
	{ 
		$mode = "a";
	} else { 
		$mode = "w";
	}

	$fh=fopen("parameters","$mode");
	if(!$fh)
	{
		$err = "Failed to open parameter file";
		return $err;
	}
	reset($_POST);
	writeParameterArray($fh,&$_POST);

	if(isset($_FILES))
	{
		writeParameterArray($fh,&$_FILES,"file.");
		fputs($fh,"mbjob.uniqid=".$uniqId);
		fclose($fh); 
	  
		// upload files
		foreach($_FILES as $file) 
		{
			if(!empty($file['name']) && !move_uploaded_file($file['tmp_name'],$workingDir."/".$file['name']))
			{
				$err = "Failed to upload files";
				return $err;
			}
		}
	}
	  
	$fh=fopen("status","w");
	if (!$fh)
	{
		$err="Failed to open status file";
		return $err;
	}
	fputs($fh,"-1");
	fclose($fh);
	  
	return false;
}

//_____________________________________________________________________________
function SubmitJob()
{
	global $uniqId;
	global $workingDir;
	global $jobDir;
	global $job;
	global $wikidir;
	global $phpdir;
	     
	$stagedIn=stageIn();
  
	if(! $stagedIn)
	{
		chdir($workingDir);
		exec("echo '-2' > status");
		exec("(perl $phpdir/script/$job.pl > stdout) >& stderr");
	} else {
		$err = "Submission has failed because: ".$stagedIn;
		exec("echo '11' > status");
		return $stagedIn;
	}
}

