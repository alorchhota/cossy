<?php 
require_once './rserve-php/settings/config.php';
require './rserve-php/Connection.php';
require './cossy-helper.php';

$cossydebug = FALSE;

if ($cossydebug) 
	require './rserve-php/segments/head.php';


$starttime = time();

//// make output
function constructCossyOutput($x, $mis, $outtype="plain") {
	$outputText = "";
	
	if($x["status"]=="ERROR"){
		$outputText = constructErrorOutput($x['error'], $outtype);
		return $outputText;
	}
	
    if($outtype=="plain"){
		$outputText = $x["status"] . "\r\n" . $x["classes"][0] . "|" . $x["classes"][1] . "\r\n" . $x["network"];
	} elseif($outtype=="xml"){
		$outputText = "<cossy><status>OK</status><classes><positive>" . $x["classes"][0]  . "</positive><negative>" . $x["classes"][1]  . "</negative></classes><network>" .  $x["network"] . "</network></cossy>";
	} else{
		$outputText = $outputText . constructErrorOutput($x, 'xml');
	}
	
	return $outputText;
}

function constructErrorOutput($errormsg, $outtype="plain"){
	$outputText = "";
	if($outtype=="plain"){
		$outputText = "CossyException: " . $errormsg;
	} elseif($outtype=="xml"){
		$outputText = "<cossy><status>ERROR</status><error>" . $errormsg . "</error></cossy>";
	}
	
	return $outputText;
	
}

//// test output
function mydump($x, $title=NULL) {
    if($title) {
        echo '<h4>'.$title.'</h4>';
    }
    echo '<pre>';
    var_dump($x);
	echo '</pre>';
}

//////////////////////// read variables //////////////////////

$inputOK = TRUE;
$errormsg= "No error.";

$profiletype = $_POST["profiletype"];
if($profiletype == NULL){
	$profiletype = "microarray";
}	
$profiletype = strtolower($profiletype);
	
$format = $_POST["format"];
if($format==NULL ||  strtolower($format) !="xml"){
	$format="plain";
}	
$format = strtolower($format);

	
if (!fileUploaded("gctfile") || !fileUploaded("clsfile") || ($profiletype=="microarray" && !fileUploaded("chipfile")))
{
	$inputOK = FALSE;
	$errormsg = "Sorry, you need to upload all the files.";
	
}


$output = "";



$chipfilepath = "";
//// TODO: more checking is required
//// check if the file sizes are within the limit
if ($inputOK)
{
	
		
	//// if files are OK, start processing.
	try {
		//// connect to Rserve
		if ($cossydebug) echo('<p>Connecting to Rserve ... ');
		
		$r = new Rserve_Connection(RSERVE_HOST, RSERVE_PORT);
		if ($cossydebug) echo (' OK</p>');
		
		//// read the uploaded filepaths
		$gctfilepath = $_FILES["gctfile"]["tmp_name"];
		$clsfilepath = $_FILES["clsfile"]["tmp_name"];
		
		$chipfilepath = "NA";
		if($profiletype=="microarray")	$chipfilepath = $_FILES["chipfile"]["tmp_name"];

		//// read inputs
		$mis = (int) $_POST["mis"]; 	// should read from POST request
		$net = $_POST["network"];
		
		
		//// perform cossy analysis
		// load icossy library
		$cmd = "library(icossy)";
		$x = $r->evalString($cmd, Rserve_Connection::PARSER_NATIVE);
		
		// run icossy function
		$formatedchipfile = ($profiletype=="microarray" ? "'". $chipfilepath."'" : $chipfilepath);
		$cmd = "icossy(gctfile = '" . $gctfilepath . "', chipfile = " . $formatedchipfile .", clsfile = '" . $clsfilepath . "', network = '" . $net . "', nmis = " . $mis . ", frank = T, qnorm = F, ztrans = F, sig.test = 'ttest')";
		$x = $r->evalString($cmd, Rserve_Connection::PARSER_NATIVE);
		
		// build output
		$output = constructCossyOutput($x, $mis, $format);
		//$output = "";
		
		if ($cossydebug) echo ( "<p>cleaning resources.... </p>");
		unlink($gctfilepath );
		unlink($clsfilepath );
		if($profiletype=="microarray")	unlink($chipfilepath );
		$r->close();
		
		
	} catch(Exception $e) {
		$errmsg = "CossyException: [" . get_class($e) . "] " .  $e->getMessage();
		$output = constructErrorOutput($errmsg, $format);
		
		try{
			// free resources
			$r->close();
			if($profiletype=="microarray")	unlink($chipfilepath );
			unlink($clsfilepath );
			unlink($gctfilepath );
		} catch(Exception $e2){
		}
	}
}
else
{
	//echo "CossyException: " . $errormsg;
	$errmsg = "CossyException: " . $errormsg;
	$output = constructErrorOutput($errmsg, $format);
}


if($format=="xml")
	header("Content-Type: text/xml");
else
	header("Content-Type: text/plain");
	
if($cossydebug)
	header("Content-Type: text/html");


echo $output;

$endtime = time();
if ($cossydebug) {
	echo "<p>Execution time: ";
	print date("i:s", $endtime-$starttime);
	echo "</p>";
}

if ($cossydebug) require 'rserve-php/segments/foot.php';


?>
