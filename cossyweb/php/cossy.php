<?php 
require_once './rserve-php/settings/config.php';
require './rserve-php/Connection.php';
require './cossy-helper.php';

$cossydebug = TRUE;

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
		$outputText = "";
		for($i=0; $i<$mis; $i++){
			$misobj = $x['results'][$i];
			for($j=0; $j<5; $j++)
				 $outputText = $outputText . $misobj['representative.genes'][$j] . "|" . $misobj['representative.expression'][$j] . "\r\n";
		}
		
	} elseif($outtype=="xml"){
		$outputText = "<cossy><status>OK</status><results>";
		
		for($i=0; $i<$mis; $i++){
			$misobj = $x['results'][$i];
			$outputText = $outputText . "<mis>";
			for($j=0; $j<5; $j++){
				$outputText = $outputText . "<gene>";
				$outputText = $outputText . "<symbol>". $misobj['representative.genes'][$j] . "</symbol><expression>" . $misobj['representative.expression'][$j] . "</expression>";
				$outputText = $outputText . "</gene>";
			}
			$outputText = $outputText . "</mis>";
		}
		$outputText = $outputText . "</results></cossy>";
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
		
		//// get the currect R directory
		$cmd = 'getwd()';
		$curdir = $r->evalString($cmd, Rserve_Connection::PARSER_NATIVE); 
		if ($cossydebug) echo ("Current directory: ". $curdir . "</p>");
		
		// //// move the uploaded files to the current R direcoty.
		// $gctfilepath = $curdir."/userdb.gct";
		// $clsfilepath = $curdir."/userdb.cls";
		
		// if($profiletype=="microarray")	$chipfilepath = $curdir."/userdb.chip";
		
		// move_uploaded_file($_FILES["gctfile"]["tmp_name"], $gctfilepath );
		// move_uploaded_file($_FILES["clsfile"]["tmp_name"], $clsfilepath );
		// if($profiletype=="microarray")	move_uploaded_file($_FILES["chipfile"]["tmp_name"], $chipfilepath );

		// //// read inputs
		// $mis = (int) $_POST["mis"]; 	// should read from POST request
		// $net = $_POST["network"];
		
		
		// //// run cossy
		// $cmd = "cossy(dataset='userdb', network='" . $net . "', T=" . $mis . ", data.type='" . $profiletype . "')";
		// //$cmd = "mis";
		// $x = $r->evalString($cmd, Rserve_Connection::PARSER_NATIVE);
		// $output = constructCossyOutput($x, $mis, $format);
		// //$output = "";
		
		// if ($cossydebug) echo ( "<p>cleaning resources.... </p>");
		// unlink($gctfilepath );
		// unlink($clsfilepath );
		// if($profiletype=="microarray")	unlink($chipfilepath );
		$r->close();
		
		
		
		
	} catch(Exception $e) {
		// echo "CossyException: [" . get_class($e) . "] " .  $e->getMessage();
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
