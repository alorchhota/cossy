<html>
<body>

	<form action="cossy-dev.php" method="post" enctype="multipart/form-data">
		<label for="gctfile">Gene Expression Profiles </label>
		<input type="file" name="gctfile" id="gctfile"> <br>
		<label for="clsfile">Select Class Label </label>
		<input type="file" name="clsfile" id="clsfile"><br>
		<label for="mapfile">Select Chip File </label>
		<input type="file" name="mapfile" id="mapfile"><br>
		<input type="hidden" name="mis" id="mis" value="5"><br>
		<input type="hidden" name="network" id="network" value="kegg"><br>
		<input type="submit" name="submit" value="Analyze & Explore">
	</form>

</body>
</html> 
