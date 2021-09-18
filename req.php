<?php
  ini_set('zlib.output_compression', 1);
  header('Content-type: application/json');
  // echo $_SERVER['QUERY_STRING'];
  $_GET['edges'] = explode(' ',$_GET['edges']);
  $J = json_encode($_GET,JSON_NUMERIC_CHECK);
  // echo $J;
  echo exec('bin/hgamweb '.'\''.$J.'\' 2>&1')."\n";
  // print_r($_GET);
  // echo json_encode($_GET);
?>
