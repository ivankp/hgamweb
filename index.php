<!DOCTYPE html>
<html lang="en-US">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>H &#8594; &gamma;&gamma; binning</title>
<link rel="icon" type="image/svg+xml" href="data:image/svg+xml;base64,
PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxMDAg
NTAiPjxwYXRoIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzgwQjVFOSIgc3Ryb2tlLXdpZHRoPSI1IiBk
PSJNNSA0NWgxMFYzNWgxMHYxMGgxMFYzNUgyNVYyNWgxMFYxNWgxMFY1aDEwdjEwaDEwdjEwaDEw
djEwSDY1djEwaDEwVjM1aDEwdjEwaDEwIi8+PC9zdmc+Cg==
">
<link rel="stylesheet" href="styles.css" type="text/css">
<script>
const vars = <?php
  $a1 = array();
  $d1 = scandir('.');
  foreach($d1 as &$f1) {
    if (preg_match('/^h\d+$/',$f1) && is_dir($f1)) {
      $a2 = array();
      $d2 = scandir($f1);
      foreach($d2 as &$f2) {
        if (preg_match('/^[^-]+/',$f2,$m) && is_file($f1.'/'.$f2)) {
          $a2[] = $m[0];
        }
      }
      $a1[$f1] = array_unique($a2,SORT_REGULAR);
    }
  }
  echo json_encode($a1);
?>;
<?php
  $f = file_get_contents('default_binning.json');
  if ($f) echo 'const default_binning = '.$f.";\n";
?>
</script>
<script src="main.js"></script>
</head>
<body>

<div id="first">
<p>H &#8594; &gamma;&gamma; binning tool</p>
<form id="q">
  <div>
    <label>Luminosity:<input name="L" type="text" size="6">ifb</label>
    <span id="data_lumi"></span>
  </div><div>
    <select name="H"></select>
    <select name="V"></select>
    <input list="edges_list" type="text" name="edges" size="30" autocomplete="off">
    <datalist id="edges_list"></datalist>
    <input type="submit" value="Rebin">
    <!-- <img id="loading" src="img/loading.gif" alt="loading"> -->
    <span id="time"></span>
  </div><div>
    <label>Background model:
      <select name="Bf">
        <option value="poly">Poly</option>
        <option value="exppoly">ExpPoly</option>
      </select>
    </label>
    <input type="number" name="Bdeg" min="0" max="6" value="2">
  </div><div>
    m&gamma;&gamma; bins/GeV:
    <label>Bkg:
      <input type="number" name="Bdiv" min="1" max="10" value="1">
    </label>
    <label>Sig:
      <input type="number" name="Sdiv" min="1" max="10" value="1">
    </label>
  </div><div>
    <label>Fine hist nbins:
      <input type="number" name="nV" min="1" max="200" value="100">
    </label>
  </div><div>
    <label>Significance method:
      <select name="method">
        <option value="bkg_rew">Bkg reweighting</option>
        <option value="sig_reg">Sig region</option>
      </select>
    </label>
  </div><div>
    <label><input name="click" type="checkbox">
      click row to show backround fit</label>
    <label><input name="unc" type="checkbox">
      show uncertainties</label>
  </div>
</form>
</div>

<div id="second"></div>

<table id="t"></table>

</body>
</html>
