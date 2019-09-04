#delim ;
version 12.0;
*
 Creation of Figure 1.
 SSC packages used:
 bspline, keyby, expgen.
 Output graphics files created:
 figseq3.eps
*;

clear all;
set scheme sj;

sysuse auto, clear;
keyby foreign make;
desc;

*
 Linear spline regression models
*;
local xmin=1500;
local xmax=5100;
local xrange=`xmax'-`xmin';
local intmin=`xrange'/12;
global tflist;
foreach nint of num 1(1)4 {;
  local knotdiff=(`xmax'-`xmin')/`nint';
  flexcurv, xvar(weight) power(1) refp(`xmin'(`knotdiff')`xmax') gene(ls`nint'_);
  desc ls`nint'_*;
  regress mpg ls`nint'_*, noconst nohead;
  * Compute predicted values *;
  preserve;
  clear;
  set obs `=`xmax'-`xmin'+1';
  gene byte nint=`nint';
  gene long predseq=_n;
  lab var predseq "Prediction sequence number";
  gene weight=`xmin'+_n-1;
  compress;
  flexcurv, xvar(weight) power(1) refp(`xmin'(`knotdiff')`xmax') gene(ls`nint'_);
  predict mpghat;
  summ predseq weight mpghat, de;
  tempfile tfcur;
  keyby nint weight;
  desc;
  save `"`tfcur'"', replace;
  global tflist `"$tflist `"`tfcur'"'"';
  restore;
};

*
 Produce plots
*;
preserve;
expgen =4, copyseq(nint);
lab var nint "Number of intervals between knots";
keyby nint foreign make;
append using $tflist;
gene byte nknot=nint+1;
lab var nknot "number of knots";
keyby nknot foreign make predseq, miss;
desc;
scatter mpg weight || line mpghat weight, sort lpattern(solid) || ,
  by(nknot, compact row(2) legend(off))
  xlab(`xmin'(`intmin')`xmax', angle(270) grid)
  ylab(10(5)45, angle(0))
  ytitle("Mileage (mpg)")
  xsize(4) ysize(4);
graph export figseq3.eps, replace;
more;
restore;

exit;
