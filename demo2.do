#delim ;
version 12.0;
*
 Creation of Figure 2.
 SSC packages used:
 bspline, keyby, expgen.
 Output graphics files created:
 figseq4.eps
*;

clear all;
set scheme sj;

sysuse auto, clear;
keyby foreign make;
desc;

*
 Spline regression models
*;
local xmin=1500;
local xmax=5100;
local xrange=`xmax'-`xmin';
local minint=`xrange'/12;
global tflist;
foreach power of num 0(1)3 {;
  if `power'==0 {;
    flexcurv, xvar(weight) power(0) refp(1500(900)4200) inc(5101) krule(interpolate) gene(sp0_);
  };
  else {;
    flexcurv, xvar(weight) power(`power') refp(1500(900)5100) gene(sp`power'_);
  };
  desc sp`power'_*;
  regress mpg sp`power'_*, noconst nohead;
  * Compute predicted values *;
  preserve;
  clear;
  set obs `=`xmax'-`xmin'+1';
  gene long power=`power';
  gene long predseq=_n;
  lab var predseq "Prediction sequence number";
  gene weight=`xmin'+_n-1;
  compress;
  if `power'==0 {;
    flexcurv, xvar(weight) power(0) refp(1500(900)4200) inc(5101) krule(interpolate) gene(sp0_);
  };
  else {;
    flexcurv, xvar(weight) power(`power') refp(1500(900)5100) gene(sp`power'_);
  };
  predict mpghat;
  summ predseq weight mpghat, de;
  tempfile tfcur;
  keyby power weight;
  desc;
  save `"`tfcur'"', replace;
  global tflist `"$tflist `"`tfcur'"'"';
  restore;  
};

*
 Produce plots
*;
preserve;
expgen =4, copyseq(power);
replace power=power-1;
lab def power 0 "Constant" 1 "Linear" 2 "Quadratic" 3 "Cubic";
lab val power power;
lab var power "degree of spline";
keyby power foreign make;
append using $tflist;
keyby power foreign make predseq, miss;
desc;
scatter mpg weight || scatter mpghat weight, msym(square) mcolor(black) msize(0.2) || ,
  by(power, compact row(2) legend(off))
  xlab(`xmin'(`minint')`xmax', angle(270) grid)
  ylab(10(5)45, angle(0))
  ytitle("Mileage (mpg)")
  xsize(4) ysize(4);
graph export figseq4.eps, replace;
more;
restore;

exit;
