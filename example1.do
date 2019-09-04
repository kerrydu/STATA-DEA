#delim cr
version 12.0
/*
 Examples from the article:
 Sensible parameters for univariate and multivariate splines
 SSC packages used:
 bspline, parmest, eclplot
 Output graphics files created:
 figseq1.eps, figseq2.eps
*/

clear all
set scheme sj

/*
 Input data
*/
sysuse auto, clear
describe

/*
 Example 1
*/
sjlog using bspline1, replace
sysuse auto
flexcurv, xvar(weight) power(3) refpts(1500(900)5100) generate(cs_)
describe cs_*
sjlog close, replace

sjlog using bspline2, replace
regress mpg cs_*, noconstant noheader
sjlog close, replace

sjlog using bspline3, replace
flexcurv, xvar(weight) power(3) refpts(1500(900)5100) base(3300) generate(bcs_)
describe bcs_*
regress mpg bcs_*, noheader
sjlog close, replace

/*
 Example 2
*/
sjlog using bspline4, replace
flexcurv, xvar(weight) power(2) refpts(2000 3000 4000) generate(qs_)
describe qs_*
regress mpg qs_*, noconstant noheader
sjlog close, replace


/*
 Create Figure 3
*/
preserve
tempfile pf1
parmest, label saving(`"`pf1'"', replace)
predict mpghat
append using `"`pf1'"', gene(param)
gene newweight=subinstr(label,"Spline at ","",1)
replace newweight=subinstr(newweight,",","",1)
replace weight=real(newweight) if param
replace mpghat=estimate if param
drop newweight
eclplot estimate min* max* weight, ///
  estopts(msize(3)) ciopts(msize(3) lwidth(0.5)) ///
  addplot(scatter mpg weight || line mpghat weight, sort lpattern(solid)) ///
  xlab(1500(500)5000) ylab(0(5)40) ///
  xtitle("Weight (lbs.)") ytitle("Mileage (mpg)")
graph export figseq1.eps, replace
more
restore

flexcurv, xvar(weight) power(2) refpts(2000 3000 4000) ///
  base(2000) generate(bqs_)
describe bqs_*
regress mpg bqs_*, noheader
ereturn list
regress mpg c.weight c.weight#c.weight, noheader
ereturn list

sjlog using bspline6, replace
regress mpg foreign bqs_*, noheader
sjlog close, replace

/*
 Example 3
*/
sjlog using bspline7, replace
flexcurv, xvar(weight) power(1) krule(interpolate) ///
  refpts(1500 2000 2500 3000 4000 5000) generate(ls_)
describe ls_*
regress mpg ls_*, noconstant noheader
sjlog close, replace

/*
 Create Figure 4
*/
preserve
tempfile pf2
parmest, label saving(`"`pf2'"', replace)
predict mpghat
append using `"`pf2'"', gene(param)
gene newweight=subinstr(label,"Spline at ","",1)
replace newweight=subinstr(newweight,",","",1)
replace weight=real(newweight) if param
replace mpghat=estimate if param
drop newweight
eclplot estimate min* max* weight, ///
  estopts(msize(3)) ciopts(msize(3) lwidth(0.5)) ///
  addplot(scatter mpg weight || line mpghat weight, sort lpattern(solid)) ///
  xlab(1500(500)5000) ylab(0(5)40) ///
  xtitle("Weight (lbs.)") ytitle("Mileage (mpg)")
graph export figseq2.eps, replace
more
restore

sjlog using bspline8, replace
flexcurv, xvar(weight) power(1) krule(interpolate) ///
  refpts(1500 2000 2500 3000 4000 5000) base(3000) generate(bls_)
describe bls_*
regress mpg bls_*, noheader
sjlog close, replace

/*
 Example 4
*/
generate odd=mod(_n,2)
lab def odd 0 "Even" 1 "Odd"
lab val odd odd
lab var odd "Oddness"
describe odd
tab odd, m

sjlog using bspline9, replace
flexcurv, xvar(weight) power(3) refpts(1760(616)4840) base(1760) ///
  generate(a_) labprefix(weight==) labfmt(%9.0g)
describe a_*
sjlog close, replace

sjlog using bspline10, replace
fvprevar ibn.odd, generate(b_)
describe b_*
sjlog close, replace

sjlog using bspline11, replace
prodvars a_*, rvarlist(b_*) generate(c_) lseparator(" & ")
describe c_*
sjlog close, replace

sjlog using bspline12, replace
regress mpg b_* c_*, noconstant noheader
sjlog close, replace

parmest, label omit ///
  list(parm label omit estimate min* max* p, noobs sepa(2))  ///
  format(label %-80s estimate min* max* %8.2f p %-8.2g)

regress mpg ibn.odd c_*, noconstant noheader
parmest, label omit ///
  list(parm label omit estimate min* max* p, noobs sepa(2))  ///
  format(label %-80s estimate min* max* %8.2f p %-8.2g)

exit
