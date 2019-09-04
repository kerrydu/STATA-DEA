****************************************************************************************************************************************************
** APPLICATION: auto.dta ***************************************************************************************************************************
****************************************************************************************************************************************************
sjlog using orderalpha1 ,replace
sysuse auto
gen gpm = 1/mpg
set seed 987654321
orderalpha weight length displacement gpm = rep78 headroom trunk if foreign
sjlog close, replace
sjlog using orderalpha2 ,replace
orderalpha weight length displacement gpm = rep78 headroom trunk if foreign, dmu(make) reps(200) table(full) nogenerate
sjlog close, replace
sjlog using orderalpha3 ,replace
orderalpha weight length displacement gpm = rep78 headroom trunk if foreign, dmu(make) ort(output) reps(200) table(full)
sjlog close, replace
sjlog using orderalpha4 ,replace
orderalpha weight length displacement gpm = rep78 headroom trunk if foreign, dmu(make) alpha(90) reps(200) table(full) generate(escore erank eref) replace
sjlog close, replace
sjlog using orderalpha5 ,replace
test _b[Audi_5000]-_b[Volvo_260]=0
sjlog close, replace
sjlog using orderalpha6 ,replace
orderm weight length displacement gpm = rep78 headroom trunk if foreign, dmu(make) m(16) d(1000) table(full) dot(2)
sjlog close, replace

****************************************************************************************************************************************************
** FIGURES: Frontiers ******************************************************************************************************************************
****************************************************************************************************************************************************
**cd I:\orderalpha\Paper
clear all
set more off
set seed 987654328
global ssize = 40
set obs $ssize
gen dmu = _n
gen y = 1
gen xtheo1 = 0.5+ 2*runiform()
gen xtheo2 = 1/xtheo1

gen x1 = xtheo1 +0.1*exp(rnormal())
gen x2 = xtheo2 +0.1*exp(rnormal())

** Outliers **
local mmfac1 = 0.5
local mmfac2 = 0.5
replace x1 = `mmfac1'*0.5          if _n == 1
replace x1 = `mmfac1'*0.8- 0.05    if _n == 2
replace x1 = `mmfac1'*1.25+ 0.2    if _n == 3
replace x1 = `mmfac1'*2            if _n == 4
replace x2 = `mmfac2'*0.5^-1       if _n == 1
replace x2 = `mmfac2'*0.8^-1       if _n == 2
replace x2 = `mmfac2'*1.25^-1 +0.1 if _n == 3
replace x2 = `mmfac2'*2^-1         if _n == 4

** Specifications **
global inp1 "x1"
global inp2 "x2"
global inp "${inp1} ${inp2}"
global out "y"
global alp = 95
global mm = 12 
global dd = 2000

** Efficiency Analyses **
orderalpha ${inp}=${out}, dmu(dmu) tab(f) replace
gen e${inp1}_fdh = _fdh_input*${inp1} 
gen e${inp2}_fdh = _fdh_input*${inp2} 

orderalpha ${inp}=${out}, alpha(${alp}) dmu(dmu) tab(f) replace
gen e${inp1}_oa = _oa_input_${alp}*${inp1} 
gen e${inp2}_oa = _oa_input_${alp}*${inp2} 

orderm ${inp}=${out}, m(${mm}) d(${dd}) dmu(dmu) tab(f) dots(2) replace
gen e${inp1}_om = _om_input_${mm}*${inp1} 
gen e${inp2}_om = _om_input_${mm}*${inp2} 

dea ${inp} = ${out} , rts(vrs) ort(i) 
mat dea_input = r(dearslt)
mat dea_input = dea_input[.,"theta"]
svmat dea_input, names(col)
gen e${inp1}_dea = theta*${inp1} 
gen e${inp2}_dea = theta*${inp2} 

** Comple DEA and FDH Hull **
sum ${inp1}
local min1 = r(min)
local max1 = r(max)
sum ${inp2}
local min2 = r(min)
local max2 = r(max)

tab dmu if theta < 1, matrow(inefdea)
local dea1 = inefdea[1,1]
local dea2 = inefdea[2,1]

tab dmu if _fdh_input < 1, matrow(ineffdh)
local fdh1 = ineffdh[1,1]
local fdh2 = ineffdh[2,1]

replace e${inp1}_dea = max(2 ,`max1') if dmu == `dea1'
replace e${inp2}_dea = `min2' if dmu == `dea1'
replace e${inp1}_dea = `min1' if dmu == `dea2'
replace e${inp2}_dea = max(2 ,`max2') if dmu == `dea2'
capture drop hulldea 
gen hulldea = 1 if dmu == `dea1' | dmu == `dea2' | theta == 1

replace e${inp1}_fdh = max(2 ,`max1') if dmu == `fdh1'
replace e${inp2}_fdh = `min2' if dmu == `fdh1'
replace e${inp1}_fdh = `min1' if dmu == `fdh2'
replace e${inp2}_fdh = max(2 ,`max2') if dmu == `fdh2'
capture drop hullfdh
gen hullfdh = 1 if dmu == `fdh1' | dmu == `fdh2' | _fdh_input == 1

** additional Sorting Variables
capture drop sort_dea
gen sort_dea = 1/e${inp1}_dea
capture drop sort_fdh
gen sort_fdh = 1/e${inp1}_fdh
capture drop sort_oa
gen sort_oa  = 1/e${inp1}_oa  
capture drop sort_om
gen sort_om  = 1/e${inp1}_om 

** Non-Parametric Frontiers **
twoway (scatter ${inp1} ${inp2} if _n>4, msize(medium) mcolor(black)) /*
*/ (scatter ${inp1} ${inp2} if _n<=4, msize(medium) mcolor(black) msymbol(diamond)) /*
*/ (function y = 1/x, range(0.4 2.2 /*e${inp2}_om*/)  lpattern(solid) lwidth(thin) lcolor(gs10)) /*
*/ (line e${inp1}_dea e${inp2}_dea if hulldea == 1, lpattern(solid) lwidth(thick) lcolor(gs12) sort(e${inp2}_dea sort_dea)) /*
*/ (line e${inp1}_fdh e${inp2}_fdh if hullfdh == 1, lpattern(dot) lwidth(thick) lcolor(gs8) connect(stairstep) sort(e${inp2}_fdh sort_fdh)) /*
*/ (line e${inp1}_oa  e${inp2}_oa, lpattern(/*dash*/ solid) lwidth(thick medthick) lcolor(gs8) connect(stairstep) sort(e${inp2}_oa sort_oa)) /*
*/ (line e${inp1}_om  e${inp2}_om, lpattern(/*dash_dot*/ solid) lwidth(thick medthick) lcolor(gs0) connect(stairstep) sort(e${inp2}_om sort_om)), /*
*/ graphregion(fcolor(gs15)) /*
*/ xscale(range(-0.05 2.1) /*lstyle(none)*/) yscale(range(-0.05 2.4) /*lstyle(none)*/) /*
*/ ytitle(input 2) ylabel(none) ytick(0) /*
*/ xtitle(input 1) xlabel(none) xtick(0) /*
*/ legend(label(3 "true frontier") label(4 "DEA") label(5 "FDH") label(6 "order-{&alpha} (${alp})") label(7 "order-m ($mm)")/*
*/ order(4 6 3 5 7) /*off size(small) */ col(3)) 

graph export illu-frontiers-new.png, width(1800) replace
graph export illu-frontiers-new.eps, logo(on) /*orientation(landscape)*/ replace

** Scatter Plot **
gen l1 = 0 if _n == 1
gen l2 = 0 if _n == 1
replace l1 = x1[5] if _n == 3
replace l2 = x2[5] if _n == 3
replace l1 = (x1[5]/x2[5])^0.5 if _n == 2
replace l2 = (x1[5]/x2[5])^-0.5  if _n == 2
capture drop mark
gen str2 mark = "O"  if _n == 1
replace  mark = "A*" if _n == 2
replace  mark = "A"  if _n == 3


twoway (scatter ${inp1} ${inp2} if _n>4, msize(medium) mcolor(black)) /*
*/ (scatter ${inp1} ${inp2} if _n<=4, msize(medium) mcolor(black) msymbol(diamond)  /*mlabel(mark) mlabposition(8) mlabcolor(black)*/) /*
*/ (function y = 1/x, range(0.4 2.2)  lpattern(solid) lwidth(medthick) lcolor(gs10)) /*
*/ (line l1 l2, lpattern(dash) lwidth(thin) lcolor(gs0) mlabel(mark) mlabposition(8) mlabcolor(black) sort)/*
*/ (scatter  l1 l2 if _n >1, /*msymbol(none)*/ msymbol(circle_hollow) msize(small) mcolor(black) mlabel(mark) mlabsize(small) mlabposition(9) mlabcolor(black))/*
*/ (scatter  l1 l2 if _n==1, /*msymbol(none)*/ msymbol(circle_hollow) msize(small) mcolor(black) mlabel(mark) mlabsize(small) mlabposition(8) mlabcolor(black)), /*
*/ graphregion(fcolor(gs15)) /*
*/ xscale(range(-0.05 2.1) /*lstyle(none)*/) yscale(range(-0.05 2.4) /*lstyle(none)*/) /*
*/ ytitle(input 2) ylabel(none) ytick(0) /*
*/ xtitle(input 1) xlabel(none) xtick(0) /*
*/ legend(label(1 "regular DMUs") label(2 "irregular DMUs (outliers)") label(3 "true frontier")/*
*/ order(1 2 3) /*off size(small) */ col(2)) 

graph export illu-scatter-new.png, width(1800) replace
graph export illu-scatter-new.eps, logo(on) /*orientation(landscape)*/ replace

