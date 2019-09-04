*! version 1.0.1  
capture program drop MTCP
program define MTCP, rclass
    version 12.1

// syntax checking and validation-----------------------------------------------
// rts - return to scale, ort - orientation
// -----------------------------------------------------------------------------
    // returns 1 if the first nonblank character of local macro `0' is a comma,
    // or if `0' is empty.
		if replay() {
    		dis as err "ivars and ovars must be inputed."
        exit 198
		}

		// get and check invarnames
    gettoken word 0 : 0, parse(" =:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
        local invars `invars' `word'
        gettoken word 0 : 0, parse("=:,")
    }
    unab invars : `invars'

	gettoken word 0 : 0, parse(" =:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
        local gopvars `gopvars' `word'
        gettoken word 0 : 0, parse(" =:,")
    }
    unab gopvars : `gopvars'
	
    syntax varlist(min=1) [if] [in], ///
	dmu(varname) time(varname) gind(varname)
	local bopvars "`varlist'"
	
	disp("-----")
	foreach j in MMCPI MCPI EFFCH TECCH CATCHUP TGR CEFF MCEFF {
               cap drop `j'
              }

	qui TCP2010 `invars' = `gopvars' : `bopvars', dmu(`dmu') time(`time')
	cap drop MMCPI MCEFF 
	rename MCPI MMCPI
	rename CEFF MCEFF
	cap drop EFFCH TECCH 
	tempvar gg
	qui egen `gg'=group(`gind')
	disp("-----")
	qui su `gg'
	local gN=r(max)
	disp("-----")
	forvalues i=1/`gN' {
	preserve
	qui keep if `gg'==`i'
	
	qui TCP2010 `invars' = `gopvars' : `bopvars', dmu(`dmu') time(`time')
	if `i'==1 {
	   qui save tcp1,replace
	}
	else {
	   qui append using tcp1
	   qui save tcp1,replace
	}
	restore
	}
	qui {
	merge 1:1 `dmu' `time' using tcp1
	drop _merge
	erase tcp1.dta
	}
	cap drop CATCHUP TGR
	qui gen CATCHUP=MMCPI/MCPI
	qui gen TGR=MCEFF/CEFF
	cap drop  MCPI
	preserve
	qui drop if missing(MMCPI)
	qui sort `dmu' `time'
	set more off
	list `dmu' Tbetween MMCPI CATCHUP EFFCH TECCH, sep(0)
	restore
	end
	
	
	*! version 1.0.1  
capture program drop TCP2010
program define TCP2010, rclass
    version 12.1

// syntax checking and validation-----------------------------------------------
// rts - return to scale, ort - orientation
// -----------------------------------------------------------------------------
    // returns 1 if the first nonblank character of local macro `0' is a comma,
    // or if `0' is empty.
		if replay() {
    		dis as err "ivars and ovars must be inputed."
        exit 198
		}

		// get and check invarnames
    gettoken word 0 : 0, parse(" =:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
        local invars `invars' `word'
        gettoken word 0 : 0, parse("=:,")
    }
    unab invars : `invars'

	gettoken word 0 : 0, parse(" =:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
        local gopvars `gopvars' `word'
        gettoken word 0 : 0, parse(" =:,")
    }
    unab gopvars : `gopvars'
	
	
	
    syntax varlist(min=1) [if] [in], ///
	dmu(varname) time(varname)
	
	
	set matsize 2000
    set more off
	
local bopvars "`varlist'"	
local ninp: word count `invars'
local ngo: word count `gopvars'
local nbo: word count `bopvars'
local nout=`ngo'+`nbo'
qui {
order `gopvars' `bopvars' `invars'
	
	tempvar tid
	sort `time' `dmu'
	*cap drop tid
    egen `tid'=group(`time')
	*tempname T
	qui su `tid'
	local T=r(max)
	count if `tid'==1
	local N1=r(N)
	local N2=_N-`N1'
	*count if `tid'==`T'
	*local N2=r(N)
	*local N2=_N-`N2'
	
	mat d11=J(_N,1,.)
	mat d12=J(`N2',1,.)
	mat d21=J(`N2',1,.)
	}
	
	disp "Computing Dt[t] & Dt+1[t+1]"
	disp "Pls wait..."
	
qui {
	local k=1
	*local i=1
	*local j=1

	forvalues i=1/`T' {
	   
	   *preserve
	   *keep if `tid'==`i'
	   *keep `gopvars' `bopvars' `invars'
	   cap mat drop Xmat dmutemp2 dmutemp1 objfun Xmat2  lamd fobj
	   mkmat `gopvars' `bopvars' `invars' if `tid'<=`i', mat(Xmat)
	   mkmat `gopvars' `bopvars' `invars' if `tid'==`i', mat(Xmat0)
	   
	   qui count if `tid'==`i'
	   local nobs=r(N)
	   *disp `nobs'
	  
	   qui count if `tid'<=`i'
	   local M0=r(N)
	   mat objfun=J(`M0',1,0)
	   *mat lamd=J(`nobs',1,1)
	   forvalues j=1/`nobs' {
	   cap mat drop dmutemp2 dmutemp1
	  * mat list Xmat
	       mat dmutemp2=Xmat0[`j',....]
		   *mat list dmutemp2
		   
	       mat dmutemp2[1,`nout'+1]=J(1,`ninp',0)
		   mat dmutemp2[1,1]=J(1,`ngo',0)
	       mat dmutemp1=Xmat0[`j',....]
		   *display "kerry"
		   mat dmutemp1[1,`ngo'+1]=J(1,`nbo',0)
	       *mat dmutemp1[1,1]=J(1,`nout',0)
		   
		   mat dmutemp2=-dmutemp2
		   *mat dmutemp2[1,`ngo'+1]=-dmutemp2[1,`ngo'+1..`nout']
		   
	       mat dmutemp1=[0 \ dmutemp1' ]
	       mat dmutemp2=[1 \ dmutemp2' ]
		  * mat Xmat2=[objfun, Xmat,lamd]
		  mat Xmat2=[objfun, Xmat]
		   mat Xmat2=Xmat2'
		   
		  preserve
		   clear
		   
		   svmat Xmat2
		   *local vnames : colfullnames Xmat2
	       svmat dmutemp2, names(theta)
	       gen rel="<="
		   replace rel="=" in 1
	       replace rel=">=" if _n<=`ngo'+1 & _n>1
	       replace rel="=" if _n>`ngo'+1 & _n<=`nout'+1
	       svmat dmutemp1, names(rhp)
		   list rel
		  lp Xmat* theta1, min rhs(rhp1)
		  mat fobj=r(lprslt)
		  mat d11[`k',1]=1/fobj[1,1]
		  local ++k
		  *list rel
		   *mat dir 
		   restore
		   }
	}
}
    cap drop CEFF1 CEFF
	svmat d11, names(CEFF)
	rename CEFF1 CEFF
	replace CEFF=1/CEFF
	disp "Computing Dt+1[t] "
	disp "Pls wait..."
	
qui {
	local k=1
	
	forvalues i=2/`T' {
	   *preserve
	   *keep if `tid'==`i'
	   *keep `gopvars' `bopvars' `invars'
	   cap mat drop Xmat dmutemp2 dmutemp1 objfun Xmat2  lamd
	   mkmat `gopvars' `bopvars' `invars' if `tid'<=`i', mat(Xmat)
	   mkmat `gopvars' `bopvars' `invars' if `tid'==`i'-1, mat(Xmat0)
	    qui count if `tid'==`i'-1
	   local nobs=r(N)
	   *disp `nobs'
	   qui count if `tid'<=`i'
	   local M0=r(N)
	   mat objfun=J(`M0',1,0)
	  
	   *mat lamd=J(`nobs'+1,1,1)
	   forvalues j=1/`nobs' {
	   cap mat drop dmutemp2 dmutemp1
	  * mat list Xmat
	       mat dmutemp2=Xmat0[`j',....]
		   *mat list dmutemp2
	       mat dmutemp2[1,`nout'+1]=J(1,`ninp',0)
		   mat dmutemp2[1,1]=J(1,`ngo',0)
	       mat dmutemp1=Xmat0[`j',....]
		   mat dmutemp1[1,`ngo'+1]=J(1,`nbo',0)
	       *mat dmutemp1[1,1]=J(1,`nout',0)
		   mat dmutemp2=-dmutemp2
		   *mat dmutemp2[1,`ngo'+1]=-dmutemp2[1,`ngo'+1..`nout']
	       mat dmutemp1=[0 \ dmutemp1']
	       mat dmutemp2=[1 \ dmutemp2']
		   *mat Xmm=[Xmat \ Xmat0[`j',....]]
		   mat Xmat2=[objfun, Xmat]
		   mat Xmat2=Xmat2'
		   
		   preserve
		   clear
		   
		   svmat Xmat2
		   *local vnames : colfullnames Xmat2
	       svmat dmutemp2, names(theta)
	       gen rel="<="
		   replace rel="=" in 1
	       replace rel=">=" if _n<=`ngo'+1 & _n>1
	       replace rel="=" if _n>`ngo'+1 & _n<=`nout'+1
	       svmat dmutemp1, names(rhp)
		   
		  lp Xmat2* theta1, min rhs(rhp1)
		  mat fobj=r(lprslt)
		  mat d21[`k',1]=1/fobj[1,1]
		  local ++k
		  *list rel
		   *mat dir 
		   restore
		   }
	  }
 }
	
    disp "Computing Dt[t+1]"
	disp "Pls wait..."
	   * local k=`N1'+1
qui {
	local k=1
	
	forvalues i=1/`=`T'-1' {
	   *preserve
	   *keep if `tid'==`i'
	   *keep `gopvars' `bopvars' `invars'
	   cap mat drop Xmat dmutemp2 dmutemp1 objfun Xmat2  lamd
	   mkmat `gopvars' `bopvars' `invars' if `tid'<=`i', mat(Xmat)
	   mkmat `gopvars' `bopvars' `invars' if `tid'==`i'+1, mat(Xmat0)
	   qui count if `tid'==`i'+1
	   local nobs=r(N)
	   qui count if `tid'<=`i'
	   local M0=r(N)
	   mat objfun=J(`M0',1,0)
	   *mat lamd=J(`nobs'+1,1,1)
	   forvalues j=1/`nobs' {
	   cap mat drop dmutemp2 dmutemp1
	  * mat list Xmat
	       mat dmutemp2=Xmat0[`j',....]
		   *mat list dmutemp2
	       mat dmutemp2[1,`nout'+1]=J(1,`ninp',0)
		   mat dmutemp2[1,1]=J(1,`ngo',0)
	       mat dmutemp1=Xmat0[`j',....]
		   mat dmutemp1[1,`ngo'+1]=J(1,`nbo',0)
	       *mat dmutemp1[1,1]=J(1,`nout',0)
		   mat dmutemp2=-dmutemp2
		   *mat dmutemp2[1,`ngo'+1]=-dmutemp2[1,`ngo'+1..`nout']
	       mat dmutemp1=[0 \ dmutemp1']
	       mat dmutemp2=[1 \ dmutemp2']
		   *mat Xmm=[Xmat \ Xmat0[`j',....]]
		   mat Xmat2=[objfun, Xmat]
		   mat Xmat2=Xmat2'
		   
		   preserve
		   clear
		   
		   svmat Xmat2
		   *local vnames : colfullnames Xmat2
	       svmat dmutemp2, names(theta)
	       gen rel="<="
		   replace rel="=" in 1
	       replace rel=">=" if _n<=`ngo'+1 & _n>1
	       replace rel="=" if _n>`ngo'+1 & _n<=`nout'+1
	       svmat dmutemp1, names(rhp)
		   
		  lp Xmat2* theta1, min rhs(rhp1)
		  mat fobj=r(lprslt)
		  mat d12[`k',1]=1/fobj[1,1]
		  local ++k
		  *list rel
		   *mat dir 
		   restore
		   }
	}
	
	cap drop d111 d211 d121 d221 
	
	*list d111 in 1/10
	mat d22=d11[`N1'+1...,1]
	mat d11=d11[1..`N2',1]
	svmat d11 
	svmat d22 
	svmat d21 
	svmat d12
	cap drop MCPI EFFCH TECCH 
	*gen MCPI=sqrt(d11*d21/d12/d22)
	*replace MCPI=d21/d22 if missing(MCPI)
	
	cap drop EFFCH TECCH MCPI 
	gen EFFCH=d11/d22
	*gen CEEF=d11
	gen TECCH=sqrt(d21/d11*d22/d12)
	replace TECCH=d21/d11 if missing(TECCH)
	gen MCPI=EFFCH*TECCH
	cap drop d11 d22 d21 d12 d111 d211 d121 d221 
	*replace mleffch=. if mleffch>=5
	*replace mltech=. if  mltech>=5
	*replace mlindex=. if missing(mleffch)|missing(mltech)	
	}
display "Computation is completed!"

*cap drop d11 d22 d21 d12
sort `dmu' `time'
tempvar tru
qui bys `dmu': gen `tru'=`time'[_n+1]
cap drop Tbetween
qui gen Tbetween=string(`time')+"~"+string(`tru')
preserve 
qui drop if missing(`tru')
qui keep `dmu' Tbetween MCPI EFFCH TECCH
order `dmu' Tbetween MCPI EFFCH TECCH
disp _n(2) "Total carbon emissions performance (MCPI) and its decomposition"
list, sep(0) 
restore
qui{
bys `dmu': gen Tbetween1=Tbetween[_n-1]
bys `dmu': gen MCPI1=MCPI[_n-1]
bys `dmu': gen EFFCH1=EFFCH[_n-1]
bys `dmu': gen TECCH1=TECCH[_n-1]
drop Tbetween MCPI EFFCH TECCH
rename Tbetween1 Tbetween
rename MCPI1 MCPI
rename EFFCH1 EFFCH
rename TECCH1 TECCH
}
dis "Results are also plasted in the data set!"
dis "Pls check it!"
dis _newline
dis "------------------------------------------"
dis "@This code is written by Kerry@"
dis "@All rights are reserved@"
end
	
*! version 1.0.0  30OCT2012
capture program drop lp
program define lp, rclass
    version 11.0

// syntax checking and validation-----------------------------------------------
// rel - relational
// rhs - right hand side
// example:
//     lp x1 x2 x3, min
//     lp x1 x2 x3, min rel(rel_var) rhs(rhs_var)
// -----------------------------------------------------------------------------
	// returns 1 if the first nonblank character of local macro `0' is a comma,
    // or if `0' is empty.
	if replay() {
        dis as err "vars required."
        exit 198
    }
	
	#del ;
    syntax varlist(min=1) [if] [in] [using/]
    [,
        REL(varname)        // default is "rel", relational
        RHS(varname)        // default is "rhs"
		MIN	                // the objective is to minimize optimizaion 
		MAX                 // the objective is to maximize optimization 
		INTVARS(varlist)    // Integer(Mixed Integer Condition) Variables
		
        TOL1(real 1e-14)    // entering or leaving value tolerance
        TOL2(real 1e-8)     // B inverse tolerance: 2.22e-12
        TRACE               // Whether or not to do the log
        SAVing(string)      // result data file name
        REPLACE             // Whether or not to replace the result data file
    ];
    #del cr

	// default rel == "rel"
	if ("`rel'" == "") local rel = "rel"
	
	// default rhs == "rel"
	if ("`rhs'" == "") local rhs = "rhs"
	
	// optimization check
	local opt = "`min'`max'"
	if (!("`opt'" == "min" || "`opt'" == "max")) {
		dis as err "optimization is must min or max, and exclusively."
        exit 198
	}
	
	if ("`using'" != "") use "`using'", clear
    if (~(`c(N)' > 0 & `c(k)' > 0)) {
        dis as err "dataset required!"
        exit 198
    }
	
// end of syntax checking and validation ---------------------------------------

	set more off
    capture log close lp_log
    log using "lp.log", replace text name(lp_log)
    preserve
	
	if ("`if'" != "" | "`in'" != "") {
        qui keep `in' `if'  // filtering : keep in range [if exp]
    }

	// -------------------------------------------------------------------------
    // LP Start
    // -------------------------------------------------------------------------
	if ("`intvars'" == "") {
		lpmain `varlist', rel(`rel') rhs(`rhs') opt(`opt') ///
    			tol1(`tol1') tol2(`tol2') `trace'
	}
	else {
		milp `varlist', rel(`rel') rhs(`rhs') opt(`opt') ///
    			intvars(`intvars') tol1(`tol1') tol2(`tol2') `trace'
	}

	tempname tableau lprslt temp_t
	matrix `tableau' = r(tableau)
	matrix `lprslt' = r(lprslt)
	local nvars = r(nvars)
	local nslacks = r(nslacks)
	local nartificials = r(nartificials)
	
	// setup lprslt colnames and rownames
	matrix `temp_t' = `tableau'[1...,1..`=colsof(`lprslt')'] 
	matrix colnames `lprslt' = `: colnames `temp_t''
	matrix rownames `lprslt' = "opt_val"

	// -------------------------------------------------------------------------
    // REPORT
    // -------------------------------------------------------------------------
	di as result _n(2) "Input Values:"
	matrix list `tableau', noblank nohalf noheader f(%9.6g)
	
	di as result _n(2) "LP Results: options(`opt')"
	matrix list `lprslt', noblank nohalf noheader f(%9.6g)
	
	di as text _n(2) ""

	return matrix tableau = `tableau'
	return matrix lprslt = `lprslt'
	return local nvars = `nvars'
	return local nslacks = `nslacks'
	return local narticials = `nartificials' 

	set more on
    restore, preserve
    log close lp_log

end

********************************************************************************
* MILP - Mixed Integer Linear Programming
********************************************************************************
program define milp, rclass
	#del ;
    syntax varlist, rel(varname) rhs(varname) opt(string) intvars(varlist)
    [
         cnt(integer 0) tol1(real 1e-14) tol2(real 1e-8) trace
    ];
    #del cr 
    
    tempname tableau lprslt baseval
    
    // #L0 
    lpmain `varlist', rel(`rel') rhs(`rhs') opt(`opt') ///
    				  tol1(`tol1') tol2(`tol2') `trace'
    matrix `tableau' = r(tableau)
    matrix `lprslt' = r(lprslt)

// for debug
di as result _n(2) "MILP L`cnt' Input Values:"
list
matrix list `tableau', noblank nohalf noheader f(%9.6g)

di as result _n(2) "MILP L`cnt' Results: options(`opt')"
matrix list `lprslt', noblank nohalf noheader f(%9.6g)
di as text _n "--------------------------------------------------"
di as text _n

	// infeasible
	if (`lprslt'[1,1] >= .) {
		return add // all results of lpmain
	}
	else {
		// check that all variables is an integer
		local max_varname = ""
		local max_mantissa = 0
		foreach varname of varlist `intvars' {
			// because tableau and lprslt are same order
			local varvalue = ///  
				round(`lprslt'[1, colnumb(`tableau',"`varname'")], `tol1')
			local mantissa = `varvalue' - floor(`varvalue')
			if (`mantissa' > `max_mantissa') {
				local max_mantissa = `mantissa'
				local max_varname = "`varname'"
				local `baseval' = `varvalue'
			}
		}

		// if all variables is an integer
		if ("`max_varname'" == "") {
			return add // all results of lpmain
		}
		// some variables is not an integer
		else {
			// #L1
			preserve
			qui { 
				set obs `=c(N)+1'
				replace `max_varname' = 1 in `c(N)'
				replace `rel' = ">=" in `c(N)'
				replace `rhs' = ceil(``baseval'') in `c(N)'
				foreach varname of varlist  `varlist' {
					if ("`max_varname'" != "`varname'") {
						replace `varname' = 0 in `c(N)'
					}
				}
			}
	
			// recursive call
			milp `varlist', rel(`rel') rhs(`rhs') opt(`opt') cnt(`=`cnt'+1') ///
    			intvars(`intvars') tol1(`tol1') tol2(`tol2') `trace'

    		matrix `tableau' = r(tableau)
			matrix `lprslt' = r(lprslt)
			local nvars = r(nvars)
			local nslacks = r(nslacks)
			local nartificials = r(nartificials)
    		
    		// #L2
			restore, preserve
			qui { 
				set obs `=c(N)+1'
				replace `max_varname' = 1 in `c(N)'
				replace `rel' = "<=" in `c(N)'
				replace `rhs' = floor(``baseval'') in `c(N)'
				foreach varname of varlist  `varlist' {
					if ("`max_varname'" != "`varname'") {
						replace `varname' = 0 in `c(N)'
					}
				}
			}

			// recursive call
			milp `varlist', rel(`rel') rhs(`rhs') opt(`opt') cnt(`=`cnt'+2') ///
    			intvars(`intvars') tol1(`tol1') tol2(`tol2') `trace'
 
    		// #L1 and #L2 are infeasible or feasible
    		// if #L1 is infeasible or #L2 > #L1 then select #L2	
   			tempname L2
   			matrix `L2' = r(lprslt)
   			
   			
   			if ("`opt'" == "max") {
				if (`lprslt'[1,1] >= . | `L2'[1,1] > `lprslt'[1,1]) {
					matrix `tableau' = r(tableau)
					matrix `lprslt' = r(lprslt)
					local nvars = r(nvars)
					local nslacks = r(nslacks)
					local nartificials = r(nartificials)
				}
			}
			else { // else if ("`opt'" == "min") {
				if (`lprslt'[1,1] >= . | `L2'[1,1] < `lprslt'[1,1]) {
					matrix `tableau' = r(tableau)
					matrix `lprslt' = r(lprslt)
					local nvars = r(nvars)
					local nslacks = r(nslacks)
					local nartificials = r(nartificials)
				}
			}
			
			restore

			// return the final results			
			return matrix tableau = `tableau'
			return matrix lprslt = `lprslt'
			return local nvars = `nvars'
			return local nslacks = `nslacks'
			return local narticials = `nartificials' 
		}
	}
    
end


********************************************************************************
* LP Main - Linear Programming Main
********************************************************************************
program define lpmain, rclass
	#del ;
    syntax varlist, rel(varname) rhs(varname) opt(string)
    [
        tol1(real 1e-14) tol2(real 1e-8) trace
    ];
    #del cr 

	tempname tableau

	// make tableau
	mktableau `varlist' `rhs', opt(`opt') rel(`rel') tableau(`tableau')
	local nvars : list sizeof varlist    // number of variables
	local nslacks = r(nslacks)           // number of slacks
	local nartificials = r(nartificials) // number of artificials
	
	// run lp phase I and II
	mata: _lp_phase("`tableau'", "`opt'", ///
					`nvars', `nslacks', `nartificials', ///
					`tol1', `tol2', "`trace'")

	// return results for lp
	return local nvars = `nvars'
	return local nslacks = `nslacks'
	return local narticials = `nartificials'
	return matrix tableau = `tableau'
	return add // r(lprslt)
end


********************************************************************************
* LP Main - Linear Programming Main
********************************************************************************
program define lpmain_1, rclass
	#del ;
    syntax varlist, rel(varname) rhs(varname) opt(string) lprslt(name)
					tableau(name) vars(real) slacks(real) artificials(real)
    [
        intvars(varlist) tol1(real 1e-14) tol2(real 1e-8) trace
    ];
    #del cr 

	mata: _lp_phase("`tableau'", "`opt'", ///
					`vars', `slacks', `artificials', ///
					`tol1', `tol2', "`trace'")

	tempname c_lprslt // current lprslt
	matrix `c_lprslt' = r(lprslt)
	matrix colnames `c_lprslt' = `: colnames(`lprslt')'
	matrix rownames `c_lprslt' = `: rownames(`lprslt')'

// FIXME
// di as result _n "lprslt:"
// matrix list `lprslt', noblank nohalf noheader f(%9.6g)
// di as result _n "c_lprslt:"
// matrix list `c_lprslt', noblank nohalf noheader f(%9.6g)

	if ("`intvars'" != "" && `c_lprslt'[1,1] < .) { // if MILP then,
		local max_varname = ""
		local max_mantissa = 0
		foreach varname of varlist `intvars' {
			local varvalue = ///
				round(`c_lprslt'[1, colnumb(`c_lprslt',"`varname'")], `tol1')
			local varvalue = `varvalue' - floor(`varvalue')
			if (`varvalue' > `max_mantissa') {
				local max_mantissa = `varvalue'
				local max_varname = "`varname'"
			}
		}

		if ("`max_varname'" != "") { // variables is not at all integer
			tempname t_tableau t_obj t_vars t_slacks t_artificials t_rhs t_st
			tempname r1_lprslt r2_lprslt temp_t

			local varvalue = `c_lprslt'[1, colnumb(`c_lprslt',"`max_varname'")]
			
			preserve
			qui { 
				set obs `=c(N)+1'
				replace `max_varname' = 1 in `c(N)'
				replace `rel' = ">=" in `c(N)'
				replace `rhs' = ceil(`varvalue') in `c(N)'
				foreach varname of varlist  `varlist' {
					if ("`max_varname'" != "`varname'") {
						replace `varname' = 0 in `c(N)'
					}
				}
			}
			
	// make tableau
	mktableau `varlist' `rhs', opt(`opt') rel(`rel') tableau(`t_tableau')
	local r1_vars = `vars'
	local r1_slacks = r(nslacks)
	local r1_artificials = r(nartificials)
	
	// make lprslt and setup lprslt colnames and rownames
	matrix `r1_lprslt' = J(1, `=(1 + `vars' + `r1_slacks')', .)
	matrix `temp_t' = `t_tableau'[1...,1..`=colsof(`r1_lprslt')'] 
	matrix colnames `r1_lprslt' = `: colnames `temp_t''
	matrix rownames `r1_lprslt' = "opt_val"
	
	// call the lp main function
	lpmain `varlist', rel(`rel') rhs(`rhs') opt(`opt') ///
			lprslt(`r1_lprslt') tableau(`t_tableau') ///
			vars(`vars') slacks(`r1_slacks') artificials(`r1_artificials') ///
			intvars(`intvars') tol1(`tol1') tol2(`tol2') `trace'

	// setup result of lprslt
	matrix `r1_lprslt' = r(lprslt)
/*
	if (`r1_lprslt'[1,1] >= .) {
		break
	}
	*/		restore, preserve
	
		}
		else { // select lprslt because all variables are integer
			if (`lprslt'[1,1] >= .) {
				matrix `lprslt' = `c_lprslt'
			}
			else if ("`opt'" == "max") {
				if (`c_lprslt'[1,1] > `lprslt'[1,1]) {
					matrix `lprslt' = `c_lprslt'
				}
			}
			else { // else if ("`opt'" == "min") {
				if (`c_lprslt'[1,1] < `lprslt'[1,1]) {
					matrix `lprslt' = `c_lprslt'
				}
			}
		}
	}
	else if (`c_lprslt'[1,1] < .) {
		matrix `lprslt' = `c_lprslt'
	}
	
// FIXME
di as result _n "final lprslt:"
matrix list `lprslt', noblank nohalf noheader f(%9.6g)

	return matrix lprslt = `lprslt'
end

// Make Tableau Matrix ---------------------------------------------------------
program define mktableau, rclass
    syntax varlist(numeric) [if] [in], opt(string) rel(varname) tableau(name)

    // make matrix	
    mkmat `varlist' `if' `in', matrix(`tableau') rownames(`rel')

	// r_vec: row vector, s_mat: slacks matrix, a_mat: artificials matrix
	tempname r_vec s_mat a_mat
	
	local s_names = ""
	local a_names = ""
	local rel_values : rownames `tableau'
	forvalues i = 2/`=rowsof(`tableau')' {
		matrix `r_vec' = J(rowsof(`tableau'), 1, 0)
		local rel_value = word("`rel_values'", `i')
		
		if ("`rel_value'" == "<" || "`rel_value'" == "<=" ) {
			// slack
			matrix `r_vec'[`i', 1] = 1
			matrix `s_mat' = nullmat(`s_mat'), `r_vec'
			local s_names = "`s_names' s`=colsof(`s_mat')'"		
		}
		else if ("`rel_value'" == ">" || "`rel_value'" == ">=" ) {
			// slcak
			matrix `r_vec'[`i', 1] = -1
			matrix `s_mat' = nullmat(`s_mat'), `r_vec'
			local s_names = "`s_names' s`=colsof(`s_mat')'"
			
			// artificial
			matrix `r_vec'[1, 1] = 1 // coefficients of aritificial
			matrix `r_vec'[`i', 1] = 1
			matrix `a_mat' = nullmat(`a_mat'), `r_vec'
			local a_names = "`a_names' a`=colsof(`a_mat')'"
		}
		else if ("`rel_value'" == "=") {
			// artificial
			matrix `r_vec'[1, 1] = 1 // coefficients of aritificial
			matrix `r_vec'[`i', 1] = 1
			matrix `a_mat' = nullmat(`a_mat'), `r_vec' 
			local a_names = "`a_names' a`=colsof(`a_mat')'"
		}
		else {
			di as err "not allowed value of relational. :[`rel_value'] "
        	exit 198 // TODO error code confirm
		}
	} // end of forvalues statements
	
	// make return values
	tempname ret_tableau
	
	// #01. init objective and variables
	matrix `r_vec' = J(rowsof(`tableau'), 1, 0)
	matrix `r_vec'[1,1] = 1
	matrix colnames `r_vec' = "z" // Objective name
	 
	matrix `ret_tableau' = `r_vec', `tableau'[1...,1..(colsof(`tableau')-1)]
	
	// #02. append slacks
	if ("`s_names'" != "") {
		matrix colnames `s_mat' = `s_names'
		matrix `ret_tableau' = `ret_tableau', `s_mat'
		return local nslacks = colsof(`s_mat') // number of slacks
	} 
	else return local nslacks = 0
	
	// #03. append artificials
	if ("`a_names'" != "") {
		matrix colnames `a_mat' = `a_names'
		matrix `ret_tableau' = `ret_tableau', `a_mat'
		return local nartificials = colsof(`a_mat') // number of artificials
	} 
	else return local nartificials = 0 
	
	// #04. append rhs
	matrix `ret_tableau' = `ret_tableau', `tableau'[1...,colsof(`tableau')]
	
	// #05. return results
    matrix `tableau' = `ret_tableau'
end

// Start of the MATA Definition Area -------------------------------------------
version 10
mata:
mata set matastrict on

void function _lp_phase (
		string scalar tableau,
		string scalar opt,
		real scalar vars,
		real scalar slacks,
		real scalar artificials,
		real scalar tol1,
		real scalar tol2,
		string scalar trace )
{
	real matrix M, VARS
	real fcols

	struct BoundCond matrix boundM
	struct LpParam scalar param
	struct LpResultStruct scalar lpresult

	// 1st. load matrix and variable indexes
	M = st_matrix(tableau)
	VARS = (0, 1..vars+slacks, -1..-artificials, 0)
	
	// 2rd. make boundary matrix
	// 0 <= weight, slacks, atrificials <= INFINITE
	boundM = J(1, cols(M), BoundCond());		
	for (i=1; i<cols(M); i++) { 
		boundM[1,i].val = 0; boundM[1,i].lower = 0; boundM[1,i].upper = .
	}
	
	// 3th. set the lp's parameters
	param.minYn			 = (opt == "min"); // 0: max, 1: min

	param.vars           = vars
	param.slacks         = slacks
	param.artificials    = artificials
	param.tol1           = tol1
	param.tol2           = tol2
	param.trace          = trace
	param.tracename      = "LP for RSM"

	lpresult = lp_phase(M, boundM, VARS, param)
	
	// -------------------------------------------------------------------------
    // final.
    // -------------------------------------------------------------------------
	if(lpresult.rc) {
		LPRSLT = J(1, 1+param.vars+param.slacks, .)
	}
	else {
		// lpresult = theta(1) + vars + slacks
		LPRSLT = J(1, param.vars+param.slacks, 0)
		for (j=1; j<=rows(lpresult.XB) ; j++) {
			if (VARS[1,j+1] > 0) LPRSLT[1, VARS[1,j+1]] = lpresult.XB[j, 1]
		}
		LPRSLT = lpresult.xVal, LPRSLT
	}

    if (param.trace == "trace") {
        msg = sprintf("%s-FINAL", param.tracename);
        // printf("\n%s: original VARS.\n", msg); orgVARS
        printf("\n%s: VARS.\n", msg); VARS
        printf("\n%s: XB.\n", msg); lpresult.XB
        printf("\n%s: LPRSLT.\n", msg); LPRSLT
    }

    st_matrix("r(lprslt)", LPRSLT)
}

/**
 * @param VARS 	- Variable Index Matrix
 *              [z, B, N, b]'s index in the original Tableau
 * @param M 	- Tableau: [z, A, S, Af, b] --> [z, B, N, b]
 * @param phase - if have artificials, then phase 1 and 2, 
 *				  otherwise only phase 2
 * @param param - parameter struct for Lp RSM 
 *
 * @return result of LP
 */
struct LpResultStruct function lp_phase ( 
	real matrix M,
	struct BoundCond matrix boundM,
	real matrix VARS,
	struct LpParam scalar param )
{
	real scalar phase, mrows, mcols, j, idx
	string scalar tracename
	real vector reorderidx, bfsidx, nonbfsidx
	real vector coef_of // coefficient of objective function
	struct LpParamStruct scalar lpParam
	struct LpResultStruct scalar lpResult
	
	// validation checking.
	if (param.minYn >= .) { // 
		displayas("err");
		_error(3351, "You have to set the minimization(1) or maximization(0) "
					+ "at the LpParam.minYn")
	}
	
	coef_of = M[1, 2..1+param.vars] // keep the objective function
	replacesubmat(M, 1, 2, J(1, param.vars, 0))
	
// initialize matrix.
if (param.trace == "trace") {
	displayas("txt")
	printf("\n\n%s: initialize matrix.\n", param.tracename); M
}

	mrows = rows(M); mcols = cols(M)
    
    // classify basic and nonbasic.
	bfsidx = J(1, mrows-1, .); nonbfsidx = J(1, 0, .)
	for (j = 2+param.vars; j <= mcols-1; j++) {
		T = M[2::mrows,j]
		if ((sum(T :!= 0) == 1) && (sum(T) == 1)) {
		    maxindex(T, 1, i, w); bfsidx[i] = j
		}
		else nonbfsidx = nonbfsidx, j
	}
	reorderidx = (1, bfsidx[1,], 2..1+param.vars, nonbfsidx[1,], mcols)
	VARS = VARS[,reorderidx]; 
	M = M[,reorderidx]; boundM = boundM[,reorderidx]

    if (param.trace == "trace") {
        displayas("txt")
        printf("\n%s: classify basic and nonbasic.\n", tracename); M; VARS
    }
    
    // set the lp's parameters
	lpParam.dmus		= param.vars
	lpParam.slacks		= param.slacks
	lpParam.artificials = param.artificials
	lpParam.tol1		= param.tol1
	lpParam.tol2		= param.tol2
	lpParam.trace		= param.trace
	
	// solve the linear programming(LP): phase I
	if (param.artificials > 0) {
		phase = 1
		lpParam.minYn = 1; // min because of phase 1
		tracename = param.tracename + "-PI"
		lpResult = lp(M, boundM, VARS, 0, phase, tracename, lpParam)

		if (lpResult.rc) return(lpResult)
	}
	
	// solve the linear programming(LP): phase II
	phase = 2
	lpParam.minYn = param.minYn // according to the optimization.
	tracename = param.tracename + "-PII"

	// set the objective function.
	mcols = cols(M)
	for (j=2; j<mcols; j++) {
		idx = VARS[1,j]
		if (0 < idx && idx <= param.vars) {
			M[1,j] = coef_of[idx] // according to variable's index
		}
	}
	lpResult = lp(M, boundM, VARS, 0, phase, tracename, lpParam)

	// return result.	
	return(lpResult)
}

end
// End of the MATA Definition Area ---------------------------------------------

	
	
	

	
	
	