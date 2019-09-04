capture program drop TCPvrs
program define TCPvrs
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
	
	foreach j in CEFF MCPI EFFCH TECCH Tbetween{
	cap drop `j'
	}
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
	local N1=_N-`N1'
	count if `tid'==`T'
	local N2=r(N)
	local N2=_N-`N2'
	
	mat d11=J(_N,1,.)
	mat d12=J(`N1',1,.)
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
	   mat lamd=J(`M0',1,1)
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
		   
	       mat dmutemp1=[0 \ dmutemp1' \1]
	       mat dmutemp2=[1 \ dmutemp2' \0]
		   mat Xmat2=[objfun, Xmat,lamd]
		 * mat Xmat2=[objfun, Xmat]
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
		   *replace rel="=" if _n==2+`ninp'+`nout'
		   replace rel="=" in `=_N'
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
	qui replace CEFF=1/CEFF
	disp "Computing Dt+1[t] "
	disp "Pls wait..."
	
qui {
	local k=1
	
	forvalues i=2/`T' {
	*local i=2
	   *preserve
	   *keep if `tid'==`i'
	   *keep `gopvars' `bopvars' `invars'
	   cap mat drop Xmat dmutemp2 dmutemp1 objfun Xmat2  lamd
	   mkmat `gopvars' `bopvars' `invars' if `tid'<=`i', mat(Xmat)
	   mkmat `gopvars' `bopvars' `invars' if `tid'==`i'-1, mat(Xmat0)
	   
	   *mat list Xmat0
	    qui count if `tid'==`i'-1
	   local nobs=r(N)
	   *disp `nobs'
	   qui count if `tid'<=`i'
	   local M0=r(N)
	   mat objfun=J(`M0',1,0)
	   mat lamd=J(`M0',1,1)
	   
	   
	  forvalues j=1/`nobs' {
	  * local j=1
	
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
		   
	       mat dmutemp1=[0 \ dmutemp1'\1]
	       mat dmutemp2=[1 \ dmutemp2'\0]
		   *mat Xmm=[Xmat \ Xmat0[`j',....]]
		   
		   mat Xmat2=[objfun, Xmat,lamd]
		   mat Xmat2=Xmat2'
		   
		   preserve
		   clear
		   
		svmat Xmat2
		   *local vnames : colfullnames Xmat2
	    svmat dmutemp2, names(theta)
	       gen rel="<="
		  qui replace rel="=" in 1
	      qui replace rel=">=" if _n<=`ngo'+1 & _n>1
	      qui replace rel="=" if _n>`ngo'+1 & _n<=`nout'+1
		  qui replace rel="="  in `=_N'
		   
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
	   mat lamd=J(`M0',1,1)
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
	       mat dmutemp1=[0 \ dmutemp1' \1]
	       mat dmutemp2=[1 \ dmutemp2' \0]
		   *mat Xmm=[Xmat \ Xmat0[`j',....]]
		   mat Xmat2=[objfun, Xmat, lamd]
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
		   replace rel="="  in `=_N'
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
	qui count if `tid'==1
	local N1=r(N)
	*qui count if `tid'==T
	*local N2=r(N)
	*local N2
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
	
