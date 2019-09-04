*! version 1.0.1  27Sep2014
capture program drop TCPI
program define TCPI
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
	
	
	
    syntax varlist(min=1),dmu(varname) time(varname)
	
	cap drop TCPI
	set matsize 2000
    set more off
	
local bopvars "`varlist'"	
local ninp: word count `invars'
local ngo: word count `gopvars'
local nbo: word count `bopvars'
local nout=`ngo'+`nbo'

	tempvar tid
	**********************sort `dmu' `time'
	sort `time' `dmu' 
	*cap drop tid
    qui egen `tid'=group(`time')
	*tempname T
	qui su `tid'
	qui local T=r(max)

qui {
order `invars' `gopvars' `bopvars'
	}
mat eff=J(_N,1,.)	
mat temp1=[1/3,J(1,`ngo',1/`ngo'/3),J(1,`nbo',1/`nbo'/3)]

    disp "Computing... ..."
	disp "Pls wait..."
	
qui {
    local k=1
	forvalues i=1/`T' {
	
	count if `tid'==`i'
	local nobs=r(N)
	mkmat `gopvars' `bopvars' `invars' if `tid'==`i', mat(Xmat)

	   cap mat drop m1 m2 obj m3 fobj temp2 temp3 temp4 XZ
	   mat obj=J(`nobs',1,0)
       mat m2=[obj,Xmat]
       mat m2=m2'
	  forvalues j=1/`nobs' {
	       mat m3=Xmat[`j',....]
		   mat m3=[0 \ m3']
		   mat temp2=Xmat[`j',....]
		   mat temp2[1,`ninp'+1]=-temp2[1,`ninp'+1..`ninp'+`ngo']
		   mat temp4=temp2[1,`ninp'..`ninp'+`ngo'+`nbo']
		   mat temp4=temp4'
		   mat m1=[temp1 \ J(`ninp'-1,`ngo'+`nbo'+1,0)\ diag(temp4)]
		   *disp("kerry3")
		   mat XZ=[m1,m2]
		  preserve
		   clear
		   
		   svmat XZ
		   svmat m3, names(rhp)
		   *local vnames : colfullnames Xmat2
	     
	       gen rel="<="
		   replace rel="=" in 1
	       *replace rel=">=" if _n<=`ninp'+1 & _n>1
	       replace rel=">=" if _n>=`ninp'+2 & _n<=`ninp'+`ngo'+1
		   replace rel="=" if _n>=`ninp'+`ngo'+2 
	      * list rel 
		   
		   *mat list m1 
		   *mat list m2
		   *mat list m3
		  
		  lp XZ*, max rhs(rhp1)
		  mat fobj=r(lprslt)
		 * mat temp4=fobj[1,2..6]
		  mat eff[`k',1]=(1-fobj[1,4])/(1+fobj[1,3])
		  local ++k
		   restore
		   }
	}
 }
	*svmat d11, names(beta)
cap drop TCPI1	
svmat eff, names(TCPI)
rename TCPI1 TCPI
display "Computation is completed!"	
dis "Results are plasted in the data set!"
dis "Pls check it!"
dis _newline
dis "------------------------------------------"
dis "@This code is written by Kerry@"
dis "@All rights are reserved@"
end
