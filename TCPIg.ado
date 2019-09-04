*! version 1.0.1  27Sep2014
*注意energy要放在第三个投入位置
capture program drop TCPIg
program define TCPIg
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
	
	
	
    syntax varlist(min=1) [if] [in]
	cap drop TCPIg
	
	set matsize 2000
    set more off
	
local bopvars "`varlist'"	
local ninp: word count `invars'
local ngo: word count `gopvars'
local nbo: word count `bopvars'
local nout=`ngo'+`nbo'
qui {
order `invars' `gopvars' `bopvars'
	}
mat eff=J(_N,1,.)	
*cap mat drop Xmat 
mkmat `invars' `gopvars' `bopvars', mat(Xmat)
*mat lamd=J(_N,1,1)
mat obj=J(_N,1,0)
mat m2=[obj,Xmat]
mat m2=m2'
mat temp1=[1/3,J(1,`ngo',1/`ngo'/3),J(1,`nbo',1/`nbo'/3)]
*disp("kerry")
	*local i=1
	*disp(_N)
	local nob=_N
    disp "Computing... ..."
	disp "Pls wait..."
	
qui {

	forvalues i=1/`nob' {
	  *disp("kerry")

	   *preserve
	   cap mat drop m1 m3 fobj temp2 temp3 temp4 XZ
	   
	  * mat list Xmat
	       mat m3=Xmat[`i',....]
		   mat m3=[0 \ m3']
		   
		   mat temp2=Xmat[`i',....]
		   *disp("kerry1")
		   *mat temp4=temp2[1,`ninp'+1..`ninp'+`ngo']
		   *mat list temp4
		   mat temp2[1,`ninp'+1]=-temp2[1,`ninp'+1..`ninp'+`ngo']
		   mat temp4=temp2[1,`ninp'..`ninp'+`ngo'+`nbo']
		   *disp("kerry2")
		   mat temp4=temp4'
		   *disp("kerry2")
		   *mat list temp1 
		   *mat temp3=diag(temp2)
		   *mat list temp3
		   *mat list temp1
		   *mat list temp4
		   mat m1=[temp1 \ J(`ninp'-1,`ngo'+`nbo'+1,0)\ diag(temp4)]
		   *disp("kerry3")
		   mat XZ=[m1,m2]
		   
		   *disp("kerry")
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
		  mat eff[`i',1]=(1-fobj[1,4])/(1+fobj[1,3])
		  *list rel
		   *mat dir 
		   restore
		  }
 }
	*svmat d11, names(beta)
cap drop TCPIg1	
svmat eff, names(TCPIg)
rename TCPIg1 TCPIg
display "Computation is completed!"	
dis "Results are plasted in the data set!"
dis "Pls check it!"
dis _newline
dis "------------------------------------------"
dis "@This code is written by Kerry@"
dis "@All rights are reserved@"
end
	

	