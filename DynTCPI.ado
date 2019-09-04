*! version 1.0.1  27Sep2014
capture program drop DynTCPI
program define DynTCPI
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
	
	
	
    syntax varlist(min=1),dmu(varname) time(varname) gind(varname)
	set matsize 2000
    set more off
	
	disp("computing...")
	disp("Pls wait... ...")
	foreach j in TCPIg TCPIi TCPIc EC BPC TGC NMMCPI{
	cap drop `j'
	}
	
	
    local bopvars "`varlist'"	
    local ninp: word count `invars'
    local ngo: word count `gopvars'
    local nbo: word count `bopvars'
    local nout=`ngo'+`nbo'
	*disp("########################")
	*qui TCPIg `invars'=`gopvars':`bopvars'
	
	*rename UEIg1 UEIg
	
	tempvar gg
	qui egen `gg'=group(`gind')
	
	qui su `gg'
	local gN=r(max)
	disp("$$$$$$$$$$$$$$$$$$$$$$$$$")
	forvalues i=1/`gN' {
	preserve
	qui keep if `gg'==`i'
	
	qui TCPIg `invars' = `gopvars' : `bopvars'
	rename TCPIg ITCPI
	*list TCPIt in 1/5
	qui TCPI `invars' = `gopvars' : `bopvars', dmu(`dmu') time(`time')
	*list TCPIi in 1/5
	rename TCPI TCPIc
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
disp("xxxxxxxxxxxxxxxxxxxxxxxxx")
	
qui TCPIg `invars'=`gopvars':`bopvars'
disp("########################")
qui {
tempvar CUEIt 	
sort `dmu' `time'
bys `dmu' : gen NMMCPI=TCPIg/TCPIg[_n-1]
bys `dmu' : gen `CUEIt'=ITCPI/ITCPI[_n-1]
rename ITCPI TCPIi
bys `dmu': gen EC=TCPIc/TCPIc[_n-1]

gen BPC=`CUEIt'/EC
gen TGC=NMMCPI/`CUEIt'
}
disp("Computation is completed!")
disp("Results are reported in the dataset.")
disp("Enjoy it!")

end

