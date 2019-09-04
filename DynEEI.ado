*! version 1.0.1  27Sep2014
capture program drop DynEEI
program define DynEEI
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
	foreach j in EEIg EEIi EEIc EC BPC TGC MMEEI {
	cap drop `j'
	}
	
	
    local bopvars "`varlist'"	
    local ninp: word count `invars'
    local ngo: word count `gopvars'
    local nbo: word count `bopvars'
    local nout=`ngo'+`nbo'
	disp("########################")
	qui EEIg `invars'=`gopvars':`bopvars'
	
	rename EEIg1 EEIg
	
	tempvar gg
	qui egen `gg'=group(`gind')
	
	qui su `gg'
	local gN=r(max)
	disp("$$$$$$$$$$$$$$$$$$$$$$$$$")
	forvalues i=1/`gN' {
	preserve
	qui keep if `gg'==`i'
	
	qui EEIg `invars' = `gopvars' : `bopvars'
	rename EEIg1 EEIi
	qui EEI `invars' = `gopvars' : `bopvars', dmu(`dmu') time(`time')
	rename EEI1 EEIc
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
qui {
tempvar CEEIt 	
sort `dmu' `time'
bys `dmu' : gen MMEEI=EEIg/EEIg[_n-1]
bys `dmu' : gen `CEEIt'=EEIi/EEIi[_n-1]

bys `dmu': gen EC=EEIc/EEIc[_n-1]

gen BPC=`CEEIt'/EC
gen TGC=MMEEI/`CEEIt'
}
disp("Computation is completed!")
disp("Results are reported in the dataset.")
disp("Enjoy it!")

end

