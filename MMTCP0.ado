*! version 1.0.1  
capture program drop MMTCP0
program define MMTCP0, rclass
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
	foreach j in MMCPI MCPI EFFCH TECCH CATCHUP TGR CEFF MCEFF SEC PTCU FCU{
               cap drop `j'
              }

	qui TCP2010 `invars' = `gopvars' : `bopvars', dmu(`dmu') time(`time')
	*cap drop MMCPI MCEFF 
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
	rename MCPI cMCPI 
	qui TCPvrs `invars' = `gopvars' : `bopvars', dmu(`dmu') time(`time')
	rename MCPI vMCPI
	*list CEFF in 1/5
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
	qui gen CATCHUP=MMCPI/cMCPI
	qui gen SEC=cMCPI/vMCPI
	qui gen TGR=MCEFF/CEFF
	qui sort `dmu' `time'
	qui bys `dmu': gen PTCU=TGR/TGR[_n-1]
	qui gen FCU=CATCHUP/PTCU
	cap drop  MCPI
	preserve
	qui drop if missing(MMCPI)
	qui sort `dmu' `time'
	set more off
	list `dmu' Tbetween MMCPI PTCU FCU EFFCH TECCH SEC, sep(0)
	restore
	
	
dis "Results are also plasted in the data set!"
dis "Pls check it!"
dis _newline
dis "------------------------------------------"
dis "@This code is written by Kerry@"
dis "@All rights are reserved@"
end
	
