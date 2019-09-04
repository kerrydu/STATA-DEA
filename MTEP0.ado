*! version 1.0.1  
capture program drop MTEP0
program define MTEP0, rclass
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
/*
	gettoken word 0 : 0, parse(" =:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
        local gopvars `gopvars' `word'
        gettoken word 0 : 0, parse(" =:,")
    }
    unab gopvars : `gopvars'
	*/
	*disp("-----")
    syntax varlist(min=1), ///
	dmu(varname) time(varname) gind(varname)
	local outvars "`varlist'"
	foreach j in Tbetween TFEE MTFEE EFFCH TECCH MMEPI TGR CATCH {
	cap drop `j'
	}
	
	
	*disp("`outvars'")
	
	qui TEP0 `invars' = `outvars', dmu(`dmu') time(`time')
	cap drop MMEPI MTFEE
	rename MEPI MMEPI
	rename TFEE MTFEE
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
	
	qui TEP0 `invars' = `outvars', dmu(`dmu') time(`time')
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
	cap drop CATCH
	qui gen CATCH=MMEPI/MEPI
	cap drop MEPI
	qui gen TGR=MTFEE/TFEE
	preserve
	qui drop if missing(MMEPI)
	qui sort `dmu' `time'
	set more off
	list `dmu' Tbetween MMEPI CATCH EFFCH TECCH, sep(0)
	restore
	end
	/*
	
qui{
bys `dmu': gen Tbetween1=Tbetween[_n-1]
bys `dmu': gen MMEPI1=MMEPI[_n-1]
bys `dmu': gen EFFCH1=EFFCH[_n-1]
bys `dmu': gen TECCH1=TECCH[_n-1]
drop Tbetween MMEPI EFFCH TECCH
rename Tbetween1 Tbetween
rename MMEPI1 MMEPI
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
	
	*/
	
	