*! version 1.0.1  1Mar2016

capture program drop sbm
program define sbm, rclass
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
	cap drop sbm
	
	set matsize 2000
    set more off
	
	local bopvars "`varlist'"	
local ninp: word count `invars'
local ngo: word count `gopvars'
local nbo: word count `bopvars'

qui {
order `invars' `gopvars' `bopvars'
	}
mat eff=J(_N,1,.)	

mkmat `invars' `gopvars' `bopvars', mat(Xmat)


local nob=_N
    disp "Computing... ..."
	disp "Pls wait..."
	
qui {

	forvalues i=1/`nob' {
	
	mat rr=xmat(i,....)
	mat rr=diag(rr)
	mat rr=vecdiag(inv(rr))
	mat xr=rr[1,1..`ninp']
	mat lr=rr[1,`ninp'+1...]
	
	mat f=(1,-1/`ninp'*xr)
	
