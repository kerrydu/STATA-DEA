// Start of the MATA Definition Area -------------------------------------------
version 10
mata:
mata clear
mata set matastrict on

/** HISTORY:
 * -----------------------------------------------------------------------------
 * 2012-09-22(SAT): Deprecate and Remove minsubscript Option
 * ---------------------------------------------------------------------------*/
 
/**
 * Declare the variable's boundary condition structure.
 */
struct BoundCond {
    real scalar val, lower, upper, free
}

/**
 * Declare the LP(RSM: Revised Simplex Method)'s parameter structure.
 * Substitute for the LpParamStruct in the feature.
 */
struct LpParam {
	real scalar minYn       // whether minimization or maximization
	
    real scalar vars        // number of variables
    real scalar slacks      // number of slacks
	real scalar artificials // number of artificials
	
	real scalar tol1        // tolerance 1
    real scalar tol2        // tolerance 2
    string scalar trace     // whether trace or not.
    string scalar tracename // trancename
}

/**
 * Declare the LP(for DEA)'s parameter structure.
 */
struct LpParamStruct {
    string scalar rts       // return to scale(CRS|VRS|IRS|DRS)
	real scalar isin        // if 1 then 'in', other then 'out'
	real scalar stagestep   // stage step. 1 or 2
	real scalar minYn       // whether minimization or maximization
	
    real scalar dmus        // number of dmus
	real scalar dmuins      // number of inputs per dmu
	real scalar dmuouts     // number of outputs per dmu
    real scalar slacks      // number of slacks
	real scalar artificials // number of artificials
	
	real scalar tol1        // tolerance 1
    real scalar tol2        // tolerance 2
    real scalar isminsubscript // whether min subscript or not.(Deprecated)
    string scalar trace     // whether trace or not.
}

/**
 * Declare the LP's result structure.
 */
struct LpResultStruct {
	real scalar xVal	// objective funtion value.
	real matrix XB		// basic feasible solution.
	real scalar rc		// return code(zero means success)
	string scalar rmsg	// return message
}

/**
 * Declare the LP's tableau structure.
 */
struct LpTableauStruct {
	pointer(real matrix) scalar CB, CNj
	pointer(real matrix) scalar  B,  Nj, b
	pointer(real matrix) scalar Bi,CBBi, rawXB, XB
}

/**
 * make frame matrix and set matrix value at the param frameMat
 * rts - return to scale, ort - orientation
 */
function _mkframemat( string scalar frameMat,
                      string scalar dmuIn,
                      string scalar dmuOut,
                      string scalar rts,
                      string scalar ort )
{
    real matrix F, DI, DO

    DI = st_matrix(dmuIn)
    DO = st_matrix(dmuOut)
	
	F = mkframemat(DI, DO, rts, ort)

    // return result
    st_matrix(frameMat, F)
}

/**
 * make simplex method frame matrix
 * rts - return to scale, ort - orientation
 */
real matrix function mkframemat( 
		real matrix DI,
		real matrix DO,
		string scalar rts,
		string scalar ort )
{
	real matrix F
    real scalar row, col, sig
    real scalar dmus, slackins, slackouts, slacks
    real scalar frows, fcols

    if (cols(DI) != cols(DO)) {
        _error(3200, "in and out count of dmu is not match!")
    }

    // basic value setting for artificial variabels
    sig = ((ort == "IN") ? -1 : 1)

    dmus = cols(DI) // or cols(DO), because cols(DI) == cols(DO)
    slackins = rows(DI); slackouts = rows(DO)
    if (rts == "CRS") {
		slacks = slackins + slackouts
	    // target coefficient\slackins\slackouts
        frows = 1 + slackins + slackouts
		// target coefficient,theta,dmus,slackins,slackouts,rhs
        fcols = 1 + 1 + dmus + slackins + slackouts + 1
    }
    else if (rts == "VRS") {
		slacks = slackins + slackouts
		// target coefficient\slackins\slackouts\sum of lamda
        frows = 1 + slackins + slackouts + 1
		// target coefficient,theta,dmus,slackins,slackouts,rhs
        fcols = 1 + 1 + dmus + slackins + slackouts + 1
    }
	else if (rts == "IRS") {
        slacks = slackins + slackouts + 1
        // target coefficient\slackins\slackouts\sum of lamda
        frows = 1 + slackins + slackouts + 1
		// target coefficient,theta,dmus,slackins,slackouts,sum of lamda,rhs
        fcols = 1 + 1 + dmus + slackins + slackouts + 1 + 1
    }
    else if (rts == "DRS") {
        slacks = slackins + slackouts + 1
        // target coefficient\slackins\slackouts\sum of lamda
        frows = 1 + slackins + slackouts + 1
		// target coefficient,theta,dmus,slackins,slackouts,sum of lamda,rhs
        fcols = 1 + 1 + dmus + slackins + slackouts + 1 + 1
    }
    else {
        _error(3498, "invalid rts optoin.")
    }

    // make frame matrix for CRS(CCR)
    F = J(frows, fcols, 0)
    F[1, 1] = 1
    replacesubmat(F, 2, 3, sig * DI)
    replacesubmat(F, 2 + slackins, 3, -sig * DO)
    replacesubmat(F, 2, 3 + dmus, sig * I(slacks))

    // adjustment
	if (rts == "VRS") {
        replacesubmat(F, frows, 3, J(1, dmus, 1))
        F[frows,fcols] = 1
    }
	else if (rts == "IRS") {
        replacesubmat(F, frows, 3, J(1, dmus, 1))
        F[frows,2 + dmus + slacks] = -1
        F[frows,fcols] = 1
    }
    else if (rts == "DRS") {
        replacesubmat(F, frows, 3, J(1, dmus, 1))
        F[frows,2 + dmus + slacks] = 1
        F[frows,fcols] = 1
    }

    // return result
	return(F)
}

/**
 * DEA Loop - Data Envelopment Analysis Loop for DMUs
 */
function _dealp ( string scalar frameMat,
                  string scalar dmuIn,
                  string scalar dmuOut,
                  string scalar rts,
                  string scalar ort,
                  real scalar stagestep,
                  real scalar tol1,
                  real scalar tol2,
                  string scalar minsubscript, // Deprecated
                  string scalar efficientVec,
                  string scalar trace, 
				  | real scalar dmui )
{
    real matrix F, DI, DO, DEALPRSLT
    real colvector effvec

    F  = st_matrix(frameMat)
    DI = st_matrix(dmuIn)
    DO = st_matrix(dmuOut)
	if (stagestep == 2) {
	    effvec = st_matrix(efficientVec)
	}
	
	DEALPRSLT = dealp(F, DI, DO, rts, ort, stagestep, tol1, tol2,
			minsubscript, effvec, trace, dmui)
	
    st_matrix("r(dealprslt)", DEALPRSLT)
}

/**
 * DEA Loop - Data Envelopment Analysis Loop for DMUs
 */
real matrix function dealp ( 
		real matrix F,
		real matrix DI,
		real matrix DO,
		string scalar rts,
		string scalar ort,
		real scalar stagestep,
		real scalar tol1,
		real scalar tol2,
		string scalar minsubscript, // Deprecated
		real colvector effvec,
		string scalar trace, 
		real scalar _dmui  )
{
    real matrix M, VARS, LPRSLT, DEALPRSLT, ARTIF
    real scalar dmus, slackins, slackouts, slacks, artificials, artificialrow
    real scalar frows, fcols, isin, i, dmui, mindmui, maxdmui
    real colvector l_effvec, skipdmu
    string scalar tracename
	
	struct BoundCond matrix boundF, boundM
	struct LpParamStruct scalar param

    if (cols(DI) != cols(DO)) {
        _error(3200, "in and out count of dmu is not match!")
    }
	if (!(rts == "CRS" || rts == "VRS" || rts == "IRS" || rts == "DRS")) {
		_error(3498, "rts must be one of CRS, VRS, IRS, DRS")
	}

    // basic value setting for artificial variabels
    isin = (ort == "IN")
    frows = rows(F); fcols = cols(F)
    dmus = cols(DI) // or cols(DO), because cols(DI) == cols(DO)
    slackins = rows(DI); slackouts = rows(DO)

    tracename = rts + "-" + ort + "-" + (stagestep == 1 ? "SI" : "SII")
	
	// -------------------------------------------------------------------------
	// define number of slacks by rts
	if (rts == "CRS" || rts == "VRS") slacks = slackins + slackouts
	else if (rts == "IRS" || rts == "DRS") slacks = slackins + slackouts + 1
	
	// define number of artificials by rts, ort, stage
	if (rts == "CRS" || rts == "DRS") {
		if (stagestep == 1) {
			if (isin) {
				artificials = slackins+slackouts; artificialrow = 2;
			}
			else artificials = 0
		}
		else {
			artificials = slackouts; artificialrow = 2+slackins;
		}
	}
	else if (rts == "VRS" || rts == "IRS") {
		if (stagestep == 1) {
			if (isin) {
				artificials = slackins+slackouts+1; artificialrow = 2;
			}
			else {
				artificials = 1; artificialrow = frows //== 2+slackins+slackouts
			}
		}
		else {
			artificials = slackouts+1; artificialrow = 2+slackins
		}
	}
	if (artificials > 0) {
		ARTIF = J(1, artificials, 1) \ J(frows-1, artificials, 0)
		replacesubmat(ARTIF, artificialrow, 1, I(artificials))
		F = F[,1..fcols-1], ARTIF, F[,fcols]
		frows = rows(F); fcols = cols(F) // revise frows, fcols
	}
	// -------------------------------------------------------------------------
	
	// constants value to right-hand side(rhs) and both sides multiplied by -1.
	if (stagestep == 2) {
		l_effvec = effvec
		skipdmu = (effvec :== .)
		if (isin) {
		    replacesubmat(F, 2, 3, -F[2..1+slackins,3::2+dmus+slackins])
		}
		else {
		    replacesubmat(F, 2+slackins, 3,
                -F[2+slackins..1+slackins+slackouts,3::2+dmus+slacks])
		}
	}
	else skipdmu = J(1, dmus, 0)
	// -------------------------------------------------------------------------
    boundF = J(1, fcols, BoundCond());
	// set the boundary for the efficiency variable(theta, eta):
	// -INFINITE <= efficiency <= INFINITE
	boundF[1,2].val = 0; boundF[1,2].lower = 0; boundF[1,2].upper = .
	
	// set boundary for the weight variable(lamda, mu):
	// 0 <= weight <= INFINITE
	for (i=3; i<dmus+3; i++) {
		boundF[1,i].val = 0; boundF[1,i].lower = 0; boundF[1,i].upper = .
	}
		
	// set boundary for the non-structural variable(slack, artificial).
	// 0 <= slacks and atrificials <= INFINITE
	for (i=dmus+3; i<fcols; i++) { 
		boundF[1,i].val = 0; boundF[1,i].lower = 0; boundF[1,i].upper = .
	}
	// liststruct(boundF); // for debug
	
	// set the lp's parameters
	param.rts            = rts
	param.isin			 = isin
	param.stagestep      = stagestep
	param.dmus           = dmus
	param.slacks         = slacks
	param.artificials    = artificials
	param.tol1           = tol1
	param.tol2           = tol2
	param.trace          = trace
	// liststruct(param); // for debug
	// -------------------------------------------------------------------------
    DEALPRSLT = J(0, 1+ dmus + slacks, 0)
	
	// Added by Brian(2012.06.30)
	if (_dmui <= 0 || _dmui >= .) {
		mindmui = 1; maxdmui = dmus;
	}
	else {
		mindmui = _dmui; maxdmui = _dmui;
		if (stagestep == 2) {
			l_effvec = J(1, dmus, effvec[1])
			skipdmu = (l_effvec :== .)
		}
	}
	
    if (isin) {
        for (dmui=mindmui; dmui<=maxdmui; dmui++) {
			if (skipdmu[dmui]) {
				LPRSLT = J(1, cols(DEALPRSLT), .)
			}
			else {
				M = F; boundM = boundF
				if (stagestep == 1) replacesubmat(M, 2, 2, DI[,dmui])
				else replacesubmat(M, 2, fcols, DI[,dmui]*l_effvec[dmui])
				replacesubmat(M, 2+slackins, fcols, DO[,dmui])

				// execute LP
				VARS   = lp_phase1(M, boundM, dmui, tracename, param)
				if (VARS[1,1] == .) {
					LPRSLT = J(1, cols(DEALPRSLT), .)
				}
				else {
					LPRSLT = lp_phase2(M, boundM, VARS, dmui, tracename, param);
				}
			}

            DEALPRSLT = DEALPRSLT \ LPRSLT
        }
    }
    else {
        for (dmui=mindmui; dmui<=maxdmui; dmui++) {
			if (skipdmu[dmui]) {
				LPRSLT = J(1, cols(DEALPRSLT), .)
			}
			else {
				M = F; boundM = boundF
				replacesubmat(M, 2, fcols, DI[,dmui])
				if (stagestep == 1) {
					if (rts == "CRS" || rts == "DRS") M[1,2] = -1
					replacesubmat(M, 2+slackins, 2, DO[,dmui])
				}
				else replacesubmat(M, 2+slackins, fcols, DO[,dmui]*l_effvec[dmui])

				// execute LP
				if (artificials == 0) { // if artificials == 0 then skip phase 1
					VARS   = (0, 2+dmus..1+dmus+slacks, 1..1+dmus, 0)
					M = M[,1], 
						M[,VARS[,2::cols(VARS)-1] :+ 1], 
						M[,cols(M)]
					boundM = boundM[,1], 
							 boundM[,VARS[,2::cols(VARS)-1] :+ 1], 
							 boundM[,cols(M)]
				}
				else {
					VARS = lp_phase1(M, boundM, dmui, tracename, param)
				}
				
				if (VARS[1,1] == .) {
					LPRSLT = J(1, cols(DEALPRSLT), .)
				}
				else {
					LPRSLT = lp_phase2(M, boundM, VARS, dmui, tracename, param);
				}
			}
            DEALPRSLT = DEALPRSLT \ LPRSLT
        }
    }

    // adjust efficiency
    if (stagestep == 2) {
        replacesubmat(DEALPRSLT, 1, 1, effvec)
    }
	return(DEALPRSLT)
}

real matrix function lp_phase1 ( real matrix M,
                                 struct BoundCond matrix boundM,
                                 real scalar dmui,
                                 string scalar aTracename,
                                 struct LpParamStruct scalar param )
{
    real matrix T, VARS
    real scalar i, j, w, mrows, mcols, phase
	real vector reorderidx, bfsidx, nonbfsidx
    string scalar tracename, msg
	struct LpResultStruct scalar lpresult

    mrows = rows(M); mcols = cols(M)
    tracename = aTracename + "-PI"

    // 1st: initialize matrix.
    if (param.trace == "trace") {
        displayas("txt")
        printf("\n\n\n----------[PHASE I]----------")
        printf("\n[DMUi=%g]%s: initialize matrix.\n",
            dmui, tracename); M
    }

    // 2nd: classify basic and nonbasic.
	VARS = (0, 1..1+param.dmus+param.slacks, -1..-param.artificials, 0)
	bfsidx = J(1, mrows-1, .); nonbfsidx = J(1, 0, .)
	for (j = 3+param.dmus; j <= mcols-1; j++) {
		/* Old Code
		T = (M[2::mrows,j] :== 1)
		if (sum(T) == 1) {
		    maxindex(T, 1, i, w); bfsidx[i] = j
		}
		else nonbfsidx = nonbfsidx, j
		*/
		
		// Modified by Brian(2012.08.25): because critical logic error.
		T = M[2::mrows,j]
		if ((sum(T :!= 0) == 1) && (sum(T) == 1)) {
		    maxindex(T, 1, i, w); bfsidx[i] = j
		}
		else nonbfsidx = nonbfsidx, j
	}
	reorderidx = (1, bfsidx[1,], 2..2+param.dmus, nonbfsidx[1,], mcols)
	VARS = VARS[,reorderidx]; 
	M = M[,reorderidx]; boundM = boundM[,reorderidx]

    if (param.trace == "trace") {
        displayas("txt")
        printf("\n[DMUi=%g]%s: classify basic and nonbasic.\n",
            dmui, tracename); M; VARS
    }

    // 3rd: solve the linear programming(LP).
	phase = 1
    lpresult = lp(M, boundM, VARS, dmui, phase, tracename, param)
	
	if(lpresult.rc) VARS[1,1] = .
    return(VARS)
}

real matrix function lp_phase2 ( real matrix M,
                                 struct BoundCond matrix boundM,
                                 real matrix VARS,
                                 real scalar dmui,
                                 string scalar aTracename,
                                 struct LpParamStruct scalar param )
{
    real matrix T, XB, orgVARS, LPRSLT
    real scalar i, j, phase, mrows, mcols, realslacks
	real vector slackidx
    string scalar tracename, msg
	struct LpResultStruct scalar lpresult

    orgVARS = VARS
    mrows = rows(M); mcols = cols(M)

    tracename = aTracename + "-PII"

	// modify target function value:
	M[1,] = J(1,mcols,0); M[1,1] = 1
	if (param.stagestep == 1) { // X = theta
		for (j=2; j<mcols; j++) {
			if (VARS[1,j] == 1) M[1,j] = 1 // because of theta index == 1
		}
	}
	else if (param.stagestep == 2) { // X = S1 + S2 + ... + Sn
	    realslacks = (param.rts == "IRS" || param.rts == "DRS") ?
		        param.slacks-1 : param.slacks;
		slackidx = (2+param.dmus..1+param.dmus+realslacks)
		for (j=2; j<mcols; j++) {
			for (i=1; i<=realslacks; i++) {
				if (VARS[1,j] == slackidx[i] && !allof(M[,j], 0)) M[1,j] = 1
			}
		}
	}

    if (param.trace == "trace") {
        displayas("txt")
        printf("\n----------[PHASE II]----------")
        printf("\n[DMUi=%g]%s: initialize matrix.\n",
            dmui, tracename); M
        printf("\n[DMUi=%g]%s: VARS.\n", dmui, tracename); VARS
    }

    phase = 2
    lpresult = lp(M, boundM, VARS, dmui, phase, tracename, param)

    // -------------------------------------------------------------------------
    // phase 2 final.
    // -------------------------------------------------------------------------
	if(lpresult.rc) {
		LPRSLT = J(1, 1+param.dmus+param.slacks, .)
	}
	else {
		// lpresult = theta(1) + dmus + slacks
		LPRSLT = J(1, 1+param.dmus+param.slacks, 0) 
		for (j=1; j<=rows(lpresult.XB) ; j++) {
			if (VARS[1,j+1] > 0) LPRSLT[1, VARS[1,j+1]] =lpresult.XB[j, 1]
		}
		if (param.stagestep == 1 && LPRSLT[1, 1] <= 0) {
			LPRSLT[1, 1] = lpresult.xVal
		}
	}

    if (param.trace == "trace") {
        msg = sprintf("[DMUi=%g]%s-FINAL", dmui, tracename);
        printf("\n%s: original VARS.\n", msg); orgVARS
        printf("\n%s: VARS.\n", msg); VARS
        printf("\n%s: XB.\n", msg); lpresult.XB
        printf("\n%s: LPRSLT.\n", msg); LPRSLT
    }

    return(LPRSLT)
}

/**
 * return 0: sucess
 * return 1: B inverse error
 * return 2: XB has negative value.
 */
real scalar function decompsition(real matrix M,
								struct BoundCond matrix boundM,
								real scalar mrows, 
								real scalar mcols,
								real scalar slacks,
								struct LpTableauStruct scalar tbl,
								struct LpParamStruct scalar param )
{
	real matrix CB, CNj
	real matrix  B,  Nj, b
	real matrix Bi,CBBi, rawXB, XB, BiNjXj
	real scalar j, Njcols, result
	
	// set the tableau.
	tbl.CB = &CB; tbl.CNj  = &CNj;  
	tbl.B  = &B;  tbl.Nj   = &Nj;   tbl.b = &b
	tbl.Bi = &Bi; tbl.CBBi = &CBBi
	tbl.rawXB = &rawXB; tbl.XB = &XB;
	
	CB  = M[1,2::1+slacks]; 		CNj = M[1,2+slacks::mcols-1]
	B   = M[2..mrows,2::1+slacks];  Nj  = M[2..mrows,2+slacks::mcols-1]   
	b   = M[2..mrows,mcols]

	Bi = lusolve(B, I(rows(B)), 1e-14)
	if (any(Bi :== .)) { // B is singular matrix.
		return (result = 1);
		// FIXME use??
		// Bi = svsolve(B, I(rows(B)), 1e-14)
		// if (any(Bi :== .)) return (result = 1);
	}
	
	CBBi = CB*Bi
	// BFS(basic feasible solution)
	Njcols = cols(Nj); BiNjXj = J(rows(Nj), 1, 0)
	for (j=1; j<=Njcols; j++) {
		BiNjXj = BiNjXj :+ (Bi*Nj[.,1]*(boundM[1,1+slacks+j].val))
	}
	rawXB = Bi*b - BiNjXj
	XB = edittozerotol(rawXB, param.tol2) // BFS(basic feasible solution)
	// BFS(basic feasible solution) must be nonnegative.
	if (any(XB :< 0)) return (result = 2);

	return (result = 0); // sucess
}

/**
 * Refactoring Target: lp_for_dea
 *
 */
struct LpResultStruct function lp ( real matrix M,
									struct BoundCond matrix boundM,
									real matrix VARS,
									real scalar dmui,
									real scalar phase,
									string scalar tracename,
									struct LpParamStruct scalar param )
{
    real matrix B, CB, Bi, b, XB, rawXB, BiNjXj, CBBi, Nj, CNj, Aj, alpha
    real matrix T, TH1, TH2, valT, lowerT, upperT, LVi, V
	real scalar i, j, w, mi, boundi, evj, lvj, Njcols, leavingCase
    real scalar mrows, mcols, enteringVar, leavingVar, calci, maxiter
    real scalar existArtificial, xVal, alphaVal
	real scalar minYn, tcols, tempVal, minVal, maxVal
	struct BoundCond matrix boundT
	struct LpResultStruct scalar lpresult
	// struct LpTableauStruct scalar tbl
	
	real colvector enterings, leavings
	real scalar    enteringi, leavingi
	
	real scalar slacks, isin, tol1, tol2
	string scalar trace, msg
	
	// -------------------------------------------------------------------------
	slacks = rows(M) - 1 // number of basic feasible solution
	tol1 = param.tol1; tol2 = param.tol2
	trace = param.trace
	if (param.minYn >= .) {
		minYn = 0
		if (param.stagestep == 1) minYn = (phase == 1) ? 1 : param.isin
		else minYn = (phase == 1)
	}
	else {
		minYn = param.minYn
	}
	
    // -------------------------------------------------------------------------
	mrows = rows(M); mcols = cols(M)
    LVi = J(slacks, 1, .) // leaving variable index matrix.
    // -------------------------------------------------------------------------
if (trace == "trace") {
    displayas("txt"); msg = "initial tableau in the LP."
    printf("\n[DMUi=%g]%s: %s\n", dmui, tracename, msg); M
}
	lpresult.rc = 0; lpresult.rmsg = ""
    existArtificial = (phase == 2 && any(VARS[,2::1+slacks] :< 0));
    maxiter = st_numscalar("c(maxiter)")
    for (calci=1 ; calci<=maxiter ; calci++) { // prevent infinite loop

if (trace == "trace") {
	printf("\n[DMUi=%g]%s-LOOP[%g] Start...\n", dmui, tracename, calci)
}
        B  = M[2..mrows,2::1+slacks];       CB  = M[1,2::1+slacks]
        Nj = M[2..mrows,2+slacks::mcols-1]; CNj = M[1,2+slacks::mcols-1]
		b  = M[2..mrows,mcols]

		Bi = lusolve(B, I(rows(B)), 1e-14)
		if (any(Bi :== .)) { // B is singular matrix.
		    Bi = svsolve(B, I(rows(B)), 1e-14)
			if (any(Bi :== .)) {
				lpresult.rc = 3498; 
				lpresult.rmsg = sprintf("%s[DMUi=%g][LOOP=%g]%s",
					"No Solution(BFS's inverse is not exist):",
					dmui, calci, tracename)
				break;
				
				/* // TODO Confirm?
				display("B:");B
				display("rank(B) : det(B)"); rank(B), det(B)
				_error(3498, "No Solution(BFS's inverse is not exist):" 
						+ "[DMUi=" + strofreal(dmui) + "]" 
						+ "[LOOP=" + strofreal(calci) + "]" 
						+ tracename)
				*/
			}
		}
		CBBi = CB*Bi
		
		// BFS(basic feasible solution)
		Njcols = cols(Nj); BiNjXj = J(rows(Nj), 1, 0)
		for (j=1; j<=Njcols; j++) {
			BiNjXj = BiNjXj :+ (Bi*Nj[.,1]*(boundM[1,1+slacks+j].val))
		}
		rawXB = Bi*b - BiNjXj
		XB = edittozerotol(rawXB, tol2) // BFS(basic feasible solution)
		
		if (any(XB :== .)) {
			lpresult.rc = 3498; 
			lpresult.rmsg = sprintf("%s[DMUi=%g][LOOP=%g]%s",
				"No Solution(XB contains missing value):",
				dmui, calci, tracename)
			break;
			
			/* // TODO Confirm?
			displayas("err"); msg = "If XB contains missing value, error."
			printf("\n[DMUi=%g]%s: %s\n", dmui, tracename, msg); B;Bi;XB;
			_error(3498, "No Solution(XB have the missing value):" 
				+ "[DMUi=" + strofreal(dmui) + "]" + tracename);
			*/
		}
		
		// BFS(basic feasible solution) must be nonnegative.
        /* // TODO Confirm ??
		if (any(XB :< 0)) {
			displayas("err"); msg = "If XB contains negative value, error."
			printf("\n[DMUi=%g]%s: %s\n", dmui, tracename, msg); XB
			_error(3498, "No Solution(XB contains negative value):" 
				+ "[DMUi=" + strofreal(dmui) + "]" + tracename);
		}
        */
		
		boundT = boundM[1,2+slacks::mcols-1]
		tcols = cols(boundT); valT = J(1,tcols,.)
		for(j=1; j<=tcols; j++) {
			valT[1,j] = boundT[1,j].val
		}
		xVal = edittozerotol((CB*rawXB + CNj*valT'), tol2) // objective funtion value
        
        T = VARS[1,2::1+rows(XB)] // rows(XB) equal to number of slacks.

if (trace == "trace") {
    printf("\n[DMUi=%g]%s-LOOP[%g]: CBBi * b = %g\n",
        dmui, tracename, calci, xVal);
    display("Entered index(if value is negative, that's artificial):"); T;
    display("XB = Bi*b - BiNjXj:"); rawXB;
}

        // ---------------------------------------------------------------------
        // loop terminated condition.
        // ---------------------------------------------------------------------
        if (phase == 1) {
            // objective function value is zero
			if (xVal == 0 ) {
				// if all artificals are out or remaining artificals are at zero, 
				// stop and go phaseII
				T = (T :< 0) // artificial remain or not?
				if (allof(T, 0) || allof(select(XB, T'), 0)) break;
				
				// If remaining artificails are not zero, No Solution.
				lpresult.rc = 3498
				lpresult.rmsg = sprintf("%s[DMUi=%g][LOOP=%g]%s",
					"No Solution(Remaining artificails are not zero):",
					dmui, calci, tracename)
				break;
				
				// TODO confirm!
				// display("[BFS index | XB]");T;XB;
				// _error(3498, "No Solution(Remaining artificails are not zero):" 
				//	+ "[DMUi=" + strofreal(dmui) + "]" + tracename)
			}
		}

        // ---------------------------------------------------------------------
        // Select entering variable.
        // ---------------------------------------------------------------------
        enteringVar = 0; tempVal = 0; minVal = 0; maxVal = 0; boundi = 0
        Njcols = cols(Nj); T = J(1, Njcols, .)
        for (j=1; j<=Njcols; j++) {
            tempVal = CBBi * Nj[,j] - CNj[1,j]
			if (abs(tempVal) < tol1) continue;

			boundi = 1 + slacks + j
			if (boundM[1,boundi].val == boundM[1,boundi].lower) { // lower bound
				T[1,j] = tempVal
			}
			else { // upper bound
				T[1,j] = -tempVal
			}
        } // end of for
		
		if (!minYn) { // maximization.
			T = T :/ (T :< 0)
			if (!allof(T, .)) {
				minindex(T, 1, enterings, w); enteringi = 1
				enteringVar = enterings[enteringi]
				evj = 1+slacks+enteringVar
			}
		}
		else { // minimization.
			T = T :/ (T :> 0)
			if (!allof(T, .)) {
				maxindex(T, 1, enterings, w); enteringi = 1
				enteringVar = enterings[enteringi]
				evj = 1+slacks+enteringVar
			}
		}

        // No more candidate for entering variable.
        if (enteringVar == 0) {
			if (trace == "trace") {
                printf("\n[DMUi=%g]%s-LOOP[%g]:", dmui, tracename, calci)
                printf("No more candidate for entering variable.\n:(CB*Bi*Nj)-Cj\n");T
            }
			if (phase == 1) {
				lpresult.rc = 3498; 
				lpresult.rmsg = sprintf("%s[DMUi=%g][LOOP=%g]%s",
					"No Solution(No more candidate for entering variable):",
					dmui, calci, tracename)
				
				// TODO Confirm?
				// _error(3498, "No Solution(No more select entering variable):"
				//	+ "[DMUi=" + strofreal(dmui) + "]" + tracename)
			}
            break
        }

if (trace == "trace") {
    displayas("txt"); msg = "Select entering variable."
    printf("\n[DMUi=%g]%s-LOOP[%g]: %s(%g:%g)\n:(CB*Bi*Nj)-Cj\n",
        dmui, tracename, calci, msg, enteringVar, T[enteringVar]); T
}

        // ---------------------------------------------------------------------
        // Select leaving variable.
        // ---------------------------------------------------------------------
        leavingVar = 0
        Aj = Nj[,enteringVar]
        alpha = edittozerotol(Bi*Aj, tol1)
        if (existArtificial) {
            T = VARS[1,2::1+slacks]; tcols = cols(T)
            for (j=1; j<=tcols; j++) {
                if (T[1,j] < 0 && alpha[j,1] != 0) {
                    leavingVar = j; lvj = 1+leavingVar
                    break;
                }
            }
        }
		
		if (leavingVar == 0) {
			boundT = boundM[1,2::1+slacks]
			tcols = cols(boundT)
			lowerT = upperT = J(1,tcols,.)
			for(j=1; j<=tcols; j++) {
				lowerT[1,j] = boundT[1,j].lower
				upperT[1,j] = boundT[1,j].upper
			}
			
			// XB=Bi*b
			leavingCase = 0
			if (boundM[1,evj].val ==  boundM[1,evj].lower) {
				minVal = .;
				// 1. alpha's positive min value
				TH1 = boundM[1,evj].lower :+ ((rawXB :- lowerT') 
											:/ (alpha :* (alpha :> 0)))
				// TH1 = edittozerotol(TH1, tol1)
				if (any(TH1 :< minVal)) {
					minindex(TH1, 1, mi, w)
					leavingVar = mi[1]; lvj = 1+leavingVar
					if (phase == 1 && w[1,2] >= 2 && VARS[1,lvj] > 0) {
						// if phase 1 and same min ratio test result,
						// artificial variable must leave first.
						tcols = w[1,2]
						for (j=2; j<=tcols; j++) {
							if (VARS[1,mi[j]+1] < 0) {
								leavingVar = mi[j]; lvj = 1+leavingVar
								break;
							}
						}
					}
					minVal = TH1[leavingVar,1]; 
					leavingCase = 1
				}
				// 2. alpha's negative min value
				TH2 = boundM[1,evj].lower :+ ((rawXB :- upperT') 
											:/ (alpha :* (alpha :< 0)))
				// TH2 = edittozerotol(TH2, tol1)
				if (any(TH2 :< minVal)) {
					minindex(TH2, 1, mi, w)
					leavingVar = mi[1]; lvj = 1+leavingVar
					minVal = TH1[leavingVar,1]; leavingCase = 2
				}
				// 3. get the enteringVar's upper value
				if (boundM[1,evj].upper < minVal) {
					minVal = boundM[1,evj].upper; leavingCase = 3
				}
				
if (trace == "trace") {
    displayas("txt"); msg = "Select leaving variable.[MinVal]"
    printf("\n[DMUi=%g]%s-LOOP[%g]: %s(%g:%g)\n",
        dmui, tracename, calci, msg, leavingVar, minVal)
    display("XB | alpha:(XB=Bi*b, alpha=Bi*Aj):");rawXB,alpha
    display("[MinVal]enteringVar's upper | theta1 | theta2")
	printf("\n[boundM[1,%g].upper:%g][leavingCase:%g]\n",
			evj, boundM[1,evj].upper, leavingCase); TH1,TH2
}

				if (leavingCase == 1) {
					boundM[1,lvj].val = boundM[1,lvj].lower
				}
				else if (leavingCase == 2) {
					boundM[1,lvj].val = boundM[1,lvj].upper
				}
				else { // if (leavingCase == 3)
					boundM[1,evj].val = boundM[1,evj].upper; continue;
				}
			}
			else { // if (boundM[1,evj].val ==  boundM[1,evj].upper)
				maxVal = 0;
				// 1. alpha's positive min value
				TH1 = boundM[1,evj].upper :+ ((rawXB :- upperT') 
											:/ (alpha :* (alpha :> 0)))
				// TH1 = edittozerotol(TH1, tol1)
				if (any(TH1 :> maxVal)) {
					maxindex(TH1, 1, mi, w)
					leavingVar = mi[1]; lvj = 1+leavingVar
					maxVal = TH1[leavingVar,1]; leavingCase = 1
				}
				// 2. alpha's negative min value
				TH2 = boundM[1,evj].upper :+ ((rawXB :- lowerT') 
											:/ (alpha :* (alpha :< 0)))
				// TH2 = edittozerotol(TH2, tol1)
				if (any(TH2 :> maxVal)) {
					maxindex(TH2, 1, mi, w)
					leavingVar = mi[1]; lvj = 1+leavingVar
					maxVal = TH1[leavingVar,1]; leavingCase = 2
				}
				// 3. get the enteringVar's lower value
				if (boundM[1,evj].lower > maxVal) {
					maxVal = boundM[1,evj].lower; leavingCase = 3
				}

if (trace == "trace") {
    displayas("txt"); msg = "Select leaving variable.[MaxVal]"
    printf("\n[DMUi=%g]%s-LOOP[%g]: %s(%g:%g)\n",
        dmui, tracename, calci, msg, leavingVar, maxVal)
    display("XB | alpha:(XB=Bi*b, alpha=Bi*Aj):");XB,alpha
    display("[MaxVal]enteringVar's lower | theta1 | theta2")
	printf("\n[boundM[1,%g].lower:%g][leavingCase:%g]\n",
			evj, boundM[1,evj].lower, leavingCase); TH1,TH2
}

			    if (leavingCase == 1) {
					boundM[1,lvj].val = boundM[1,lvj].upper
				}
				else if (leavingCase == 2) {
					boundM[1,lvj].val = boundM[1,lvj].lower
				}
				else  { // if (leavingCase == 3)
					boundM[1,evj].val = boundM[1,evj].lower; continue;
				}
			}
        }

        // If no leaving variable exits
        if (leavingVar == 0) {
            if (trace == "trace")
                display("Break: No more candidate for leaving variable.")
			
			lpresult.rc = 3498; 
			lpresult.rmsg = sprintf("%s[DMUi=%g][LOOP=%g]%s",
				"No Solution(No more candidate for leaving variable):",
				dmui, calci, tracename)
					
            break
        }
        // When theta is leaving at phase 2, break! // FIXME: is correct ?
        if (phase == 2 && VARS[,lvj] == 1) {
            if (trace == "trace") display("Break: theta(еш) is not leaving.")
            break
        }

        // ---------------------------------------------------------------------
        // reply calculatation result.
        // ---------------------------------------------------------------------
		LVi[leavingVar,1] = VARS[,evj]
		_swapcols(M, lvj, evj)
		_swapcols(boundM, lvj, evj)
		_swapcols(VARS, lvj, evj)

        // Clear artificial variable
		if (VARS[,evj] < 0) {
			T = J(1, cols(VARS), 1); T[1,evj] = 0
			
			VARS   = select(VARS, T)
			M      = select(M, T)
			boundM = select(boundM, T)
			mcols  = cols(M)
			
			// if exist artificial(phase II)
			if (existArtificial) {
				existArtificial = any(VARS[,2::1+slacks] :< 0)
			}
		}

if (trace == "trace") {
    printf("\n[DMUi=%g]%s-LOOP[%g]: updated tableau.[%g(%g) <--> %g(%g)]\n",
        dmui, tracename, calci, lvj, leavingVar, evj, enteringVar); M
    display("LVi: Entered variable's VARS index value."); LVi
    display("VARS: Variable's index."); VARS
}

    } //end of main for

    // return lpresult
	if (calci > maxiter) {
		lpresult.rc = 3498; 
		lpresult.rmsg = sprintf("%s[DMUi=%g][LOOP=%g]%s",
			"No Solution(LOOP greater than maxiter):",
			dmui, calci, tracename)
	}
	if(lpresult.rc) display(lpresult.rmsg)
	lpresult.xVal = xVal
	lpresult.XB = XB
    return(lpresult)
}

/* A[.,lvj] <--> A[.,evj] */
function _swapcols( transmorphic matrix A, 
					real scalar lvj, 
					real scalar evj )
{
		transmorphic colvector  v

		v = A[., lvj]
		A[., lvj] = A[., evj]
		A[., evj] = v
}

function replacesubmat ( transmorphic matrix M,
                         real scalar row,
                         real scalar col,
                         transmorphic matrix T )
{
    M[|row,col\row + rows(T) - 1, col + cols(T) - 1|] = T
}

function _setup_dearslt_names(string scalar dearsltmat,
                              string scalar dmuinmat,
                              string scalar dmuoutmat )
{
    string matrix DMU_CS     // dmu in matrix column stripes
    string matrix DEARSLT_CS // dea result matrix column stripes
    string matrix DEARSLT_RS // dea result matrix row stripes
    real matrix M
    real scalar mcols, cnt, i

    M = st_matrix(dearsltmat)
    mcols = cols(M)

    // TODO replace the chars. ex) if ([a-z][A-Z][0-9][_]) is not. replace '_'
    DMU_CS = st_matrixcolstripe(dmuinmat)
    for (i = 1; i <= rows(DMU_CS); i++) {
        DMU_CS[i, 1] = "ref"
    }

    DEARSLT_CS = ("","rank"\"","theta")\DMU_CS\ // column join
        st_matrixrowstripe(dmuinmat)\st_matrixrowstripe(dmuoutmat)
    if (mcols - rows(DEARSLT_CS) > 0) {
        cnt = 0
        for (i = rows(DEARSLT_CS)+1 ; i <= mcols ; i++) {
            DEARSLT_CS = DEARSLT_CS \ ("slack", "slack_" + strofreal(++cnt))
        }
    }

    DEARSLT_RS = st_matrixcolstripe(dmuinmat)

    // name the row and column of dea result matrix
    st_matrixrowstripe(dearsltmat, DEARSLT_RS)
    st_matrixcolstripe(dearsltmat, DEARSLT_CS)
}

/**
 * deamat - dmucount x ( 1(theta) + dmu count + slcak(in, out) count)
 */
function _dmurank( string scalar deamat,
                   real scalar dmuincount,
                   real scalar dmuoutcount,
                   real scalar minrank,
                   real scalar tol )
{
    real matrix M
    real rowvector v, vv, retvec, slcaksum
    real scalar m, mm, row, i, ii, w, ww

    M = st_matrix(deamat)
    v = round(M[,1], tol)
    if (minrank) minindex(v, rows(v), i, w)
    else maxindex(v, rows(v), i, w)

    retvec = J(rows(v), 1, .)
    if (allof(w[,2], 1)) {
        retvec[i[1::rows(v)]] = (1::rows(v))
    }
    else {
        // rank correction for ties
        slcaksum = rowsum(M[|1,cols(M) - (dmuincount + dmuoutcount - 1)\.,.|])
        for (m = 1; m <= rows(w); m++) {
            if (w[m,2] >= 2) {
                vv = i[w[m,1]::(w[m,1] + w[m,2] - 1)]
				minindex(slcaksum[vv], w[m,2], ii, ww)
                for (mm = 1; mm <= rows(ww); mm++) {
                    for (row = ww[mm,1]; row < ww[mm,1] + ww[mm,2]; row++) {
                        retvec[vv[ii[row]]] = w[m,1] + ww[mm,1] - 1
                    }
                }
            }
            else {
                retvec[i[w[m,1]]] = w[m,1] // row = w[m,1]
            }
        }
    }
    st_matrix("r(rank)", retvec)
}

function maxvecindex( string scalar vecname )
{
    real matrix A
	real scalar i, w

    A = st_matrix(vecname)
    maxindex(A, 1, i, w)

    st_numscalar("r(maxval)", A[i[1]])
    st_numscalar("r(maxindex)", i[1])
    st_matrix("r(maxindexes)", i)
}

function minvecindex( string scalar vecname )
{
    real matrix A
	real scalar i, w

    A = st_matrix(vecname)
    if (sum(A :< .) > 0) {
        minindex(A, 1, i, w)
        st_numscalar("r(minval)", A[i[1]])
        st_numscalar("r(minindex)", i[1])
        st_matrix("r(minindexes)", i)
    }
    // if overall missing value.
    else {
        st_numscalar("r(minval)", .)
        st_numscalar("r(minindex)", 0)
        st_matrix("r(minindexes)", 0)
    }
}

function _roundmat( string scalar matname, real scalar tol )
{
    real matrix A
    A = round(st_matrix(matname), tol)
    st_matrix(matname, A)
}

function _uniqrowmat( string scalar matname, string scalar varname )
{
    st_matrix(matname, sort(uniqrows(st_data(., varname)), 1))
}

function _file_exists( string scalar fn )
{
    st_numscalar("r(fileexists)", fileexists(fn))
}

mata mlib create ldea, replace
mata mlib add ldea *()
mata mlib index

end
// End of the MATA Definition Area ---------------------------------------------
