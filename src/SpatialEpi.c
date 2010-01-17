#include <R.h> 
#include <Rmath.h> 
#include <Rinternals.h>


/*
 ------------------------------------------------------------------------
 
 Utility Functions
 
 ------------------------------------------------------------------------
 */


/*
	getListElement:  Get list element named "str"
*/
static SEXP getListElement(SEXP list, char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < length(list); i++)
	if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	}
    return elmt;
}


/*
	vecSum:  Computes the sum of the elements of a vector
*/
SEXP vecSum(SEXP Rvec){
	// SEXP coercion of R objects --------------------------	
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);

	int i, n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP value;
	PROTECT(value = allocVector(REALSXP,1) );
	REAL(value)[0] = 0.0;

	// Compute sum
	for (i = 0; i < n; i++) 
		REAL(value)[0] = REAL(value)[0] + vec[i];

	// Output results
	UNPROTECT(2);
	return value;
}

/*
	vecSumSq:  Computes the sum of the squares elements of a vector
*/
SEXP vecSumSq(SEXP Rvec){
	// SEXP coercion of R objects --------------------------	
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);

	int i, n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP value;
	PROTECT(value = allocVector(REALSXP,1) );
	REAL(value)[0] = 0.0;

	// Compute sum
	for (i = 0; i < n; i++) 
		REAL(value)[0] = REAL(value)[0] + vec[i]*vec[i];

	// Output results
	UNPROTECT(2);
	return value;
}


/*
	max:	Computes the max of the elements of a vector
*/
SEXP max(SEXP Rvec){
	// SEXP coercion of R objects --------------------------	
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);
	
	int i, n = length(Rvec);
	
	// Store result here + Initialize to first element
	SEXP maxValue;
	PROTECT(maxValue = allocVector(REALSXP,1) );
	REAL(maxValue)[0] = vec[0];
	
	// Compute Max
	for (i = 1; i < n; i++){
		if(REAL(maxValue)[0] < vec[i]){
			REAL(maxValue)[0] = vec[i];
		}
	}

	
	// Output results
	UNPROTECT(2);	
	return maxValue;
}


/*
	whichMax:	Computes the index of the max of the vector
				If there are ties, it returns the smallest index
*/
SEXP whichMax(SEXP Rvec){
	// SEXP coercion of R objects --------------------------	
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);

	int i, n = length(Rvec);
	
	// Store result here + Initialize to first element
	SEXP whichMaxValue, maxValue;
	PROTECT(whichMaxValue = allocVector(INTSXP,1) );
	INTEGER(whichMaxValue)[0] = 0;
	
	PROTECT(maxValue = allocVector(REALSXP,1) );	
	REAL(maxValue)[0] = vec[0];

	// Compute Index of Max
	for (i = 1; i < n; i++){
		if(REAL(maxValue)[0] < vec[i]){
			REAL(maxValue)[0] = vec[i];
			INTEGER(whichMaxValue)[0] = i;
		}
	}

	// Output results
	UNPROTECT(3);		
	return whichMaxValue;
}


/*
	min:	Computes the min of the elements of a vector
*/
SEXP min(SEXP Rvec){
	// SEXP coercion of R objects --------------------------	
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);
	
	int i, n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP minValue;
	PROTECT(minValue = allocVector(REALSXP,1) );
	REAL(minValue)[0] = vec[0];
	
	// Compute Min
	for (i = 1; i < n; i++){
		if(REAL(minValue)[0] > vec[i]){
			REAL(minValue)[0] = vec[i];
		}
	}

	// Output results
	UNPROTECT(2);	
	return minValue;
}


/*
	whichMin:	Computes the index of the min of the vector
				If there are ties, it returns the smallest index
*/
SEXP whichMin(SEXP Rvec){
	// SEXP coercion of R objects --------------------------	
	PROTECT( Rvec = coerceVector(Rvec, REALSXP) );
	double *vec = REAL(Rvec);

	int i,n = length(Rvec);

	// Store result here + Initialize to first element
	SEXP whichMinValue, minValue;
	PROTECT(whichMinValue = allocVector(INTSXP,1) );
	INTEGER(whichMinValue)[0] = 0;

	PROTECT(minValue = allocVector(REALSXP,1) );
	REAL(minValue)[0] = vec[0];

	// Compute Index of Min
	for (i = 1; i < n; i++){

		if(REAL(minValue)[0] > vec[i]){
			REAL(minValue)[0] = vec[i];
			INTEGER(whichMinValue)[0] = i;
		}
	}

	// Output results
	UNPROTECT(3);		
	return whichMinValue;
}


/*
	compute log of density of multinomial 
 */
SEXP ldmultinom(SEXP x_INPUT, SEXP n_INPUT, SEXP prob_INPUT){
	// SEXP coercion of R objects --------------------------	
	PROTECT( x_INPUT = coerceVector(x_INPUT, REALSXP) );	
	PROTECT( n_INPUT = coerceVector(n_INPUT, REALSXP) );	
	PROTECT( prob_INPUT = coerceVector(prob_INPUT, REALSXP) );	
	
	double n = REAL(n_INPUT)[0], *x = REAL(x_INPUT);
	int i;
	R_len_t len = length(x_INPUT);
	
	// Store result here + Initialize to first element
	SEXP p, value, sum_p;
	PROTECT(p = allocVector(REALSXP,len) );
	PROTECT(value = allocVector(REALSXP,1) );
	PROTECT(sum_p = allocVector(REALSXP,1) );
	
	sum_p = vecSum(prob_INPUT);
	
	// normalize probabilities if needed
	for(i=0; i<len; i++){
		REAL(p)[i] = REAL(prob_INPUT)[i]/REAL(sum_p)[0];
	}
	
	REAL(value)[0] = lgamma(n + 1);
	for(i=0; i<len; i++){
		if(REAL(p)[i]!=0){
			REAL(value)[0] -= lgamma(x[i] + 1);
			REAL(value)[0] += x[i]*log(REAL(p)[i]);
		}
	}
	
	
	// Output results
	UNPROTECT(6);		
	return value;
}













/*
------------------------------------------------------------------------

Main Functions

------------------------------------------------------------------------
*/

/*
	countZones:  compute the total number of zones defined by the list of vectors
	nearestNeighbors_INPUT
	
	INPUT:
	-nearestNeighbors_INPUT:  a list of length n, where each element corresponds
	 to the nearest neighbors of a particular area, listed in order of distance
	 
	OUTPUT:
	-nZones:  the total number of zones
*/

SEXP countZones( SEXP nearestNeighbors_INPUT ){
	// SEXP coercion of R objects --------------------------	
	PROTECT( nearestNeighbors_INPUT = coerceVector(nearestNeighbors_INPUT, VECSXP) );

	int  i, j, n = length(nearestNeighbors_INPUT);
	SEXP nZones;
	PROTECT( nZones = allocVector(INTSXP,1) );
	INTEGER(nZones)[0] = 0; 

	for(i = 0; i < n; i++){
			for(j = 0; j < length(VECTOR_ELT(nearestNeighbors_INPUT,i)); j++){
				INTEGER(nZones)[0]++;
			}
	}

	// Output results
	UNPROTECT(2);
	return(nZones);
}


/*
	Computes the internally standardized expected number of counts
	
	IMPORTANT:  We assume that the population and cases vectors are balanced:
				all counts are sorted by area first, and then within each area the 
				counts for all strata are listed (even if 0 count) in the same order.	
				That way we can loop over indices easily.
 
	INPUT:
	-population_INPUT:  population for each strata in each area
	-cases_INPUT:  cases for each strata in each ares
	-nStrata_INPUT:  number of strata considered
 
	OUTPUT:
	-E:  expected counts
	-q:  strata-specific rates
*/
SEXP expected( SEXP population_INPUT, SEXP cases_INPUT, SEXP nStrata_INPUT ){
	// SEXP coercion of R objects --------------------------	
	PROTECT( population_INPUT = coerceVector(population_INPUT, REALSXP) );
	PROTECT( cases_INPUT = coerceVector(cases_INPUT, REALSXP) );
	PROTECT( nStrata_INPUT = coerceVector(nStrata_INPUT, INTSXP) );

	double	*population = REAL(population_INPUT), 
			*cases = REAL(cases_INPUT);
	int		nStrata = INTEGER(nStrata_INPUT)[0],
			n = length(population_INPUT)/nStrata;
	int i, j;


	// Creation of SEXP objects --------------------------	
	SEXP qNum, qDenom, q, E;	
	PROTECT(q = allocVector(REALSXP,nStrata) );
	PROTECT(E = allocVector(REALSXP,n) );
	PROTECT(qNum = allocVector(REALSXP,nStrata) );
	PROTECT(qDenom = allocVector(REALSXP,nStrata) );


	// Initialize vectors to 0
	for(i = 0; i < nStrata; i++){
		REAL(qNum)[i] = 0.0; REAL(qDenom)[i] = 0.0; REAL(q)[i]=0.0;
	}
	for(i = 0; i < n; i++){
		REAL(E)[i] = 0.0;
	}


	// Compute q:  strata-specific rates.  We do the numerator and denominator separately
	for(i = 0; i < nStrata; i++){
		for(j = 0; j < n; j++){
			REAL(qNum)[i] = REAL(qNum)[i] +  cases[j*nStrata + i];
			REAL(qDenom)[i] = REAL(qDenom)[i] +  population[j*nStrata + i];
		}		
	}
	for(i = 0; i < nStrata; i++){
		REAL(q)[i] = REAL(qNum)[i]/REAL(qDenom)[i];
	}
	// unprotect num and denominator
	UNPROTECT(2);


	// Compute E expected counts
	for(i=0; i < n; i++)
		for(j=0; j < nStrata; j++){
			REAL(E)[i] = REAL(E)[i] + REAL(q)[j] * population[i*nStrata + j];
		}


	// output results as R list
	SEXP list, list_names;
	PROTECT(list_names = allocVector(STRSXP, 2));
	SET_STRING_ELT(list_names, 0,  mkChar("E"));
	SET_STRING_ELT(list_names, 1,  mkChar("q"));
	
	PROTECT(list = allocVector(VECSXP, 2));	// Creating list w/2 vector elements 
	SET_VECTOR_ELT(list, 0, E);				// attaching vector to list
	SET_VECTOR_ELT(list, 1, q);      

	setAttrib(list, R_NamesSymbol, list_names); //and attaching the vector names


	// Output results
	UNPROTECT(7);
	return list;
}










/*
------------------------------------------------------------------------

Kulldorff Functions

------------------------------------------------------------------------
*/

/*
	Compute binomial log-likelihood statistic as described by Kulldorff

	NO EXPECTED COUNTS (STRATA LEVEL DATA) HERE
	
	cz: count of cases in Zone z
	nz: population (or expected count) in Zone z
	C: Total number of cases in study area
	N: Total population (total expected count) in study area
	
*/
SEXP binomialLogLkhd(SEXP cz_INPUT, SEXP nz_INPUT, SEXP N_INPUT, SEXP C_INPUT){
	// SEXP coercion of R objects --------------------------	
	PROTECT( cz_INPUT = coerceVector(cz_INPUT, REALSXP) );
	PROTECT( nz_INPUT = coerceVector(nz_INPUT, REALSXP) );
	PROTECT( N_INPUT = coerceVector(N_INPUT, REALSXP) );
	PROTECT( C_INPUT = coerceVector(C_INPUT, REALSXP) );
	
	double	cz = REAL(cz_INPUT)[0], nz = REAL(nz_INPUT)[0],
			N = REAL(N_INPUT)[0], C = REAL(C_INPUT)[0];


	// Creation of SEXP objects --------------------------	
	SEXP logLkhd;
	PROTECT( logLkhd = allocVector(REALSXP,1) );
	double *xlogLkhd; 
	xlogLkhd = REAL(logLkhd);
	
	// Compute Binomial log-likelihood
	if(cz / nz <= (C - cz)/(N - nz)){
		xlogLkhd[0] = 0;
	}else{
		xlogLkhd[0]  =	N * ( log( N-nz-C+cz ) - log( N-C )  + log( N ) - log( N-nz ) ) +
						nz * ( log( nz-cz ) - log( N-nz-C+cz )  + log( N-nz ) - log( nz ) )+
						cz * ( log( cz ) - log( nz-cz ) + log( N-nz-C+cz ) - log( C-cz ) )+ 
						C * ( log( C-cz ) - log( C ) + log( N-C ) - log( N-nz-C+cz ) );								
	}
	

	// Output results
	UNPROTECT(5);
	return(logLkhd);
}


/*
	Compute Poisson log-likelihood statistic as described by Kulldorff
	
	cz: count of cases in Zone z
	nz: population (or expected count) in Zone z
	C: Total number of cases in study area
	N: Total population (total expected count) in study area
*/
SEXP poissonLogLkhd(SEXP cz_INPUT, SEXP nz_INPUT, SEXP N_INPUT, SEXP C_INPUT){
	// SEXP coercion of R objects --------------------------	
	PROTECT( cz_INPUT = coerceVector(cz_INPUT, REALSXP) );
	PROTECT( nz_INPUT = coerceVector(nz_INPUT, REALSXP) );
	PROTECT( N_INPUT = coerceVector(N_INPUT, REALSXP) );
	PROTECT( C_INPUT = coerceVector(C_INPUT, REALSXP) );
	
	double	cz = REAL(cz_INPUT)[0], nz = REAL(nz_INPUT)[0],
			N = REAL(N_INPUT)[0], C = REAL(C_INPUT)[0];
	
	// Creation of SEXP objects --------------------------	
	SEXP logLkhd;
	PROTECT( logLkhd = allocVector(REALSXP,1) );
	double *xlogLkhd; 
	xlogLkhd = REAL(logLkhd);
	

	// Compute Poisson log-likelihood
	if(cz / nz <= (C - cz)/(N - nz)){
		xlogLkhd[0] = 0;
	}else{
		xlogLkhd[0]  =	cz * log(  (cz / nz) )  +
						cz * log(  ( (N - nz)/( C - cz) ) )  +
						C * log(  ( (C-cz)/(N-nz) )  ) + 
						C * log(  ( N/ C )  );
	}
	

	// Output results
	UNPROTECT(5);
	return(logLkhd);
}




/*
	computeAllLogLkhd:  computes the log-likelihood ratio statistic, for all possible zones
	as defined in the list nearestNeighbors_INPUT.
	
	INPUT:	
	-observedCases_INPUT:	the number of cases in each of the n areas
	-expectedCases_INPUT:	the expected number of cases in each of the n areas
	-nearestNeighbors_INPUT:  a list of length n, each element being that area's nearest neighbors
	-logLkhdType_INPUT:  the type of log-likelihood (either poisson or binomial)
	
	OUTPUT:
	-logLkhd:  a vector of the log-likelihood of each zone, in order they are listed in nearestNeighbors_INPUT
*/
SEXP computeAllLogLkhd( SEXP observedCases_INPUT, SEXP expectedCases_INPUT, SEXP nearestNeighbors_INPUT, SEXP logLkhdType_INPUT ){
	// SEXP coercion of R objects --------------------------	
	PROTECT( observedCases_INPUT = coerceVector(observedCases_INPUT, REALSXP) );
	PROTECT( expectedCases_INPUT = coerceVector(expectedCases_INPUT, REALSXP) );
	PROTECT( nearestNeighbors_INPUT = coerceVector(nearestNeighbors_INPUT, VECSXP) );
	PROTECT( logLkhdType_INPUT = coerceVector(logLkhdType_INPUT, STRSXP) );

	// Define variables
	double *countyCases = REAL(observedCases_INPUT), *expectedCases = REAL(expectedCases_INPUT);
	R_len_t n = length(nearestNeighbors_INPUT);						// number of areas
	int	counter = 0, i, j, lkhd_type, 
		nZones = INTEGER(countZones(nearestNeighbors_INPUT))[0];	// number of zones
	
	// Creation of SEXP objects --------------------------	
	SEXP logLkhd, cz, nz, N, C;
	
	PROTECT(logLkhd = allocVector(REALSXP,nZones) );	// store results here
	PROTECT(cz = allocVector(REALSXP,1) );				// number of cases in a particular zone
	PROTECT(nz = allocVector(REALSXP,1) );				// population or expected cases in a particular zone
	PROTECT(C = allocVector(REALSXP,1) );				// total number of cases in study region
	C = vecSum(observedCases_INPUT);
	PROTECT(N = allocVector(REALSXP,1) );				// total population or total expected cases of study region
	N = vecSum(expectedCases_INPUT);
	
	
	// Determine likelihood type
	if( strcmp( CHAR(STRING_ELT(logLkhdType_INPUT,0)),"poisson") == 0 ){
		lkhd_type = 0;
	}else if( strcmp( CHAR(STRING_ELT(logLkhdType_INPUT,0)),"binomial") == 0 ){
		lkhd_type = 1;
	}



	/*
		Main loop:  For all possible zones, compute the log-likelihood ratio statistic
	*/
	// For all areas
	for(i = 0; i < n; i++){

		REAL(cz)[0] = 0;  REAL(nz)[0] = 0;
	
		// For each area's nearest neighbours
		for(j = 0; j < length(VECTOR_ELT(nearestNeighbors_INPUT,i)); j++){
		
			// Watch off by 1 vector indexing in C as opposed to R
			REAL(cz)[0] = REAL(cz)[0] + countyCases[ INTEGER(VECTOR_ELT(nearestNeighbors_INPUT,i))[j] - 1 ];
			REAL(nz)[0] = REAL(nz)[0] + expectedCases[ INTEGER(VECTOR_ELT(nearestNeighbors_INPUT,i))[j] - 1  ];			
				
			// Compute likelihood statistic
			switch (lkhd_type){
				case 0:
					REAL(logLkhd)[counter] = REAL(poissonLogLkhd(cz, nz, N, C))[0];
				break;

				case 1:
					REAL(logLkhd)[counter] = REAL(binomialLogLkhd(cz, nz, N, C))[0];
				break;
			}
			counter++;
		}
	}

	// Output results
	UNPROTECT(9);
	return(logLkhd);
}


/*
	kulldorffMC:  does the Monte Carlo simulation described by Kulldorff:  based on permutations of 
		the total number of cases under the null hypothesis, compute the most likely cluster
	
	INPUT:	
	-permutedCases_INPUT:	a matrix of simulated numbers of cases under the null hypothesis
	-expectedCases_INPUT:	the expected number of cases in each of the n areas
	-nearestNeighbors_INPUT:  a list of length n, each element being that area's nearest neighbors
	-logLkhdType_INPUT:  the type of log-likelihood (either poisson or binomial)
	
	OUTPUT:
	-maxLogLkhd:  a vector of the max log-likelihood for each simulation
*/
SEXP kulldorffMC( SEXP permutedCases_INPUT, SEXP expectedCases_INPUT, SEXP nearestNeighbors_INPUT, SEXP logLkhdType_INPUT ){
	// SEXP coercion of R objects --------------------------	
	PROTECT( permutedCases_INPUT = coerceVector(permutedCases_INPUT, REALSXP) );
	PROTECT( expectedCases_INPUT = coerceVector(expectedCases_INPUT, REALSXP) );
	PROTECT( nearestNeighbors_INPUT = coerceVector(nearestNeighbors_INPUT, VECSXP) );
	PROTECT( logLkhdType_INPUT = coerceVector(logLkhdType_INPUT, STRSXP) );

	// Define variables
	R_len_t n = length(expectedCases_INPUT);
	int	i, j, nSimulations = length(permutedCases_INPUT)/n, nZones = INTEGER(countZones(nearestNeighbors_INPUT))[0];

	// Creation of SEXP objects --------------------------	
	SEXP maxLogLkhd, permutedCases, logLkhd;
	PROTECT( maxLogLkhd = allocVector(REALSXP, nSimulations) );		// store max log-likelihood for each simulation here
	PROTECT( logLkhd = allocVector(REALSXP, nZones) );				// store all log-likelihoods for a particular simulation
	PROTECT( permutedCases = allocVector(REALSXP, n) );				// used to pick out cases from a particular simulation


	/*
		Main loop:  For each simulation, compute all the log-likelihood ratio statistics for each zone and then return the max
	*/
	for( j=0; j < nSimulations; j++){
		// Load simulatedCases
		for( i = 0; i < n; i++){
			REAL(permutedCases)[i] = REAL(permutedCases_INPUT)[i + n*j];
		}
		
		// compute max log likelihood
		logLkhd =  computeAllLogLkhd( permutedCases, expectedCases_INPUT, nearestNeighbors_INPUT, logLkhdType_INPUT );
		REAL(maxLogLkhd)[j] =  REAL(max(logLkhd))[0];
	}


	// Output results	
	UNPROTECT(7);
	return maxLogLkhd;
}





/*
------------------------------------------------------------------------

Besag Newell Functions

------------------------------------------------------------------------
*/


/*
	besagNewell:  does the Monte Carlo simulation described by Kulldorff:  based on permutations of 
		the total number of cases under the null hypothesis, compute the most likely cluster
	
	Note:  this method differs from the original Besag-Newell method in that:
	-we are looking for k cases, not k OTHER cases
	-we base p-values on the observed number of cases in the circle that contains at least k cases,
	 not k cases

	INPUT:	
	-observedCases_INPUT:	the number of cases in each of the n areas
	-expectedCases_INPUT:	the expected number of cases in each of the n areas
	-nearestNeighbors_INPUT:  a list of length n, each element being that area's nearest neighbors
	-k_INPUT:  the number of cases we wish to include in our circles
		
	OUTPUT:
	A list with 3 elements
	-observedPValues:  the p-values for each area
	-observedMValues:  the observed number of areas we need to consider to observe k cases
	-observedKValues:  the observed number of cases in each circle
	
*/
SEXP besagNewell( SEXP observedCases_INPUT, SEXP expectedCases_INPUT, SEXP nearestNeighbors_INPUT, SEXP k_INPUT){
	// SEXP coercion of R objects --------------------------	
	PROTECT( observedCases_INPUT = coerceVector(observedCases_INPUT, REALSXP) );
	PROTECT( expectedCases_INPUT = coerceVector(expectedCases_INPUT, REALSXP) );
	PROTECT( nearestNeighbors_INPUT = coerceVector(nearestNeighbors_INPUT, VECSXP) );
	PROTECT( k_INPUT = coerceVector(k_INPUT, INTSXP));

	// Define variables
	double *observedCases = REAL(observedCases_INPUT), *expectedCases = REAL(expectedCases_INPUT);
	double sumObserved, sumExpected;
	R_len_t n = length(nearestNeighbors_INPUT);
	int i, j, k = INTEGER(k_INPUT)[0];

	// Creation of SEXP objects --------------------------	
	SEXP observedPValues, observedMValues, observedKValues;
	PROTECT(observedPValues = allocVector(REALSXP,n) );		// store p-values for each area here
	PROTECT(observedMValues = allocVector(INTSXP,n) );		// store number of neighbors needed to have k cases
	PROTECT(observedKValues = allocVector(INTSXP,n) );		// store number of observed cases s.t. we are considering k or more cases
	
	
	
	/*
		Main loop:	For all areas, consider neighbors until k cases are observed, then base p-values on
					observed cases, not k.
	*/
	for(i=0; i<n; i++){
		// For all nearest neighbours.  
		sumObserved = 0;  sumExpected = 0;

		for(j = 0; j < length(VECTOR_ELT(nearestNeighbors_INPUT,i)); j++){
				// Watch off by 1 vector indexing in C as opposed to R
				sumObserved = sumObserved + observedCases[ INTEGER(VECTOR_ELT(nearestNeighbors_INPUT,i))[j] - 1 ];
				sumExpected = sumExpected + expectedCases[ INTEGER(VECTOR_ELT(nearestNeighbors_INPUT,i))[j] - 1 ];
			
				// Difference from original method:  stop when we see k cases, not k OTHER cases
				if( sumObserved >= k ){
					INTEGER(observedKValues)[i] = sumObserved;
					break;
				}
		}
	
		// Watch off by 1 vector indexing in C as opposed to R
		INTEGER(observedMValues)[i] = j + 1;

		// Difference from original method:  use actual observed k, not specified k
		REAL(observedPValues)[i] = ppois(sumObserved-1, sumExpected, 1, 0);
		REAL(observedPValues)[i] = 1 - REAL(observedPValues)[i];	
	}

	
	// output results as R list
	SEXP list, list_names;
	PROTECT(list_names = allocVector(STRSXP, 3));
	SET_STRING_ELT(list_names, 0,  mkChar("observed.p.values"));
	SET_STRING_ELT(list_names, 1,  mkChar("observed.m.values"));
	SET_STRING_ELT(list_names, 2,  mkChar("observed.k.values"));

	PROTECT(list = allocVector(VECSXP, 3));			// Creating list with 3 vector elements 
	SET_VECTOR_ELT(list, 0, observedPValues);		// attaching vector to list
	SET_VECTOR_ELT(list, 1, observedMValues);      
	SET_VECTOR_ELT(list, 2, observedKValues);      

	setAttrib(list, R_NamesSymbol, list_names); //and attaching the vector names


	// Output results
	UNPROTECT(9);
	return(list);
}





