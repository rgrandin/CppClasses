/**
 * @file RootFinding.cpp
 * @author Robert Grandin
 * @brief Implementation of RootFinding class.
 */


#include "RootFinding.h"    /* Added for syntax hilighting */


/* ==================================================================
 * ================
 * ================    PRIVATE FUNCTIONS
 * ================
 */


//	NONE




/* ==================================================================
 * ================
 * ================    PUBLIC FUNCTIONS
 * ================
 */
template <class T>
RootFinding<T>::RootFinding()
{
    warnenabled = true;
}


template <class T>
RootFinding<T>::RootFinding(const RootFinding<T> &a) : RootFinding()
{
    RootFindingSwap(*this, a);
}


#ifdef CXX11
template <class T>
RootFinding<T>::RootFinding(RootFinding<T> &&a) : RootFinding<T>()
{
    RootFindingSwap(*this, a);
}
#endif


template <class T>
RootFinding<T>::~RootFinding()
{
	; // NOTHING TO BE DONE
}


template <class T>
RootFinding<T>& RootFinding<T>::operator=(RootFinding<T> a)
{
    RootFindingSwap(*this, a);
    return *this;
}


template <class T>
T RootFinding<T>::f(Array1D<T> &varval)
{
	; // NO CODE.  THIS FUNCTION EXISTS FOR DEMONSTRATION PURPOSES ONLY!
}


// BISECTION METHOD
template <class T>
T RootFinding<T>::Bisection(Array1D<T> &varval, const int rvar, const T a, const T b, const T tol,
                            const int maxiter, T &tolachieved, int &itersused, int &fevals, T (*f)(Array1D<T> &),
                            const bool warn)
{
    T retval = (T)0.0e0;
    T currval = (T)1.0e5;   /* Defaults large to avoid incorrectly recognizing convergence */
    T leftval = (T)1.0e5;
    T rightval = (T)1.0e5;
    T aa = a;
    T bb = b;
    T midpoint = (T)0.5e0*(aa + bb);
    itersused = 0;
    fevals = 0;
    tolachieved = 1.0e0;

    currval = f(varval);
    fevals++;

    /* Check if function satisfies convergence condition */
    if(fabs(currval) < tol){
        retval = varval(rvar);
    } else {
        /* Check that a < b.  If not, swap their values. */
        if(aa > bb){
            T tmp = bb;
            bb = aa;
            aa = tmp;
        }

        /* Loop through routine while either the convergence tolerance has not been met or the
           maximum number of iterations has not been exceeded.  For both cases, the return value
           must not be equal to NaN.  NaN is checked by comparing retval to itself since NaN != NaN.
        */
        while(fabs(currval) >= tol && itersused <= maxiter && retval == retval){
            /* Evaluate function at section endpoints */
            varval(rvar) = aa;
            leftval = f(varval);
            fevals++;
            varval(rvar) = bb;
            rightval = f(varval);
            fevals++;

            /* Check that signs of leftval and rightval differ */
            if((leftval > 0 && rightval > 0) || (leftval < 0 && rightval < 0)){
                if(warn == true){
                    std::cerr << "ERROR: Specified region does not contain a zero-point!" << std::endl;
                }
                retval = std::numeric_limits<T>::quiet_NaN();
            }

            midpoint = (T)0.5e0*(aa + bb);
            varval(rvar) = midpoint;
            currval = f(varval);
            fevals++;

            /* Debugging Output *//*
            std::cout << "f(" << a << ") = " << leftval << "     f(" << midpoint << ") = " << currval <<
                         "     f(" << b << ") = " << rightval << std::endl;
            */

            if((leftval > 0 && currval < 0) || (leftval < 0 && currval > 0)){
                /* Solution is within [a,midpoint] */
                bb = midpoint;
            } else {
                /* Solution is within [midpoint,b] */
                aa = midpoint;
            }

            /* Increment iteration counter */
            itersused++;
        }

        /* Set return value to current section midpoint location if it is not NaN. */
        if(retval == retval){
            retval = midpoint;
        }
    }

    tolachieved = fabs(currval);
    return retval;

}


// NEWTON'S METHOD
template <class T>
T RootFinding<T>::Newton(Array1D<T> &varval, const int rvar, const T p, const T tol, const int maxiter,
        T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T> &), const bool warn)
{
	/*
	 * THIS FUNCTION ACCEPTS A VECTOR OF INPUT VARIABLES WITH THE NUMBER OF
	 * VARIABLES INICATED BY nvar AND THE VARIABLE TO BE ITERATED ON INDICATED
	 * BY rvar.  ITERATION CONTINUES UNTIL EITHER A TOLERANCE tol IS REACHED OR
	 * THE MAXIMUM NUMBER OF ALLOWED ITERATIONS iterlim IS REACHED.  THE FUNCTION
	 * TO BE ZEROED IS PASSED-IN AS THE FINAL ARGUMENT AND IS EXPECTED TO ACCEPT
	 * A VECTOR OF VARIABLE VALUES AND THE NUMBER OF VARIABLES IN THAT FUNCTION.
	 */

	// ALLOCATE TEMPORARY VECTOR ARRAYS
	int nvar = varval.GetDim();
    Array1D<T> tmpvars(nvar,(T)0.0e0);
	T xlower, xupper, xnew;
	T f1, f2, fcurrent, fp;
    T retval = (T)0.0e0;
    itersused = 0;
    fevals = 0;
    tolachieved = 1.0e0;

	fcurrent = f(varval);	// GET CURRENT FUNCTION VALUE
    fevals++;

    while(fabs(fcurrent) >= tol && itersused <= maxiter){

		// ESTIMATE DERIVATIVE WITH CENTERED DIFFERENCE FORMULA

		// CREATE TEMPORARY VECTOR FOR FIRST POINT
		for(int i=0; i<nvar; i++){
			if(i != rvar){
				tmpvars(i) = varval(i);
			}
			if(i == rvar){
                tmpvars(i) = varval(i)*((T)1.0e0 - p);
                if(abs(varval(i)) < (T)1.0e-4){
                    tmpvars(i) = varval(i) -= (T)0.01e0;
				}
				xlower = tmpvars(i);
			}
		}

		// FIRST POINT FOR DERIVATIVE
		f1 = f(tmpvars);

		// CREATE TEMPORARY VECTOR FOR SECOND POINT
		for(int i=0; i<nvar; i++){
			if(i != rvar){
				tmpvars(i) = varval(i);
			}
			if(i == rvar){
                tmpvars(i) = varval(i)*((T)1.0e0 + p);
                if(abs(varval(i)) < (T)1.0e-4){
                    tmpvars(i) = varval(i) += (T)0.01e0;
				}
				xupper = tmpvars(i);
			}
		}

		// SECOND POINT FOR DERIVATIVE
		f2 = f(tmpvars);
        fevals++;

		// DERIVATIVE ESTIMATE
		fp = (f2 - f1)/(xupper - xlower);

		// ROOT ESTIMATE
		xnew = varval(rvar) - fcurrent/fp;
		retval = xnew;

		// ASSIGN ROOT ESTIMATE TO APPROPRIATE VARIABLE IN THE VARIABLE VECTOR
		varval(rvar) = xnew;

		fcurrent = f(varval);	// GET UPDATED FUNCTION VALUE
        fevals++;
        itersused += 1;	// UPDATE ITERATION COUNT
	}

    tolachieved = fabs(fcurrent);
	return retval;
}


// SECANT METHOD
template <class T>
T RootFinding<T>::Secant(Array1D<T> &varval, const int rvar, const T p, const T tol, const int maxiter,
        T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T> &), const bool warn)
{
	/*
	 * THIS FUNCTION ACCEPTS A VECTOR OF INPUT VARIABLES WITH THE NUMBER OF
	 * VARIABLES INICATED BY nvar AND THE VARIABLE TO BE ITERATED ON INDICATED
	 * BY rvar.  ITERATION CONTINUES UNTIL EITHER A TOLERANCE tol IS REACHED OR
	 * THE MAXIMUM NUMBER OF ALLOWED ITERATIONS iterlim IS REACHED.  THE FUNCTION
	 * TO BE ZEROED IS PASSED-IN AS THE FINAL ARGUMENT AND IS EXPECTED TO ACCEPT
	 * A VECTOR OF VARIABLE VALUES AND THE NUMBER OF VARIABLES IN THAT FUNCTION.
	 */

    T xold, xoldtmp;
    T f1, f2, fcurrent;
    T retval = (T)0.0e0;
    itersused = 0;
    fevals = 0;
    tolachieved = 1.0e0;
	T xnew;


	// START-UP: GET TWO FUNCTION VALUES

	// GET CURRENT FUNCTION VALUE
	fcurrent = f(varval);
    fevals++;
	f1 = fcurrent;
	xold = varval(rvar);

    // INCREASE BY p FOR SECOND FUNCTION VALUE
    if(varval(rvar) < tol){
        varval(rvar) = p;   /* Handle case where x ~ 0 */
    } else {
        varval(rvar) *= (1.0e0 + p);
    }
	xnew = varval(rvar);
    f2 = f(varval);
    fevals++;

    while(fabs(fcurrent) >= tol && itersused <= maxiter){

        //std::cout << "i: " << itersused << "    f1: " << f1 << "   f2: " << f2 <<
        //             "    xnew: " << xnew << std::endl;
		// ROOT ESTIMATE
        xoldtmp = xnew;
		xnew = xnew - f2*(xnew - xold)/(f2 - f1);
		retval = xnew;
        xold = xoldtmp;


		// ASSIGN ROOT ESTIMATE TO APPROPRIATE VARIABLE IN THE VARIABLE VECTOR
		varval(rvar) = xnew;

		f1 = f2;
		f2 = f(varval);
        fevals++;
        fcurrent = f2;	// GET UPDATED FUNCTION VALUE
        itersused += 1;	// UPDATE ITERATION COUNT
	}

    tolachieved = fabs(fcurrent);
	return retval;
}


// FIXED POINT ITERATION
template <class T>
T RootFinding<T>::FixedPointIteration(Array1D<T> &varval, const int rvar, const T p, const T tol,
        int maxiter, T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T> &), const bool warn)
{
	/*
	 * THIS FUNCTION ACCEPTS A VECTOR OF INPUT VARIABLES WITH THE NUMBER OF
	 * VARIABLES INICATED BY nvar AND THE VARIABLE TO BE ITERATED ON INDICATED
	 * BY rvar.  ITERATION CONTINUES UNTIL EITHER A TOLERANCE tol IS REACHED OR
	 * THE MAXIMUM NUMBER OF ALLOWED ITERATIONS iterlim IS REACHED.  THE FUNCTION
	 * TO BE ZEROED IS PASSED-IN AS THE FINAL ARGUMENT AND IS EXPECTED TO ACCEPT
	 * A VECTOR OF VARIABLE VALUES AND THE NUMBER OF VARIABLES IN THAT FUNCTION.
	 */

    T retval = (T)0.0e0;
	T xnew;
    itersused = 0;
    fevals = 0;

	// GET CURRENT FUNCTION VALUE
	xnew = f(varval);

    while(fabs(xnew) >= tol && itersused <= maxiter){

		varval(rvar) = xnew;

		// GET UPDATED FUNCTION VALUE
		xnew = f(varval);
        fevals++;

		retval = xnew;

        itersused += 1;	// UPDATE ITERATION COUNT
	}

    tolachieved = fabs(xnew);
	return retval;
}


// INVERSE PARABOLIC
template <class T>
T RootFinding<T>::InverseParabolic(Array1D<T> &varval, const int rvar, const T p, const T tol,
        int maxiter, T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T> &), const bool warn)
{
	/*
	 * THIS FUNCTION ACCEPTS A VECTOR OF INPUT VARIABLES WITH THE NUMBER OF
	 * VARIABLES INICATED BY nvar AND THE VARIABLE TO BE ITERATED ON INDICATED
	 * BY rvar.  ITERATION CONTINUES UNTIL EITHER A TOLERANCE tol IS REACHED OR
	 * THE MAXIMUM NUMBER OF ALLOWED ITERATIONS iterlim IS REACHED.  THE FUNCTION
	 * TO BE ZEROED IS PASSED-IN AS THE FINAL ARGUMENT AND IS EXPECTED TO ACCEPT
	 * A VECTOR OF VARIABLE VALUES AND THE NUMBER OF VARIABLES IN THAT FUNCTION.
	 */

	T f1, f2, f3, x1, x2, x3;
    T retval = (T)0.0e0;
    itersused = 0;
    fevals = 0;
	T xnew;

	f1 = f(varval);
    fevals++;
	x1 = varval(rvar);
    T xtmp = x1;

    if(varval(rvar) < tol){
        varval(rvar) = p;   /* Handle case where x ~ 0 */
    } else {
        varval(rvar) = xtmp*(1.0e0 + p);
    }
	f2 = f(varval);
    fevals++;
	x2 = varval(rvar);

    if(varval(rvar) < tol){
        varval(rvar) = -p;   /* Handle case where x ~ 0 */
    } else {
        varval(rvar) = xtmp*(1.0e0 - p);
    }
	f3 = f(varval);
    fevals++;
	x3 = varval(rvar);


    while(fabs(f3) >= tol && itersused <= maxiter){

		xnew = x1*(f2)*(f3)/((f1-f2)*(f1-f3)) + x2*(f1)*(f3)/((f2-f1)*(f2-f3)) +
				x3*(f1)*(f2)/((f3-f1)*(f3-f2));

		varval(rvar) = xnew;
		retval = xnew;

		x1 = x2;
		x2 = x3;
		x3 = xnew;
		f1 = f2;
		f2 = f3;
		f3 = f(varval);
        fevals++;

        itersused += 1;						// UPDATE ITERATION COUNT
	}

    tolachieved = fabs(f3);
	return retval;
}


template <class T>
T RootFinding<T>::Muller(Array1D<T> &varval, const int rvar, const T p, const T tol, const int maxiter,
        T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T> &), const bool warn)
{
	; // NO CODE FOR METHOD YET
}

