/**
 * @file Stats.cpp
 * @author Robert Grandin
 * @brief Implementation of Stats class.
 */


#include "Stats.h"

template <class T>
Stats<T>::Stats()
{
    PI_Stats = (T)3.0e0;
    if(typeid(T) == typeid(float)){
        PI_Stats = 3.1415926e0;
    }
    if(typeid(T) == typeid(double)){
        PI_Stats = 3.141592653589793e0;
    }
    if(typeid(T) == typeid(long double)){
        PI_Stats = 3.141592653589793e0;
    }
}


template <class T>
Stats<T>::Stats(const Stats<T> &a) :
    PI_Stats(a.PI_Stats)
{
}


#ifdef CXX11
template <class T>
Stats<T>::Stats(Stats<T> &&a) : Stats()
{
    StatsSwap(*this, a);
}
#endif


template <class T>
Stats<T>::~Stats()
{
}


template <class T>
Stats<T>& Stats<T>::operator=(Stats<T> a)
{
    StatsSwap(*this, a);
    return *this;
}


template <class T>
T Stats<T>::BinomialCoeffSingle(const int n, const int k) const
{
    /*
      This code from the Stat580 HW1 assignment didn't want to work, so the work-around
      used involves calling the function to calculate the entire row of Pascal's triangle
      and then returning the desired entry from the row rather than the entire row itself.

    int np = n + 0;
    int nterms = 0;
    if(k > (np-k)){
        nterms = np - k;
    } else {
        nterms = k;
    }

    T result = (T)1.0e0;
    for(int i=0; i<nterms; i++){
        result *= ((T)np - (T)i)/((T)nterms - (T)i);
    }

    return result;
    */
    Array1D<T> tmp;
    Stats<T>::BinomialCoeffTriangleRow(n,tmp);
    if(k < tmp.GetDim()){
        return tmp(k);
    } else {
        return (T)0.0e0;
    }
}


template <class T>
void Stats<T>::BinomialCoeffTriangleRow(const int n, Array1D<T> &array) const
{
    /* Increase n by 1 since the first row corresponds to n=0 */
    int np = n + 1;

    /* Reset size of array.  Initialied to 1 to create boundary conditions. */
    array.ResetSize(np,(T)1.0e0);

    /* Loop through rows of triangle until desired row is reached. */
    for(int row=0; row < np; row++){
        if(row < np){
            array(row) = (T)1.0e0;
        }
        /* Recurrence relation. */
        for(int k=row-1; k>0; k--){
            array(k) = array(k-1) + array(k);
        }
    }
}


template <class T>
void Stats<T>::BinomialCoeffFullTriangle(const int n, Array2D<T> &array) const
{
    /* Increase n by 1 since the first row corresponds to n=0 */
    int np = n + 1;

    /* Reset size of array.  Initialized to 0 */
    array.ResetSize(np,np,(T)0.0e0);

    /* Set boundary conditions in first two rows */
    array(0,0) = (T)1.0e0;
    array(1,0) = (T)1.0e0;
    array(1,1) = (T)1.0e0;

    /* Loop through array rows. */
    for(int row=2; row<np; row++){   /* Start on 3rd row since rows 1 and 2 are already known */
        array(row,0) = (T)1.0e0;
        array(row,row) = (T)1.0e0;

        /* Recurrence relation */
        for(int k=1; k<row; k++){ /* Loop through columns, skipPI_Statsng first and last */
            array(row,k) = array(row-1,k-1) + array(row-1,k);
        }
    }
}


template <class T>
T Stats<T>::PolyVal(const Array1D<T> &array, const T x) const
{
    T result = (T)0.0e0;
    int n = array.GetDim() - 1; /* "-1" used since size of array is n+1 */

    /* Evaluate using the Horner scheme */
    for(int i=n; i>=0; i--){
        result = x*result + array(i);
    }

    return result;
}


template <class T>
T Stats<T>::PolyDerivVal(const Array1D<T> &array, const int m, const T x)
{
    int n = array.GetDim() - 1;     /* "-1" used since size of array is n+1 */
    int Nm = n - m;

    /*
       The appropriate coefficient for each term of the polynomial is calculated by first
       calculating the leading term.  This term is then corrected for each lower-power term
       by multiplying it by a correction fraction.  These coefficients, when combined
       with the original coefficients of the polynomial, are then the correct coefficients
       for the derivative of the polynomial.
    */

    T multterm = 1.0e0;             /* Build factorial term: n!/(n-m)! */
    for(int i=n; i>Nm; i--){
        multterm *= (T)i;
    }

    T result = array(n)*multterm;   /* Determine first coefficient */

    T num = (T)Nm;                  /* Define numerator for correction of lower-power terms */
    T den = (T)n;                   /* Define denominator for correction of lower-power terms */

    for(int i=n; i>m; i--){
        multterm *= (num/den);      /* Correct factorial term for lower-power terms */
        num -= 1.0e0;               /* Update numerator value */
        den -= 1.0e0;               /* Update denominator value */

        result = x*result + multterm*array(i-1);    /* Compute polynomial */
    }

    return result;
}


template <class T>
T Stats<T>::HermiteEval(const int k, const T x) const
{
    T retval = (T)0.0e0;

    /* Check if k=0, or k=1.  Otherwise, use recurrence relation. */
    switch(k){
    case 0:
        retval = (T)1.0e0;
        break;

    case 1:
        retval = (T)2.0e0*x;
        break;

    default:
        T hkm1 = (T)2.0e0*x;
        T hkm2 = (T)1.0e0;
        for(int kk=2; kk<=k; kk++){
            retval = (T)2.0e0*x*hkm1 - (T)2.0e0*((T)kk - (T)1.0e0)*hkm2;
            hkm2 = hkm1;
            hkm1 = retval;
        }
        break;
    }

    return retval;
}


template <class T>
T Stats<T>::HermiteProbEval(const int k, const T x) const
{
    T retval = (T)0.0e0;

    /* Check if k=0, or k=1.  Otherwise, use recurrence relation. */
    switch(k){
    case 0:
        retval = (T)1.0e0;
        break;

    case 1:
        retval = x;
        break;

    default:
        T hkm1 = x;
        T hkm2 = (T)1.0e0;
        for(int kk=2; kk<=k; kk++){
            retval = x*hkm1 - ((T)kk - (T)1.0e0)*hkm2;
            hkm2 = hkm1;
            hkm1 = retval;
        }
        break;
    }

    return retval;
}


template <class T>
T Stats<T>::LnGamma(const T x) const
{
    T retval = (T)0.0e0;
    if(x <= 0){
        std::cerr << "ERROR: Argument to LnGamma() must be positive" << std::endl;
        retval = std::numeric_limits<T>::quiet_NaN();
    } else {
        double xx = (T)0.0e0;
        double yy = (T)0.0e0;
        double tmp = (T)0.0e0;
        double ser = (T)0.0e0;
        double coeffs[6] = {76.18009172947146e0, -86.50532032941677e0, 24.01409824083091e0,
                       -1.231739572450155e0, 0.1208650973866179e-2, -0.5395239384953e-5};

        yy = (double)x;
        xx = (double)x;
        tmp = (double)x + 5.5e0;
        tmp -= ((double)x + 0.5e0)*log(tmp);
        ser = 1.000000000190015e0;
        for(int j=0; j<6; j++){
            ser += coeffs[j]/++yy;
        }

        retval = (T)(-tmp + log(2.5066282746310005*ser/xx));
    }

    return retval;
}


template <class T>
T Stats<T>::GammaHW(const T x) const
{
    T xx = x - 1.0e0;

    T denom = 12.0e0*xx + 2.0e0/(5.0e0*xx + 53.0e0/42.0e0*xx);
    return pow(xx,xx)*sqrt(2.0e0*PI_Stats*xx)*exp(1.0e0/denom - xx);
}


template <class T>
void Stats<T>::ComputeMoments(const int n, const int k, const T tol, const Array1D<T> &params,
                              Array1D<T> &w_array, T (*f)(const Array1D<T>&, const int), int &termsreqd) const
{
    /* Array to store coefficients of Taylor series expansion of MGF.  Oversized due to unknown
       number of terms needed.
    */
    Array1D<T> tmp(100,(T)1.0e0);
    int nterms = 0;
    T coeffval = (T)1.0e0;
    while(coeffval >= tol && nterms < tmp.GetDim()){
        coeffval = f(params,nterms+1);
        tmp(nterms) = coeffval;

        /* Change coeffval if it's too small and nterms = 0 to require at least one good term */
        if(coeffval < tol && nterms < 5){
            coeffval = 1.0e0;
        }
        nterms++;
    }

    /* Copy coefficients into new array of the proper size.  This array is passed to the derivative
       function to evaluate each moment of the distribution.
    */
    Array1D<T> coeffs(nterms,(T)0.0e0);
    for(int i=0; i<nterms; i++){
        coeffs(i) = tmp(i);
    }

    /* Save the number of terms required */
    termsreqd = nterms;

    /* Resize tmp array to free memory (negligible value here, but good practice */
    tmp.ResetSize(1);

    /* Calculate moments for distribution.  The first k moments are calculated since they are the ones
       to be combined for multiple distributions.  Derivative is evaluated at 0 since this based on a
       moment generating function.
    */
    Stats<T> s;
    Array1D<T> distmoments(k+1,(T)0.0e0);
    for(int kk=1; kk<=k; kk++){
        distmoments(kk) = s.PolyDerivVal(coeffs,kk,(T)0.0e0);
    }

    /* Use the moments for the distribution to determine the moments for the sum of n iid distributions */
    w_array.ResetSize(k+1,(T)1.0e0);
    Array1D<T> ichoosej(k+1,(T)1.0e0);
    for(int i=1; i<=k; i++){
        w_array(i) = (T)0.0e0;
        for(int j=i-1; j>= 0; j--){
            if(i > 1 && j > 0){
                ichoosej(j) = ichoosej(j) + ichoosej(j-1);
            } else {
                ichoosej(j) = (T)1.0e0;
            }
            w_array(i) += ichoosej(j)*(T)(n*(i-j) - j)*distmoments(i-j)*w_array(j);
        }
        w_array(i) /= (T)i;
        ichoosej(i) = (T)1.0e0;
    }
}


template <class T>
void Stats<T>::ComputeMoments(const int n, const int k, const Array1D<T> &m_array, Array1D<T> &w_array) const
{
    w_array.ResetSize(k+1,(T)1.0e0);
    Array1D<T> ichoosej(k+1,(T)1.0e0);
    for(int i=1; i<=k; i++){
        w_array(i) = (T)0.0e0;
        for(int j=i-1; j>= 0; j--){
            if(i > 1 && j > 0){
                ichoosej(j) = ichoosej(j) + ichoosej(j-1);
            } else {
                ichoosej(j) = (T)1.0e0;
            }
            w_array(i) += ichoosej(j)*(T)(n*(i-j) - j)*m_array(i-j)*w_array(j);
        }
        w_array(i) /= (T)i;
        ichoosej(i) = (T)1.0e0;
    }
}


template <class T>
void Stats<T>::ConvertMomentsCumulants(const Array1D<T> &input, Array1D<T> &output, const int flag) const
{
    Stats<T> s; /* Provides access to Gamma function to undo factorials and binomials */

    /* Resize output to match input */
    output.ResetSize(input.GetDim());
    output.ResetVal(0.0e0);

    /* Check if converting cumulants to moments */
    if(flag < 0){
        /* Convert cumulants to moments */
        output(0) = 1.0e0;
        for(int k=1; k<output.GetDim(); k++){
            for(int j=0; j<k; j++){
                output(k) += s.BinomialCoeffSingle((k-1),j)*input(k-j)*output(j);
            }
        }
    } else {
        /* Convert moments to cumulants */
        output(0) = 0.0e0;
        for(int k=1; k<output.GetDim(); k++){
            output(k) = input(k);
            for(int j=1; j<k; j++){
                output(k) -= s.BinomialCoeffSingle((k-j),j)*output(k-j)*input(j);
            }
        }
    }
}


template <class T>
T Stats<T>::GammaCDF(const T k, const T theta, const T x, const T tol, const int maxiter, T &tolachieved, int &itersused) const
{
    T retval = 0.0e0;

    T z = x/theta;

    /* Initialize first term */
    T gammaval = exp(LnGamma(k + 1.0e0));
    retval = 1.0e0/gammaval;

    tolachieved = 1.0e0;
    itersused = 1;

    /* Iterate sum until either tolerance is achieved or maximum number of iterations is encountered */
    while(fabs(tolachieved) >= tol && itersused <= maxiter){
        gammaval *= (k + (T)itersused);                 /* Use relation Gamma(a+1) = a*Gamma(a) */
        tolachieved = pow(z,(T)itersused)/gammaval;
        retval += tolachieved;
        itersused++;

    }

    retval *= exp(-z)*pow(z,k);

    return retval;
}


template <class T>
T Stats<T>::GammaCDF_CFE(const T k, const T theta, const T x, const T tol, const int maxiter, T &tolachieved, int &itersused) const
{
    /*
      This function calculates the CDF of a Gamma distribution using continued fraction expansion.  This is implemented by
      calculating the upper-incomplete Gamma distribution using the CFE coefficients found on Wikipedia
      (http://en.wikipedia.org/wiki/Incomplete_gamma_function).  The returned value is the lower-incomplete Gamma function
      which is simply the upper-incomplete Gamma function subtracted from 1.  The CFE coefficients are implemented using
      the Wallis Algorithm as described in Stat 580.
    */

    T retval = 0.0e0;
    tolachieved = 1.0e0;
    itersused = 0;

    /* Initialize [A,B] coefficients in preparation of the recursion relation (assuming first n in the continued
       fraction is 1
    */
    T A_nm2 = 1.0e0;
    T B_nm2 = 0.0e0;
    T A_nm1 = 0.0e0;
    T B_nm1 = 1.0e0;

    /* Scale x */
    T z = x/theta;

    /* Calculate f_0 */
    retval = (A_nm2*z + A_nm1)/(B_nm2*z + B_nm1);
    itersused++;

    T acoeff = 0.0e0;
    T bcoeff = 0.0e0;
    T mm = 0.0e0;
    T m = 0.0e0;
    T oldval = retval;
    T A_n = 0.0e0;
    T B_n = 0.0e0;

    /* Add terms to continued fraction expansion until convergence is achieved */
    while(fabs(tolachieved) > tol && itersused <= maxiter){
        /* Calculate a_n and b_n terms of CFE.
           Coefficients found at http://en.wikipedia.org/wiki/Incomplete_gamma_function
        */
        m = (T)itersused;
        if(itersused % 2 == 1){
            mm = (m - 1.0e0)*0.5e0;
            if(mm < 1.0e-4){    /* mm ~= 0 */
                //acoeff = pow(z,k)*exp(-z);    /* Factored in later to avoid overflow for large z and k */
                acoeff = 1.0e0;
            } else {
                acoeff = mm*z;
            }
        } else {
            mm = 0.5e0*m;
            acoeff = -(k + mm - 1.0e0)*z;
        }
        bcoeff = k + m - 1.0e0;

        /* Calculate A_n and B_n */
        A_n = bcoeff*A_nm1 + acoeff*A_nm2;
        B_n = bcoeff*B_nm1 + acoeff*B_nm2;

        retval = (A_nm1*z + A_n)/(B_nm1*z + B_n);

        /* Update A and B coefficients */
        A_nm2 = A_nm1;
        B_nm2 = B_nm1;
        A_nm1 = A_n;
        B_nm1 = B_n;

        /* Calculate tolerance as change in 'retval' and update iteration counter */
        tolachieved = retval - oldval;
        oldval = retval;
        itersused++;

    }

    /*
      Current value of 'retval' reflects the continued fraction expansion of the incomplete gamma function
      with '1' in the highest numerator rather than 'z^k * exp(-z)'.  This was done to avoid overflow for
      large z and k.  To return the Gamma CDF, the current value of 'retval' must be multiplied by
      'z^k * exp(-z) / Gamma(k)'.  This term is calculated in the log-domain to help avoid overflow issues
      when calculating the terms which comprise it (the end result is of reasonable magnitude).
    */
    retval = retval*exp(k*log(z) - z - LnGamma(k));

    /* Check for unusual results and report to user if necessary.  Corrections performed if necessary. */
    if(retval != retval){
//        std::cerr << "WARNING: NaN created in GammaCDF_CFE()" << std::endl;
    }
    if(retval < 0.1e0 && retval > 0.0e0 && z > k){
        /*
          For very large values of z, with respect to shape parameter k, the calculated CDF value is
          near-0.  When z is greater than k, the CDF value should be nearly
          equal to 1 and is manually overridden as such.
        */
//        std::cerr << "WARNING: manual override from " << retval << " to 1 in GammaCDF_CFE()" << std::endl;
        retval = 1.0e0;
    }
    if(retval < 0.0e0){
        /*
          As discussed in the above if-block, very large values of z and k can produce unexpected results.
          If the result is negative, it appears as if adding '1' to the result produces an approximately-
          correct answer.  The user-defined accuracy is not met, but several decimal places may be accurate.
        */
//        std::cerr << "WARNING: negative result " << retval
//                  << " added to 1 to produce answer.  Accuracy may suffer!" << std::endl;
        retval += 1.0e0;
    }
    return retval;
}


template <class T>
T Stats<T>::GammaPDF(const T k, const T theta, const T x) const
{
    return exp((k - 1.0e0)*log(x) - (x/theta) - k*log(theta) - LnGamma(k));
}


template <class T>
T Stats<T>::ChiSquaredCDF(const T n, const T twolambda, const T x, const T tol, const int maxiter, T &tolachieved, int &itersused) const
{
    T lambda = 0.5e0*twolambda;

    T retval = 0.0e0;

    /* If lambda is greater than 0 (tol used to allow a little buffer), calculate the non-central chi-squared CDF value */
    if(fabs(lambda) > tol){
        T loctolachieved = 1.0e0;
        int locitersused = 0;

        tolachieved = 1.0e0;

        T denominator = 1.0e0;
        T lambdapower = 1.0e0;

        T z = 0.5e0*x;
        T s = 0.5e0*n;

        T F;        /* Current F term */
        T gammaval; /* Gamma term in denominator */

        /* Calculate first term to start the series */
        F = GammaCDF(s,2.0e0,x,tol,maxiter,loctolachieved,locitersused);
        retval = F;
        gammaval = exp(LnGamma(1.0e0 + 0.5e0*n));
        itersused = 1;

        /* Iterate sum while either tolerance is too large or iteration limit has not been reached */
        while(fabs(tolachieved) >= tol && itersused <= maxiter){
            denominator *= (T)itersused;    /* k! term in denominator */
            lambdapower *= lambda;          /* lambda^k term in numerator */

            s = 0.5e0*n + (T)itersused;
            F -= pow(z,(s-1.0e0))*exp(-z)/gammaval;

            tolachieved = lambdapower*F/denominator;
            retval += tolachieved;
            itersused++;

            gammaval *= s;                  /* Use Gamma(a+1) = a*Gamma(a) */
        }
        retval *= exp(-lambda);
    } else {
        /* If lambda is equal to 0 (within 'tol' of zero), calculate the central chi-squared CDF value */
        retval = GammaCDF(0.5e0*n,2.0e0,x,tol,maxiter,tolachieved,itersused);
    }

    return retval;
}


template <class T>
T Stats<T>::ChiSquaredPDF(const T n, const T twolambda, const T x, const T tol, const int maxiter, T &tolachieved, int &itersused) const
{
    T lambda = 0.5e0*twolambda;

    tolachieved = 1.0e0;
    itersused = 0;

    T denominator = 1.0e0;
    T lambdapower = 1.0e0;

    T retval = 0.0e0;

    /* If lambda is greater than 0 (tol used to allow a little buffer), calculate the non-central chi-squared PDF value */
    if(fabs(lambda) > tol){
        /* Iterate sum while either tolerance is too large or iteration limit has not been reached */
        while(fabs(tolachieved) >= tol && itersused <= maxiter){
            denominator *= (T)itersused;    /* k! term in denominator */
            lambdapower *= lambda;          /* lambda^k term in numerator */
            if(itersused == 0){
                denominator = 1.0e0;
                lambdapower = 1.0e0;
            }
            tolachieved = lambdapower*GammaPDF(0.5e0*(n+2.0e0*itersused),2.0e0,x)/denominator;
            retval += tolachieved;
            itersused++;
        }
        retval *= exp(-lambda);
    } else {
        /* If lambda is equal to 0 (within 'tol' of zero), calculate the central chi-squared PDF value */
        retval = GammaPDF(0.5e0*n,2.0e0,x);
    }

    return retval;
}


template <class T>
T Stats<T>::BetaPDF(const T x, const T alpha, const T beta) const
{
    return (pow(x,(alpha-1.0e0))*pow((1.0e0-x),(beta-1.0e0)))/(exp(LnGamma(alpha)+LnGamma(beta)-LnGamma(alpha+beta)));
}


template <class T>
T Stats<T>::BetaPDFRoot(const T x, const T alpha, const T beta, const T lambda) const
{
    return (pow(x,(alpha-1.0e0))*pow((1.0e0-x),(beta-1.0e0)))/(exp(LnGamma(alpha)+LnGamma(beta)-LnGamma(alpha+beta))) - lambda;
}


template <class T>
T Stats<T>::BetaCDF(const T x, const T alpha, const T beta, const T tol, const int maxiter, T &tolachieved, int &itersused) const
{
    /*
      Per the note in Numerical Recipie's, this approach converges quickly for x < (a + 1)/(a + b + 2).  Using
      symmetry, the alpha and beta parameters can be switched for x greater than this threshold, and x is replaced
      with (1 - x).  This is done to provide rapid convergence for x greater than this threshold.
    */
    T BT = exp(LnGamma(alpha+beta) - LnGamma(alpha) - LnGamma(beta) + alpha*log(x) + beta*log(1.0e0-x));
    T threshold = (alpha + 1.0e0)/(alpha + beta + 2.0e0);
    if(x < threshold){
        return BT*BetaCFE(x,alpha,beta,tol,maxiter,tolachieved,itersused)/alpha;
    } else {
        return 1.0e0 - BT*BetaCFE(1.0-x,beta,alpha,tol,maxiter,tolachieved,itersused)/beta;
    }
}


template <class T>
T Stats<T>::BetaCDF(const T x, const T alpha, const T beta) const
{
    T tolreqd = std::numeric_limits<T>::min();
    T tolachieved = 1.0e0;
    int maxiter = 1000;
    int itersused = 0;

    return BetaCDF(x,alpha,beta,tolreqd,maxiter,tolachieved,itersused);
}


template <class T>
void Stats<T>::BetaCI(const T alpha, const T beta, const T conf, const T tol, T &tolachieved, T &xleft, T &xright, int &fevals_cdf, int &fevals_pdf)
{
    T f1, f2, f3, x1, x2, x3, xnew;
    T tolreqd = tol;
    T p = 0.01e0;
    T xrightinit = xright;
    int maxiter = 1000;
    int itersused = 0;
    int fevals_cdf_tmp = 0;
    int fevals_pdf_tmp = 0;

    fevals_cdf = 0;
    fevals_pdf = 0;

    x1 = 4.5e0; /* Tempermental.  DO NOT CHANGE */
    x2 = x1;
    x3 = x1;
    f1 = BetaCIRoot(alpha,beta,x1,conf,xleft,xright,tolreqd,maxiter,tolachieved,itersused,fevals_cdf_tmp,fevals_pdf_tmp);
    fevals_cdf += fevals_cdf_tmp;
    fevals_pdf += fevals_pdf_tmp;

    if(fabs(xleft-xright) < tol){
        xright = xrightinit;
    }

    if(x2 < tol){
        x2 = p;   /* Handle case where x ~ 0 */
    } else {
        x2 *= ((T)1.0e0 + p);
    }
    f2 = BetaCIRoot(alpha,beta,x2,conf,xleft,xright,tolreqd,maxiter,tolachieved,itersused,fevals_cdf_tmp,fevals_pdf_tmp);
    fevals_cdf += fevals_cdf_tmp;
    fevals_pdf += fevals_pdf_tmp;

    if(fabs(xleft-xright) < tol){
        xright = xrightinit;
    }

    if(x3 < tol){
        x3 = -p;   /* Handle case where x ~ 0 */
    } else {
        x3 *= ((T)1.0e0 - p);
    }
    f3 = BetaCIRoot(alpha,beta,x3,conf,xleft,xright,tolreqd,maxiter,tolachieved,itersused,fevals_cdf_tmp,fevals_pdf_tmp);
    fevals_cdf += fevals_cdf_tmp;
    fevals_pdf += fevals_pdf_tmp;

    if(fabs(xleft-xright) < tol){
        xright = xrightinit;
    }


    while(fabs(f3) >= tol && itersused <= maxiter){

        xnew = x1*(f2)*(f3)/((f1-f2)*(f1-f3)) + x2*(f1)*(f3)/((f2-f1)*(f2-f3)) +
                x3*(f1)*(f2)/((f3-f1)*(f3-f2));

        x1 = x2;
        x2 = x3;
        x3 = xnew;
        f1 = f2;
        f2 = f3;
        f3 = BetaCIRoot(alpha,beta,x3,conf,xleft,xright,tolreqd,maxiter,tolachieved,itersused,fevals_cdf_tmp,fevals_pdf_tmp);
        fevals_cdf += fevals_cdf_tmp;
        fevals_pdf += fevals_pdf_tmp;

        if(fabs(xleft-xright) < tol){
            xright = xrightinit;
        }

        itersused += 1;
    }

    tolachieved = fabs(f3);
}


template <class T>
T Stats<T>::BetaCIRoot(const T alpha, const T beta, const T lambda, const T conf, T &xleft, T &xright,
                       const T tol, const int maxiter, T &tolachieved, int &itersused, int &fevals_cdf, int &fevals_pdf) const
{
    tolachieved = 1.0e0;
    itersused = 0;
    fevals_cdf = 0;
    fevals_pdf = 0;

    /* Find the lower-bound corresponding to the desired PDF value 'lambda' */
    T p = 0.01e0;
    T x1 = xleft;
    T x2 = xleft;
    T x3 = xleft;
    T xnew = xleft;
    T f1 = 1.0e0;
    T f2 = 1.0e0;
    T f3 = 1.0e0;
    f1 = BetaPDFRoot(x1,alpha,beta,lambda);
    fevals_pdf++;

    if(x2 < tol){
        x2 = p;   /* Handle case where x ~ 0 */
    } else {
        x2 *= ((T)1.0e0 + p);
    }
    f2 = BetaPDFRoot(x2,alpha,beta,lambda);
    fevals_pdf++;

    if(x3 < tol){
        x3 = -p;   /* Handle case where x ~ 0 */
    } else {
        x3 *= ((T)1.0e0 - p);
    }
    f3 = BetaPDFRoot(x3,alpha,beta,lambda);
    fevals_pdf++;


    while(fabs(f3) >= tol && itersused <= maxiter){

        xnew = x1*(f2)*(f3)/((f1-f2)*(f1-f3)) + x2*(f1)*(f3)/((f2-f1)*(f2-f3)) +
                x3*(f1)*(f2)/((f3-f1)*(f3-f2));

        if(xnew < 0.0e0){
            xnew = 0.0e0;
        }
        if(xnew > 1.0e0){
            xnew = 1.0e0;
        }

        x1 = x2;
        x2 = x3;
        x3 = xnew;
        f1 = f2;
        f2 = f3;
        f3 = BetaPDFRoot(x3,alpha,beta,lambda);
        fevals_pdf++;

        itersused += 1;
    }
    xleft = xnew;
    tolachieved = fabs(f3);





    /* Find the upper-bound corresponding to the desired PDF value 'lambda' */
    x1 = xright;
    x2 = xright;
    x3 = xright;
    xnew = xright;
    f1 = 1.0e0;
    f2 = 1.0e0;
    f3 = 1.0e0;
    f1 = BetaPDFRoot(x1,alpha,beta,lambda);
    fevals_pdf++;

    if(x2 < tol){
        x2 = p;   /* Handle case where x ~ 0 */
    } else {
        x2 *= ((T)1.0e0 + p);
    }
    f2 = BetaPDFRoot(x2,alpha,beta,lambda);
    fevals_pdf++;

    if(x3 < tol){
        x3 = -p;   /* Handle case where x ~ 0 */
    } else {
        x3 *= ((T)1.0e0 - p);
    }
    f3 = BetaPDFRoot(x3,alpha,beta,lambda);
    fevals_pdf++;


    while(fabs(f3) >= tol && itersused <= maxiter){

        xnew = x1*(f2)*(f3)/((f1-f2)*(f1-f3)) + x2*(f1)*(f3)/((f2-f1)*(f2-f3)) +
                x3*(f1)*(f2)/((f3-f1)*(f3-f2));

        if(xnew < 0.0e0){
            xnew = 0.0e0;
        }
        if(xnew > 1.0e0){
            xnew = 1.0e0;
        }

        x1 = x2;
        x2 = x3;
        x3 = xnew;
        f1 = f2;
        f2 = f3;
        f3 = BetaPDFRoot(x3,alpha,beta,lambda);
        fevals_pdf++;

        itersused += 1;
    }
    xright = xnew;
    tolachieved = fabs(f3);




    /* Compute integral of PDF within [xleft,xright] */
    T CIsize = BetaCDF(xright,alpha,beta) - BetaCDF(xleft,alpha,beta);
    fevals_cdf += 2;

    /* Return difference between integral size and desired CI size */
    return (CIsize - conf);
}


template <class T>
T Stats<T>::BetaCFE(const T x, const T alpha, const T beta, const T tol=std::numeric_limits<T>::min(),
                    const int maxiter=1000, T &tolachieved=(T)NULL, int &itersused=(T)NULL) const
{
    T retval = 0.0e0;
    tolachieved = 1.0e0;
    itersused = 0;

    /* Initialize [A,B] coefficients in preparation of the recursion relation (assuming first n in the continued
       fraction is 1
    */
    T A_nm2 = 1.0e0;
    T B_nm2 = 0.0e0;
    T A_nm1 = 1.0e0;
    T B_nm1 = 1.0e0;

    /* Calculate f_0 */
    retval = (A_nm2*x + A_nm1)/(B_nm2*x + B_nm1);
    itersused++;

    T acoeff = 0.0e0;
    T bcoeff = 0.0e0;
    T oldval = retval;
    T A_n = 0.0e0;
    T B_n = 0.0e0;
    T m = 0.0e0;
    T mm = 0.0e0;

    /* Add terms to continued fraction expansion until convergence is achieved */
    while(fabs(tolachieved) > tol && itersused <= maxiter){
        /* Calculate a_n and b_n terms of CFE.
           Coefficients found at http://apps.nrbook.com/c/index.html, page 227
        */
        m = (T)itersused;
        if(itersused % 2 == 1){
            mm = (m - 1.0e0)*0.5e0;
            acoeff = -(alpha + mm)*(alpha + beta + mm)*x/((alpha + 2.0e0*mm)*(alpha + 2.0e0*mm + 1));
        } else {
            mm = m*0.5e0;
            acoeff = mm*(beta - mm)*x/((alpha + 2.0e0*mm - 1.0e0)*(alpha + 2.0e0*mm));
        }
        bcoeff = 1.0e0;

        /* Calculate A_n and B_n */
        A_n = bcoeff*A_nm1 + acoeff*A_nm2;
        B_n = bcoeff*B_nm1 + acoeff*B_nm2;

        retval = (A_nm1*x + A_n)/(B_nm1*x + B_n);

        /* Update A and B coefficients */
        A_nm2 = A_nm1;
        B_nm2 = B_nm1;
        A_nm1 = A_n;
        B_nm1 = B_n;

        /* Calculate tolerance as change in 'retval' and update iteration counter */
        tolachieved = retval - oldval;
        oldval = retval;
        itersused++;
    }

    return 1.0e0/retval;    /* Return inverse as b_0 = 0, and a_0 = 1 */
}


template <class T>
void Stats<T>::Histogram(const Array1D<T> &input, const int nbins, Array2D<float> &result) const
{
    /* Reset 'result' array based on the specified number of bins */
    result.ResetSize(nbins,2,(float)0.0e0);

    /* Calculate the value range for each bin, assuming uniform bin size */
    size_t di;                                                 /* Dummy variables */
    float fmin = (float)input.MinVal(di);                   /* Minimum value in array */
    float fmax = (float)input.MaxVal(di);                   /* Maximum value in array */
    float binsize = (fmax - fmin)/((float)nbins - 1.0e0);   /* Range of values to be stored in each bin */
    float firstbin = fmin + binsize/2.0e0;                  /* Centered-value of first bin */

    /* Fill first column of 'result' with center values of each bin */
    result(0,0) = firstbin;
    for(int i=1; i<nbins; i++){
        result(i,0) = (T)i*binsize + firstbin;
    }

    /* Loop through array and count the number of elements in each bin */
    int idx = 0;
    for(int i=0; i<input.GetDim(); i++){
        idx = (int)(((float)input(i) - fmin)/binsize);
        result(idx,1) += 1.0e0;
    }
}


template <class T>
void Stats<T>::Histogram(const Array2D<T> &input, const int nbins, Array2D<float> &result) const
{
    /* Reset 'result' array based on the specified number of bins */
    result.ResetSize(nbins,2,(float)0.0e0);

    /* Calculate the value range for each bin, assuming uniform bin size */
    size_t di,dj;                                              /* Dummy variables */
    float fmin = (float)input.MinVal(di,dj);                /* Minimum value in array */
    float fmax = (float)input.MaxVal(di,dj);                /* Maximum value in array */
    float binsize = (fmax - fmin)/((float)nbins - 1.0e0);   /* Range of values to be stored in each bin */
    float firstbin = fmin + binsize/2.0e0;                  /* Centered-value of first bin */

    /* Fill first column of 'result' with center values of each bin */
    result(0,0) = firstbin;
    for(int i=1; i<nbins; i++){
        result(i,0) = (T)i*binsize + firstbin;
    }

    /* Loop through array and count the number of elements in each bin */
    int idx = 0;
    for(int i=0; i<input.GetDim(1); i++){
        for(int j=0; j<input.GetDim(2); j++){
            idx = (int)(((float)input(i,j) - fmin)/binsize);
            result(idx,1) += 1.0e0;
        }
    }
}


template <class T>
void Stats<T>::Histogram(const Array3D<T> &input, const int nbins, Array2D<float> &result) const
{
    /* Reset 'result' array based on the specified number of bins */
    result.ResetSize(nbins,2,(float)0.0e0);

    /* Calculate the value range for each bin, assuming uniform bin size */
    size_t di,dj,dk;                                           /* Dummy variables */
    float fmin = (float)input.MinVal(di,dj,dk);             /* Minimum value in array */
    float fmax = (float)input.MaxVal(di,dj,dk);             /* Maximum value in array */
    float binsize = (fmax - fmin)/((float)nbins - 1.0e0);   /* Range of values to be stored in each bin */
    float firstbin = fmin + binsize/2.0e0;                  /* Centered-value of first bin */

    /* Fill first column of 'result' with center values of each bin */
    result(0,0) = firstbin;
    for(int i=1; i<nbins; i++){
        result(i,0) = (T)i*binsize + firstbin;
    }

    /* Loop through array and count the number of elements in each bin */
    int idx = 0;
    for(int k=0; k<input.GetDim(3); k++){
        for(int i=0; i<input.GetDim(1); i++){
            for(int j=0; j<input.GetDim(2); j++){
                idx = (int)(((float)input(i,j,k) - fmin)/binsize);
                result(idx,1) += 1.0e0;
            }
        }
    }
}


template <class T>
void Stats<T>::Histogram(const Array4D<T> &input, const int nbins, Array2D<float> &result) const
{
    /* Reset 'result' array based on the specified number of bins */
    result.ResetSize(nbins,2,(float)0.0e0);

    /* Calculate the value range for each bin, assuming uniform bin size */
    size_t di,dj,dk,dl;                                     /* Dummy variables */
    float fmin = (float)input.MinVal(di,dj,dk,dl);          /* Minimum value in array */
    float fmax = (float)input.MaxVal(di,dj,dk,dl);          /* Maximum value in array */
    float binsize = (fmax - fmin)/((float)nbins - 1.0e0);   /* Range of values to be stored in each bin */
    float firstbin = fmin + binsize/2.0e0;                  /* Centered-value of first bin */

    /* Fill first column of 'result' with center values of each bin */
    result(0,0) = firstbin;
    for(int i=1; i<nbins; i++){
        result(i,0) = (T)i*binsize + firstbin;
    }

    /* Loop through array and count the number of elements in each bin */
    int idx = 0;
    for(int l=0; l<input.GetDim(4); l++){
        for(int k=0; k<input.GetDim(3); k++){
            for(int i=0; i<input.GetDim(1); i++){
                for(int j=0; j<input.GetDim(2); j++){
                    idx = (int)(((float)input(i,j,k,l) - fmin)/binsize);
                    result(idx,1) += 1.0e0;
                }
            }
        }
    }
}


//template <class T>
//void Stats<T>::OptimizeTreeTransProb_SteepestAscent(Array1D<T> &params, Tree2<T> &tree, T(*f)(Array1D<T>&, Tree2<T>&),
//                                                    const T hmin, const T tolreqd, const int maxiter, T &tolachieved,
//                                                    int &itersused, int &fevals, const T derivfrac) const
//{
//    T h = 1.0e0;
//    tolachieved = 1.0e0;
//    T tol0 = 1.0e0;
//    itersused = 0;
//    fevals = 0;
//    T tmpparam = 0.0e0;
//    T tmpval1 = 0.0e0;
//    T tmpval2 = 0.0e0;
//    T fval = f(params,tree);
//    T tmpfval = fval - 1.0e0;
//    int nparam = params.GetDim();
//    Array1D<T> dparams(nparam);
//    Array1D<T> tmpparams(nparam);

//    while(h > hmin && tolachieved > tolreqd && itersused < maxiter){
//        /* Numerically approximate derivative of parameter vector using central difference */
//        for(int i=0; i<nparam; i++){
//            tmpparam = params(i);

//            if(tmpparam > derivfrac){
//                params(i) = tmpparam*(1.0e0 + derivfrac);
//            } else {
//                params(i) = tmpparam + derivfrac;
//            }
//            tmpval1 = f(params,tree);
//            fevals++;

//            if(tmpparam > derivfrac){
//                params(i) = tmpparam*(1.0e0 - derivfrac);
//            } else {
//                params(i) = tmpparam - derivfrac;
//            }
//            tmpval2 = f(params,tree);
//            fevals++;

//            if(params(i) > derivfrac){
//                dparams(i) = (tmpval2 - tmpval1)/(2.0e0*derivfrac*tmpparam);
//            } else {
//                dparams(i) = (tmpval2 - tmpval1)/(2.0e0*derivfrac);
//            }

//            params(i) = tmpparam;
//        }

//        while(h > hmin && tmpfval < fval){
//            /* Calculate new parameter estimate */
//            for(int i=0; i<nparam; i++){
//                tmpparams(i) = params(i) - h*dparams(i);
//            }

//            /* Evaluate function with updated parameters */
//            tmpfval = f(tmpparams,tree);

//            /* Ensure ascent */
//            if(tmpfval < fval){
//                h *= 0.5e0;
//            }

//        }

//        /* Calculate tolerance measure */
//        tolachieved = 0.0e0;
//        for(int i=0; i<nparam; i++){
//            tolachieved += dparams(i)*dparams(i);
//        }
//        if(itersused == 0){
//            tol0 = tolachieved;
//        }
//        tolachieved /= tol0;

//        /* Update iteration count */
//        itersused++;
//    }
//}

