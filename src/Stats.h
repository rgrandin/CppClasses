/**
 * @file Stats.h
 * @author Robert Grandin
 * @date 7 Feb 2012
 * @brief Definition of Stats class.
 *
 * @section Class Description & Notes
 *
 * This class contains statistical functions developed during my Statistics 580
 * coursework.  Many of these functions were provided to the class by the instructor,
 * Dr. Karin Dorman, in C-form.  They have been reimplemented into C++ by me.
 * For these functions, their original author will be noted in the function
 * documentation.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 *
 * @section Revisions
 *
 * @date 7 February 2012
 *	- Creation date.
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2012, Robert Grandin
 * All rights reserved.
 *
 * Redistribution and use of this file is permitted provided that the following
 * conditions are met:
 * 	-# 	Redistributions must produce the above copyright notice, this list of
 * 		conditions, and the following disclaimer in the documentation and/or
 * 		other materials provided with the distribution.
 * 	-#	Neither the name of the organization nor the names of its contributors
 * 		may be used to endorse or promote products derived from this software
 * 		without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
 * EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING BUT NOT
 * LIMITING TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 */




#ifndef Stats_
#define Stats_

#include <math.h>
#include <limits>

#include <Array1D.h>
#include <Array2D.h>
#include <Array3D.h>
#include <Array4D.h>
#include <Tree2.h>

/**
  @brief Class to contain statistical functions.
  @warning C++11 features, such as move-constructor and move-assignment, require the symbol
   "CXX11" to be defined.
  */
template <class T>
class Stats
{
public:
    /**
      @brief Constructor.  This class has no member variables to initialize.
      @pre None.
      @post Stats object created.  Member functions can now be accessed.
      @return None.
      */
    Stats();


    /**
     * @brief Copy constructor.
     * @param a Reference to existing Stats object to be copied.
     */
    Stats(Stats<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing Stats object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    Stats(Stats<T> &&a);
#endif


    /**
      @brief Destructor.
      @pre Stats object exists.
      @post Stats object destroyed.
      @return None.
      */
    virtual ~Stats();


    /**
     * @brief Copy-assignment operator.
     * @param a Stats object being assigned.
     * @return Reference to instance of Stats.
     */
    Stats& operator=(Stats<T> a);


#ifdef CXX11
    /**
     * @brief Move-assignment operator (C++11).
     * @param a Rvalue to Stats object being assigned.
     * @return Reference to instance of Stats.
     * @warning This function requires C++11 compiler support.
     */
    Stats& operator=(Stats<T> &&a);
#endif


    /**
      @brief Compute a single binomial coefficient @f$ n \choose k @f$.

      This is
        makes use of BinomialCoeffTriangleRow() and allows for calculating a single
        coefficient without the user needing to create a temporary array (as required
        in the calculation of a full row of Pascal's Triangle).
      @pre Stats object exists.
      @param n Order of polynomial: @f$ (1 + x)^n @f$.
      @param k @f$ x^k @f$ term for which coefficient is desired.
      @post Coefficient calculated.
      @return Calculated coefficient of a type matching the type of the Stats object.
      */
    T BinomialCoeffSingle(const int n, const int k) const;


    /**
      @brief Compute the entire row of Pascal's Triangle of binomial coefficients for
        a specified polynomial power, n.

      This is reimplemented from Dr. Dorman's
        'pascal_triangle()' function provided as part of homework 1 in Stat 580.
      @pre Stats object exists.
      @param n Order of polynomial: @f$ (1 + x)^n @f$.
      @param array Array1D in which triangle-row values are to be stored.  Array will be
        properly sized within this function.  Entry k will contain the
        coefficient on the @f$ x^k @f$ term in the polynomial.
      @post Full triangle-row calculated and stored in array.
      @return None.  Coefficients stored in array.
      */
    void BinomialCoeffTriangleRow(const int n, Array1D<T> &array) const;


    /**
      @brief Compute the entire Pascal's Triangle of binomial coefficients through
        a specified polynomial power, n.

      This is an extension of Dr. Dorman's
        'pascal_triangle()' function provided as part of homework 1 in Stat 580.
      @pre Stats object exists.
      @param n Order of polynomial: @f$ (1 + x)^n @f$.
      @param array Array2D in which triangle values are to be stored.  Array will be
        properly sized within this function.  Rows can be thought of as corresponding to
        n.  Columns can be thought of as corresponding to k.  Column k will contain the
        coefficient on the @f$ x^k @f$ term in the polynomial.
      @post Full triangle, up to order n, calculated and stored in array.  Elements which
        are not part of the triangle are set to 0.
      @return None.  Coefficients stored in array.
      */
    void BinomialCoeffFullTriangle(const int n, Array2D<T> &array) const;


    /**
      @brief Evaluate a polynomial with user-specified coefficients.
      @pre Stats object exists.
      @param array Array1D of n+1 coefficients.  The 0-based index of the coefficient within the
        array identifies the power of x with which it is associated (i.e., element 0 is
        the constant term and element n is associated with @f$ x^n @f$).
      @param x Point at which polynomial is to be evaluated.
      @post Polynomial value is calculated.
      @return Value of polynomial at specified point.
     */
    T PolyVal(const Array1D<T> &array, const T x) const;


    /**
      @brief Evaluate the derivative of a polynomial with user-specified coefficients.
      @pre Stats object exists.
      @param array Array1D object containing the n+1 coefficients for the original polynomial.  The
        0-based index of the coefficient within the array identifies the power of x with which it
        is associated (i.e., element 0 is the constant term in the original polynomial and element
        n is associated with @f$ x^n @f$ in the original polynomial).
      @param m Derivative to be calculated.
      @param x Point at which derivative is to be evaluated.
      @post Derivative value calculated.
      @return Value of the derivative at the specified point.
     */
    T PolyDerivVal(const Array1D<T> &array, const int m, const T x);


    /**
      @brief Evaluate the Hermite polynomial @f$ H_k(x) @f$.

      This is based on State 580 course
        notes provided by Dr. Dorman and the leading coefficient is <b>not</b> scaled to 1.
      @pre Stats object exists.
      @param k Order of polynomial.
      @param x Point at which polynomial is evaluated.
      @post Polynomial evaluated.
      @return Value of @f$ H_k(x) @f$.
     */
    T HermiteEval(const int k, const T x) const;


    /**
      @brief Evaluate probablists' Hermite polynomial @f$ He_k(x) @f$.

      The leading coefficient is
        scaled to 1.  Algorithm is based on material found on Wikipedia.
       @pre Stats object exists.
       @param k Order of polynomial.
       @param x Point at which polynomial is evaluted.
       @post Polynomial evaluated.
       @return Value of @f$ He_k(x) @f$.
      */
    T HermiteProbEval(const int k, const T x) const;


    /**
      @brief Calculate the logarithm of the Gamma function, @f$ \Gamma(x) @f$.

      Algorithm based on
        Numerical Recipies in C (found online at http://apps.nrbook.com/c/index.html page 214).
        Calculations are done in double-precision, and the returned value is of template parameter T.
      @pre Stats object exists.
      @param x Point at which @f$ ln \left(\Gamma(x) \right) @f$ is to be calculated.  This must be greater
        than 0.
      @post Value calculated.
      @return @f$ ln \left(\Gamma(x)\right) @f$.  Returns NaN if x <= 0.
     */
    T LnGamma(const T x) const;


    /**
      @brief Approximate the Gamma function, @f$ \Gamma(x) @f$, using the approximation provided in
        homework 2 for Stat 580.
      @pre Stats object exists.
      @param x Point at which @f$ ln \left(\Gamma(x) \right) @f$ is to be calculated.  This must be greater
        than 1.
      @post No changes to object.
      @return @f$ ln \left(\Gamma(x)\right) @f$
     */
    T GammaHW(const T x) const;


    /**
      @brief Compute moments for a sum of independent and identically distributed random variables.
      @pre Stats object exists.
      @param n Number of iid random variables to be summed.
      @param k Number of moments to be calculated.
      @param tol Tolerance required for Taylor expansion of moment-generating function.  This controls
        the number of terms used in the polynomial expansion.
      @param params Array1D object containing the parameters for the distribution.
      @param w_array Array1D object containing the moments for the sum of iid random variables.
      @param f Function which computes the mth moment for the desired random variable.  This function
        has two arguments: an Array1D object containing the parameters of the distribution and an
        integer specifying which coefficient in the polynomial expression of the MGF should be returned.
        A general polynomial form is used to allow greater flexibility in the MGF calculations as
        compared to if a Taylor series was assumed.
      @param termsreqd Number of terms used when calculating the MGF.  This value is determined by
        ComputeMoments() and returned via this referenced variable.
      @post Moments for sum of iid random variables computed and stored in w_array.
      @return None.
     */
    void ComputeMoments(const int n, const int k, const T tol, const Array1D<T> &params,
                        Array1D<T> &w_array, T (*f)(const Array1D<T>&, const int), int &termsreqd) const;


    /**
      @brief Compute moments for a sum of independent and identically distributed random variables using
        a provided set of moments for a single random variable of the desired distribution.
      @pre Stats object exists.
      @param n Number of iid random variables to be summed.
      @param k Number of moments to be calculated.
      @param m_array Array1D object containing the moments for a single random variable which has
        the desired distribution.
      @param w_array Array1D object containing the moments for the sum of iid random variables.
      @post Moments for sum of iid random variables computed and stored in w_array.
      @return None.
     */
    void ComputeMoments(const int n, const int k, const Array1D<T> &m_array, Array1D<T> &w_array) const;


    /**
      @brief Convert between moments and cumulants.
      @pre Stats object exists.
      @param input Array1D of input values.
      @param output Array1D of output values.
      @param flag Identifies if input is comprised of moment or cumulant values.
        - Positive or 0: input values are moments
        - Negative: input values are cumulants
      @post Conversion performed and result stored in output.  If input is moments, output is cumulants.
        If input is cumulants, output is moments.
      @return None.
      @warning Values for moments do not include the factorial term from the Taylor series expansion of the MGF.
        That is, the moment values here (both input and output) are @f$ m_k@f$ in the equation
        @f$ M(s) = \sum_{k=1}^\infty \frac{m_k}{k!} s^k @f$.
     */
    void ConvertMomentsCumulants(const Array1D<T> &input, Array1D<T> &output, const int flag) const;


    /**
      @brief Evaluate the Gamma distribution using a sum as described in Stat 580 course notes.

      It is defined as
        @f$ P(a,b) = \frac{1}{\Gamma(a)} \int_0^b t^{a-1} e^{-t} dt @f$.  Using this definition, @f$ a = k @f$
        and @f$ b = \frac{x}{\theta} @f$.  This is the cumultative distribution function, CDF, of
        the gamma distribution, and the integral term is known as the lower incomplete Gamma function.
      @pre Stats object exists.
      @param k Shape parameter.
      @param theta Scale parameter.
      @param x Point at which CDF value is to be calculated.
      @param tol Desired tolerance required for convergence.
      @param maxiter Maximum number of iterations allowed for calculations.
      @param tolachieved Tolerance actually achieved.
      @param itersused Iterations used for computation.
      @post No changes to object.
      @return Value of integral with specified shape and scale parameters.
     */
    T GammaCDF(const T k, const T theta, const T x, const T tol, const int maxiter,
               T &tolachieved, int &itersused) const;


    /**
      @brief Evaluate the Gamma distribution using continued fraction expansion.

      It is defined as
        @f$ P(a,b) = \frac{1}{\Gamma(a)} \int_0^b t^{a-1} e^{-t} dt @f$.  Using this definition, @f$ a = k @f$
        and @f$ b = \frac{x}{\theta} @f$.  This is the cumultative distribution function, CDF, of
        the gamma distribution, and the integral term is known as the lower incomplete Gamma function.  This
        implementation converges faster than GammaCDF(), which calculates the same quantity using a summation.
      @pre Stats object exists.
      @param k Shape parameter.
      @param theta Scale parameter.
      @param x Point at which CDF value is to be calculated.
      @param tol Desired tolerance required for convergence.
      @param maxiter Maximum number of iterations allowed for calculations.
      @param tolachieved Tolerance actually achieved.
      @param itersused Iterations used for computation.
      @post No changes to object.
      @return Value of integral with specified shape and scale parameters.
     */
    T GammaCDF_CFE(const T k, const T theta, const T x, const T tol, const int maxiter,
                   T &tolachieved, int &itersused) const;


    /**
      @brief Evaluate the Gamma probability density function, defined as
        @f$ p(x;k,\theta) = \frac{x^{k-1} e^{-\frac{x}{\theta}}}{\theta^k \Gamma(k)} @f$.
      @pre Stats object exists.
      @param k Shape parameter.
      @param theta Scale parameter.
      @param x Point at which density is to be calculated.
      @post No changes to object.
      @return Value of probability density function.
     */
    T GammaPDF(const T k, const T theta, const T x) const;


    /**
      @brief Evaluate the @f$ \chi^2 @f$ distribution CDF.

      This makes use of the relation
        @f$ F_n(x) = P(\frac{n}{2},\frac{x}{2}) @f$ where F is the distribution function of the
        @f$ \chi^2 @f$ distribution and P is the distribution function, GammaCDF().
      @pre Stats object exists.
      @param n Degrees of freedom.
      @param twolambda Noncentrality parameter.  @f$ \lambda = 0 @f$ is for central @f$ \chi^2 @f$ distribution.
      @param x Value for which CDF is desired.
      @param tol Desired accuracy for calculation.
      @param maxiter Maximum number of iterations to use in calculation.
      @param tolachieved Actual accuracy achieved in calculation.
      @param itersused Actual number of iterations used in calculation.
      @post No changes to object.
      @return Value of @f$ \chi^2 @f$ distribution CDF at x.
     */
    T ChiSquaredCDF(const T n, const T twolambda, const T x, const T tol, const int maxiter,
                    T &tolachieved, int &itersused) const;


    /**
      @brief Evaluate the @f$ \chi^2 @f$ distribution PDF.
      @pre Stats object exists.
      @param n Degrees of freedom.
      @param twolambda Noncentrality parameter.  @f$ \lambda = 0 @f$ is for central @f$ \chi^2 @f$ distribution.
      @param x Point at which density is to be calculated.
      @param tol Desired accuracy for calculation.
      @param maxiter Maximum number of iterations to use in calculation.
      @param tolachieved Actual accuracy achieved in calculation.
      @param itersused Actual number of iterations used in calculation.
      @post No changes to object.
      @return Value of @f$ \chi^2 @f$ distribution PDF at x.
     */
    T ChiSquaredPDF(const T n, const T twolambda, const T x, const T tol, const int maxiter,
                    T &tolachieved, int &itersused) const;


    /**
      @brief Evaluate the @f$ \beta @f$ distribution PDF.
      @pre Stats object exists.
      @param x Point at which density is to be calculated.
      @param alpha First shape parameter.
      @param beta Second shape parameter.
      @post No changes to object.
      @return Value of @f$ \beta @f$ distribution PDF at x.
    */
    T BetaPDF(const T x, const T alpha, const T beta) const;


    /**
      @brief Evaluate the @f$ \beta @f$ distribution CDF.
      @pre Stats object exists.
      @param x Point at which density is to be calculated.
      @param alpha First shape parameter.
      @param beta Second shape parameter.
      @param tol Desired accuracy for calculation.
      @param maxiter Maximum number of iterations to use in calculation.
      @param tolachieved Actual accuracy achieved in calculation.
      @param itersused Actual number of iterations used in calculation.
      @post No changes to object.
      @return Value of @f$ \beta @f$ distribution CDF at x.
    */
    T BetaCDF(const T x, const T alpha, const T beta, const T tol, const int maxiter,
              T &tolachieved, int &itersused) const;


    /**
      @brief Evaluate the @f$ \beta @f$ distribution CDF.
      @pre Stats object exists.
      @param x Point at which density is to be calculated.
      @param alpha First shape parameter.
      @param beta Second shape parameter.
      @post No changes to object.
      @return Value of @f$ \beta @f$ distribution CDF at x.
    */
    T BetaCDF(const T x, const T alpha, const T beta) const;


    /**
      @brief Determine the bounds of the minimum-length interval of the @f$ \beta @f$ distribution for
        a user-specified confidence value.
      @pre Stats object exists.
      @param alpha First shape parameter of the distribution.
      @param beta Second shape parameter of the distribution.
      @param conf Confidence level.  A confidence level of 95 percent should have 'conf = 0.05'.
      @param tol Tolerance required for convergence.
      @param tolachieved Actual tolerance achieved with respect to the determination of the interval size.
      @param xleft Left-bound of interval.
      @param xright Right-bound of interval.
      @param fevals_cdf Number of times the cumulative distribution function was evaluated.
      @param fevals_pdf Number of times the probability density function was evaluated.
      @post No changes to object.
      @return None.
    */
    void BetaCI(const T alpha, const T beta, const T conf, const T tol, T &tolachieved, T &xleft, T &xright,
                int &fevals_cdf, int &fevals_pdf);


    /**
      @brief Calculate histogram of Array1D object and place result into 2D array.

      First column
        of result array contains the center-value of each bin.  Second column contains the
        number of array elements with values in each bin.  Each bin is assumed to be of
        uniform size.
      @pre Stats object exists.
      @param input Reference to 1D array containing input data.
      @param nbins Number of bins to use.
      @param result Reference to 2D array containing histogram data.  This array is of type
        'float' in order to handle non-integer bin-center values (in the case of integer arrays)
        and large bin counts ('float' will allow ~1e30).
      @post No changes to object.
      @return None.
     */
    void Histogram(const Array1D<T> &input, const int nbins, Array2D<float> &result) const;


    /**
      @brief Calculate histogram of Array3D object and place result into 2D array.

      First column
        of result array contains the center-value of each bin.  Second column contains the
        number of array elements with values in each bin.  Each bin is assumed to be of
        uniform size.
      @pre Stats object exists.
      @param input Reference to 2D array containing input data.
      @param nbins Number of bins to use.
      @param result Reference to 2D array containing histogram data.  This array is of type
        'float' in order to handle non-integer bin-center values (in the case of integer arrays)
        and large bin counts ('float' will allow ~1e30).
      @post No changes to object.
      @return None.
     */
    void Histogram(const Array2D<T> &input, const int nbins, Array2D<float> &result) const;


    /**
      @brief Calculate histogram of Array3D object and place result into 2D array.

      First column
        of result array contains the center-value of each bin.  Second column contains the
        number of array elements with values in each bin.  Each bin is assumed to be of
        uniform size.
      @pre Stats object exists.
      @param input Reference to 3D array containing input data.
      @param nbins Number of bins to use.
      @param result Reference to 2D array containing histogram data.  This array is of type
        'float' in order to handle non-integer bin-center values (in the case of integer arrays)
        and large bin counts ('float' will allow ~1e30).
      @post No changes to object.
      @return None.
     */
    void Histogram(const Array3D<T> &input, const int nbins, Array2D<float> &result) const;


    /**
      @brief Calculate histogram of Array4D object and place result into 2D array.

      First column
        of result array contains the center-value of each bin.  Second column contains the
        number of array elements with values in each bin.  Each bin is assumed to be of
        uniform size.
      @pre Stats object exists.
      @param input Reference to 4D array containing input data.
      @param nbins Number of bins to use.
      @param result Reference to 2D array containing histogram data.  This array is of type
        'float' in order to handle non-integer bin-center values (in the case of integer arrays)
        and large bin counts ('float' will allow ~1e30).
      @post No changes to object.
      @return None.
     */
    void Histogram(const Array4D<T> &input, const int nbins, Array2D<float> &result) const;



protected:
    /** @brief Mathematical constant PI.  The value is hard-coded to avoid calculation
      round-off errors. */
    T PI_Stats;


    /**
      @brief Evaluate the continued fraction for the @f$ \beta @f$ distribution.
      @pre Stats object exists.
      @param x Point at which density is to be calculated.
      @param alpha First shape parameter.
      @param beta Second shape parameter.
      @param tol Desired accuracy for calculation.
      @param maxiter Maximum number of iterations to use in calculation.
      @param tolachieved Actual accuracy achieved in calculation.
      @param itersused Actual number of iterations used in calculation.
      @post No changes to object.
      @return None.
    */
    T BetaCFE(const T x, const T alpha, const T beta, const T tol, const int maxiter,
              T &tolachieved, int &itersused) const;


    /**
      @brief Function used to determine value of x corresponding to a specified PDF value.
      @pre Stats object exists.
      @param x Point at which density is to be calculated.
      @param alpha First shape parameter.
      @param beta Second shape parameter.
      @param lambda Desired PDF value.
      @post No changes to object.
      @return Value of @f$ \beta @f$ distribution PDF at x subtracted from the desired value, 'lambda'.
    */
    T BetaPDFRoot(const T x, const T alpha, const T beta, const T lambda) const;


    /**
      @brief Function used to determine the minimum-length confidence interval for a beta distribution.
      @pre Stats object exists.
      @param alpha First shape parameter.
      @param beta Second shape parameter.
      @param lambda Desired value of PDF at interval endpoints.
      @param conf Desired size of confidence interval.
      @param xleft Current left-bound of interval.
      @param xright Current right-bound of interval.
      @param tol Tolerance required for convergence.
      @param maxiter Maximum number of allowed iterations.
      @param tolachieved Actual tolerance achieved with respect to the determination of the
        xleft and xright locations for the specified 'lambda'.
      @param itersused Actual number of iterations used when finding the xleft and xright locations.
      @param fevals_cdf Number of times the cumulative distribution function was evaluated.
      @param fevals_pdf Number of times the probability density function was evaluated.
      @post No changes to object.
      @return Confidence interval size subtracted from the CDF evaluated at each endpoint.
    */
    T BetaCIRoot(const T alpha, const T beta, const T lambda, const T conf, T &xleft, T &xright,
                 const T tol, const int maxiter, T &tolachieved, int &itersused,
                 int &fevals_cdf, int &fevals_pdf) const;


private:
    /**
     * @brief StatsSwap swaps member information between two Stats objects.
     * @param first First Stats object.
     * @param second Second Stats object.
     */
    friend void StatsSwap(Stats<T> &first, Stats<T> &second)
    {
        std::swap(first.PI_Stats, second.PI_Stats);
    }



};

#include "Stats.cpp"

#endif /* Stats_ */
