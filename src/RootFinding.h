/**
 * @file RootFinding.h
 * @author Robert Grandin
 * @date 28 November 2010
 * @brief Definition of RootFinding class.
 *
 * @section Class Description & Notes
 *
 * This class contains routines to be used for finding the zero-points of
 * functions.  All calculations are done using the precision of the class as-
 * indicated by template parameter T.  As such, unexpected behavior may result
 * if this class is used for integer calculations.
 *
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 28 November 2010
 *	- Creation date.
 * @date 14 February 2012
 *  - Added code for bisection method
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2010-2012 Robert Grandin
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

#ifndef RootFinding_
#define RootFinding_


#include <cmath>
#include <limits>

#include <Array1D.h>
#include <Array2D.h>


/**
 * @brief Class for routines used to find roots of equations.
 */
template <class T> class RootFinding{

public:

	/**
	 * @brief Constructor.
	 * @pre None.
	 * @post RootFinding object created.
	 * @return None.
	 */
	RootFinding();


	/**
	 * @brief Deconstructor.
	 * @pre RootFinding object exists.
	 * @post RootFinding object destroyed.
	 * @return None.
	 */
    virtual ~RootFinding();


	/**
	 * @brief Example of format for defining the function for which zero is to
	 * 			be found.
	 * @pre None.
	 * @param varval Array1D containing the values for all variables used in the
	 * 			function.
	 * @post None.
	 * @return Value of the function with the specified variable values.
	 */
	T f(Array1D<T> &varval);


	/**
	 * @brief Bisection method for finding roots.
	 * @pre RootFinding object exists.
     * @param varval Array1D containing the values for all variables used in the
     *  function for which zeros are found.
     * @param rvar Index of varval which contains the variable to be modified
     *  when searching for the function's root.
     * @param a Starting point of section to be searched.
     * @param b Ending point of section to be searched.
	 * @param tol Tolerance required to declare convergence.
	 * @param maxiter Maximum number of iterations allowed.
     * @param tolachieved Actual tolerance achieved.  This is the function value at
     *  the determined root location.
     * @param itersused Number of iterations used.
     * @param fevals Number of function evaluations required.
	 * @param f Function for which the root is to be found.
     * @param warn Enable/disable warning message output from function.
	 * @return Variable value for which the function is equal to zero.  Note
     * 			that this is the variable identified by rvar.
     * @warning The function must have a zero-value within the range [a,b].
     * @warning It is expected that a < b.  If this is not true, their values will
     *      be swapped to satisfy this condition.
	 */
    T Bisection(Array1D<T> &varval, const int rvar, const T a, const T b, const T tol, const int maxiter,
            T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T>&), const bool warn);



	/**
	 * @brief Newton's method for finding roots.
	 * @pre RootFinding object exists.
     * @param varval Array1D containing the values for all variables used in the
     *  function for which zeros are found.
     * @param rvar Index of varval which contains the variable to be modified
     *  when searching for the function's root.
	 * @param p Perterbation amount as a fraction of the variable's value.  A
	 * 			typical value is in the neighborhood of 0.01 (1%).
	 * @param tol Tolerance required to declare convergence.
	 * @param maxiter Maximum number of iterations allowed.
     * @param tolachieved Actual tolerance achieved.  This is the function value at
     *  the determined root location.
     * @param itersused Number of iterations used.
     * @param fevals Number of function evaluations required.
	 * @param f Function for which the root is to be found.
     * @param warn Enable/disable warning message output from function.
	 * @post Number of iterations used saved to public attribute
	 * 			"iterations_used".
	 * @return Variable value for which the function is equal to zero.  Note
	 * 			that this is the variable identified by rvar.
	 */
    T Newton(Array1D<T> &varval, const int rvar, const T p, const T tol, const int maxiter,
            T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T>&), const bool warn);


	/**
	 * @brief Secant method for finding roots.
	 * @pre RootFinding object exists.
     * @param varval Array1D containing the values for all variables used in the
     *  function for which zeros are found.
     * @param rvar Index of varval which contains the variable to be modified
     *  when searching for the function's root.
     * @param p Perterbation amount as a fraction of the variable's value.  A
     *  typical value is in the neighborhood of 0.01 (1%).
	 * @param tol Tolerance required to declare convergence.
	 * @param maxiter Maximum number of iterations allowed.
     * @param tolachieved Actual tolerance achieved.  This is the function value at
     *  the determined root location.
     * @param itersused Number of iterations used.
     * @param fevals Number of function evaluations required.
	 * @param f Function for which the root is to be found.
     * @param warn Enable/disable warning message output from function.
	 * @post Number of iterations used saved to public attribute
	 * 			"iterations_used".
	 * @return Variable value for which the function is equal to zero.  Note
	 * 			that this is the variable identified by rvar.
	 */
    T Secant(Array1D<T> &varval, const int rvar, const T p, const T tol, const int maxiter,
            T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T>&), const bool warn);


	/**
	 * @brief Fixed-Point-Iteration method for finding roots.
	 * @pre RootFinding object exists.
     * @param varval Array1D containing the values for all variables used in the
     *  function for which zeros are found.
     * @param rvar Index of varval which contains the variable to be modified
     *  when searching for the function's root.
     * @param p Perterbation amount as a fraction of the variable's value.  A
     *  typical value is in the neighborhood of 0.01 (1%).
	 * @param tol Tolerance required to declare convergence.
	 * @param maxiter Maximum number of iterations allowed.
     * @param tolachieved Actual tolerance achieved.  This is the function value at
     *  the determined root location.
     * @param itersused Number of iterations used.
     * @param fevals Number of function evaluations required.
	 * @param f Function for which the root is to be found.
     * @param warn Enable/disable warning message output from function.
	 * @post Number of iterations used saved to public attribute
	 * 			"iterations_used".
	 * @return Variable value for which the function is equal to zero.  Note
	 * 			that this is the variable identified by rvar.
	 */
    T FixedPointIteration(Array1D<T> &varval, const int rvar, const T p, const T tol, const int maxiter,
            T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T>&), const bool warn);


	/**
	 * @brief Inverse-Parabolic method for finding roots.
	 * @pre RootFinding object exists.
     * @param varval Array1D containing the values for all variables used in the
     *  function for which zeros are found.
     * @param rvar Index of varval which contains the variable to be modified
     *  when searching for the function's root.
     * @param p Perterbation amount as a fraction of the variable's value.  A
     *  typical value is in the neighborhood of 0.01 (1%).
	 * @param tol Tolerance required to declare convergence.
	 * @param maxiter Maximum number of iterations allowed.
     * @param tolachieved Actual tolerance achieved.  This is the function value at
     *  the determined root location.
     * @param itersused Number of iterations used.
     * @param fevals Number of function evaluations required.
	 * @param f Function for which the root is to be found.
     * @param warn Enable/disable warning message output from function.
	 * @post Number of iterations used saved to public attribute
	 * 			"iterations_used".
	 * @return Variable value for which the function is equal to zero.  Note
	 * 			that this is the variable identified by rvar.
	 */
    T InverseParabolic(Array1D<T> &varval, const int rvar, const T p, const T tol, const int maxiter,
            T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T>&), const bool warn);


	/**
	 * @brief Muller's method for finding roots.
	 * @pre RootFinding object exists.
     * @param varval Array1D containing the values for all variables used in the
     *  function for which zeros are found.
     * @param rvar Index of varval which contains the variable to be modified
     *  when searching for the function's root.
     * @param p Perterbation amount as a fraction of the variable's value.  A
     *  typical value is in the neighborhood of 0.01 (1%).
	 * @param tol Tolerance required to declare convergence.
	 * @param maxiter Maximum number of iterations allowed.
     * @param tolachieved Actual tolerance achieved.  This is the function value at
     *  the determined root location.
     * @param itersused Number of iterations used.
     * @param fevals Number of function evaluations required.
	 * @param f Function for which the root is to be found.
     * @param warn Enable/disable warning message output from function.
	 * @post Number of iterations used saved to public attribute
	 * 			"iterations_used".
	 * @return Variable value for which the function is equal to zero.  Note
	 * 			that this is the variable identified by rvar.
	 */
    T Muller(Array1D<T> &varval, const int rvar, const T tol, const T p, const int maxiter,
            T &tolachieved, int &itersused, int &fevals, T(*f)(Array1D<T>&), const bool warn);



private:
    /** @brief Track if warnings are to be output (true) or suppressed (false). */
    bool warnenabled;


};
 
#include "RootFinding.cpp"

#endif /* RootFinding_ */
