/**
 * @file ODE.h
 * @author Robert Grandin
 * @date 27 November 2010
 * @brief Definition of ODE class.
 *
 * @section Class Description & Notes
 *
 * This class contains routines for solving ordinary differential equations.
 *
 * The derivative function required by each of these routines is expected to
 * have the same form as function "f" within this class.  The function "f"
 * within this class is only intended to demonstrate the necessary calling
 * structure and thus contains no working code.
 *
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 27 November 2010
 *	- Creation date.
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2010, Robert Grandin
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

#ifndef ODE_
#define ODE_

#include <cmath>
#include <Array1D.h>
#include <Array2D.h>
#include <Grid2D.h>
#include <MPI_Custom.h>


/**
 * @brief Class for routines used to solve ordinary differential equations.
 * @warning C++11 features, such as move-constructor and move-assignment, require the symbol
 *  "CXX11" to be defined.
 */
template <class T> class ODE{

public:
	/**
	 * @brief Constructor.
	 * @pre None.
	 * @post Object created.
	 * @return None.
	 */
	ODE();


    /**
     * @brief Copy constructor.
     * @param a Reference to existing ODE object to be copied.
     */
    ODE(ODE<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing ODE object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    ODE(ODE<T> &&a);
#endif

	/**
	 * @brief Deconstructor.
	 * @pre ODE object exists.
	 * @post Object destroyed.
	 * @return None.
	 */
    virtual ~ODE();


    /**
     * @brief Copy-assignment operator.
     * @param a Reference to ODE object being assigned.
     * @return Reference to instance of ODE.
     */
    ODE& operator=(ODE<T> a);


#ifdef CXX11
    /**
     * @brief Move-assignment operator (C++11).
     * @param a Rvalue to ODE object being assigned.
     * @return Reference to instance of ODE.
     * @warning This function requires C++11 compiler support.
     */
    ODE& operator=(ODE<T> &&a);
#endif


	/**
	 * @brief Sample function showing the calling format for the function which
	 * 			evaluates the derivatives while solving an ordinary differential
     * 			equation.
     *
     * Sample calculation within the derivative function:
     * 	"yout(0) = -0.2*yin(0)".  This specific function ("f") has no
     * 	working code and is only intended to demonstrate the proper form.
	 * @pre None.
	 * @param yin Array1D object of the input state vector.
	 * @param yout Array1D object of the output derivatives which correspond to
	 * 				the elements of the input state vector.
	 * @param x Value of the independent variable (i.e., 't').
	 * @post Derivatives calculated and saved in yout.
	 * @return None.
	 */
	int f(Array1D<T> &yin,	Array1D<T> &yout,T x);


	/**
	 * @brief Adams-Bashforth 2nd Order.
	 * @pre ODE object exists and sufficient memory is available for the storage
	 * 		of temporary copies of the state vector.  This routine requires one
	 * 		temporary copy of the state vector.
	 * @param h Step size.
	 * @param x Current value of independent variable (i.e., 't').
	 * @param y State vector (Array1D object).  This contains the post-step
	 * 			values of the state vector.
	 * @param derivs Derivatives at previous time-steps (Array2D object).  Array
	 * 			is indexed (nelements x nderivs).
	 * @param f Pointer to function defining the derivatives.
	 * @post State vector, y, is updated and contains the new values.
	 * @return None.
	 */
	void AB2(T h, T x, Array1D<T> &y, Array2D<T> &derivs,
			int(*f)(Array1D<T>&, Array1D<T>&, Array1D<T>&));


	/**
	 * @brief Adams-Bashforth 4th Order.
	 * @pre ODE object exists and sufficient memory is available for the storage
	 * 		of temporary copies of the state vector.  This routine requires one
	 * 		temporary copy of the state vector.
	 * @param h Step size.
	 * @param x Current value of independent variable (i.e., 't').
	 * @param y State vector (Array1D object).  This contains the post-step
	 * 			values of the state vector.
	 * @param derivs Derivatives at previous time-steps (Array2D object).  Array
	 * 			is indexed (nelements x nderivs).
	 * @param f Pointer to function defining the derivatives.
	 * @post State vector, y, is updated and contains the new values.
	 * @return None.
	 */
	void AB4(T h, T x, Array1D<T> &y, Array2D<T> &derivs,
			int(*f)(Array1D<T>&, Array1D<T>&, Array1D<T>&));


	/**
	 * @brief Explicit Euler method (also known as Runge-Kutta 1st Order).
	 * @pre ODE object exists and sufficient memory is available for the storage
	 * 		of temporary copies of the state vector.  This routine requires one
	 * 		temporary copy of the state vector.
	 * @param h Step size.
	 * @param x Current value of independent variable (i.e., 't').
	 * @param y State vector (Array1D object).  This contains the post-step
	 * 			values of the state vector.
	 * @param f Pointer to function defining the derivatives.
	 * @post State vector, y, is updated and contains the new values.
	 * @return None.
	 */
	void Euler(T h, T x, Array1D<T> &y,
			int(*f)(Array1D<T>&, Array1D<T>&, T));


    /**
	 * @brief Runge-Kutta 2nd Order.
	 * @pre ODE object exists and sufficient memory is available for the storage
	 * 		of temporary copies of the state vector.  This routine requires two
	 * 		temporary copies of the state vector.
	 * @param h Step Size.
	 * @param x Current value of independent variable (i.e., 't').
	 * @param y State vector (Array1D object).  This contains the post-step
	 * 			values of the state vector.
	 * @param f Pointer to function defining the derivatives.
	 * @post State vector, y, is updated and contains the new values.
	 * @return None.
	 */
    void RK2(T h, T x, Array1D<T> &y, int(*f)(Array1D<T>&, Array1D<T>&, T));


    /**
     * @brief Runge-Kutta 2nd Order for use with my CFD2 functions.
     * @pre ODE object exists and sufficient memory is available for the storage
     * 		of temporary copies of the state vector.  This routine requires two
     * 		temporary copies of the state vector.
     * @param h Step Size.
     * @param x Current value of independent variable (i.e., 't').
     * @param grid Grid2D object containing the current quantity values at each point.
     * @param f Pointer to function defining the derivatives.
     * @post Grid2D object, grid, is updated and contains the new values.
     * @return Magnitude of the derivative vector after the 2nd step.
     */
    T RK2(T h, T x, Grid2D<T> &grid, void(*f)(Grid2D<T>&, Grid2D<T>&, T));


    /**
     * @brief Runge-Kutta 2nd Order for use with my CFD2 functions, modified to
     *  allow passing of extra needed parameters and work with the viscous CFD2 functions.
     * @pre ODE object exists and sufficient memory is available for the storage
     * 		of temporary copies of the state vector.  This routine requires two
     * 		temporary copies of the state vector.
     * @param h Step Size.
     * @param x Current value of independent variable (i.e., 't').
     * @param grid Grid2D object containing the current quantity values at each point.
     * @param f Pointer to function defining the derivatives.
     * @param dependency Array2D containing MPI-dependency information.
     * @param nrecv Number of MPI RECV() calls required by local MPI process.
     * @param mu Dynamic viscosity.
     * @param kthermal Thermal conductivity.
     * @post Grid2D object, grid, is updated and contains the new values.
     * @return Magnitude of the derivative vector after the 2nd step.
     */
    T RK2(T h, T x, Grid2D<T> &grid, void(*f)(Grid2D<T>&, Array3D<T>&, const T,
                                              Array2D<int>&, const int,
                                              const T, const T),
          Array2D<int> &dependency, const int nrecv,
          const T mu, const T kthermal);


	/**
	 * @brief Runge-Kutta 4th Order.
	 * @pre ODE object exists and sufficient memory is available for the storage
	 * 		of temporary copies of the state vector.  This routine requires six
	 * 		temporary copies of the state vector.
	 * @param h Step size.
	 * @param x Current value of independent variable (i.e., 't').
	 * @param y State vector (Array1D object).  This contains the post-step
	 * 			values of the state vector.
	 * @param *f Pointer to function defining the derivatives.
	 * @post State vector, y, is updated and contains the new values.
	 * @return None.
	 */
	void RK4(T h, T x, Array1D<T> &y,
			int(*f)(Array1D<T>&, Array1D<T>&, T));


	/**
	 * @brief Runge-Kutta 4th/5th Order.
	 * @pre ODE object exists and sufficient memory is available for the storage
	 * 		of temporary copies of the state vector.  This routine requires eight
	 * 		temporary copies of the state vector.
	 * @param h Step size.
	 * @param x Current value of independent variable (i.e., 't').
	 * @param y State vector (Array1D object).  This contains the 5th Order
	 * 			post-step values of the state vector.
	 * @param y4 State vector (Array1D object).  This contains the 4th Order
	 * 			post-step values of the state vector.
	 * @param f Pointer to function defining the derivatives.
	 * @post State vector, y, is updated and contains the new values.
	 * @return None.
	 */
	void RKF45(T h, T x, Array1D<T> &y, Array1D<T> &y4,
			int(*f)(Array1D<T>&, Array1D<T>&, T));


	/**
	 * @brief Lax-Wendroff method.
	 * @pre ODE object exists and sufficient memory is available for the storage
	 * 		of temporary copies of the state vector.  This routine requires three
	 * 		temporary copies of the state vector.
	 * @param h Step size.
	 * @param x Current value of independent variable (i.e., 't').
	 * @param y State vector (Array1D object).  This contains the post-step
	 * 			values of the state vector.
	 * @param c C.F.L. number.
	 * @param bc Boundary condition to be used.
	 * 			- 1:Fixed
	 * 			- 2:Periodic (for 1,2,...,n points, the n+1 point is
	 *			  mapped to 1 and point 0 is mapped to n)
	 * @param var Variable number which is being integrated.  This is needed
	 * 		calling the flux function, g, passed to this routine.
	 * @param flux Reference to the array containing the flux values.
	 * @post State vector, y, is updated and contains the new values.
	 * @return None.
	 */
	void LaxWendroff(T h, T x, Array1D<T> &y, T c, int bc, int var,
			Array2D<T> &flux);


	/**
	 * @brief MacCormack method.
	 * @pre ODE object exists and sufficient memory is available for the storage
	 * 		of temporary copies of the state vector.  This routine requires two
	 * 		temporary copies of the state vector.
	 * @param h Step size.
	 * @param x Current value of independent variable (i.e., 't').
	 * @param y State vector (Array1D object).  This contains the post-step
	 * 			values of the state vector.
	 * @param c C.F.L. number.
	 * @param bc Boundary condition to be used.
	 * 			- 1:Fixed
	 * 			- 2:Periodic (for 1,2,...,n points, the n+1 point is
	 *			  mapped to 1 and point 0 is mapped to n)
	 * @param var Variable number which is being integrated.  This is needed
	 * 		calling the flux function, g, passed to this routine.
	 * @param flux Reference to the array containing the flux values.
	 * @param disp Value for numerical dissipation to be used.
	 * @post State vector, y, is updated and contains the new values.
	 * @return None.
	 */
	void MacCormack(T h, T x, Array1D<T> &y, T c, int bc, int var,
			Array2D<T> &flux, T disp);



protected:



private:

	/** @brief Track if the Runge-Kutta 4th Order Butcher Table values have been
	 * 			set */
	bool rk4butcherset;

	/** @brief Track if the Runge-Kutta-Fehlberg 4th/5th Order Butcher Table
	 * values have been	set */
	bool rkf45butcherset;

	/** @brief Alpha values in Butcher Table for RK4 routine */
	Array1D<T> alpha_rk4;

	/** @brief Beta values in Butcher Table for RK4 routine */
	Array2D<T> beta_rk4;

	/** @brief Gamma values in Butcher Table for RK4 routine */
	Array1D<T> gamma_rk4;

	/** @brief Alpha values in Butcher Table for RKF45 routine */
	Array1D<T> alpha_rkf45;

	/** @brief Beta values in Butcher Table for RKF45 routine */
	Array2D<T> beta_rkf45;

	/** @brief Gamma values in Butcher Table for RKF45 routine */
	Array1D<T> gamma4_rkf45;

	/** @brief Gamma values in Butcher Table for RKF45 routine */
	Array1D<T> gamma5_rkf45;

	/**
	 * @brief Assign weight values to Butcher Table variables for RK4 routine.
	 * @pre ODE object exists.
	 * @post Butcher Table values defined.
	 * @return None.
	 */
	void AssignRK4Weights();

	/**
	 * @brief Assign weight values to Butcher Table variables for RKF45 routine.
	 * @pre ODE object exists.
	 * @post Butcher Table values defined.
	 * @return None.
	 */
	void AssignRKF45Weights();



    /**
     * @brief ODESwap swaps member information between two ODE objects.
     * @param first First ODE object.
     * @param second Second ODE object.
     */
    friend void ODESwap(ODE<T> &first, ODE<T> &second)
    {
        std::swap(first.rk4butcherset, second.rk4butcherset);
        std::swap(first.rkf45butcherset, second.rkf45butcherset);
        std::swap(first.alpha_rk4, second.alpha_rk4);
        std::swap(first.beta_rk4, second.beta_rk4);
        std::swap(first.gamma_rk4, second.gamma_rk4);
        std::swap(first.alpha_rkf45, second.alpha_rkf45);
        std::swap(first.beta_rkf45, second.beta_rkf45);
        std::swap(first.gamma4_rkf45, second.gamma4_rkf45);
        std::swap(first.gamma5_rkf45, second.gamma5_rkf45);
    }

};

#include "ODE.cpp"

#endif /* ODE_ */
