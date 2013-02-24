/**
 * @file ODE.cpp
 * @author Robert Grandin
 * @brief Implementation of ODE class.
 */


#include "ODE.h"

/* ==================================================================
 * ================
 * ================    PRIVATE FUNCTIONS
 * ================
 */

// CREATE AND ASSIGN VALUES FOR RK4 WEIGHT FACTORS
template <class T>
void ODE<T>::AssignRK4Weights()
{
	/*
	* THIS ROUTINE ASSIGNS THE NECESSARY WEIGHT FACTORS FOR THE RUNGE-KUTTA 4TH
	* ORDER SCHEME.  THE VALUES BELOW FILL IN THE BUTCHER TABLE FOR THIS METHOD.
	*/

	alpha_rk4.ResetSize(4,0.0);
	beta_rk4.ResetSize(4,3,0.0);
	gamma_rk4.ResetSize(4,0.0);
	rk4butcherset = true;

    alpha_rk4(0) = 0.0e0;
    alpha_rk4(1) = 0.5e0;
    alpha_rk4(2) = 0.5e0;
    alpha_rk4(3) = 1.0e0;

    beta_rk4(1,0) = 0.5e0;
    beta_rk4(2,0) = 0.0e0;
    beta_rk4(2,1) = 0.5e0;
    beta_rk4(3,0) = 0.0e0;
    beta_rk4(3,1) = 0.0e0;
    beta_rk4(3,2) = 1.0e0;

    gamma_rk4(0) = 1.0e0/6.0e0;
    gamma_rk4(1) = 1.0e0/3.0e0;
    gamma_rk4(2) = 1.0e0/3.0e0;
    gamma_rk4(3) = 1.0e0/6.0e0;
}


// CREATE AND ASSIGN VALUES FOR RK45 WEIGHT FACTORS
template <class T>
void ODE<T>::AssignRKF45Weights()
{
	/*
	* THIS ROUTINE ASSIGNS THE NECESSARY WEIGHT FACTORS FOR THE RUNGE-KUTTA-
	* FEHLBERG 4TH/5TH ORDER SCHEME.  THE WEIGHT FACTORS BELOW FILL IN THE
	* BUTCHER TABLE.
	*/

    alpha_rk4.ResetSize(6,0.0e0);
    beta_rk4.ResetSize(6,5,0.0e0);
    gamma_rk4.ResetSize(6,0.0e0);
	rkf45butcherset = true;

    alpha_rkf45(0) = 0.0e0;
    alpha_rkf45(1) = 0.25e0;
    alpha_rkf45(2) = 0.375e0;
    alpha_rkf45(3) = 12.0e0/13.0e0;
    alpha_rkf45(4) = 1.0e0;
    alpha_rkf45(5) = 0.5e0;

    beta_rkf45(1,0) = 0.25e0;
    beta_rkf45(2,0) = 3.0e0/32.0e0;
    beta_rkf45(2,1) = 9.0e0/32.0e0;
    beta_rkf45(3,0) = 1932.0e0/2197.0e0;
    beta_rkf45(3,1) = -7200.0e0/2197.0e0;
    beta_rkf45(3,2) = 7296.0e0/2197.0e0;
    beta_rkf45(4,0) = 439.0e0/216.0e0;
    beta_rkf45(4,1) = -8.0e0;
    beta_rkf45(4,2) = 3680.0e0/513.0e0;
    beta_rkf45(4,3) = -845.0e0/4104.0e0;
    beta_rkf45(5,0) = -8.0e0/27.0e0;
    beta_rkf45(5,1) = 2.0e0;
    beta_rkf45(5,2) = -3544.0e0/2565.0e0;
    beta_rkf45(5,3) = 1859.0e0/4104.0e0;
    beta_rkf45(5,4) = -11.0e0/40.0e0;

    gamma4_rkf45(0) = 25.0e0/216.0e0;
    gamma4_rkf45(1) = 0.0e0;
    gamma4_rkf45(2) = 1408.0e0/2565.0e0;
    gamma4_rkf45(3) = 2197.0e0/4104.0e0;
    gamma4_rkf45(4) = -0.2e0;
    gamma4_rkf45(5) = 0.0e0;

    gamma5_rkf45(0) = 16.0e0/135.0e0;
    gamma5_rkf45(1) = 0.0e0;
    gamma5_rkf45(2) = 6656.0e0/12825.0e0;
    gamma5_rkf45(3) = 28561.0e0/56430.0e0;
    gamma5_rkf45(4) = -9.0e0/50.0e0;
    gamma5_rkf45(5) = 2.0e0/55.0e0;
}





/* ==================================================================
 * ================
 * ================    PUBLIC FUNCTIONS
 * ================
 */
template <class T>
ODE<T>::ODE()
{
	rk4butcherset = false;
	rkf45butcherset = false;
}


template <class T>
ODE<T>::ODE(const ODE<T> &a) : ODE()
{
    ODESwap(*this, a);
}


#ifdef CXX11
template <class T>
ODE<T>::ODE(ODE<T> &&a) : ODE<T>()
{
    ODESwap(*this, a);
}
#endif


template <class T>
ODE<T>::~ODE()
{
	;	// NOTHING TO DO IN DECONSTRUCTOR
}


template <class T>
ODE<T>& ODE<T>::operator=(ODE<T> a)
{
    ODESwap(*this, a);
    return *this;
}


template <class T>
int ODE<T>::f(Array1D<T> &yin,	Array1D<T> &yout,T x)
{
 	/*
 	* THIS FUNCTION CALCULATES THE FIRST-ORDER DERIVATIVES USED WHEN EVALUATING
 	* ORDINARY DIFFERENTIAL EQUATIONS.  THE CURRENT STATE VECTOR IS PASSED-IN
	* AS A REFERENCE TO AN Array1D OBJECT AND THE DERIVATIVES ARE PASSED-OUT
 	* THROUGH THE REFERENCE TO A SECOND Array1D OBJECT.  THE ABILITY OF THE
 	* Array1D OBJECTS TO STORE THEIR SIZE REMOVES THE NEED TO PASS THE ARRAY
 	* SIZE TO THIS ROUTINE.
 	*/

 	// CALCULATE OUTPUT DERIVATIVE VECTOR
 	//yout(0) = -0.2*yin(0);
	return 0;
}


// ADAMS-BASHFORTH, 2ND ORDER
template <class T>
void ODE<T>::AB2(T h, T x, Array1D<T> &y, Array2D<T> &derivs,
		int(*f)(Array1D<T>&, Array1D<T>&, Array1D<T>&))
{
	/*
	* THIS ROUTINE SOLVES ORDINARY DIFFERENTIAL EQUATIONS USING THE 2ND ORDER
	* ADAMS-BASHFORTH METHOD.  THE PREVIOUS DERIVATIVE VALUES ARE PASSED-IN AND
	* THE CURRENT-TIME-STEP DERIVATIVES ARE CALCULATED DURING THE EXECUTION OF
	* THIS ROUTINE.
	*/

	// GET DERIVATIVES FOR CURRENT TIME STEP
	int n = y.GetDim();
	Array1D<T> currderivs(n,0.0);
	f(y,currderivs,x);

	// CALCULATE NEXT VALUE OF STATE VECTOR
	for(int i=0; i<n; i++){
		y(i) += h*(1.5*currderivs(i) - 0.5*derivs(i,0));

		// UPDATE DERIVATIVES
		derivs(i,0) = currderivs(i);
	}
}


// ADAMS-BASHFORTH, 4TH ORDER
template <class T>
void ODE<T>::AB4(T h, T x, Array1D<T> &y, Array2D<T> &derivs,
		int(*f)(Array1D<T>&, Array1D<T>&, Array1D<T>&))
{
	/*
	 * THIS ROUTINE SOLVES ORDINARY DIFFERENTIAL EQUATIONS USING THE 4TH ORDER
	 * ADAMS-BASHFORTH METHOD.  THE PREVIOUS DERIVATIVE VALUES ARE PASSED-IN AND
	 * THE CURRENT-TIME-STEP DERIVATIVES ARE CALCULATED DURING THE EXECUTION OF
	 * THIS ROUTINE.
	 */

	// GET DERIVATIVES FOR CURRENT TIME STEP
	int n = y.GetDim();
	Array1D<T> currderivs(n,0.0);
	f(y,currderivs,x);

	// CALCULATE NEXT VALUE OF STATE VECTOR
	for(int i=0; i<n; i++){
		y(i) += (h/24.0)*(55.0*currderivs(i) - 59.0*derivs(i,0) +
				37.0*derivs(i,1) - 9.0*derivs(i,2));

		// UPDATE DERIVATIVES
		derivs(i,2) = derivs(i,1);
		derivs(i,1) = derivs(i,0);
		derivs(i,0) = currderivs(i);
	}
}


// EXPLICIT EULER METHOD
template <class T>
void ODE<T>::Euler(T h, T x, Array1D<T> &y,
		int(*f)(Array1D<T>&, Array1D<T>&, T))
{
	/*
	 * THIS ROUTINE SOLVES ORDINARY DIFFERENTIAL EQUATIONS USING THE FORWARD
	 * EULER EXPLICIT METHOD.
	 */

	int n = y.GetDim();
	Array1D<T> currderivs(n,0.0);
	f(y,currderivs,x);

	for(int i=0; i<n; i++){
		y(i) += h*currderivs(i);
	}
}



// RUNGE-KUTTA 2ND-ORDER
template <class T>
void ODE<T>::RK2(T h, T x, Array1D<T> &y, int(*f)(Array1D<T>&, Array1D<T>&, T))
{
	int n = y.GetDim();
	int neq = n;
	Array1D<T> k1(n,0.0);
	Array1D<T> yinit(n,0.0);
	T xt = x + 0.5*h;

	for(int i=0; i<n; i++){
		yinit(i) = y(i);
	}
	f(yinit,y,x);

	for(int i=0; i<neq; i++){
		k1(i) = yinit(i) + 0.5*h*y(i);
	}

	f(k1,y,xt);
	for(int i=0; i<neq; i++){
		y(i) = yinit(i) + h*y(i);
	}

}


template <class T>
T ODE<T>::RK2(T h, T x, Grid2D<T> &grid, void (*f)(Grid2D<T> &, Grid2D<T> &, T))
{
    int nx = grid.GetSize(2);
    int ny = grid.GetSize(1);
    T xt = x + 0.5e0*h;

    Grid2D<T> k1;
    k1.CreateGenericHGrid(nx,ny,0.0e0,1.0e0,0.0e0,1.0e0);
    k1.CopyNodalCoordinates(grid);
    k1.AddScalarQuantity("density");
    k1.AddScalarQuantity("energy");
    k1.AddScalarQuantity("pressure");
    k1.AddVectorQuantity("velocity");
    k1.AddVectorQuantity("momentum");

    Grid2D<T> gridinit;
    gridinit.CreateGenericHGrid(nx,ny,0.0e0,1.0e0,0.0e0,1.0e0);
    gridinit.CopyNodalCoordinates(grid);
    gridinit.AddScalarQuantity("density");
    gridinit.AddScalarQuantity("energy");
    gridinit.AddScalarQuantity("pressure");
    gridinit.AddVectorQuantity("velocity");
    gridinit.AddVectorQuantity("momentum");

    /* Note: Only density, energy, and momentum are integrated through time.  Pressure
             and velocity are calculated from these values as-needed.
    */

    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            gridinit(i,j,0) = grid(i,j,0);
            gridinit(i,j,1) = grid(i,j,1);
            gridinit(i,j,0,1) = grid(i,j,0,1);
            gridinit(i,j,1,1) = grid(i,j,1,1);
        }
    }

    f(gridinit,grid,x);

    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            k1(i,j,0) = gridinit(i,j,0) + 0.5e0*h*grid(i,j,0);
            k1(i,j,1) = gridinit(i,j,1) + 0.5e0*h*grid(i,j,1);
            k1(i,j,0,1) = gridinit(i,j,0,1) + 0.5e0*h*grid(i,j,0,1);
            k1(i,j,1,1) = gridinit(i,j,1,1) + 0.5e0*h*grid(i,j,1,1);
        }
    }

    f(k1,grid,xt);

    /* Calculate L2 residual of the most-recent derivative grid and return the total residual */
    T residual = 0.0e0;
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            residual += grid(i,j,0)*grid(i,j,0);
            residual += grid(i,j,1)*grid(i,j,1);
            residual += grid(i,j,0,1)*grid(i,j,0,1);
            residual += grid(i,j,1,1)*grid(i,j,1,1);
        }
    }
    residual = sqrt(residual);

    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            grid(i,j,0) = gridinit(i,j,0) + h*grid(i,j,0);
            grid(i,j,1) = gridinit(i,j,1) + h*grid(i,j,1);
            grid(i,j,0,1) = gridinit(i,j,0,1) + h*grid(i,j,0,1);
            grid(i,j,1,1) = gridinit(i,j,1,1) + h*grid(i,j,1,1);
        }
    }

    return residual;
}


template <class T>
T ODE<T>::RK2(T h, T x, Grid2D<T> &grid, void(*f)(Grid2D<T> &, Array3D<T> &, const T, Array2D<int> &, const int, const T, const T),
              Array2D<int> &dependency, const int nrecv, const T mu, const T kthermal)
{
    int nx = grid.GetSize(2);
    int ny = grid.GetSize(1);
    T xt = x + 0.5e0*h;

    Grid2D<T> k1;
    k1.CreateGenericHGrid(nx,ny,0.0e0,1.0e0,0.0e0,1.0e0);
    k1.CopyNodalCoordinates(grid);
    k1.AddScalarQuantity("density");
    k1.AddScalarQuantity("xmomentum");
    k1.AddScalarQuantity("ymomentum");
    k1.AddScalarQuantity("energy");

    Array3D<T> dvals(ny,nx,4,0.0e0);


    /* Note: Only density, energy, and momentum are integrated through time.  Pressure
             and velocity are calculated from these values as-needed.
    */
    int qs = -1;
    int qe = -1;
    int nw = nx*ny;
    int rank = 0;
    int np = 1;
    int i = -1;
    int j = -1;
#ifdef USEMPI
    rank = MPI::COMM_WORLD.Get_rank();
    np = MPI::COMM_WORLD.Get_size();
#endif

    MPI_Custom::LocalWorkUnits(rank,np,nw,qs,qe);

    f(grid,dvals,x,dependency,nrecv,mu,kthermal);

    for(int q=qs; q<=qe; q++){
        MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
        k1(i,j,0) = grid(i,j,0) + 0.5e0*h*dvals(i,j,0);
        k1(i,j,1) = grid(i,j,1) + 0.5e0*h*dvals(i,j,1);
        k1(i,j,2) = grid(i,j,2) + 0.5e0*h*dvals(i,j,2);
        k1(i,j,3) = grid(i,j,3) + 0.5e0*h*dvals(i,j,3);
#ifndef RELEASE
        T g0 = grid(i,j,0);
        T g1 = grid(i,j,1);
        T g2 = grid(i,j,2);
        T g3 = grid(i,j,3);
        T d0 = 0.5e0*h*dvals(i,j,0);
        T d1 = 0.5e0*h*dvals(i,j,1);
        T d2 = 0.5e0*h*dvals(i,j,2);
        T d3 = 0.5e0*h*dvals(i,j,3);
        T k0 = k1(i,j,0);
        T k11 = k1(i,j,1);
        T k2 = k1(i,j,2);
        T k3 = k1(i,j,3);
        g0 = grid(i,j,0);
#endif
    }

    f(k1,dvals,xt,dependency,nrecv,mu,kthermal);

    for(int q=qs; q<=qe; q++){
        MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
        grid(i,j,0) = grid(i,j,0) + h*dvals(i,j,0);
        grid(i,j,1) = grid(i,j,1) + h*dvals(i,j,1);
        grid(i,j,2) = grid(i,j,2) + h*dvals(i,j,2);
        grid(i,j,3) = grid(i,j,3) + h*dvals(i,j,3);
#ifndef RELEASE
        T g0 = grid(i,j,0);
        T g1 = grid(i,j,1);
        T g2 = grid(i,j,2);
        T g3 = grid(i,j,3);
        T d0 = h*dvals(i,j,0);
        T d1 = h*dvals(i,j,1);
        T d2 = h*dvals(i,j,2);
        T d3 = h*dvals(i,j,3);
        g0 = grid(i,j,0);
#endif
    }

    /* Calculate L2 residual of the most-recent derivative grid and return the total residual */
    T residual = 0.0e0;
    for(int q=qs; q<=qe; q++){
        MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
        residual += dvals(i,j,0)*dvals(i,j,0);
        residual += dvals(i,j,1)*dvals(i,j,1);
        residual += dvals(i,j,2)*dvals(i,j,2);
        residual += dvals(i,j,3)*dvals(i,j,3);
    }
    residual = sqrt(residual);

    return residual;
}


// RK4
template <class T>
void ODE<T>::RK4(T h, T x, Array1D<T> &y,
		int(*f)(Array1D<T>&, Array1D<T>&, T))
{
	/*
	 * THIS ROUTINE SOLVES ORDINARY DIFFERENTIAL EQUATIONS USING A RUNGE-KUTTA
	 * 4TH ORDER SCHEME.
	 */
	if(rk4butcherset == false){ODE<T>::AssignRK4Weights();}

	int n = y.GetDim();
	Array1D<T> k1(n,0.0);
	Array1D<T> k2(n,0.0);
	Array1D<T> k3(n,0.0);
	Array1D<T> k4(n,0.0);
	Array1D<T> yinit(n,0.0);
	Array1D<T> ytmp(n,0.0);
	T t;

	// ASSIGN yinit TO EQUAL y
	for(int i=0; i<n; i++){
		yinit(i) = y(i);
	}

	// FIRST STEP
	t = x + alpha_rk4(0)*h;
	f(yinit,y,t);
	for(int i=0; i<n; i++){
		k1(i) = h*y(i);
	}

	// SECOND STEP
	t = x + alpha_rk4(1)*h;
	for(int i=0; i<n; i++){
		ytmp(i) = yinit(i) + beta_rk4(1,0)*k1(i);
	}
	f(ytmp,y,t);
	for(int i=0; i<n; i++){
		k2(i) = h*y(i);
	}

	// THIRD STEP
	t = x + alpha_rk4(2)*h;
	for(int i=0; i<n; i++){
		ytmp(i) = yinit(i) + beta_rk4(2,0)*k1(i) + beta_rk4(2,1)*k2(i);
	}
	f(ytmp,y,t);
	for(int i=0; i<n; i++){
		k3(i) = h*y(i);
	}

	// FOURTH STEP
	t = x + alpha_rk4[3]*h;
	for(int i=0; i<n; i++){
		ytmp(i) = yinit(i) + beta_rk4(3,0)*k1(i) + beta_rk4(3,1)*k2(i) +
				beta_rk4(3,2)*k3(i);
	}
	f(ytmp,y,t);
	for(int i=0; i<n; i++){
		k4(i) = h*y(i);
	}

	// WEIGHTED AVERAGE
	for(int i=0; i<n; i++){
		y(i) = yinit(i) + gamma_rk4(0)*k1(i) + gamma_rk4(1)*k2(i) +
				gamma_rk4(2)*k3(i) + gamma_rk4(3)*k4(i);
	}
}


// RK45
template <class T>
void ODE<T>::RKF45(T h, T x, Array1D<T> &y, Array1D<T> &y4,
		int(*f)(Array1D<T>&, Array1D<T>&, T))
{
	/*
	 * THIS ROUTINE SOLVES ORDINARY DIFFERENTIAL EQUATIONS USING THE RUNGE-KUTTA-
	 * FEHLBERG COMBINED 4TH/5TH ORDER SCHEME.  SOLUTION ESTIMATES FOR BOTH ORDERS
	 * ARE DETERMINED AND BOTH SETS OF VALUES ARE RETURNED FOR USE IN ERROR
	 * ESTIMATION.
	 */
	if(rkf45butcherset == false){ODE<T>::AssignRKF45Weights();}

	int n = y.GetDim();
	Array1D<T> k1(n,0.0);
	Array1D<T> k2(n,0.0);
	Array1D<T> k3(n,0.0);
	Array1D<T> k4(n,0.0);
	Array1D<T> k5(n,0.0);
	Array1D<T> k6(n,0.0);
	Array1D<T> yinit(n,0.0);
	Array1D<T> ytmp(n,0.0);
	T t;

	// ASSIGN yinit TO EQUAL y
	for(int i=0; i<n; i++){
		yinit(i) = y(i);
	}

	// FIRST STEP
	t = x + alpha_rkf45(0)*h;
	odefcn(yinit,y,t);
	for(int i=0; i<n; i++){
		k1(i) = h*y(i);
	}

	// SECOND STEP
	t = x + alpha_rkf45(1)*h;
	for(int i=0; i<n; i++){
		ytmp(i) = yinit(i) + beta_rkf45(1,0)*k1(i);
	}
	odefcn(ytmp,y,t);
	for(int i=0; i<n; i++){
		k2(i) = h*y(i);
	}

	// THIRD STEP
	t = x + alpha_rkf45(2)*h;
	for(int i=0; i<n; i++){
		ytmp(i) = yinit(i) + beta_rkf45(2,0)*k1(i) + beta_rkf45(2,1)*k2(i);
	}
	odefcn(ytmp,y,t);
	for(int i=0; i<n; i++){
		k3(i) = h*y(i);
	}

	// FOURTH STEP
	t = x + alpha_rkf45(3)*h;
	for(int i=0; i<n; i++){
		ytmp(i) = yinit(i) + beta_rkf45(3,0)*k1(i) + beta_rkf45(3,1)*k2(i) +
				beta_rkf45(3,2)*k3(i);
	}
	odefcn(ytmp,y,t);
	for(int i=0; i<n; i++){
		k4(i) = h*y(i);
	}

	// FIFTH STEP
	t = x + alpha_rkf45(4)*h;
	for(int i=0; i<n; i++){
		ytmp(i) = yinit(i) + beta_rkf45(4,0)*k1(i) + beta_rkf45(4,1)*k2(i) +
				beta_rkf45(4,2)*k3(i) + beta_rkf45(4,3)*k4(i);
	}
	odefcn(ytmp,y,t);
	for(int i=0; i<n; i++){
		k5(i) = h*y(i);
	}

	// SIXTH STEP
	t = x + alpha_rkf45(5)*h;
	for(int i=0; i<n; i++){
		ytmp(i) = yinit(i) + beta_rkf45(5,0)*k1(i) + beta_rkf45(5,1)*k2(i) +
				beta_rkf45(5,2)*k3(i) + beta_rkf45(5,3)*k4(i) + beta_rkf45(5,4)*k5(i);
	}
	odefcn(ytmp,y,t);
	for(int i=0; i<n; i++){
		k6(i) = h*y(i);
	}

	// WEIGHTED AVERAGE
	for(int i=0; i<n; i++){
		// 4th ORDER SOLUTION
		y4(i) = yinit(i) + gamma4_rkf45(0)*k1(i) + gamma4_rkf45(1)*k2(i) +
				gamma4_rkf45(2)*k3(i) +	gamma4_rkf45(3)*k4(i) +
				gamma4_rkf45(4)*k5(i) + gamma4_rkf45(5)*k6(i);

		// 5th ORDER SOLUTION
		y(i) = yinit(i) + gamma5_rkf45(0)*k1(i) + gamma5_rkf45(1)*k2(i) +
				gamma5_rkf45(2)*k3(i) +	gamma5_rkf45(3)*k4(i) +
				gamma5_rkf45(4)*k5(i) + gamma5_rkf45(5)*k6(i);

	}
}


// LAX-WENDROFF METHOD
template <class T>
void ODE<T>::LaxWendroff(T h, T x, Array1D<T> &y, T c, int bc, int var,
		Array2D<T> &flux)
{
	/*
	 * THIS ROUTINE SOLVES ORDINARY DIFFERENTIAL EQUATIONS USING THE LAX-
	 * WENDROFF METHOD.
	 */

	int n = y.GetDim();
	Array1D<T> um(n,0.0);
	Array1D<T> up(n,0.0);
	Array1D<T> yinit(n,0.0);
	for(int i=0; i<n; i++){
		yinit(i) = y(i);
	}
	T pm1, pp1, p;

	if(bc == 1){ 	// FIXED END POINTS
		for(int i=1; i<n-1; i++){
			pm1 = flux(i-1,var);
			pp1 = flux(i+1,var);
			p = flux(i,var);

			up(i) = 0.5*(p + pp1) - 0.5*c*(pp1 - p);
			um(i) = 0.5*(p + pm1) - 0.5*c*(p - pm1);
		}
		for(int i=0; i<n; i++){
			y(i) = yinit(i) - c*(up(i) - um(i));
		}
	}
	if(bc == 2){	// PERIODIC END POINTS
		for(int i=0; i<n; i++){
			if(i == 0){
				pm1 = flux(n-1,var);
				pp1 = flux(i+1,var);
				p = flux(i,var);
			}
			if(i == n-1){
				pp1 = flux(0,var);
				pm1 = flux(i-1,var);
				p = flux(i,var);
			}
			if(i > 0 && i < n-1){
				pm1 = flux(i-1,var);
				pp1 = flux(i+1,var);
				p = flux(i,var);
			}
			up(i) = 0.5*(p + pp1) - 0.5*c*(pp1 - p);
			um(i) = 0.5*(p + pm1) - 0.5*c*(p - pm1);
		}
		for(int i=0; i<n; i++){
			y(i) = yinit(i) - c*(up(i) - um(i));
		}
	}
}


// MACCORMIC METHOD
template <class T>
void ODE<T>::MacCormack(T h, T x, Array1D<T> &y, T c, int bc, int var,
		Array2D<T> &flux, T disp)
{
	/*
	 * THIS ROUTINE SOLVES ORDINARY DIFFERENTIAL EQUATIONS USING THE MACCORMICK
	 * METHOD.
	 */

	int n = y.GetDim();
	Array1D<T> us(n,0.0);
	Array1D<T> yinit(n,0.0);
	for(int i=0; i<n; i++){
		yinit(i) = y(i);
	}
	T pm1, pp1, p;

	if(bc == 1){ 	// FIXED END POINTS
		for(int i=1; i<n-1; i++){
			pm1 = flux(i-1,var);
			pp1 = flux(i+1,var);
			p = flux(i,var);
			us(i) = yinit(i) - c*(pp1 - p);
		}
		for(int i=1; i<n-1; i++){
			y(i) = 0.5*(us(i) + yinit(i)) - 0.5*c*(us(i) - us(i-1));
		}
	}
	if(bc == 2){	// PERIODIC END POINTS
		for(int i=0; i<n; i++){
			if(i == 0){
				pm1 = flux(n-1,var);
				pp1 = flux(i+1,var);
				p = flux(i,var);
			}
			if(i == n-1){
				pp1 = flux(0,var);
				pm1 = flux(i-1,var);
				p = flux(i,var);
			}
			if(i > 0 && i < n-1){
				pm1 = flux(i-1,var);
				pp1 = flux(i+1,var);
				p = flux(i,var);
			}
			us(i) = yinit(i) - c*(pp1 - p);
		}
		for(int i=0; i<n; i++){
			if(i == 0){ pm1 = us(n-1); }
			if(i > 0 && i < n-1){
				pm1 = us(i-1);
			}
			y(i) = 0.5*(us(i) + yinit(i)) - 0.5*c*(us(i) - pm1);
		}
    }
}



