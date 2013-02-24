/**
 * @file Grid2D.cpp
 * @author Robert Grandin
 * @brief Implementation of Grid2D class.
 */



#include <Grid2D.h>


/* ==================================================================
 * ================
 * ================    PRIVATE FUNCTIONS
 * ================
 */







/* ==================================================================
 * ================
 * ================    PUBLIC FUNCTIONS
 * ================
 */

template <class T>
Grid2D<T>::Grid2D()
{
	isogrid = false;
	ishgrid = false;
    pi = 4.0e0*atan(1.0e0);
    nscalars = 0;
    nvectors = 0;
    gridname = "Gridname_Not_Set";
    qtysize =256;
    inc_snapshots = -1;
    count_snapshots = 0;
    name_snapshots = "snapshot";
}


template <class T>
Grid2D<T>::Grid2D(const Grid2D<T> &a) : Grid2D()
{
    Grid2DSwap(*this, a);
}


#ifdef CXX11
template <class T>
Grid2D<T>::Grid2D(Grid2D<T> &&a) : Grid2D<T>()
{
    Grid2DSwap(*this, a);
}
#endif


template <class T>
Grid2D<T>::~Grid2D()
{
    /*
      Delete arrays of 2D and 3D quantities not required here due to
      the destructor of PArray1D doing this.
     */

}


template <class T>
void Grid2D<T>::CreateOGridFromFunctions(int nrad, int ntheta, T(*fx)(T, T, T, T),
                                         T(*fy)(T, T, T, T), T rotouter,
                                         T rotinner, T rotoffset, T ao,
                                         T bo, T ai, T bi, T Pp, T Qq)
        {
    /*
     * DIMENSION 2D ARRAYS FOR SOLUTION VALUES, DERIVATIVE METRICS, AND
     * SPATIAL COORDINATES OF COMPUTATIONAL NODES.
     */
    dzidx.ResetSize(nrad,ntheta);
    dzidy.ResetSize(nrad,ntheta);
    detadx.ResetSize(nrad,ntheta);
    detady.ResetSize(nrad,ntheta);
    xcoords.ResetSize(nrad,ntheta);
    ycoords.ResetSize(nrad,ntheta);
    iblank.ResetSize(nrad,ntheta,1.0);
    nx = ntheta;
    ny = nrad;
    isogrid = true;


    /*
     * DETERMINE THE X AND Y COORDINATES OF THE INNER AND OUTER SURFACES.
     */
    T thetamin = 0.0e0;
    T thetamax = 2.0e0*pi;
    T dtheta = (thetamax - thetamin)/((T)ntheta - 1.0e0);

    for(int i=0; i<ntheta; i++){
        T theta = thetamin + dtheta*(T)i;
        T tmpval;

        tmpval = fx(theta+rotoffset,rotinner,ai,bi);
        xcoords.SetValue(0,i,tmpval);
        tmpval = fx(theta,rotouter,ao,bo);
        xcoords.SetValue(nrad-1,i,tmpval);
        tmpval = fy(theta+rotoffset,rotinner,ai,bi);
        ycoords.SetValue(0,i,tmpval);
        tmpval = fy(theta,rotouter,ao,bo);
        ycoords.SetValue(nrad-1,i,tmpval);
    }


    /*
     * USE HERMITE INTERPOLATION TO ENFORCE ORTHOGONALITY AT BOUNDARY SURFACES
     * WHEN DETERMINING THE (X,Y) COORDINATES OF THE GRID.
     */
    for(int j=0; j<ntheta; j++){

        // DETERMINE UNIT NORMAL AT BOUNDARY SURFACES
        T dxin;
        T dxout;
        T dyin;
        T dyout;
        if(j > 0 && j < nx-1){
            dxin = xcoords.GetValue(0,j+1) - xcoords.GetValue(0,j-1);
            dxout = xcoords.GetValue(ny-1,j+1) - xcoords.GetValue(ny-1,j-1);
            dyin = ycoords.GetValue(0,j+1) - ycoords.GetValue(0,j-1);
            dyout = ycoords.GetValue(ny-1,j+1) - ycoords.GetValue(ny-1,j-1);
        }
        if(j == 0){
            dxin = xcoords.GetValue(0,j+1) - xcoords.GetValue(0,j);
            dxout = xcoords.GetValue(ny-1,j+1) - xcoords.GetValue(ny-1,j);
            dyin = ycoords.GetValue(0,j+1) - ycoords.GetValue(0,j);
            dyout = ycoords.GetValue(ny-1,j+1) - ycoords.GetValue(ny-1,j);
        }
        if(j == nx-1){
            dxin = xcoords.GetValue(0,j) - xcoords.GetValue(0,j-1);
            dxout = xcoords.GetValue(ny-1,j) - xcoords.GetValue(ny-1,j-1);
            dyin = ycoords.GetValue(0,j) - ycoords.GetValue(0,j-1);
            dyout = ycoords.GetValue(ny-1,j) - ycoords.GetValue(ny-1,j-1);
        }
        T ynin = dxin/(sqrt(dxin*dxin + dyin*dyin));
        T xnin = -dyin/(sqrt(dxin*dxin + dyin*dyin));
        T ynout = dxout/(sqrt(dxout*dxout + dyout*dyout));
        T xnout = -dyout/(sqrt(dxout*dxout + dyout*dyout));
        T xinval = xcoords.GetValue(0,j);
        T xoutval = xcoords.GetValue(ny-1,j);
        T yinval = ycoords.GetValue(0,j);
        T youtval = ycoords.GetValue(ny-1,j);

        for(int i=0; i<nrad; i++){
            T s = (T)i/((T)nrad - 1.0);
            T H1 = 1.0 - 3.0*s*s + 2.0*s*s*s;
            T H2 = 3.0*s*s - 2.0*s*s*s;
            T H3 = (s*s*s - 2.0*s*s + s)*Pp;
            T H4 = (s*s*s - s*s)*Qq;
            T xval = xinval*H1 + xoutval*H2 + xnin*H3 + xnout*H4;
            T yval = yinval*H1 + youtval*H2 + ynin*H3 + ynout*H4;
            xcoords.SetValue(i,j,xval);
            ycoords.SetValue(i,j,yval);
        }
    } /* END OF COORDINATE-DETERMINATION */

    // COMPUTE METRICS
    Grid2D<T>::ComputeMetrics();
}


template <class T>
void Grid2D<T>::CreateHGridFromFunctions(int nxx, int nyy, T(*fxin)(T),
        T(*fxout)(T), T(*fyin)(T), T(*fyout)(T), const bool hermite, T Pp, T Qq,
        T xxmin, T xxmax, T yymin, T yymax, const bool stretch1, const T alpha, const T beta)
        {
    /*
     * DIMENSION 2D ARRAYS FOR SOLUTION VALUES, DERIVATIVE METRICS, AND
     * SPATIAL COORDINATES OF COMPUTATIONAL NODES.
     */
    dzidx.ResetSize(nyy,nxx);
    dzidy.ResetSize(nyy,nxx);
    detadx.ResetSize(nyy,nxx);
    detady.ResetSize(nyy,nxx);
    xcoords.ResetSize(nyy,nxx);
    ycoords.ResetSize(nyy,nxx);
    iblank.ResetSize(nyy,nxx,1.0);
    nx = nxx;
    ny = nyy;
    ishgrid = true;
    xmin = xxmin;
    xmax = xxmax;
    ymin = yymin;
    ymax = yymax;


    /*
     * DETERMINE THE X AND Y COORDINATES OF THE INNER AND OUTER SURFACES.
     */
    T dx = (xmax - xmin)/((T)nx - 1.0);
    T dy = (ymax - ymin)/((T)ny - 1.0e0);

    for(int i=0; i<nx; i++){
        T x = xmin + dx*(T)i;
        T tmpval;

        tmpval = fxin(x);
        xcoords.SetValue(0,i,tmpval);
        tmpval = fxout(x);
        xcoords.SetValue(ny-1,i,tmpval);
        tmpval = fyin(x);
        ycoords.SetValue(0,i,tmpval);
        tmpval = fyout(x);
        ycoords.SetValue(ny-1,i,tmpval);
//        if(ycoords(0,i) < ymin){
//            ymin = ycoords(0,i);
//        }
//        if(ycoords(ny-1,i) < ymin){
//            ymin = ycoords(0,i);
//        }
//        if(ycoords(0,i) < ymax){
//            ymax = ycoords(0,i);
//        }
//        if(ycoords(ny-1,i) < ymax){
//            ymax = ycoords(0,i);
//        }
    }


    if(hermite == true){
        /*
         * USE HERMITE INTERPOLATION TO ENFORCE ORTHOGONALITY AT BOUNDARY SURFACES
         * WHEN DETERMINING THE (X,Y) COORDINATES OF THE GRID.
         */
        for(int j=0; j<nx; j++){

            // DETERMINE UNIT NORMAL AT BOUNDARY SURFACES
            T dxin = 0.0e0;
            T dxout = 0.0e0;
            T dyin = 0.0e0;
            T dyout = 0.0e0;
            if(j > 0 && j < nx-1){
                dxin = xcoords.GetValue(0,j+1) - xcoords.GetValue(0,j-1);
                dxout = xcoords.GetValue(ny-1,j+1) - xcoords.GetValue(ny-1,j-1);
                dyin = ycoords.GetValue(0,j+1) - ycoords.GetValue(0,j-1);
                dyout = ycoords.GetValue(ny-1,j+1) - ycoords.GetValue(ny-1,j-1);
            }
            if(j == 0){
                dxin = xcoords.GetValue(0,j+1) - xcoords.GetValue(0,j);
                dxout = xcoords.GetValue(ny-1,j+1) - xcoords.GetValue(ny-1,j);
                dyin = ycoords.GetValue(0,j+1) - ycoords.GetValue(0,j);
                dyout = ycoords.GetValue(ny-1,j+1) - ycoords.GetValue(ny-1,j);
            }
            if(j == nx-1){
                dxin = xcoords.GetValue(0,j) - xcoords.GetValue(0,j-1);
                dxout = xcoords.GetValue(ny-1,j) - xcoords.GetValue(ny-1,j-1);
                dyin = ycoords.GetValue(0,j) - ycoords.GetValue(0,j-1);
                dyout = ycoords.GetValue(ny-1,j) - ycoords.GetValue(ny-1,j-1);
            }
            T ynin = dxin/(sqrt(dxin*dxin + dyin*dyin));
            T xnin = -dyin/(sqrt(dxin*dxin + dyin*dyin));
            T ynout = dxout/(sqrt(dxout*dxout + dyout*dyout));
            T xnout = -dyout/(sqrt(dxout*dxout + dyout*dyout));
            T xinval = xcoords.GetValue(0,j);
            T xoutval = xcoords.GetValue(ny-1,j);
            T yinval = ycoords.GetValue(0,j);
            T youtval = ycoords.GetValue(ny-1,j);

            for(int i=0; i<ny; i++){
                T s = (T)i/((T)ny - 1.0);
                T H1 = 1.0 - 3.0*s*s + 2.0*s*s*s;
                T H2 = 3.0*s*s - 2.0*s*s*s;
                T H3 = (s*s*s - 2.0*s*s + s)*Pp;
                T H4 = (s*s*s - s*s)*Qq;
                T xval = xinval*H1 + xoutval*H2 + xnin*H3 + xnout*H4;
                T yval = yinval*H1 + youtval*H2 + ynin*H3 + ynout*H4;
                xcoords(i,j) = xval;
                ycoords(i,j) = yval;
            }
        } /* END OF COORDINATE-DETERMINATION */
    }

    if(stretch1 == true){
        T ytmp = 0.0e0;
        T base = (beta + 1.0e0)/(beta - 1.0e0);
        T exponent = 0.0e0;
        for(int i=0; i<ny; i++){
            for(int j=0; j<nx; j++){
                xcoords(i,j) = xmin + dx*(T)j;
                ytmp = ymin + dy*(T)i;
                exponent = (ytmp - alpha)/(1.0e0 - alpha);
                ycoords(i,j) = (ycoords(ny-1,j) - ycoords(0,j))*((beta + 2.0e0*alpha)*pow(base,exponent) - beta + 2.0e0*alpha)/
                        ((2.0e0*alpha + 1.0e0)*(1.0e0 + pow(base,exponent)));
            }
        }
    }

    // COMPUTE METRICS
    Grid2D<T>::ComputeMetrics();
}


template <class T>
void Grid2D<T>::CreateGenericHGrid(const int nnx, const int nny, const T xxmin, const T xxmax, const T yymin, const T yymax)
{
    xmin = xxmin;
    xmax = xxmax;
    ymin = yymin;
    ymax = yymax;
    ny = nny;
    nx = nnx;

    dzidx.ResetSize(ny,nx);
    dzidy.ResetSize(ny,nx);
    detadx.ResetSize(ny,nx);
    detady.ResetSize(ny,nx);
    iblank.ResetSize(ny,nx,1.0e0);

    ishgrid = true;

    T dx = (xmax - xmin)/((T)nx - 1.0e0);
    T dy = (ymax - ymin)/((T)ny - 1.0e0);

    xcoords.ResetSize(ny,nx,0.0e0);
    ycoords.ResetSize(ny,nx,0.0e0);

    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            xcoords(i,j) = xmin + dx*(T)j;
            ycoords(i,j) = ymin + dy*(T)i;
        }
    }

    Grid2D<T>::ComputeMetrics();
}



template <class T>
void Grid2D<T>::EllipticSmoothing(const int niterations, const T tol, int &iterationsused, T &tolachieved,
                                  const bool fixbound1, const bool fixbound2, const bool fixbound3, const bool fixbound4)
{
    iterationsused = 0;
    tolachieved = 1.0e0;

    int iterlimit = niterations;
    T tollimit = tol;

    /* Check for valid termination controls.  If both controls disabled, set default value. */
    if(niterations < 0 && tol < 0){
        iterlimit = 100000;
        tollimit = 1.0e-4;
    } else {
        if(niterations < 0){
            /* Set to large, unlikely-to-be-encountered number. */
            iterlimit = 100000;
        }

        /* Negative 'tol' value will not prevent iterations and will therefore be left negative. */
    }

    /* Set loop limits in grid based on boolean input values */
    int lowlimx = 1;
    int highlimx = nx-1;
    int lowlimy = 1;
    int highlimy = ny-1;
    if(fixbound1 == false){
        lowlimy = 0;
    }
    if(fixbound2 == false){
        highlimy = ny;
    }
    if(fixbound3 == false){
        lowlimx = 0;
    }
    if(fixbound4 == false){
        highlimx = nx;
    }

    /* Iterate while tolerance and iteration limits are satisfied */
    T xnew = 0.0e0;
    T ynew = 0.0e0;
    T dx = 1.0e0;
    T dy = 1.0e0;
    T residual = 1.0e0;
    T residual0 = 1.0e0;
    T alpha = 0.0e0;
    T beta = 0.0e0;
    T gamma = 0.0e0;
    T dxdzi = 0.0e0;
    T dydzi = 0.0e0;
    T dxdeta = 0.0e0;
    T dydeta = 0.0e0;
    int im1 = 0;
    int ip1 = 0;
    int jm1 = 0;
    int jp1 = 0;

    while(fabs(tolachieved) >= tollimit && iterationsused <= iterlimit){

        residual = 0.0e0;

        /* Loop through all interior grid points*/
        for(int i=lowlimy; i<highlimy; i++){
            for(int j=lowlimx; j<highlimx; j++){
                /* Approximate spatial derivatives with respect to computational grid.  These are 2nd-order accurate */
                /* and do not assume any boundaries are periodic (regardless of the iteration bounds).               */
                im1 = i - 1;
                ip1 = i + 1;
                jm1 = j - 1;
                jp1 = j + 1;
                if(i == 0){
                    im1 = ny-2;
                }
                if(i == ny-1){
                    ip1 = 1;
                }
                if(j == 0){
                    jm1 = nx-2;
                }
                if(j == nx-1){
                    jp1 = 1;
                }

                dxdzi = 0.5e0*(xcoords(i,jp1) - xcoords(i,jm1));
                dydzi = 0.5e0*(ycoords(i,jp1) - ycoords(i,jm1));
                dxdeta = 0.5e0*(xcoords(ip1,j) - xcoords(im1,j));
                dydeta = 0.5e0*(ycoords(ip1,j) - ycoords(im1,j));




                /* Calculate metric terms, assumed to be constant */
                alpha = dxdeta*dxdeta + dydeta*dydeta;
                gamma = dxdzi*dxdzi + dydzi*dydzi;
                beta = dxdzi*dxdeta + dydzi*dydeta;

                /* Calculate change in x-coordinate using Gauss-Seidel iteration */

                xnew = (-alpha*(xcoords(i,jp1) + xcoords(i,jm1)) +
                        0.5e0*beta*(xcoords(ip1,jp1) - xcoords(ip1,jm1) - xcoords(im1,jp1) + xcoords(im1,jm1)) -
                        gamma*(xcoords(ip1,j) + xcoords(im1,j)))/(-2.0e0*(alpha+gamma));
                dx = xnew - xcoords(i,j);
                xcoords(i,j) = xnew;

                /* Calculate change in y-coordinate using Gauss-Seidel iteration */
                ynew = (-alpha*(ycoords(i,jp1) + ycoords(i,jm1)) +
                        0.5e0*beta*(ycoords(ip1,jp1) - ycoords(ip1,jm1) - ycoords(im1,jp1) + ycoords(im1,jm1)) -
                        gamma*(ycoords(ip1,j) + ycoords(im1,j)))/(-2.0e0*(alpha+gamma));
                dy = ynew - ycoords(i,j);
                ycoords(i,j) = ynew;

                /* Calculate residual to be the magnitude-squared of the change-in-position of the point */
                residual += dx*dx + dy*dy;
            }
        }

        /* Calculate the normalized residual for this iteration, normalized to the error obtained on the first iteration. */
        if(iterationsused == 0){
            residual0 = residual;
        }
        tolachieved = residual/residual0;

        /* Increment iteration counter */
        iterationsused++;
    }

}



template <class T>
void Grid2D<T>::WriteGrid(std::string file)
{
	/*
	 * VTK STRUCTURED GRID FILE CREATED TO SHOW THE GRID.  ALL GRID POINTS ARE
	 * ASSOCIATED WITH A SCALAR VALUE OF ZERO.  THE INTENT IS FOR THE GRID TO
	 * BE VIEWED IN PARAVIEW.
	 */


    std::fstream ssfile;
    ssfile.open(file.c_str(),std::ios::out);
    ssfile << "# vtk DataFile Version 2.0" << std::endl;
    ssfile << gridname << std::endl;
    ssfile << "ASCII" << std::endl;
    ssfile << "DATASET STRUCTURED_GRID" << std::endl;
    ssfile << "DIMENSIONS " << nx << " " << ny << " 1" << std::endl;
    ssfile << "POINTS " << nx*ny << " float" << std::endl;
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            ssfile << xcoords(i,j) << " " << ycoords(i,j) << " 0.0" << std::endl;
        }
    }
    ssfile << "POINT_DATA " << nx*ny << std::endl;
    ssfile << "SCALARS grid float 1" << std::endl;
    ssfile << "LOOKUP_TABLE default" << std::endl;
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            ssfile << "0" << std::endl;
        }
    }
    ssfile.close();
}


template <class T>
void Grid2D<T>::WriteSolution(std::string file)
{
    /*
      Write VTK Structured Grid file to show the solution.  Each scalar and
      vector quantity will be included in the file.
     */

    std::fstream ssfile;
    ssfile.open(file.c_str(),std::ios::out);
    ssfile << "# vtk DataFile Version 2.0" << std::endl;
    ssfile << gridname << std::endl;
    ssfile << "ASCII" << std::endl;
    ssfile << "DATASET STRUCTURED_GRID" << std::endl;
    ssfile << "DIMENSIONS " << nx << " " << ny << " 1" << std::endl;
    ssfile << "POINTS " << nx*ny << " float" << std::endl;
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            ssfile << xcoords(i,j) << " " << ycoords(i,j) << " 0.0" << std::endl;
        }
    }
    ssfile << "POINT_DATA " << nx*ny << std::endl;
    for(int s=0; s<nscalars; s++){
        std::string name(scalar_names(s)->substr());
        ssfile << "SCALARS " << name << " float 1" << std::endl;
        ssfile << "LOOKUP_TABLE default" << std::endl;
        for(int i=0; i<ny; i++){
            for(int j=0; j<nx; j++){
                ssfile << pscalars(s)->operator ()(i,j) << std::endl;
            }
        }
    }
    for(int v=0; v<nvectors; v++){
        std::string name(vector_names(v)->substr());
        ssfile << "VECTORS " << name << " float" << std::endl;
        for(int i=0; i<ny; i++){
            for(int j=0; j<nx; j++){
                ssfile << pvectors(v)->operator ()(i,j,0) << " " << pvectors(v)->operator ()(i,j,1) << " " <<
                       pvectors(v)->operator ()(i,j,2) << std::endl;
            }
        }
    }
    ssfile.close();

}


template <class T>
int Grid2D<T>::GetSize(int dim)
{
	int retval = -1;
	if(dim == 1){retval = ny;}
	if(dim == 2){retval = nx;}
	return retval;
}


template <class T>
void Grid2D<T>::ComputeMetrics()
{
	/*
	 * COMPUTE METRICS RELATING SPATIAL-DOMAIN DERIVATIVES TO COMPUTATIONAL-
	 * DOMAIN NODES.
	 */

	Array2D<T> dxdeta(ny,nx,0.0);
	Array2D<T> dydeta(ny,nx,0.0);
	Array2D<T> dxdzi(ny,nx,0.0);
	Array2D<T> dydzi(ny,nx,0.0);

	// Y-DIRECTION DIFFERENCES
	for(int i=1; i<ny-1; i++){
		for(int j=0; j<nx; j++){
			dxdeta(i,j) = 0.5*(xcoords(i+1,j) - xcoords(i-1,j));
			dydeta(i,j) = 0.5*(ycoords(i+1,j) - ycoords(i-1,j));
		}
	}
	for(int j=0; j<nx; j++){
		dxdeta(0,j) = 0.5*(-xcoords(2,j) + 4.0*xcoords(1,j) - 3.0*xcoords(0,j));
		dydeta(0,j) = 0.5*(-ycoords(2,j) + 4.0*ycoords(1,j) - 3.0*ycoords(0,j));
	}
	for(int j=0; j<nx; j++){
		dxdeta(ny-1,j) = 0.5*(3.0*xcoords(ny-1,j) - 4.0*xcoords(ny-2,j)
				+ xcoords(ny-3,j));
		dydeta(ny-1,j) = 0.5*(3.0*ycoords(ny-1,j) - 4.0*ycoords(ny-2,j)
				+ ycoords(ny-3,j));
	}

	// X-DIRECTION DIFFERENCES
	for(int j=1; j<nx-1; j++){
		for(int i=0; i<ny; i++){
			dxdzi(i,j) = 0.5*(xcoords(i,j+1) - xcoords(i,j-1));
			dydzi(i,j) = 0.5*(ycoords(i,j+1) - ycoords(i,j-1));
		}
	}
	for(int i=0; i<ny; i++){
		dxdzi(i,0) = 0.5*(-xcoords(i,2) + 4.0*xcoords(i,1) - 3.0*xcoords(i,0));
		dydzi(i,0) = 0.5*(-ycoords(i,2) + 4.0*ycoords(i,1) - 3.0*ycoords(i,0));
	}
	for(int i=0; i<ny; i++){
		dxdzi(i,nx-1) = 0.5*(3.0*xcoords(i,nx-1) - 4.0*xcoords(i,nx-2)
				+ xcoords(i,nx-3));
		dydzi(i,nx-1) = 0.5*(3.0*ycoords(i,nx-1) - 4.0*ycoords(i,nx-2)
				+ ycoords(i,nx-3));
	}

	// COMPUTE JACOBIAN AND METRICS
	for(int j=0; j<nx; j++){
		for(int i=0; i<ny; i++){
			T jacobian = dxdzi(i,j)*dydeta(i,j) - dxdeta(i,j)*dydzi(i,j);
			dzidx(i,j) = dydeta(i,j)/jacobian;
			detadx(i,j) = -dydzi(i,j)/jacobian;
			dzidy(i,j) = -dxdeta(i,j)/jacobian;
			detady(i,j) = dxdzi(i,j)/jacobian;
		}
	}
}


template <class T>
void Grid2D<T>::AddScalarQuantity(const std::string name)
{
    // Require that a grid has been generated
    if(isogrid == true || ishgrid == true){
        // Increase count of number of scalar quantities
        nscalars++;

        /*
          Increase the size of the PArray1D object containing pointers to
          Array2D objects containing data.  The data must be copied into temporary
          arrays since the current Array2D objects are destroyed when the
          PArray1D object is resized.  The data is then copied back into the
          resized PArray1D object.

          During the copy process, the source object is deleted immediately after
          its data has been copied.  This is to reduce memory requirements so that
          only one extra object's worth of memory is required during the
          copy process.  Note that setting the resulting pointer to NULL is
          needed to avoid freeing the memory more than once upon deletion of
          the temporary PArray1D object.

          After the copy is complete, add a new Array2D object to the final
          entry in the resized PArray1D object, set its size equal to the grid
          size, and all points initialized to 0.0.
         */
        PArray1D<Array2D<T>*> tmparray(nscalars);
        for(int i=0; i<nscalars-1; i++){
            tmparray(i) = new Array2D<T>;
            tmparray(i)->ResetSize(ny,nx);
            for(int ii=0; ii<ny; ii++){
                for(int jj=0; jj<nx; jj++){
                    tmparray(i)->operator ()(ii,jj) = pscalars(i)->operator ()(ii,jj);
                }
            }
            delete pscalars(i);
            pscalars(i) = NULL;
        }
        pscalars.ResetSize(nscalars);
        for(int i=0; i<nscalars-1; i++){
            pscalars(i) = new Array2D<T>;
            pscalars(i)->ResetSize(ny,nx);
            for(int ii=0; ii<ny; ii++){
                for(int jj=0; jj<nx; jj++){
                    pscalars(i)->operator ()(ii,jj) = tmparray(i)->operator ()(ii,jj);
                }
            }
            delete tmparray(i);
            tmparray(i) = NULL;
        }
        pscalars(nscalars-1) = new Array2D<T>;
        pscalars(nscalars-1)->ResetSize(ny,nx);
        pscalars(nscalars-1)->ResetVal((T)0.0e0);


        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nscalars);
        for(int i=0; i<nscalars-1; i++){
            tmparray2(i) = new std::string;
            tmparray2(i)->reserve(qtysize);
            tmparray2(i)->assign(scalar_names(i)->substr());
            delete scalar_names(i);
            scalar_names(i) = NULL;
        }
        scalar_names.ResetSize(nscalars);
        for(int i=0; i<nscalars-1; i++){
            scalar_names(i) = new std::string;
            scalar_names(i)->reserve(qtysize);
            scalar_names(i)->assign(tmparray2(i)->substr());
            delete tmparray2(i);
            tmparray2(i) = NULL;
        }
        scalar_names(nscalars-1) = new std::string;
        scalar_names(nscalars-1)->reserve(qtysize);
        scalar_names(nscalars-1)->assign(name);
    }
}


template <class T>
void Grid2D<T>::RemoveScalarQuantity(const int qty)
{
    // Require that a grid has been generated
    if((isogrid == true || ishgrid == true) && qty < nscalars){
        // Decrease count of number of scalar quantities
        nscalars--;

        /*
          Decrease the size of the PArray1D object containing pointers to
          Array2D objects containing data.  The data must be copied into temporary
          arrays since the current Array2D objects are destroyed when the
          PArray1D object is resized.  The data is then copied back into the
          resized PArray1D object.

          During the copy process, the source object is deleted immediately after
          its data has been copied.  This is to reduce memory requirements so that
          only one extra object's worth of memory is required during the
          copy process.  Note that setting the resulting pointer to NULL is
          needed to avoid freeing the memory more than once upon deletion of
          the temporary PArray1D object.

          After the copy is complete, add a new Array2D object to the final
          entry in the resized PArray1D object, set its size equal to the grid
          size, and all points initialized to 0.0.
         */
        PArray1D<Array2D<T>*> tmparray(nscalars+1);
        for(int i=0; i<nscalars+1; i++){
            tmparray(i) = new Array2D<T>;
            tmparray(i)->ResetSize(ny,nx);
            for(int ii=0; ii<ny; ii++){
                for(int jj=0; jj<nx; jj++){
                    tmparray(i)->operator ()(ii,jj) = pscalars(i)->operator ()(ii,jj);
                }
            }
            delete pscalars(i);
            pscalars(i) = NULL;
        }
        pscalars.ResetSize(nscalars);
        for(int i=0; i<nscalars+1; i++){
            int itmp = i;
            if(i != qty){
                if(i > qty){
                    itmp = i - 1;   // Reduce index by 1 since loop is past removed quantity
                }
                pscalars(itmp) = new Array2D<T>;
                pscalars(itmp)->ResetSize(ny,nx);
                for(int ii=0; ii<ny; ii++){
                    for(int jj=0; jj<nx; jj++){
                        pscalars(itmp)->operator ()(ii,jj) = tmparray(i)->operator ()(ii,jj);
                    }
                }
            }
            delete tmparray(i);
            tmparray(i) = NULL;
        }


        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nscalars+1);
        for(int i=0; i<nscalars+1; i++){
            tmparray2(i) = new std::string;
            tmparray2(i)->reserve(qtysize);
            tmparray2(i)->assign(scalar_names(i)->substr());
            delete scalar_names(i);
            scalar_names(i) = NULL;
        }
        scalar_names.ResetSize(nscalars);
        for(int i=0; i<nscalars+1; i++){
            int itmp = i;
            if(i != qty){
                if(i > qty){
                    itmp = i - 1;
                }
                scalar_names(itmp) = new std::string;
                scalar_names(itmp)->reserve(qtysize);
                scalar_names(itmp)->assign(tmparray2(i)->substr());
                delete tmparray2(i);
                tmparray2(i) = NULL;
            }
        }
    }
}


template <class T>
void Grid2D<T>::AddVectorQuantity(const std::string name)
{
    // Require that a grid has been generated
    if(isogrid == true || ishgrid == true){
        // Increase count of number of vector quantities
        nvectors++;

        /*
          Increase the size of the PArray1D object containing pointers to
          Array3D objects containing data.  The data must be copied into temporary
          arrays since the current Array3D objects are destroyed when the
          PArray1D object is resized.  The data is then copied back into the
          resized PArray1D object.

          During the copy process, the source object is deleted immediately after
          its data has been copied.  This is to reduce memory requirements so that
          only one extra object's worth of memory is required during the
          copy process.  Note that setting the resulting pointer to NULL is
          needed to avoid freeing the memory more than once upon deletion of
          the temporary PArray1D object.

          After the copy is complete, add a new Array3D object to the final
          entry in the resized PArray1D object, set its size equal to the grid
          size, and all points initialized to 0.0.
         */
        PArray1D<Array3D<T>*> tmparray(nvectors);
        for(int i=0; i<nvectors-1; i++){
            tmparray(i) = new Array3D<T>;
            tmparray(i)->ResetSize(ny,nx,3);
            for(int kk=00; kk<3; kk++){
                for(int ii=0; ii<ny; ii++){
                    for(int jj=0; jj<nx; jj++){
                        tmparray(i)->operator ()(ii,jj,kk) = pvectors(i)->operator ()(ii,jj,kk);
                    }
                }
            }
            delete pvectors(i);
            pvectors(i) = NULL;
        }
        pvectors.ResetSize(nvectors);
        for(int i=0; i<nvectors-1; i++){
            pvectors(i) = new Array3D<T>;
            pvectors(i)->ResetSize(ny,nx,3);
            for(int kk=0; kk<3; kk++){
                for(int ii=0; ii<ny; ii++){
                    for(int jj=0; jj<nx; jj++){
                        pvectors(i)->operator ()(ii,jj,kk) = tmparray(i)->operator ()(ii,jj,kk);
                    }
                }
            }
            delete tmparray(i);
            tmparray(i) = NULL;
        }
        pvectors(nvectors-1) = new Array3D<T>;
        pvectors(nvectors-1)->ResetSize(ny,nx,3);
        pvectors(nvectors-1)->ResetVal((T)0.0e0);


        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nvectors);
        for(int i=0; i<nvectors-1; i++){
            tmparray2(i) = new std::string;
            tmparray2(i)->reserve(qtysize);
            tmparray2(i)->assign(vector_names(i)->substr());
            delete vector_names(i);
            vector_names(i) = NULL;
        }
        vector_names.ResetSize(nvectors);
        for(int i=0; i<nvectors-1; i++){
            vector_names(i) = new std::string;
            vector_names(i)->reserve(qtysize);
            vector_names(i)->assign(tmparray2(i)->substr());
            delete tmparray2(i);
            tmparray2(i) = NULL;
        }
        vector_names(nvectors-1) = new std::string;
        vector_names(nvectors-1)->reserve(qtysize);
        vector_names(nvectors-1)->assign(name);
    }
}


template <class T>
void Grid2D<T>::RemoveVectorQuantity(const int qty)
{
    // Require that a grid has been generated
    if((isogrid == true || ishgrid == true) && qty < nvectors){
        // Decrease count of number of scalar quantities
        nvectors--;

        /*
          Decrease the size of the PArray1D object containing pointers to
          Array3D objects containing data.  The data must be copied into temporary
          arrays since the current Array3D objects are destroyed when the
          PArray1D object is resized.  The data is then copied back into the
          resized PArray1D object.

          During the copy process, the source object is deleted immediately after
          its data has been copied.  This is to reduce memory requirements so that
          only one extra object's worth of memory is required during the
          copy process.  Note that setting the resulting pointer to NULL is
          needed to avoid freeing the memory more than once upon deletion of
          the temporary PArray1D object.

          After the copy is complete, add a new Array3D object to the final
          entry in the resized PArray1D object, set its size equal to the grid
          size, and all points initialized to 0.0.
         */
        PArray1D<Array3D<T>*> tmparray(nvectors+1);
        for(int i=0; i<nvectors+1; i++){
            tmparray(i) = new Array3D<T>;
            tmparray(i)->ResetSize(ny,nx,3);
            for(int kk=0; kk<3; kk++){
                for(int ii=0; ii<ny; ii++){
                    for(int jj=0; jj<nx; jj++){
                        tmparray(i)->operator ()(ii,jj,kk) = pvectors(i)->operator ()(ii,jj,kk);
                    }
                }
            }
            delete pvectors(i);
            pvectors(i) = NULL;
        }
        pvectors.ResetSize(nvectors);
        for(int i=0; i<nvectors+1; i++){
            int itmp = i;
            if(i != qty){
                if(i > qty){
                    itmp = i - 1;   // Reduce index by 1 since loop is past removed quantity
                }
                pvectors(itmp) = new Array3D<T>;
                pvectors(itmp)->ResetSize(ny,nx,3);
                for(int kk=0; kk<3; kk++){
                    for(int ii=0; ii<ny; ii++){
                        for(int jj=0; jj<nx; jj++){
                            pvectors(itmp)->operator ()(ii,jj,kk) = tmparray(i)->operator ()(ii,jj,kk);
                        }
                    }
                }
            }
            delete tmparray(i);
            tmparray(i) = NULL;
        }


        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nvectors+1);
        for(int i=0; i<nvectors+1; i++){
            tmparray2(i) = new std::string;
            tmparray2(i)->reserve(qtysize);
            tmparray2(i)->assign(vector_names(i)->substr());
            delete vector_names(i);
            vector_names(i) = NULL;
        }
        vector_names.ResetSize(nvectors);
        for(int i=0; i<nvectors+1; i++){
            int itmp = i;
            if(i != qty){
                if(i > qty){
                    itmp = i - 1;
                }
                vector_names(itmp) = new std::string;
                vector_names(itmp)->reserve(qtysize);
                vector_names(itmp)->assign(tmparray2(i)->substr());
                delete tmparray2(i);
                tmparray2(i) = NULL;
            }
        }
    }
}


template <class T>
int Grid2D<T>::GetNumScalars() const
{
    return nscalars;
}


template <class T>
int Grid2D<T>::GetNumVectors() const
{
    return nvectors;
}



// () OPERATOR
template < class T > inline
T& Grid2D<T>::operator()(int ind1, int ind2, int qty)
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 >= 0 && ind1 < ny);
        assert(ind2 >= 0 && ind2 < nx);
        return pscalars(qty)->operator()(ind1,ind2);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pscalars(qty)->operator()(ind1,ind2);
    #endif
}


template < class T > inline
const T& Grid2D<T>::operator()(int ind1, int ind2, int qty) const
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 >= 0 && ind1 < ny);
        assert(ind2 >= 0 && ind2 < nx);
        return pscalars(qty)->operator ()(ind1,ind2);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pscalars(qty)->operator()(ind1,ind2);
    #endif
}


template < class T > inline
T& Grid2D<T>::operator()(int ind1, int ind2, int ind3, int qty)
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 >= 0 && ind1 < ny);
        assert(ind2 >= 0 && ind2 < nx);
        assert(ind3 >= 0 && ind3 < 3);
        return pvectors(qty)->operator ()(ind1,ind2,ind3);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pvectors(qty)->operator ()(ind1,ind2,ind3);
    #endif
}


template < class T > inline
const T& Grid2D<T>::operator()(int ind1, int ind2, int ind3, int qty) const
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 >= 0 && ind1 < ny);
        assert(ind2 >= 0 && ind2 < nx);
        assert(ind3 >= 0 && ind3 < 3);
        return pvectors(qty)->operator ()(ind1,ind2,ind3);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pvectors(qty)->operator ()(ind1,ind2,ind3);
    #endif
}


template <class T>
Grid2D<T>& Grid2D<T>::operator=(const Grid2D<T> a)
{
    Grid2DSwap(*this, a);
    return *this;
}


template <class T>
void Grid2D<T>::setGridName(std::string name)
{
    gridname = name;
}


template <class T>
std::string Grid2D<T>::GridName() const
{
    return gridname;
}


template <class T>
void Grid2D<T>::setScalarQuantityName(const int qty, const std::string name)
{
    scalar_names(qty)->assign(name);
}


template <class T>
std::string Grid2D<T>::ScalarQuantityName(const int qty)
{
    return scalar_names(qty)->substr();
}


template <class T>
void Grid2D<T>::setVectorQuantityName(const int qty, const std::string name)
{
    vector_names(qty)->assign(name);
}


template <class T>
std::string Grid2D<T>::VectorQuantityName(const int qty)
{
    return vector_names(qty)->substr();
}


template <class T>
void Grid2D<T>::CopyNodalCoordinates(Grid2D<T> &sgrid)
{
    /* Require that sgrid has the same dimensions as this grid */
    assert(sgrid.GetSize(2) == nx);
    assert(sgrid.GetSize(1) == ny);

    /* Loop through grid points and copy values */
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            xcoords(i,j) = sgrid.xcoords(i,j);
            ycoords(i,j) = sgrid.ycoords(i,j);
        }
    }

    /* Compute metrics for this grid */
    Grid2D<T>::ComputeMetrics();
}


template <class T>
double Grid2D<T>::EstimateMemoryUsage() const
{
    double retval = 0.0e0;

    /* Tally memory required for data stored at grid points */
    for(int i=0; i<nscalars; i++){
        retval += pscalars(i)->GetMemoryUsage();
    }
    for(int i=0; i<nvectors; i++){
        retval += pvectors(i)->GetMemoryUsage();
    }

    /* Add memory required for member arrays */
    retval += dzidx.GetMemoryUsage();
    retval += dzidy.GetMemoryUsage();
    retval += detadx.GetMemoryUsage();
    retval += detady.GetMemoryUsage();
    retval += xcoords.GetMemoryUsage();
    retval += ycoords.GetMemoryUsage();
    retval += iblank.GetMemoryUsage();
    retval += converge.GetMemoryUsage();

    /* Add memory required for pointers to data stored at grid points */
    retval += pscalars.GetMemoryUsage();
    retval += pvectors.GetMemoryUsage();

    /* Add memory required for member variables */
    retval += 8*sizeof(int);    /* nx, ny, itersused, nscalars, nvectors, qty_size, inc_snapshots, count_snapshots */
    retval += 2*sizeof(bool);   /* isogrid, ishgrid */
    retval += 5*sizeof(T);      /* min/max x, min/max y, pi */

    return retval;
}


template <class T>
void Grid2D<T>::MPI_BCAST(const int src)
{
    /* Only perform MPI operations if USEMPI is defined.  Otherwise, perform no actions */
#ifdef USEMPI
    /* Set MPI datatype */
    MPI::Datatype mpidtype;
    if(typeid(T) == typeid(char)){
        mpidtype = MPI::CHAR;
    }
    if(typeid(T) == typeid(short int)){
        mpidtype = MPI::SHORT;
    }
    if(typeid(T) == typeid(int)){
        mpidtype = MPI::INT;
    }
    if(typeid(T) == typeid(float)){
        mpidtype = MPI::FLOAT;
    }
    if(typeid(T) == typeid(double)){
        mpidtype = MPI::DOUBLE;
    }
    if(typeid(T) == typeid(long double)){
        mpidtype = MPI::LONG_DOUBLE;
    }
    /* Broadcast scalar quantities, ignoring names */
    for(int s=0; s<nscalars; s++){
        MPI::COMM_WORLD.Bcast(&pscalars(s)->operator ()(0,0),nx*ny,mpidtype,src);
    }

    /* Broadcast vector quantities, ignoring names */
    for(int v=0; v<nvectors; v++){
        MPI::COMM_WORLD.Bcast(&pvectors(v)->operator ()(0,0,0),nx*ny*3,mpidtype,src);
    }
#endif
}


template <class T>
void Grid2D<T>::MPI_BCAST_GRID(const int src)
{
    /* Only perform MPI operations if USEMPI is defined.  Otherwise, perform no actions */
#ifdef USEMPI
    /* Set MPI datatype */
    MPI::Datatype mpidtype;
    if(typeid(T) == typeid(char)){
        mpidtype = MPI::CHAR;
    }
    if(typeid(T) == typeid(short int)){
        mpidtype = MPI::SHORT;
    }
    if(typeid(T) == typeid(int)){
        mpidtype = MPI::INT;
    }
    if(typeid(T) == typeid(float)){
        mpidtype = MPI::FLOAT;
    }
    if(typeid(T) == typeid(double)){
        mpidtype = MPI::DOUBLE;
    }
    if(typeid(T) == typeid(long double)){
        mpidtype = MPI::LONG_DOUBLE;
    }

    /* Broadcast nodal coordinates and metrics */
    MPI::COMM_WORLD.Bcast(&dzidx(0,0),nx*ny,mpidtype,src);
    MPI::COMM_WORLD.Bcast(&detadx(0,0),nx*ny,mpidtype,src);
    MPI::COMM_WORLD.Bcast(&dzidy(0,0),nx*ny,mpidtype,src);
    MPI::COMM_WORLD.Bcast(&detady(0,0),nx*ny,mpidtype,src);
    MPI::COMM_WORLD.Bcast(&xcoords(0,0),nx*ny,mpidtype,src);
    MPI::COMM_WORLD.Bcast(&ycoords(0,0),nx*ny,mpidtype,src);
    MPI::COMM_WORLD.Bcast(&iblank(0,0),nx*ny,mpidtype,src);

    /* Broadcast scalar quantities, ignoring names */
    for(int s=0; s<nscalars; s++){
        MPI::COMM_WORLD.Bcast(&pscalars(s)->operator ()(0,0),nx*ny,mpidtype,src);
    }

    /* Broadcast vector quantities, ignoring names */
    for(int v=0; v<nvectors; v++){
        MPI::COMM_WORLD.Bcast(&pvectors(v)->operator ()(0,0,0),nx*ny*3,mpidtype,src);
    }
#endif
}


template <class T>
void Grid2D<T>::MPI_GATHER(const int dest, const int rank, const int np)
{
#ifdef USEMPI
    /* Set MPI datatype */
    MPI::Datatype mpidtype;
    if(typeid(T) == typeid(char)){
        mpidtype = MPI::CHAR;
    }
    if(typeid(T) == typeid(short int)){
        mpidtype = MPI::SHORT;
    }
    if(typeid(T) == typeid(int)){
        mpidtype = MPI::INT;
    }
    if(typeid(T) == typeid(float)){
        mpidtype = MPI::FLOAT;
    }
    if(typeid(T) == typeid(double)){
        mpidtype = MPI::DOUBLE;
    }
    if(typeid(T) == typeid(long double)){
        mpidtype = MPI::LONG_DOUBLE;
    }

    int ls = 0;    int le = 0;  /* Starting and ending values for single-index */
    int i = 0;     int j = 0;
    int idx = 0;

    /* Send scalar quantities to destination process */
    if(nscalars > 0){
        if(rank != dest){
            MPI_Custom::LocalWorkUnits(rank,np,nx*ny,ls,le);
            int nwork = le - ls + 1;
            int npts = nwork*nscalars;
            Array1D<T> tmparray(npts,0.0e0);

            for(int s=0; s<nscalars; s++){
                for(int q=ls; q<=le; q++){
                    MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
                    idx = s*nwork + q - ls;
                    tmparray(idx) = pscalars(s)->operator ()(i,j);
                }
            }

            MPI::COMM_WORLD.Send(&tmparray(0),npts,mpidtype,dest,0);

        } else {
            for(int p=0; p<np; p++){
                if(p != rank){
                    MPI_Custom::LocalWorkUnits(p,np,nx*ny,ls,le);
                    int nwork = le - ls + 1;
                    int npts = nwork*nscalars;
                    Array1D<T> tmparray(npts,0.0e0);

                    MPI::COMM_WORLD.Recv(&tmparray(0),npts,mpidtype,p,0);

                    for(int s=0; s<nscalars; s++){
                        for(int q=ls; q<=le; q++){
                            MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
                            idx = s*nwork + q - ls;
                            pscalars(s)->operator ()(i,j) = tmparray(idx);
                        }
                    }
                }
            }
        }
    }


    /* Send vector quantities to destination process */
    if(nvectors > 0){
        if(rank != dest){
            MPI_Custom::LocalWorkUnits(rank,np,nx*ny,ls,le);
            int nwork = le - ls + 1;
            int npts = nwork*nvectors*3;
            Array1D<T> tmparray(npts,0.0e0);

            for(int v=0; v<nvectors; v++){
                for(int q=ls; q<=le; q++){
                    MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
                    for(int k=0; k<3; k++){
                        idx = v*nwork + (q-ls)*3 + k;
                        tmparray(idx) = pvectors(v)->operator ()(i,j,k);
                    }
                }
            }

            MPI::COMM_WORLD.Send(&tmparray(0),npts,mpidtype,dest,0);

        } else {
            for(int p=0; p<np; p++){
                if(p != rank){
                    MPI_Custom::LocalWorkUnits(p,np,nx*ny,ls,le);
                    int nwork = le - ls + 1;
                    int npts = nwork*nvectors*3;
                    Array1D<T> tmparray(npts,0.0e0);

                    MPI::COMM_WORLD.Recv(&tmparray(0),npts,mpidtype,p,0);

                    for(int v=0; v<nvectors; v++){
                        for(int q=ls; q<=le; q++){
                            MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
                            for(int k=0; k<3; k++){
                                idx = v*nwork + (q-ls)*3 + k;
                                pvectors(v)->operator ()(i,j,k) = tmparray(idx);
                            }
                        }
                    }
                }
            }
        }
    }


#endif
}
