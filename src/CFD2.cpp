/**
 * @file CFD2.cpp
 * @author Robert Grandin
 * @brief Implementation of CFD2 namespace functions.
 */


#include "CFD2.h"

template <class T>
void CFD2::dQdt_FirstOrder(Grid2D<T> &qtygrid, Grid2D<T> &dgrid, const T t)
{
    /* Require that the grids are the same size */
    int nx = qtygrid.GetSize(2);
    int ny = qtygrid.GetSize(1);
    assert(dgrid.GetSize(2) == nx);
    assert(dgrid.GetSize(1) == ny);


    /* Reset dgrid values to 0 */
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            dgrid(i,j,0) = 0.0e0;       /* Density */
            dgrid(i,j,1) = 0.0e0;       /* Total energy */
            dgrid(i,j,2) = 0.0e0;       /* Pressure */
            dgrid(i,j,0,0) = 0.0e0;     /* Velocity, x-direction */
            dgrid(i,j,1,0) = 0.0e0;     /* Velocity, y-direction */
            dgrid(i,j,0,1) = 0.0e0;     /* Momentum, x-direction */
            dgrid(i,j,1,1) = 0.0e0;     /* Momentum, y-direction */
        }
    }

    /* Calculate velocity and pressure from equation quantities */
    T rho = 0.0e0;
    T u = 0.0e0;
    T v = 0.0e0;
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            rho = qtygrid(i,j,0);
            qtygrid(i,j,0,0) = qtygrid(i,j,0,1)/rho;
            qtygrid(i,j,1,0) = qtygrid(i,j,1,1)/rho;
            u = qtygrid(i,j,0,0);
            v = qtygrid(i,j,1,0);
            qtygrid(i,j,2) = 0.4e0*(qtygrid(i,j,1) - 0.5e0*rho*(u*u + v*v));
        }
    }


    /*
      Calculate time-derivatives of cell quantities in parallel.  Inner (across columns) loop
      of the double loop required to traverse the entire grid is placed within dQdt_Parallel().

      This is allowed since the calculation of dgrid(i,j) values only depend upon the calculations
      performed while current cell is qtygrid(i,j).  Each OpenMP thread will have access to both
      qtygrid and dgrid, and race conditions are avoided due to the independence of dgrid values
      (i.e., calculating dgrid(i,j) is independent of other dgrid values).
    */
//    #pragma omp parallel for
//    for(int i=0; i<ny; i++){
//        CFD2::dQdt_Parallel(1,i,qtygrid,dgrid,t);
//    }

} /* End of function */


template <class T>
void CFD2::dQdt_SecondOrder(Grid2D<T> &qtygrid, Grid2D<T> &dgrid, const T t)
{
    /* Require that the grids are the same size */
    int nx = qtygrid.GetSize(2);
    int ny = qtygrid.GetSize(1);
    assert(dgrid.GetSize(2) == nx);
    assert(dgrid.GetSize(1) == ny);

    /* Discover rank number of current MPI process and number of processes */
    int rank = 0;
    int np = 1;
#ifdef USEMPI
    rank = MPI::COMM_WORLD.Get_rank();
    np = MPI::COMM_WORLD.Get_size();

    /* Broadcast current quantity grid to all MPI processes */
    qtygrid.MPI_BCAST(0);
#endif

    /* Reset dgrid values to 0 */
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            dgrid(i,j,0) = 0.0e0;       /* Density */
            dgrid(i,j,1) = 0.0e0;       /* Total energy */
            dgrid(i,j,2) = 0.0e0;       /* Pressure */
            dgrid(i,j,0,0) = 0.0e0;     /* Velocity, x-direction */
            dgrid(i,j,1,0) = 0.0e0;     /* Velocity, y-direction */
            dgrid(i,j,0,1) = 0.0e0;     /* Momentum, x-direction */
            dgrid(i,j,1,1) = 0.0e0;     /* Momentum, y-direction */
        }
    }

    /* Calculate velocity and pressure from equation quantities */
    T rho = 0.0e0;
    T u = 0.0e0;
    T v = 0.0e0;
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            rho = qtygrid(i,j,0);
            qtygrid(i,j,0,0) = qtygrid(i,j,0,1)/rho;
            qtygrid(i,j,1,0) = qtygrid(i,j,1,1)/rho;
            u = qtygrid(i,j,0,0);
            v = qtygrid(i,j,1,0);
            qtygrid(i,j,2) = 0.4e0*(qtygrid(i,j,1) - 0.5e0*rho*(u*u + v*v));
        }
    }


    /*
      Calculate time-derivatives of cell quantities in parallel.  Inner (across columns) loop
      of the double loop required to traverse the entire grid is placed within dQdt_Parallel().

      This is allowed since the calculation of dgrid(i,j) values only depend upon the calculations
      performed while current cell is qtygrid(i,j).  Each OpenMP thread will have access to both
      qtygrid and dgrid, and race conditions are avoided due to the independence of dgrid values
      (i.e., calculating dgrid(i,j) is independent of other dgrid values).
    */
    int qs = 0; int qe = 0;
    MPI_Custom::LocalWorkUnits(rank,np,nx*ny,qs,qe);

    int i = 0;  int j = 0;
    for(int q=qs; q<=qe; q++){
        MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
        if(j > 0 && j < nx-1){
            /* If-block required for my implementation of boundary conditions at inlet and exit */
            CFD2::dQdt_Parallel<T>(2,i,j,qtygrid,dgrid,t);
        }
    }

    /* Gather dgrid values to root process */
#ifdef USEMPI
    dgrid.MPI_GATHER(0,rank,np);
#endif
}


template <class T>
void CFD2::dQdt_SecondOrder_Viscous_MPI(Grid2D<T> &qtygrid, Array3D<T> &dvals, const T t,
                                        Array2D<int> &dependency, const int nrecv,
                                        const T mu, const T kthermal)
{
    /* Send updated values of conserved quantities to other MPI processes which need them */
    int nx = qtygrid.GetSize(2);
    int ny = qtygrid.GetSize(1);
    int qs = -1;
    int qe = -1;
    int nw = nx*ny;
    int rank = 0;
    int np = 1;
    int i = -1;
    int j = -1;

#ifdef USEMPI       /* Only include MPI commands if MPI is to be used */
    rank = MPI::COMM_WORLD.Get_rank();
    np = MPI::COMM_WORLD.Get_size();

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

    T qsend[4] = {0.0e0, 0.0e0, 0.0e0, 0.0e0};
    T qrecv[4] = {0.0e0, 0.0e0, 0.0e0, 0.0e0};

    MPI_Custom::LocalWorkUnits(rank,np,nw,qs,qe);
    for(int q=qs; q<=qe; q++){
        MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
        for(int p=0; p<np; p++){
            if(dependency(q,p) == 1 && p != rank){
                qsend[0] = qtygrid(i,j,0);
                qsend[1] = qtygrid(i,j,1);
                qsend[2] = qtygrid(i,j,2);
                qsend[3] = qtygrid(i,j,3);
                MPI::COMM_WORLD.Send(qsend,4,mpidtype,p,q); /* Tag used to identify work unit */
            }
        }
    }


    /* Receive values of conserved quantities needed by this MPI process */
    MPI::Status recvstatus;
    int q = -1;
    for(int r=0; r<nrecv; r++){
        MPI::COMM_WORLD.Recv(qrecv,4,mpidtype,MPI::ANY_SOURCE,MPI::ANY_TAG,recvstatus); /* Receive from any process */
        q = recvstatus.Get_tag();                   /* Tag used to identify work unit */
        MPI_Custom::Idx_Transfer_2D(q,nx,i,j);      /* Determine (i,j) indices */
        qtygrid(i,j,0) = qrecv[0];
        qtygrid(i,j,1) = qrecv[1];
        qtygrid(i,j,2) = qrecv[2];
        qtygrid(i,j,3) = qrecv[3];
    }

    MPI::COMM_WORLD.Barrier();
#endif


    /* Calculate dQ/dt for conserved quantities */
    dvals.ResetVal(0.0e0);

    if(np == 1){
        /* Use OpenMP to parallelize for multiple cores since only 1 MPI process exists */
        #pragma omp parallel for private(i,j)
        for(int q=qs; q<=qe; q++){
            MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
            if(j > 0 && j < nx-1){
                CFD2::dQdt_Parallel_Viscous(2,i,j,qtygrid,dvals,t,mu,kthermal);
            }
        }

    } else {
        /* No OpenMP in order to avoid over-saturating processors */
        for(int q=qs; q<=qe; q++){
            MPI_Custom::Idx_Transfer_2D(q,nx,i,j);
            if(j > 0 && j < nx-1){
                /* If-block required for my implementation of boundary conditions at inlet and exit */
                CFD2::dQdt_Parallel_Viscous(2,i,j,qtygrid,dvals,t,mu,kthermal);
            }
        }
    }
}



template <class T>
void CFD2::dQdt_ThirdOrder(Grid2D<T> &qtygrid, Grid2D<T> &dgrid, const T t)
{
    /* Require that the grids are the same size */
    int nx = qtygrid.GetSize(2);
    int ny = qtygrid.GetSize(1);
    assert(dgrid.GetSize(2) == nx);
    assert(dgrid.GetSize(1) == ny);


    /* Reset dgrid values to 0 */
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            dgrid(i,j,0) = 0.0e0;       /* Density */
            dgrid(i,j,1) = 0.0e0;       /* Total energy */
            dgrid(i,j,2) = 0.0e0;       /* Pressure */
            dgrid(i,j,0,0) = 0.0e0;     /* Velocity, x-direction */
            dgrid(i,j,1,0) = 0.0e0;     /* Velocity, y-direction */
            dgrid(i,j,0,1) = 0.0e0;     /* Momentum, x-direction */
            dgrid(i,j,1,1) = 0.0e0;     /* Momentum, y-direction */
        }
    }

    /* Calculate velocity and pressure from equation quantities */
    T rho = 0.0e0;
    T u = 0.0e0;
    T v = 0.0e0;
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            rho = qtygrid(i,j,0);
            qtygrid(i,j,0,0) = qtygrid(i,j,0,1)/rho;
            qtygrid(i,j,1,0) = qtygrid(i,j,1,1)/rho;
            u = qtygrid(i,j,0,0);
            v = qtygrid(i,j,1,0);
            qtygrid(i,j,2) = 0.4e0*(qtygrid(i,j,1) - 0.5e0*rho*(u*u + v*v));
        }
    }


    /*
      Calculate time-derivatives of cell quantities in parallel.  Inner (across columns) loop
      of the double loop required to traverse the entire grid is placed within dQdt_Parallel().

      This is allowed since the calculation of dgrid(i,j) values only depend upon the calculations
      performed while current cell is qtygrid(i,j).  Each OpenMP thread will have access to both
      qtygrid and dgrid, and race conditions are avoided due to the independence of dgrid values
      (i.e., calculating dgrid(i,j) is independent of other dgrid values).
    */
//    #pragma omp parallel for
//    for(int i=0; i<ny; i++){
//        CFD2::dQdt_Parallel<T>(3,i,qtygrid,dgrid,t);
//    }
}


template <class T>
void CFD2::RusanovFlux(Array1D<T> &ql, Array1D<T> &qr, Array1D<T> &nvec, Array1D<T> &flux)
{
    /*
      Calculate flux values at interface from quantity values

      Calculating  flux = L (dot) n    where L = (E,F)

      L (dot) n = [ rho*u       ]         + [ rho*v       ]
                  [ rho*u*u + p ] * nx    + [ rho*u*v     ] * ny
                  [ rho*u*v     ]         + [ rho*v*v + p ]
                  [ (Et + p)*u  ]         + [ (Et + p)*v  ]

                = [ rho*vn          ]
                  [ rho*u*vn + p*nx ]   ,     vn = u*nx + v*ny
                  [ rho*v*vn + p*ny ]
                  [ (Et + p)*vn     ]
    */
        T rho = ql(0);
        T u = ql(1)/ql(0);
        T v = ql(2)/ql(0);
        T p = 0.4e0*(ql(3) - 0.5e0*rho*(u*u + v*v));
        T Et = ql(3);
        T vn = u*nvec(0) + v*nvec(1);
        T fl0 = rho*vn;
        T fl1 = rho*u*vn + p*nvec(0);
        T fl2 = rho*v*vn + p*nvec(1);
        T fl3 = vn*(Et + p);

        rho = qr(0);
        u = qr(1)/qr(0);
        v = qr(2)/qr(0);
        p = 0.4e0*(qr(3) - 0.5e0*rho*(u*u + v*v));
        Et = qr(3);
        vn = u*nvec(0) + v*nvec(1);
        T fr0 = rho*vn;
        T fr1 = rho*u*vn + p*nvec(0);
        T fr2 = rho*v*vn + p*nvec(1);
        T fr3 = vn*(Et + p);

        /* Calculate maximum wave speed (flow velocity magnitude + speed of sound) */
        rho = ql(0);
        u = ql(1)/ql(0);
        v = ql(2)/ql(0);
        p = 0.4e0*(ql(3) - 0.5e0*rho*(ql(1)*ql(1) + ql(2)*ql(2)));
        Et = ql(3);
        T c = sqrt(1.4e0*p/rho);            /* Uses current cell's values */
        T lambdamax = sqrt(u*u+v*v) + c;

        /* Calculate flux value */
        flux(0) = 0.5e0*(fl0 + fr0 - lambdamax*(qr(0) - ql(0)));
        flux(1) = 0.5e0*(fl1 + fr1 - lambdamax*(qr(1) - ql(1)));
        flux(2) = 0.5e0*(fl2 + fr2 - lambdamax*(qr(2) - ql(2)));
        flux(3) = 0.5e0*(fl3 + fr3 - lambdamax*(qr(3) - ql(3)));

    /* Dummy variables to assure flux values are not 'nan' (not checked for infinite values, though) */
    #ifndef RELEASE
        T ql0 = ql(0);
        T ql1 = ql(1);
        T ql2 = ql(2);
        T ql3 = ql(3);
        T qr0 = qr(0);
        T qr1 = qr(1);
        T qr2 = qr(2);
        T qr3 = qr(3);
        T f0 = flux(0);
        T f1 = flux(1);
        T f2 = flux(2);
        T f3 = flux(3);
        assert(f0 == f0 && f1 == f1 && f2 == f2 && f3 == f3);
    #endif
}


template <class T>
void CFD2::RusanovFlux_Viscous(Array1D<T> &ql, Array1D<T> &qr, Array1D<T> &nvec, Array1D<T> &flux)
{
    /*
      Calculate flux values at interface from quantity values

      Calculating  flux = L (dot) n    where L = (E,F)

      L (dot) n = [ rho*u       ]         + [ rho*v       ]
                  [ rho*u*u + p ] * nx    + [ rho*u*v     ] * ny
                  [ rho*u*v     ]         + [ rho*v*v + p ]
                  [ (Et + p)*u  ]         + [ (Et + p)*v  ]

                = [ rho*vn          ]
                  [ rho*u*vn + p*nx ]   ,     vn = u*nx + v*ny
                  [ rho*v*vn + p*ny ]
                  [ (Et + p)*vn     ]
    */
        T rho = ql(0);
        T u = ql(1)/ql(0);
        T v = ql(2)/ql(0);
        T p = 0.4e0*(ql(3) - 0.5e0*rho*(u*u + v*v));
        T Et = ql(3);
        T vn = u*nvec(0) + v*nvec(1);
        T fl0 = rho*vn;
        T fl1 = rho*u*vn + p*nvec(0);
        T fl2 = rho*v*vn + p*nvec(1);
        T fl3 = vn*(Et + p);

        rho = qr(0);
        u = qr(1)/qr(0);
        v = qr(2)/qr(0);
        p = 0.4e0*(qr(3) - 0.5e0*rho*(u*u + v*v));
        Et = qr(3);
        vn = u*nvec(0) + v*nvec(1);
        T fr0 = rho*vn;
        T fr1 = rho*u*vn + p*nvec(0);
        T fr2 = rho*v*vn + p*nvec(1);
        T fr3 = vn*(Et + p);

        /* Calculate maximum wave speed (flow velocity magnitude + speed of sound) */
        rho = ql(0);
        u = ql(1)/ql(0);
        v = ql(2)/ql(0);
        vn = u*nvec(0) + v*nvec(1);
        p = 0.4e0*(ql(3) - 0.5e0*rho*(ql(1)*ql(1) + ql(2)*ql(2)));
        Et = ql(3);
        T c = sqrt(1.4e0*p/rho);            /* Uses current cell's values */
        T lambdamax = fabs(vn) + c;

        /* Calculate flux value */
        flux(0) = 0.5e0*(fl0 + fr0 - lambdamax*(qr(0) - ql(0)));
        flux(1) = 0.5e0*(fl1 + fr1 - lambdamax*(qr(1) - ql(1)));
        flux(2) = 0.5e0*(fl2 + fr2 - lambdamax*(qr(2) - ql(2)));
        flux(3) = 0.5e0*(fl3 + fr3 - lambdamax*(qr(3) - ql(3)));

    /* Dummy variables to assure flux values are not 'nan' (not checked for infinite values, though) */
    #ifndef RELEASE
        T ql0 = ql(0);
        T ql1 = ql(1);
        T ql2 = ql(2);
        T ql3 = ql(3);
        T qr0 = qr(0);
        T qr1 = qr(1);
        T qr2 = qr(2);
        T qr3 = qr(3);
        T f0 = flux(0);
        T f1 = flux(1);
        T f2 = flux(2);
        T f3 = flux(3);
        assert(f0 == f0 && f1 == f1 && f2 == f2 && f3 == f3);
    #endif
}


template <class T>
void CFD2::Reconstruction(int order, Grid2D<T> &qtygrid, int i, int j,
                          T di, T dj, Array1D<T> &nvec, Array1D<T> &qvec)
{
    /* Require that only 0 or one interface direction is specified.  Values of -0.5, 0, or 0.5 are expected for each. */
    assert((fabs(di) > 0.4e0 && fabs(dj) < 0.4e0) || (fabs(di) < 0.4e0 && fabs(dj) > 0.4e0) ||
           (fabs(di) < 0.4e0 && fabs(dj) < 0.4e0));

    /* Require that qvec be of size 4.  Values are only assigned, not retrieved, so initializing to 0 isn't required */
    assert(qvec.GetDim() == 4);

    /* Get grid size */
    int nx = qtygrid.GetSize(2);
    int ny = qtygrid.GetSize(1);

    /* Determine sign of direction to interface */
    int isigni = CFD2::sgn<T>(di);
    int isignj = CFD2::sgn<T>(dj);

    /* Over-ride if indices are outside the grid */
    int iorig = i;
    int jorig = j;
    if(i <= 0){
        i = 0;
    }
    if(i >= ny){
        i = ny-1;
    }
    if(j <= 0){
        j = 0;
    }
    if(j >= nx){
        j = nx-1;
    }

    /* Determine (i,j) indices for the surrounding points */
    int im = i - abs(isigni);
    int jm = j - abs(isignj);
    int im2 = i - 2*abs(isigni);
    int jm2 = j - 2*abs(isignj);
    int ip = i + abs(isigni);
    int jp = j + abs(isignj);

    /* Override if indices are outside the grid */
    if(im <= 0){
        im = 0;
    }
    if(im2 <= 0){
        im2 = 0;
    }
    if(ip >= ny){
        ip = ny-1;
    }
    if(jm <= 0){
        jm = 0;
    }
    if(jm2 <= 0){
        jm2 = 0;
    }
    if(jp >= nx){
        jp = nx-1;
    }

    /* Declare variables needed for reconstruction */
    T rho = 0.0e0;       /* For "central" point of expansion */
    T Et = 0.0e0;
    T xmom = 0.0e0;
    T ymom = 0.0e0;
    T rhom = 0.0e0;      /* For "i-1" or "j-1" point of expansion ("m" for minus) */
    T Etm = 0.0e0;
    T xmomm = 0.0e0;
    T ymomm = 0.0e0;
    T rhom2 = 0.0e0;     /* For "i-2" or "j-2" point of expansion ("m2" for minus two) */
    T Etm2 = 0.0e0;
    T xmomm2 = 0.0e0;
    T ymomm2 = 0.0e0;
    T rhop = 0.0e0;      /* For "i+1" or "j+1" point of expansion ("p" for plus) */
    T Etp = 0.0e0;
    T xmomp = 0.0e0;
    T ymomp = 0.0e0;

    /* Assign values to primitive variables surrounding desired point */
    rho = qtygrid(i,j,0);
    Et = qtygrid(i,j,1);
    xmom = qtygrid(i,j,0,1);
    ymom = qtygrid(i,j,1,1);
    rhom = qtygrid(im,jm,0);
    Etm = qtygrid(im,jm,1);
    xmomm = qtygrid(im,jm,0,1);
    ymomm = qtygrid(im,jm,1,1);
    rhom2 = qtygrid(im2,jm2,0);
    Etm2 = qtygrid(im2,jm2,1);
    xmomm2 = qtygrid(im2,jm2,0,1);
    ymomm2 = qtygrid(im2,jm2,1,1);
    rhop = qtygrid(ip,jp,0);
    Etp = qtygrid(ip,jp,1);
    xmomp = qtygrid(ip,jp,0,1);
    ymomp = qtygrid(ip,jp,1,1);

    /* Correct values as-needed due to boundaries */
    if(i == 0 && isigni != 0){     /* Assumed symmetry axis */
        rhom = rho;
        Etm = Et;
        xmomm = xmom;
        ymomm = -ymom;
        rhom2 = rho;
        Etm2 = Etm;
        xmomm2 = xmomm;
        ymomm2 = ymomm;
    }
    if(i == ny-1 && isigni != 0){   /* Assumed impermeable slip-wall */
        rhop = rho;
        Etp = Et;
        T vn = (xmom*nvec(0) + ymom*nvec(1))/rho;
        xmomp = xmom - 2.0e0*vn*nvec(0)*rho;
        ymomp = ymom - 2.0e0*vn*nvec(1)*rho;
    }
    if(j == 1 && isignj != 0){      /* Inlet, use one-sided difference */
        rhom = rho;
        Etm = Et;
        xmomm = xmom;
        ymomm = ymom;
        rhom2 = rhom;
        Etm2 = Etm;
        xmomm2 = xmomm;
        ymomm2 = ymomm;
    }
    if(j == 0 && isignj != 0){      /* Inlet, use analytic solution stored in first column of grid */
        rhom = rho;
        Etm = Et;
        xmomm = xmom;
        ymomm = ymom;
        rhom2 = rhom;
        Etm2 = Etm;
        xmomm2 = xmomm;
        ymomm2 = ymomm;
        rhop = rho;
        Etp = Et;
        xmomp = xmom;
        ymomp = ymom;
    }
    if(j == nx-2 && isignj != 0){   /* Exit, use one-sided difference */
        rhop = rho;
        Etp = Et;
        xmomp = xmom;
        ymomp = ymom;
    }
    if(j == nx-1 && isignj != 0){   /* Exit, use analytic solution stored in last column of grid */
        rhop = rho;
        Etp = Et;
        xmomp = xmom;
        ymomp = ymom;
        rhom = rho;
        Etm = Et;
        xmomm = xmom;
        ymomm = ymom;
        rhom2 = rho;
        Etm2 = Et;
        xmomm2 = xmom;
        ymomm2 = ymom;
    }
    if(iorig >= ny){    /* Outside assumed wall */
        rho = rhop;
        Et = Etp;
        xmom = xmomp;
        ymom = ymomp;
    }
    if(iorig < 0){      /* Outside assumed symmetry boundary */
        rho = rhom;
        Et = Etm;
        xmom = xmomm;
        ymom = ymomm;
    }
    if(jorig >= nx){    /* Exit conditions */
        rho = rhop;
        Et = Etp;
        xmom = xmomp;
        ymom = ymomp;
    }
    if(jorig < 0){      /* Inlet conditions */
        rho = rhom;
        Et = Etm;
        xmom = xmomm;
        ymom = ymomm;
    }

    /*
      At this point, the flow variables are known in the current cell, as well as the cells on either side
      in the desired reconstruction direction.  The MUSCL scheme requires these values, so now MUSCL
      reconstruction can be used for either direction and at any grid point.

      The MinMod() flux limiter is multiplied by the 2nd and 3rd terms in the MUSCL reconstruction.  The
      flux limiter has no effect for first-order reconstructions as these two terms in the MUSCL scheme
      are both 0.
    */

    /* Set MUSCL parameters for reconstruction */
    T kappa = 0.0e0;    /* kappa term in MUSCL reconstruction */
    T coeff = 0.0e0;    /* Set to zero and combine with kappa=0 for first-order */
    switch(order){
    case 1:     /* First-order reconstruction */
        kappa = 0.0e0;
        coeff = 0.0e0;
        break;
    case 2:     /* Second-order reconstruction */
        kappa = 0.0e0;
        coeff = 1.0e0;
        break;
    case 3:     /* Third-order reconstruction */
        kappa = 1.0e0/3.0e0;
        coeff = 1.0e0;
        if(j == 1 || j == nx-2){    /* Override kappa value to force linear approximation of flux */
            kappa = 0.0e0;
        }
        break;
    default:
        std::cerr << "WARNING: Undefined order specified in CFD2::Reconstruction<T>()" << std::endl;
        break;
    }



    /* Calculate limiter value */
    T limiter = 1.0e0;
    if(order > 1){
        /* Only calculate for 2nd and 3rd order since limiter has no effect on 1st order */

        T mmx = (rho - rhom)*(rho - rhom) + (xmom - xmomm)*(xmom - xmomm) + (ymom - ymomm)*(ymom-ymomm) + (Et - Etm)*(Et - Etm);
        T mmy = (rhom - rhom2)*(rhom - rhom2) + (xmomm - xmomm2)*(xmomm - xmomm2) + (ymomm - ymomm2)*(ymomm-ymomm2) + (Etm - Etm2)*(Etm - Etm2);
        mmx = sqrt(mmx);
        mmy = sqrt(mmy);
        limiter = CFD2::MinMod<T>(mmx,mmy);
    }


    /* Perform the reconstruction. */
    T deltax = 1.0e0;   /* Unit grid spacing.  Interpolation is in computational domain, not physical. */
    T dx = 0.0e0;       /* This should be -0.5, or +0.5, depending on the interface at which we are reconstructing. */
    if(isigni != 0){
        dx = di;
    } else {
        dx = dj;
    }

    /* Density */
    T q = rho;
    T qm = rhom;
    T qp = rhop;
    T qm2 = rhom2;
//    T mmx = q - qm;
//    T mmy = qm - qm2;
//    qvec(0) = q + CFD2::MinMod<T>(mmx,mmy)*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
//            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));
    qvec(0) = q + limiter*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));

    /* X-momentum */
    q = xmom;
    qm = xmomm;
    qp = xmomp;
    qm2 = xmomm2;
//    mmx = q - qm;
//    mmy = qm - qm2;
//    qvec(1) = q + CFD2::MinMod<T>(mmx,mmy)*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
//            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));
    qvec(1) = q + limiter*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));

    /* Y-momentum */
    q = ymom;
    qm = ymomm;
    qp = ymomp;
    qm2 = ymomm2;
//    mmx = q - qm;
//    mmy = qm - qm2;
//    qvec(2) = q + CFD2::MinMod<T>(mmx,mmy)*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
//            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));
    qvec(2) = q + limiter*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));

    /* Energy */
    q = Et;
    qm = Etm;
    qp = Etp;
    qm2 = Etm2;
//    mmx = q - qm;
//    mmy = qm - qm2;
//    qvec(3) = q + CFD2::MinMod<T>(mmx,mmy)*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
//            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));
    qvec(3) = q + limiter*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));


    /* Override values for end points */
    if((j == 0 || j == nx-1) && isignj != 0){
        qvec(0) = qtygrid(i,j,0);
        qvec(1) = qtygrid(i,j,0,1);
        qvec(2) = qtygrid(i,j,1,1);
        qvec(3) = qtygrid(i,j,1);
    }
    if(j == 1 && isignj != 0){
        qvec(0) = qtygrid(i,j,0) - 0.5e0*(qtygrid(i,j+1,0) - qtygrid(i,j,0));
        qvec(1) = qtygrid(i,j,0,1) - 0.5e0*(qtygrid(i,j+1,0,1) - qtygrid(i,j,0,1));
        qvec(2) = qtygrid(i,j,1,1) - 0.5e0*(qtygrid(i,j+1,1,1) - qtygrid(i,j,1,1));
        qvec(3) = qtygrid(i,j,1) - 0.5e0*(qtygrid(i,j+1,1) - qtygrid(i,j,1));
    }
    if(j == nx-2 && isignj != 0){
        qvec(0) = qtygrid(i,j,0) + 0.5e0*(qtygrid(i,j,0) - qtygrid(i,j-1,0));
        qvec(1) = qtygrid(i,j,0,1) + 0.5e0*(qtygrid(i,j,0,1) - qtygrid(i,j-1,0,1));
        qvec(2) = qtygrid(i,j,1,1) + 0.5e0*(qtygrid(i,j,1,1) - qtygrid(i,j-1,1,1));
        qvec(3) = qtygrid(i,j,1) + 0.5e0*(qtygrid(i,j,1) - qtygrid(i,j-1,1));
    }

    #ifndef RELEASE
        T q0 = qvec(0);
        T q1 = qvec(1);
        T q2 = qvec(2);
        T q3 = qvec(3);
        assert(q0 == q0 && q1 == q1 && q2 == q2 && q3 == q3);
    #endif
}


template <class T>
void CFD2::Reconstruction_v2(int order, Grid2D<T> &qtygrid, int i, int j,
                          T di, T dj, Array1D<T> &nvec, Array1D<T> &qvec)
{
    /* Require that only 0 or one interface direction is specified.  Values of -0.5, 0, or 0.5 are expected for each. */
//    assert((fabs(di) > 0.4e0 && fabs(dj) < 0.4e0) || (fabs(di) < 0.4e0 && fabs(dj) > 0.4e0) ||
//           (fabs(di) < 0.4e0 && fabs(dj) < 0.4e0));

    /* Require that qvec be of size 4.  Values are only assigned, not retrieved, so initializing to 0 isn't required */
    assert(qvec.GetDim() == 4);

    /* Get grid size */
    int nx = qtygrid.GetSize(2);
    int ny = qtygrid.GetSize(1);

    /* Determine sign of direction to interface */
    int isigni = CFD2::sgn<T>(di);
    int isignj = CFD2::sgn<T>(dj);

    /* Over-ride if indices are outside the grid */
    int iorig = i;
    int jorig = j;
    if(i <= 0){
        i = 0;
    }
    if(i >= ny){
        i = ny-1;
    }
    if(j <= 0){
        j = 0;
    }
    if(j >= nx){
        j = nx-1;
    }

    /* Determine (i,j) indices for the surrounding points */
    int im = i - abs(isigni);
    int jm = j - abs(isignj);
    int im2 = i - 2*abs(isigni);
    int jm2 = j - 2*abs(isignj);
    int ip = i + abs(isigni);
    int jp = j + abs(isignj);

    /* Override if indices are outside the grid */
    if(im <= 0){
        im = 0;
    }
    if(im2 <= 0){
        im2 = 0;
    }
    if(ip >= ny){
        ip = ny-1;
    }
    if(jm <= 0){
        jm = 0;
    }
    if(jm2 <= 0){
        jm2 = 0;
    }
    if(jp >= nx){
        jp = nx-1;
    }

    /* Declare variables needed for reconstruction */
    T rho = 0.0e0;       /* For "central" point of expansion */
    T Et = 0.0e0;
    T xmom = 0.0e0;
    T ymom = 0.0e0;
    T rhom = 0.0e0;      /* For "i-1" or "j-1" point of expansion ("m" for minus) */
    T Etm = 0.0e0;
    T xmomm = 0.0e0;
    T ymomm = 0.0e0;
    T rhom2 = 0.0e0;     /* For "i-2" or "j-2" point of expansion ("m2" for minus two) */
    T Etm2 = 0.0e0;
    T xmomm2 = 0.0e0;
    T ymomm2 = 0.0e0;
    T rhop = 0.0e0;      /* For "i+1" or "j+1" point of expansion ("p" for plus) */
    T Etp = 0.0e0;
    T xmomp = 0.0e0;
    T ymomp = 0.0e0;

    /* Assign values to primitive variables surrounding desired point */
    rho = qtygrid(i,j,0);
    Et = qtygrid(i,j,3);
    xmom = qtygrid(i,j,1);
    ymom = qtygrid(i,j,2);
    rhom = qtygrid(im,jm,0);
    Etm = qtygrid(im,jm,3);
    xmomm = qtygrid(im,jm,1);
    ymomm = qtygrid(im,jm,2);
    rhom2 = qtygrid(im2,jm2,0);
    Etm2 = qtygrid(im2,jm2,3);
    xmomm2 = qtygrid(im2,jm2,1);
    ymomm2 = qtygrid(im2,jm2,2);
    rhop = qtygrid(ip,jp,0);
    Etp = qtygrid(ip,jp,3);
    xmomp = qtygrid(ip,jp,1);
    ymomp = qtygrid(ip,jp,2);

    /* Correct values as-needed due to boundaries */
    if(i == 0 && isigni != 0){     /* Assumed symmetry axis */
        rhom = rho;
        Etm = Et;
        xmomm = xmom;
        ymomm = -ymom;
        rhom2 = rho;
        Etm2 = Etm;
        xmomm2 = xmomm;
        ymomm2 = ymomm;
    }
    if(i == ny-1 && isigni != 0){   /* Assumed impermeable nonslip-wall */
        rhop = rho;
        Etp = Et;
        xmomp = -xmom;
        ymomp = -ymom;
    }
    if(j == 1 && isignj != 0){      /* Inlet, use one-sided difference */
        rhom = rho;
        Etm = Et;
        xmomm = xmom;
        ymomm = ymom;
        rhom2 = rhom;
        Etm2 = Etm;
        xmomm2 = xmomm;
        ymomm2 = ymomm;
    }
    if(j == 0 && isignj != 0){      /* Inlet, use analytic solution stored in first column of grid */
        rhom = rho;
        Etm = Et;
        xmomm = xmom;
        ymomm = ymom;
        rhom2 = rhom;
        Etm2 = Etm;
        xmomm2 = xmomm;
        ymomm2 = ymomm;
        rhop = rho;
        Etp = Et;
        xmomp = xmom;
        ymomp = ymom;
    }
    if(j == nx-2 && isignj != 0){   /* Exit, use one-sided difference */
        rhop = rho;
        Etp = Et;
        xmomp = xmom;
        ymomp = ymom;
    }
    if(j == nx-1 && isignj != 0){   /* Exit, use analytic solution stored in last column of grid */
        rhop = rho;
        Etp = Et;
        xmomp = xmom;
        ymomp = ymom;
        rhom = rho;
        Etm = Et;
        xmomm = xmom;
        ymomm = ymom;
        rhom2 = rho;
        Etm2 = Et;
        xmomm2 = xmom;
        ymomm2 = ymom;
    }
    if(iorig >= ny){    /* Outside assumed wall */
        rho = rhop;
        Et = Etp;
        xmom = xmomp;
        ymom = ymomp;
    }
    if(iorig < 0){      /* Outside assumed symmetry boundary */
        rho = rhom;
        Et = Etm;
        xmom = xmomm;
        ymom = ymomm;
    }
    if(jorig >= nx){    /* Exit conditions */
        rho = rhop;
        Et = Etp;
        xmom = xmomp;
        ymom = ymomp;
    }
    if(jorig < 0){      /* Inlet conditions */
        rho = rhom;
        Et = Etm;
        xmom = xmomm;
        ymom = ymomm;
    }

    /*
      At this point, the flow variables are known in the current cell, as well as the cells on either side
      in the desired reconstruction direction.  The MUSCL scheme requires these values, so now MUSCL
      reconstruction can be used for either direction and at any grid point.

      The MinMod() flux limiter is multiplied by the 2nd and 3rd terms in the MUSCL reconstruction.  The
      flux limiter has no effect for first-order reconstructions as these two terms in the MUSCL scheme
      are both 0.
    */

    /* Set MUSCL parameters for reconstruction */
    T kappa = 0.0e0;    /* kappa term in MUSCL reconstruction */
    T coeff = 0.0e0;    /* Set to zero and combine with kappa=0 for first-order */
    switch(order){
    case 1:     /* First-order reconstruction */
        kappa = 0.0e0;
        coeff = 0.0e0;
        break;
    case 2:     /* Second-order reconstruction */
        kappa = 0.0e0;
        coeff = 1.0e0;
        break;
    case 3:     /* Third-order reconstruction */
        kappa = 1.0e0/3.0e0;
        coeff = 1.0e0;
        if(j == 1 || j == nx-2){    /* Override kappa value to force linear approximation of flux */
            kappa = 0.0e0;
        }
        break;
    default:
        std::cerr << "WARNING: Undefined order specified in CFD2::Reconstruction<T>()" << std::endl;
        break;
    }



    /* Calculate limiter value */
    T limiter = 1.0e0;
//    if(order > 1){
//        /* Only calculate for 2nd and 3rd order since limiter has no effect on 1st order */

//        T mmx = (rho - rhom)*(rho - rhom) + (xmom - xmomm)*(xmom - xmomm) + (ymom - ymomm)*(ymom-ymomm) + (Et - Etm)*(Et - Etm);
//        T mmy = (rhom - rhom2)*(rhom - rhom2) + (xmomm - xmomm2)*(xmomm - xmomm2) + (ymomm - ymomm2)*(ymomm-ymomm2) + (Etm - Etm2)*(Etm - Etm2);
//        mmx = sqrt(mmx);
//        mmy = sqrt(mmy);
//        limiter = CFD2::MinMod<T>(mmx,mmy);
//    }


    /* Perform the reconstruction. */
    T deltax = 1.0e0;   /* Unit grid spacing.  Interpolation is in computational domain, not physical. */
    T dx = 0.0e0;       /* This should be -0.5, or +0.5, depending on the interface at which we are reconstructing. */
    if(isigni != 0){
        dx = di;
    } else {
        dx = dj;
    }

    /* Density */
    T q = rho;
    T qm = rhom;
    T qp = rhop;
    qvec(0) = q + limiter*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));

    /* X-momentum */
    q = xmom;
    qm = xmomm;
    qp = xmomp;
    qvec(1) = q + limiter*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));

    /* Y-momentum */
    q = ymom;
    qm = ymomm;
    qp = ymomp;
    qvec(2) = q + limiter*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));

    /* Energy */
    q = Et;
    qm = Etm;
    qp = Etp;
    qvec(3) = q + limiter*(coeff/(2.0e0*deltax)*(qp - qm)*dx +
            3.0e0*kappa/(2.0e0*deltax*deltax)*(qp - 2.0e0*q + qm)*(dx*dx - deltax*deltax/12.0e0));


    /* Override values for end points */
    if((j == 0 || j == nx-1) && isignj != 0){
        qvec(0) = qtygrid(i,j,0);
        qvec(1) = qtygrid(i,j,1);
        qvec(2) = qtygrid(i,j,2);
        qvec(3) = qtygrid(i,j,3);
    }
    if(j == 1 && isignj != 0){
        qvec(0) = qtygrid(i,j,0) - 0.5e0*(qtygrid(i,j+1,0) - qtygrid(i,j,0));
        qvec(1) = qtygrid(i,j,1) - 0.5e0*(qtygrid(i,j+1,1) - qtygrid(i,j,1));
        qvec(2) = qtygrid(i,j,2) - 0.5e0*(qtygrid(i,j+1,2) - qtygrid(i,j,2));
        qvec(3) = qtygrid(i,j,3) - 0.5e0*(qtygrid(i,j+1,3) - qtygrid(i,j,3));
    }
    if(j == nx-2 && isignj != 0){
        qvec(0) = qtygrid(i,j,0) + 0.5e0*(qtygrid(i,j,0) - qtygrid(i,j-1,0));
        qvec(1) = qtygrid(i,j,1) + 0.5e0*(qtygrid(i,j,1) - qtygrid(i,j-1,1));
        qvec(2) = qtygrid(i,j,2) + 0.5e0*(qtygrid(i,j,2) - qtygrid(i,j-1,2));
        qvec(3) = qtygrid(i,j,3) + 0.5e0*(qtygrid(i,j,3) - qtygrid(i,j-1,3));
    }

    #ifndef RELEASE
        T q0 = qvec(0);
        T q1 = qvec(1);
        T q2 = qvec(2);
        T q3 = qvec(3);
        assert(q0 == q0 && q1 == q1 && q2 == q2 && q3 == q3);
    #endif
}


template <class T>
T CFD2::MinMod(T x, T y)
{
    T retval = 0.0e0;
    if(x*y >= 0.0e0){
        T  minval = fabs(x);    /* This is faster than using fmin(fabs(x),fabs(y)) */
        if(fabs(y) < minval){
            minval = fabs(y);
        }
        retval = CFD2::sgn<T>(x)*minval;
    }
    return retval;
}


template <class T>
int CFD2::sgn(T val) {
    /* Code source: http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c/4609795#4609795 */
    return (val > T(0)) - (val < T(0));
}


template <class T>
void CFD2::dQdt_Parallel(const int order, const int i, const int j, Grid2D<T> &qtygrid, Grid2D<T> &dgrid, const T t)
{
    /* These lines done to avoid "unused variables" compiler warnings */
    T tloc = t + 1.0e0;
    tloc *= 2.0e0;

    /*
      Loop through grid points to calculate derivative values.  As stated in the
      header documentation, the first and last columns in the grid are not modified
      so that they function as boundary conditions.
    */
    int nx = qtygrid.GetSize(2);
    int ny = qtygrid.GetSize(1);
    T volume = 0.0e0;
    T ax = 0.0e0; T ay = 0.0e0;
    T bx = 0.0e0; T by = 0.0e0;
    T cx = 0.0e0; T cy = 0.0e0;
    T dx = 0.0e0; T dy = 0.0e0;
    T Sab = 0.0e0;
    T Sbc = 0.0e0;
    T Scd = 0.0e0;
    T Sda = 0.0e0;
    T nnx = 0.0e0;
    T nny = 0.0e0;
    T len = 0.0e0;
    T di = 0.0e0;
    T dj = 0.0e0;
    T drho = 0.0e0;
    T dEt = 0.0e0;
    T dxmom = 0.0e0;
    T dymom = 0.0e0;
    Array1D<T> ql(4,0.0e0);
    Array1D<T> qr(4,0.0e0);
    Array1D<T> nvec(2,0.0e0);
    Array1D<T> flux(4,0.0e0);

    /* Debugging variables */
    #ifndef RELEASE
        T rhonew = 0.0e0;
        T Etnew = 0.0e0;
        T xmomnew = 0.0e0;
        T ymomnew = 0.0e0;
        T f0 = 0.0e0;
        T f1 = 0.0e0;
        T f2 = 0.0e0;
        T f3 = 0.0e0;
        T ql0 = 0.0e0;
        T ql1 = 0.0e0;
        T ql2 = 0.0e0;
        T ql3 = 0.0e0;
        T qr0 = 0.0e0;
        T qr1 = 0.0e0;
        T qr2 = 0.0e0;
        T qr3 = 0.0e0;
    #endif

    /*
       Calculate volume of cell by finding the magnitude of the cross-product
       of the vectors which span the opposing corners of the quadrilateral element.

       Cell corners are found using bilinear interpolation with the surrounding
       grid points.

       Cell corners are lettered a-d, starting in the upper-right corner and
       moving CCW.  The gridpoint is at position 'N' in the sketch below.

            b ------- a
            |         |
            |    N    |
            |         |
            c ------- d

        Normal vectors to cell sides are assumed to be outward-positive.

    */

    CFD2::CalculateCellVertices<T>(qtygrid,i,j,ax,ay,bx,by,cx,cy,dx,dy);

    volume = fabs((ax-cx)*(dy-by) - (dx-bx)*(ay-cy));

    /* Calculate side lengths of cell */
    Sab = sqrt((ax-bx)*(ax-bx) + (ay-by)*(ay-by));
    Sbc = sqrt((bx-cx)*(bx-cx) + (by-cy)*(by-cy));
    Scd = sqrt((cx-dx)*(cx-dx) + (cy-dy)*(cy-dy));
    Sda = sqrt((dx-ax)*(dx-ax) + (dy-ay)*(dy-ay));


    /*
      For each face, calculate flux through the face and add the flux to the
      time-derivative value.
    */

    /* Face AB */
    di = 0.5e0;
    dj = 0.0e0;
    nnx = by - ay;
    nny = ax - bx;
    len = sqrt(nnx*nnx + nny*nny);
    nnx /= len;
    nny /= len;
    nvec(0) = nnx;
    nvec(1) = nny;
    CFD2::Reconstruction<T>(order,qtygrid,i,j,di,dj,nvec,ql);
    di = -0.5e0;
    CFD2::Reconstruction<T>(order,qtygrid,i+1,j,di,dj,nvec,qr);

    #ifndef RELEASE
        ql0 = ql(0);
        ql1 = ql(1);
        ql2 = ql(2);
        ql3 = ql(3);
        qr0 = qr(0);
        qr1 = qr(1);
        qr2 = qr(2);
        qr3 = qr(3);
        assert(ql0 == ql0 && ql1 == ql1 && ql2 == ql2 && ql3 == ql3);
        assert(qr0 == qr0 && qr1 == qr1 && qr2 == qr2 && qr3 == qr3);
    #endif

    CFD2::RusanovFlux<T>(ql,qr,nvec,flux);

    drho = flux(0);
    dEt = flux(3);
    dxmom = flux(1);
    dymom = flux(2);

    #ifndef RELEASE
        f0 = flux(0);
        f1 = flux(1);
        f2 = flux(2);
        f3 = flux(3);
        assert(f0 == f0 && f1 == f1 && f2 == f2 && f3 == f3);
        rhonew = qtygrid(i,j,0) + drho;
        Etnew = qtygrid(i,j,1) + dEt;
        xmomnew = qtygrid(i,j,0,1) + dxmom;
        ymomnew = qtygrid(i,j,1,1) + dymom;
        assert(drho == drho && dxmom == dxmom && dymom == dymom && dEt == dEt);
    #endif

    dgrid(i,j,0) -= drho*Sab;
    dgrid(i,j,1) -= dEt*Sab;
    dgrid(i,j,0,1) -= dxmom*Sab;
    dgrid(i,j,1,1) -= dymom*Sab;





    /* Face BC */
    di = 0.0e0;
    dj = -0.5e0;
    nnx = cy - by;
    nny = bx - cx;
    len = sqrt(nnx*nnx + nny*nny);
    nnx /= len;
    nny /= len;
    nvec(0) = nnx;
    nvec(1) = nny;
    CFD2::Reconstruction<T>(order,qtygrid,i,j,di,dj,nvec,ql);

    dj = 0.5e0;
    CFD2::Reconstruction<T>(order,qtygrid,i,j-1,di,dj,nvec,qr);

    #ifndef RELEASE
        ql0 = ql(0);
        ql1 = ql(1);
        ql2 = ql(2);
        ql3 = ql(3);
        qr0 = qr(0);
        qr1 = qr(1);
        qr2 = qr(2);
        qr3 = qr(3);
        assert(ql0 == ql0 && ql1 == ql1 && ql2 == ql2 && ql3 == ql3);
        assert(qr0 == qr0 && qr1 == qr1 && qr2 == qr2 && qr3 == qr3);
    #endif

    CFD2::RusanovFlux<T>(ql,qr,nvec,flux);

    drho = flux(0);
    dEt = flux(3);
    dxmom = flux(1);
    dymom = flux(2);

    #ifndef RELEASE
        f0 = flux(0);
        f1 = flux(1);
        f2 = flux(2);
        f3 = flux(3);
        assert(f0 == f0 && f1 == f1 && f2 == f2 && f3 == f3);
        rhonew = qtygrid(i,j,0) + drho;
        Etnew = qtygrid(i,j,1) + dEt;
        xmomnew = qtygrid(i,j,0,1) + dxmom;
        ymomnew = qtygrid(i,j,1,1) + dymom;
        assert(drho == drho && dxmom == dxmom && dymom == dymom && dEt == dEt);
    #endif

    dgrid(i,j,0) -= drho*Sbc;
    dgrid(i,j,1) -= dEt*Sbc;
    dgrid(i,j,0,1) -= dxmom*Sbc;
    dgrid(i,j,1,1) -= dymom*Sbc;




    /* Face CD */
    di = -0.5e0;
    dj = 0.0e0;
    nnx = dy - cy;
    nny = cx - dx;
    len = sqrt(nnx*nnx + nny*nny);
    nnx /= len;
    nny /= len;
    nvec(0) = nnx;
    nvec(1) = nny;
    CFD2::Reconstruction<T>(order,qtygrid,i,j,di,dj,nvec,ql);

    di = 0.5e0;
    CFD2::Reconstruction<T>(order,qtygrid,i-1,j,di,dj,nvec,qr);

    #ifndef RELEASE
        ql0 = ql(0);
        ql1 = ql(1);
        ql2 = ql(2);
        ql3 = ql(3);
        qr0 = qr(0);
        qr1 = qr(1);
        qr2 = qr(2);
        qr3 = qr(3);
        assert(ql0 == ql0 && ql1 == ql1 && ql2 == ql2 && ql3 == ql3);
        assert(qr0 == qr0 && qr1 == qr1 && qr2 == qr2 && qr3 == qr3);
    #endif

    CFD2::RusanovFlux<T>(ql,qr,nvec,flux);

    drho = flux(0);
    dEt = flux(3);
    dxmom = flux(1);
    dymom = flux(2);

    #ifndef RELEASE
        f0 = flux(0);
        f1 = flux(1);
        f2 = flux(2);
        f3 = flux(3);
        assert(f0 == f0 && f1 == f1 && f2 == f2 && f3 == f3);
        rhonew = qtygrid(i,j,0) + drho;
        Etnew = qtygrid(i,j,1) + dEt;
        xmomnew = qtygrid(i,j,0,1) + dxmom;
        ymomnew = qtygrid(i,j,1,1) + dymom;
        assert(drho == drho && dxmom == dxmom && dymom == dymom && dEt == dEt);
    #endif

    dgrid(i,j,0) -= drho*Scd;
    dgrid(i,j,1) -= dEt*Scd;
    dgrid(i,j,0,1) -= dxmom*Scd;
    dgrid(i,j,1,1) -= dymom*Scd;




    /* Face DA */
    di = 0.0e0;
    dj = 0.5e0;
    nnx = ay - dy;
    nny = dx - ax;
    len = sqrt(nnx*nnx + nny*nny);
    nnx /= len;
    nny /= len;
    nvec(0) = nnx;
    nvec(1) = nny;
    CFD2::Reconstruction<T>(order,qtygrid,i,j,di,dj,nvec,ql);

    dj = -0.5e0;
    CFD2::Reconstruction<T>(order,qtygrid,i,j+1,di,dj,nvec,qr);

    #ifndef RELEASE
        ql0 = ql(0);
        ql1 = ql(1);
        ql2 = ql(2);
        ql3 = ql(3);
        qr0 = qr(0);
        qr1 = qr(1);
        qr2 = qr(2);
        qr3 = qr(3);
        assert(ql0 == ql0 && ql1 == ql1 && ql2 == ql2 && ql3 == ql3);
        assert(qr0 == qr0 && qr1 == qr1 && qr2 == qr2 && qr3 == qr3);
    #endif

    CFD2::RusanovFlux<T>(ql,qr,nvec,flux);

    drho = flux(0);
    dEt = flux(3);
    dxmom = flux(1);
    dymom = flux(2);

    #ifndef RELEASE
        f0 = flux(0);
        f1 = flux(1);
        f2 = flux(2);
        f3 = flux(3);
        assert(f0 == f0 && f1 == f1 && f2 == f2 && f3 == f3);
        rhonew = qtygrid(i,j,0) + drho;
        Etnew = qtygrid(i,j,1) + dEt;
        xmomnew = qtygrid(i,j,0,1) + dxmom;
        ymomnew = qtygrid(i,j,1,1) + dymom;
        assert(drho == drho && dxmom == dxmom && dymom == dymom && dEt == dEt);
    #endif

    dgrid(i,j,0) -= drho*Sda;
    dgrid(i,j,1) -= dEt*Sda;
    dgrid(i,j,0,1) -= dxmom*Sda;
    dgrid(i,j,1,1) -= dymom*Sda;






    /* Divide by cell volume */
    dgrid(i,j,0) /= volume;
    dgrid(i,j,1) /= volume;
    dgrid(i,j,0,1) /= volume;
    dgrid(i,j,1,1) /= volume;


} /* End of dQdt_Parallel() */


template <class T>
void CFD2::dQdt_Parallel_Viscous(const int order, const int i, const int j, Grid2D<T> &qtygrid,
                                 Array3D<T> &dvals, const T t, const T mu, const T kthermal)
{
    /* These lines done to avoid "unused variables" compiler warnings */
    T tloc = t + 1.0e0;
    tloc *= 2.0e0;

    /*
      Loop through grid points to calculate derivative values.  As stated in the
      header documentation, the first and last columns in the grid are not modified
      so that they function as boundary conditions.
    */
    T volume = 0.0e0;
    T ax = 0.0e0; T ay = 0.0e0;
    T bx = 0.0e0; T by = 0.0e0;
    T cx = 0.0e0; T cy = 0.0e0;
    T dx = 0.0e0; T dy = 0.0e0;
    T Sab = 0.0e0;
    T Sbc = 0.0e0;
    T Scd = 0.0e0;
    T Sda = 0.0e0;
    T nnx = 0.0e0;
    T nny = 0.0e0;
    T len = 0.0e0;
    T di = 0.0e0;
    T dj = 0.0e0;
    T drho = 0.0e0;
    T dEt = 0.0e0;
    T dxmom = 0.0e0;
    T dymom = 0.0e0;
    Array1D<T> nab(2,0.0e0);        /* Outward-normal vectors for each cell face */
    Array1D<T> nbc(2,0.0e0);
    Array1D<T> ncd(2,0.0e0);
    Array1D<T> nda(2,0.0e0);
    Array1D<T> flux(4,0.0e0);       /* General vector for fluxes */
    Array1D<T> qab(4,0.0e0);        /* Face values of conserved quantities (current cell) */
    Array1D<T> qbc(4,0.0e0);
    Array1D<T> qcd(4,0.0e0);
    Array1D<T> qda(4,0.0e0);
    Array1D<T> qrab(4,0.0e0);       /* Face values of conserved quantities (adjacent cells) */
    Array1D<T> qrbc(4,0.0e0);
    Array1D<T> qrcd(4,0.0e0);
    Array1D<T> qrda(4,0.0e0);

    /* Debugging variables */
    #ifndef RELEASE
        T rhonew = 0.0e0;
        T Etnew = 0.0e0;
        T xmomnew = 0.0e0;
        T ymomnew = 0.0e0;
        T f0 = 0.0e0;
        T f1 = 0.0e0;
        T f2 = 0.0e0;
        T f3 = 0.0e0;
        T ql0 = 0.0e0;
        T ql1 = 0.0e0;
        T ql2 = 0.0e0;
        T ql3 = 0.0e0;
        T qr0 = 0.0e0;
        T qr1 = 0.0e0;
        T qr2 = 0.0e0;
        T qr3 = 0.0e0;
    #endif

    /*
       Calculate volume of cell by finding the magnitude of the cross-product
       of the vectors which span the opposing corners of the quadrilateral element.

       Cell corners are found using bilinear interpolation with the surrounding
       grid points.

       Cell corners are lettered a-d, starting in the upper-right corner and
       moving CCW.  The gridpoint is at position 'N' in the sketch below.

            b ------- a
            |         |
            |    N    |
            |         |
            c ------- d

        Normal vectors to cell sides are assumed to be outward-positive.

    */

    CFD2::CalculateCellVertices<T>(qtygrid,i,j,ax,ay,bx,by,cx,cy,dx,dy);

    volume = fabs((ax-cx)*(dy-by) - (dx-bx)*(ay-cy));

    /* Calculate side lengths of cell */
    Sab = sqrt((ax-bx)*(ax-bx) + (ay-by)*(ay-by));
    Sbc = sqrt((bx-cx)*(bx-cx) + (by-cy)*(by-cy));
    Scd = sqrt((cx-dx)*(cx-dx) + (cy-dy)*(cy-dy));
    Sda = sqrt((dx-ax)*(dx-ax) + (dy-ay)*(dy-ay));


    /*
      Calculate inviscid flux terms.  These are done first since the reconstruction values are
      needed by the viscous flux calculations.
    */

    /* Face AB */
    di = 0.5e0;
    dj = 0.0e0;
    nnx = by - ay;
    nny = ax - bx;
    len = sqrt(nnx*nnx + nny*nny);
    nnx /= len;
    nny /= len;
    nab(0) = nnx;
    nab(1) = nny;
    CFD2::Reconstruction_v2<T>(order,qtygrid,i,j,di,dj,nab,qab);
    di = -0.5e0;
    CFD2::Reconstruction_v2<T>(order,qtygrid,i+1,j,di,dj,nab,qrab);

//    #ifndef RELEASE
//        ql0 = ql(0);
//        ql1 = ql(1);
//        ql2 = ql(2);
//        ql3 = ql(3);
//        qr0 = qr(0);
//        qr1 = qr(1);
//        qr2 = qr(2);
//        qr3 = qr(3);
//        assert(ql0 == ql0 && ql1 == ql1 && ql2 == ql2 && ql3 == ql3);
//        assert(qr0 == qr0 && qr1 == qr1 && qr2 == qr2 && qr3 == qr3);
//    #endif

    CFD2::RusanovFlux_Viscous<T>(qab,qrab,nab,flux);

    drho = flux(0);
    dEt = flux(3);
    dxmom = flux(1);
    dymom = flux(2);

    #ifndef RELEASE
        f0 = flux(0);
        f1 = flux(1);
        f2 = flux(2);
        f3 = flux(3);
        assert(f0 == f0 && f1 == f1 && f2 == f2 && f3 == f3);
        rhonew = qtygrid(i,j,0) + drho;
        Etnew = qtygrid(i,j,3) + dEt;
        xmomnew = qtygrid(i,j,1) + dxmom;
        ymomnew = qtygrid(i,j,2) + dymom;
        assert(drho == drho && dxmom == dxmom && dymom == dymom && dEt == dEt);
    #endif

    dvals(i,j,0) -= drho*Sab;
    dvals(i,j,1) -= dxmom*Sab;
    dvals(i,j,2) -= dymom*Sab;
    dvals(i,j,3) -= dEt*Sab;






    /* Face BC */
    di = 0.0e0;
    dj = -0.5e0;
    nnx = cy - by;
    nny = bx - cx;
    len = sqrt(nnx*nnx + nny*nny);
    nnx /= len;
    nny /= len;
    nbc(0) = nnx;
    nbc(1) = nny;
    CFD2::Reconstruction_v2<T>(order,qtygrid,i,j,di,dj,nbc,qbc);

    dj = 0.5e0;
    CFD2::Reconstruction_v2<T>(order,qtygrid,i,j-1,di,dj,nbc,qrbc);

//    #ifndef RELEASE
//        ql0 = ql(0);
//        ql1 = ql(1);
//        ql2 = ql(2);
//        ql3 = ql(3);
//        qr0 = qr(0);
//        qr1 = qr(1);
//        qr2 = qr(2);
//        qr3 = qr(3);
//        assert(ql0 == ql0 && ql1 == ql1 && ql2 == ql2 && ql3 == ql3);
//        assert(qr0 == qr0 && qr1 == qr1 && qr2 == qr2 && qr3 == qr3);
//    #endif

    CFD2::RusanovFlux_Viscous<T>(qbc,qrbc,nbc,flux);

    drho = flux(0);
    dEt = flux(3);
    dxmom = flux(1);
    dymom = flux(2);

    #ifndef RELEASE
        f0 = flux(0);
        f1 = flux(1);
        f2 = flux(2);
        f3 = flux(3);
        assert(f0 == f0 && f1 == f1 && f2 == f2 && f3 == f3);
        rhonew = qtygrid(i,j,0) + drho;
        Etnew = qtygrid(i,j,3) + dEt;
        xmomnew = qtygrid(i,j,1) + dxmom;
        ymomnew = qtygrid(i,j,2) + dymom;
        assert(drho == drho && dxmom == dxmom && dymom == dymom && dEt == dEt);
    #endif

    dvals(i,j,0) -= drho*Sbc;
    dvals(i,j,1) -= dxmom*Sbc;
    dvals(i,j,2) -= dymom*Sbc;
    dvals(i,j,3) -= dEt*Sbc;




    /* Face CD */
    di = -0.5e0;
    dj = 0.0e0;
    nnx = dy - cy;
    nny = cx - dx;
    len = sqrt(nnx*nnx + nny*nny);
    nnx /= len;
    nny /= len;
    ncd(0) = nnx;
    ncd(1) = nny;
    CFD2::Reconstruction_v2<T>(order,qtygrid,i,j,di,dj,ncd,qcd);

    di = 0.5e0;
    CFD2::Reconstruction_v2<T>(order,qtygrid,i-1,j,di,dj,ncd,qrcd);

//    #ifndef RELEASE
//        ql0 = ql(0);
//        ql1 = ql(1);
//        ql2 = ql(2);
//        ql3 = ql(3);
//        qr0 = qr(0);
//        qr1 = qr(1);
//        qr2 = qr(2);
//        qr3 = qr(3);
//        assert(ql0 == ql0 && ql1 == ql1 && ql2 == ql2 && ql3 == ql3);
//        assert(qr0 == qr0 && qr1 == qr1 && qr2 == qr2 && qr3 == qr3);
//    #endif

    CFD2::RusanovFlux_Viscous<T>(qcd,qrcd,ncd,flux);

    drho = flux(0);
    dEt = flux(3);
    dxmom = flux(1);
    dymom = flux(2);

    #ifndef RELEASE
        f0 = flux(0);
        f1 = flux(1);
        f2 = flux(2);
        f3 = flux(3);
        assert(f0 == f0 && f1 == f1 && f2 == f2 && f3 == f3);
        rhonew = qtygrid(i,j,0) + drho;
        Etnew = qtygrid(i,j,3) + dEt;
        xmomnew = qtygrid(i,j,1) + dxmom;
        ymomnew = qtygrid(i,j,2) + dymom;
        assert(drho == drho && dxmom == dxmom && dymom == dymom && dEt == dEt);
    #endif

    dvals(i,j,0) -= drho*Scd;
    dvals(i,j,1) -= dxmom*Scd;
    dvals(i,j,2) -= dymom*Scd;
    dvals(i,j,3) -= dEt*Scd;




    /* Face DA */
    di = 0.0e0;
    dj = 0.5e0;
    nnx = ay - dy;
    nny = dx - ax;
    len = sqrt(nnx*nnx + nny*nny);
    nnx /= len;
    nny /= len;
    nda(0) = nnx;
    nda(1) = nny;
    CFD2::Reconstruction_v2<T>(order,qtygrid,i,j,di,dj,nda,qda);

    dj = -0.5e0;
    CFD2::Reconstruction_v2<T>(order,qtygrid,i,j+1,di,dj,nda,qrda);

//    #ifndef RELEASE
//        ql0 = ql(0);
//        ql1 = ql(1);
//        ql2 = ql(2);
//        ql3 = ql(3);
//        qr0 = qr(0);
//        qr1 = qr(1);
//        qr2 = qr(2);
//        qr3 = qr(3);
//        assert(ql0 == ql0 && ql1 == ql1 && ql2 == ql2 && ql3 == ql3);
//        assert(qr0 == qr0 && qr1 == qr1 && qr2 == qr2 && qr3 == qr3);
//    #endif

    CFD2::RusanovFlux_Viscous<T>(qda,qrda,nda,flux);

    drho = flux(0);
    dEt = flux(3);
    dxmom = flux(1);
    dymom = flux(2);

    #ifndef RELEASE
        f0 = flux(0);
        f1 = flux(1);
        f2 = flux(2);
        f3 = flux(3);
        assert(f0 == f0 && f1 == f1 && f2 == f2 && f3 == f3);
        rhonew = qtygrid(i,j,0) + drho;
        Etnew = qtygrid(i,j,3) + dEt;
        xmomnew = qtygrid(i,j,1) + dxmom;
        ymomnew = qtygrid(i,j,2) + dymom;
        assert(drho == drho && dxmom == dxmom && dymom == dymom && dEt == dEt);
    #endif

    dvals(i,j,0) -= drho*Sda;
    dvals(i,j,1) -= dxmom*Sda;
    dvals(i,j,2) -= dymom*Sda;
    dvals(i,j,3) -= dEt*Sda;



    /* Calculate fiscous flux terms for each face */
//    ViscousFlux(Array1D<T> &q1, Array1D<T> &q2, Array1D<T> &qfl, Array1D<T> &qfr,
//                     Array1D<T> &nvec, Array1D<T> &flux, const T x1, const T y1, const T x2, const T y2,
//                     const int i, const int j, const int ir, const int jr, const Grid2D<T> &qtygrid,
//                     const T mu, const T k);


    /* Side AB */
    T x1 = 0.5e0*(bx + cx);
    T y1 = 0.5e0*(by + cy);
    T x2 = 0.5e0*(dx + ax);
    T y2 = 0.5e0*(dy + ay);
    int jr = j;
    int ir = i + 1;
    CFD2::ViscousFlux(qbc,qda,qab,qrab,nab,flux,x1,y1,x2,y2,i,j,ir,jr,qtygrid,mu,kthermal);

    for(int k=0; k<4; k++){
        dvals(i,j,k) += flux(k)*Sab;
    }



    /* Side BC */
    x1 = 0.5e0*(ax + bx);
    y1 = 0.5e0*(ay + by);
    x2 = 0.5e0*(cx + dx);
    y2 = 0.5e0*(cy + dy);
    jr = j - 1;
    ir = i;
    CFD2::ViscousFlux(qab,qcd,qbc,qrbc,nbc,flux,x1,y1,x2,y2,i,j,ir,jr,qtygrid,mu,kthermal);

    for(int k=0; k<4; k++){
        dvals(i,j,k) += flux(k)*Sbc;
    }




    /* Side CD */
    x1 = 0.5e0*(cx + bx);
    y1 = 0.5e0*(cy + by);
    x2 = 0.5e0*(ax + dx);
    y2 = 0.5e0*(ay + dy);
    jr = j;
    ir = i - 1;
    CFD2::ViscousFlux(qbc,qda,qcd,qrcd,ncd,flux,x1,y1,x2,y2,i,j,ir,jr,qtygrid,mu,kthermal);

    for(int k=0; k<4; k++){
        dvals(i,j,k) += flux(k)*Scd;
    }



    /* Side DA */
    x1 = 0.5e0*(ax + bx);
    y1 = 0.5e0*(ay + by);
    x2 = 0.5e0*(cx + dx);
    y2 = 0.5e0*(cy + dy);
    jr = j + 1;
    ir = i;
    CFD2::ViscousFlux(qab,qcd,qda,qrda,nda,flux,x1,y1,x2,y2,i,j,ir,jr,qtygrid,mu,kthermal);

    for(int k=0; k<4; k++){
        dvals(i,j,k) += flux(k)*Sda;
    }




    /* Divide by cell volume */
    dvals(i,j,0) /= volume;
    dvals(i,j,1) /= volume;
    dvals(i,j,2) /= volume;
    dvals(i,j,3) /= volume;

#ifndef RELEASE
    T d0 = dvals(i,j,0);
    T d1 = dvals(i,j,1);
    T d2 = dvals(i,j,2);
    T d3 = dvals(i,j,3);
#endif
}


template <class T>
void CFD2::CalculateCellVertices(Grid2D<T> &grid, const int i, const int j, T &ax, T &ay,
                           T &bx, T &by, T &cx, T &cy, T &dx, T &dy)
{
    int nx = grid.GetSize(2);
    int ny = grid.GetSize(1);

    /*
       Calculate volume of cell by finding the magnitude of the cross-product
       of the vectors which span the opposing corners of the quadrilateral element.

       Cell corners are found using bilinear interpolation with the surrounding
       grid points.

       Cell corners are lettered a-d, starting in the upper-right corner and
       moving CCW.  The gridpoint is at position 'N' in the sketch below.

            b ------- a
            |         |
            |    N    |
            |         |
            c ------- d

        Normal vectors to cell sides are assumed to be outward-positive.

    */
    if(i == 0){
        if(j == 0){
            ax = 0.25e0*(grid.xcoords(i,j)       + grid.xcoords(i,j+1) +
                         grid.xcoords(i+1,j+1)   + grid.xcoords(i+1,j));
            ay = 0.25e0*(grid.ycoords(i,j)       + grid.ycoords(i,j+1) +
                         grid.ycoords(i+1,j+1)   + grid.ycoords(i+1,j));
            bx = 0.5e0*(grid.xcoords(i,j)        + grid.xcoords(i+1,j));
            by = 0.5e0*(grid.ycoords(i,j)        + grid.ycoords(i+1,j));
            cx = grid.xcoords(i,j);
            cy = grid.ycoords(i,j);
            dx = 0.5e0*(grid.xcoords(i,j+1)     + grid.xcoords(i,j));
            dy = 0.5e0*(grid.ycoords(i,j+1)     + grid.ycoords(i,j));
        }
        if(j == nx-1){
            ax = 0.5e0*(grid.xcoords(i,j)       + grid.xcoords(i+1,j));
            ay = 0.5e0*(grid.ycoords(i,j)       +  grid.ycoords(i+1,j));
            bx = 0.25e0*(grid.xcoords(i,j-1)     + grid.xcoords(i,j) +
                         grid.xcoords(i+1,j)     + grid.xcoords(i+1,j-1));
            by = 0.25e0*(grid.ycoords(i,j-1)     + grid.ycoords(i,j) +
                         grid.ycoords(i+1,j)     + grid.ycoords(i+1,j-1));
            cx = 0.5e0*(grid.xcoords(i,j)       + grid.xcoords(i,j-1));
            cy = 0.5e0*(grid.ycoords(i,j)       + grid.ycoords(i,j-1));
            dx = grid.xcoords(i,j);
            dy = grid.ycoords(i,j);
        }
        if(j > 0 && j < nx-1){
            ax = 0.25e0*(grid.xcoords(i,j)       + grid.xcoords(i,j+1) +
                         grid.xcoords(i+1,j+1)   + grid.xcoords(i+1,j));
            ay = 0.25e0*(grid.ycoords(i,j)       + grid.ycoords(i,j+1) +
                         grid.ycoords(i+1,j+1)   + grid.ycoords(i+1,j));
            bx = 0.25e0*(grid.xcoords(i,j-1)     + grid.xcoords(i,j) +
                         grid.xcoords(i+1,j)     + grid.xcoords(i+1,j-1));
            by = 0.25e0*(grid.ycoords(i,j-1)     + grid.ycoords(i,j) +
                         grid.ycoords(i+1,j)     + grid.ycoords(i+1,j-1));
            cx = 0.5e0*(grid.xcoords(i,j)       + grid.xcoords(i,j-1));
            cy = 0.5e0*(grid.ycoords(i,j)       + grid.ycoords(i,j-1));
            dx = 0.5e0*(grid.xcoords(i,j+1)     + grid.xcoords(i,j));
            dy = 0.5e0*(grid.ycoords(i,j+1)     + grid.ycoords(i,j));
        }
    }
    if(i == ny-1){
        if(j == 0){
            ax = 0.5e0*(grid.xcoords(i,j)       + grid.xcoords(i,j+1));
            ay = 0.5e0*(grid.ycoords(i,j)       + grid.ycoords(i,j+1));
            bx = grid.xcoords(i,j);
            by = grid.ycoords(i,j);
            cx = 0.5e0*(grid.xcoords(i-1,j)      + grid.xcoords(i,j));
            cy = 0.5e0*(grid.ycoords(i-1,j)      + grid.ycoords(i,j));
            dx = 0.25e0*(grid.xcoords(i-1,j)     + grid.xcoords(i-1,j+1) +
                         grid.xcoords(i,j+1)     + grid.xcoords(i,j));
            dy = 0.25e0*(grid.ycoords(i-1,j)     + grid.ycoords(i-1,j+1) +
                         grid.ycoords(i,j+1)     + grid.ycoords(i,j));
        }
        if(j == nx-1){
            ax = grid.xcoords(i,j);
            ay = grid.ycoords(i,j);
            bx = 0.5e0*(grid.xcoords(i,j-1)     + grid.xcoords(i,j));
            by = 0.5e0*(grid.ycoords(i,j-1)     + grid.ycoords(i,j));
            cx = 0.25e0*(grid.xcoords(i-1,j-1)   + grid.xcoords(i-1,j) +
                         grid.xcoords(i,j)       + grid.xcoords(i,j-1));
            cy = 0.25e0*(grid.ycoords(i-1,j-1)   + grid.ycoords(i-1,j) +
                         grid.ycoords(i,j)       + grid.ycoords(i,j-1));
            dx = 0.5e0*(grid.xcoords(i-1,j)     + grid.xcoords(i,j));
            dy = 0.5e0*(grid.ycoords(i-1,j)     +  grid.ycoords(i,j));
        }
        if(j > 0 && j < nx-1){
            ax = 0.5e0*(grid.xcoords(i,j)       + grid.xcoords(i,j+1));
            ay = 0.5e0*(grid.ycoords(i,j)       + grid.ycoords(i,j+1));
            bx = 0.5e0*(grid.xcoords(i,j-1)     + grid.xcoords(i,j));
            by = 0.5e0*(grid.ycoords(i,j-1)     + grid.ycoords(i,j));
            cx = 0.25e0*(grid.xcoords(i-1,j-1)   + grid.xcoords(i-1,j) +
                         grid.xcoords(i,j)       + grid.xcoords(i,j-1));
            cy = 0.25e0*(grid.ycoords(i-1,j-1)   + grid.ycoords(i-1,j) +
                         grid.ycoords(i,j)       + grid.ycoords(i,j-1));
            dx = 0.25e0*(grid.xcoords(i-1,j)     + grid.xcoords(i-1,j+1) +
                         grid.xcoords(i,j+1)     + grid.xcoords(i,j));
            dy = 0.25e0*(grid.ycoords(i-1,j)     + grid.ycoords(i-1,j+1) +
                         grid.ycoords(i,j+1)     + grid.ycoords(i,j));
        }
    }
    if(i > 0 && i < ny-1){
        if(j == 0){
            ax = 0.25e0*(grid.xcoords(i,j)       + grid.xcoords(i,j+1) +
                         grid.xcoords(i+1,j+1)   + grid.xcoords(i+1,j));
            ay = 0.25e0*(grid.ycoords(i,j)       + grid.ycoords(i,j+1) +
                         grid.ycoords(i+1,j+1)   + grid.ycoords(i+1,j));
            bx = 0.5e0*(grid.xcoords(i,j)        + grid.xcoords(i+1,j));
            by = 0.5e0*(grid.ycoords(i,j)        + grid.ycoords(i+1,j));
            cx = 0.5e0*(grid.xcoords(i-1,j)      + grid.xcoords(i,j));
            cy = 0.5e0*(grid.ycoords(i-1,j)      + grid.ycoords(i,j));
            dx = 0.25e0*(grid.xcoords(i-1,j)     + grid.xcoords(i-1,j+1) +
                         grid.xcoords(i,j+1)     + grid.xcoords(i,j));
            dy = 0.25e0*(grid.ycoords(i-1,j)     + grid.ycoords(i-1,j+1) +
                         grid.ycoords(i,j+1)     + grid.ycoords(i,j));
        }
        if(j == nx-1){
            ax = 0.5e0*(grid.xcoords(i,j)       + grid.xcoords(i+1,j));
            ay = 0.5e0*(grid.ycoords(i,j)       + grid.ycoords(i+1,j));
            bx = 0.25e0*(grid.xcoords(i,j-1)     + grid.xcoords(i,j) +
                         grid.xcoords(i+1,j)     + grid.xcoords(i+1,j-1));
            by = 0.25e0*(grid.ycoords(i,j-1)     + grid.ycoords(i,j) +
                         grid.ycoords(i+1,j)     + grid.ycoords(i+1,j-1));
            cx = 0.25e0*(grid.xcoords(i-1,j-1)   + grid.xcoords(i-1,j) +
                         grid.xcoords(i,j)       + grid.xcoords(i,j-1));
            cy = 0.25e0*(grid.ycoords(i-1,j-1)   + grid.ycoords(i-1,j) +
                         grid.ycoords(i,j)       + grid.ycoords(i,j-1));
            dx = 0.5e0*(grid.xcoords(i-1,j)     + grid.xcoords(i,j));
            dy = 0.5e0*(grid.ycoords(i-1,j)     + grid.ycoords(i,j));
        }
        if(j > 0 && j < nx-1){
            ax = 0.25e0*(grid.xcoords(i,j)       + grid.xcoords(i,j+1) +
                         grid.xcoords(i+1,j+1)   + grid.xcoords(i+1,j));
            ay = 0.25e0*(grid.ycoords(i,j)       + grid.ycoords(i,j+1) +
                         grid.ycoords(i+1,j+1)   + grid.ycoords(i+1,j));
            bx = 0.25e0*(grid.xcoords(i,j-1)     + grid.xcoords(i,j) +
                         grid.xcoords(i+1,j)     + grid.xcoords(i+1,j-1));
            by = 0.25e0*(grid.ycoords(i,j-1)     + grid.ycoords(i,j) +
                         grid.ycoords(i+1,j)     + grid.ycoords(i+1,j-1));
            cx = 0.25e0*(grid.xcoords(i-1,j-1)   + grid.xcoords(i-1,j) +
                         grid.xcoords(i,j)       + grid.xcoords(i,j-1));
            cy = 0.25e0*(grid.ycoords(i-1,j-1)   + grid.ycoords(i-1,j) +
                         grid.ycoords(i,j)       + grid.ycoords(i,j-1));
            dx = 0.25e0*(grid.xcoords(i-1,j)     + grid.xcoords(i-1,j+1) +
                         grid.xcoords(i,j+1)     + grid.xcoords(i,j));
            dy = 0.25e0*(grid.ycoords(i-1,j)     + grid.ycoords(i-1,j+1) +
                         grid.ycoords(i,j+1)     + grid.ycoords(i,j));
        }
    }
}


template <class T>
void CFD2::ViscousFlux(Array1D<T> &q1, Array1D<T> &q2, Array1D<T> &qfl, Array1D<T> &qfr,
                       Array1D<T> &nvec, Array1D<T> &flux, const T x1, const T y1, const T x2, const T y2,
                       const int i, const int j, const int ir, const int jr, Grid2D<T> &qtygrid,
                       const T mu, const T k)
{
    /* Require proper array sizes */
    assert(q1.GetDim() == 4);
    assert(q2.GetDim() == 4);
    assert(qfl.GetDim() == 4);
    assert(qfr.GetDim() == 4);
    assert(flux.GetDim() == 4);
    assert(nvec.GetDim() == 2);


    /* Calculate values of "right" cell, handling boundaries */
    int ny = qtygrid.GetSize(1);
    assert(jr >= 0 && jr < qtygrid.GetSize(2));  /* Only handle i-direction boundaries */
    T xr = 0.0e0;
    T yr = 0.0e0;
    T rhor = 0.0e0;
    T xmomr = 0.0e0;
    T ymomr = 0.0e0;
    T Etr = 0.0e0;
    T ur = 0.0e0;
    T vr = 0.0e0;
    if(ir < 0){
        xr = qtygrid.xcoords(0,jr);
        yr = -qtygrid.ycoords(1,jr);
        rhor = qtygrid(0,jr,0);
        xmomr = qtygrid(0,jr,1);
        ymomr = -qtygrid(0,jr,2);       /* Symmetric about nozzle centerline */
        Etr = qtygrid(0,jr,3);
        ur = xmomr/rhor;
        vr = ymomr/rhor;
    }
    if(ir >= ny){
        xr = qtygrid.xcoords(ny-1,jr);
        yr = qtygrid.ycoords(ny-1,jr) + qtygrid.ycoords(ny-1,jr) - qtygrid.ycoords(ny-2,jr);
        rhor = qtygrid(ny-1,jr,0);
        xmomr = -qtygrid(ny-1,jr,1);    /* Non-slip wall */
        ymomr = -qtygrid(ny-1,jr,2);    /* Non-slip wall */
        Etr = qtygrid(ny-1,jr,3);
        ur = xmomr/rhor;
        vr = ymomr/rhor;
    }
    if(ir >= 0 && ir < ny){
        xr = qtygrid.xcoords(ir,jr);
        yr = qtygrid.ycoords(ir,jr);
        rhor = qtygrid(ir,jr,0);
        xmomr = qtygrid(ir,jr,1);
        ymomr = qtygrid(ir,jr,2);
        Etr = qtygrid(ir,jr,3);
        ur = xmomr/rhor;
        vr = ymomr/rhor;
    }


    /* Calculate path from point 1 to right grid cell */
    T lx1 = xr - x1;
    T ly1=  yr - y1;
    T l1 = sqrt(lx1*lx1 + ly1*ly1);
    lx1 /= l1;  /* Make unit vector component */
    ly1 /= l1;  /* Make unit vector component */


    /* Calculate path from point 2 to right grid cell */
    T lx2 = xr - x2;
    T ly2 = yr - y2;
    T l2 = sqrt(lx2*lx2 + ly2*ly2);
    lx2 /= l2;  /* Make unit vector component */
    ly2 /= l2;  /* Make unit vector component */


    /* Calculate gradient of x-momentum terms */
    T xgradmag1 = (ur - q1(1)/q1(0))/l1;
    T xgradmag2 = (ur - q2(1)/q1(0))/l2;


    /* Calculate gradient of y-momentum terms */
    T ygradmag1 = (vr - q1(2)/q1(0))/l1;
    T ygradmag2 = (vr - q2(2)/q1(0))/l2;


    /* Calculate velocity gradients from momentum gradients */
    T det = 1.0e0/(ly2*lx1 - lx2*ly1);
    T dudx = det*(ly2*xgradmag1 - ly1*xgradmag2);
    T dudy = det*(lx1*xgradmag2 - lx2*xgradmag1);
    T dvdx = det*(ly2*ygradmag2 - ly1*ygradmag1);
    T dvdy = det*(lx1*ygradmag1 - lx2*ygradmag2);


    /* Calculate shear stresses from velocity gradients */
    T tauxx = 2.0e0/3.0e0*mu*(2.0e0*dudx - dvdy);
    T tauyy = 2.0e0/3.0e0*mu*(2.0e0*dvdy - dudx);
    T tauxy = mu*(dudy + dvdx);



    /* Calculate temperature gradient between cells */
    T lx = xr - qtygrid.xcoords(i,j);
    T ly = yr - qtygrid.ycoords(i,j);
    T l = sqrt(lx*lx + ly*ly);

    T u = xmomr/rhor;
    T v = ymomr/rhor;
    T pright = 0.4e0*(Etr - 0.5e0*rhor*(u*u + v*v));
    u = qtygrid(i,j,1)/qtygrid(i,j,0);
    v = qtygrid(i,j,2)/qtygrid(i,j,0);
    T rholeft = qtygrid(i,j,0);
    T pleft = 0.4e0*(qtygrid(i,j,3) - 0.5e0*rholeft*(u*u + v*v));

    T Tleft = pleft/rholeft;        /* R = 1 */
    T Tright = pright/rhor;         /* R = 1 */

    T dTdn = (Tright - Tleft)/l;


    /* Calculate viscous energy flux on each side of interface */
    T uleft = qfl(1)/qfl(0);
    T vleft = qfl(2)/qfl(0);
    T uright = qfr(1)/qfr(0);
    T vright = qfr(2)/qfr(0);
    T fluxEleft = uleft*(nvec(0)*tauxx + nvec(1)*tauxy) + vleft*(nvec(0)*tauxy + nvec(1)*tauyy) + k*dTdn;
    T fluxEright = uright*(nvec(0)*tauxx + nvec(1)*tauxy) + vright*(nvec(0)*tauxy + nvec(1)*tauyy) + k*dTdn;


    /* Calculate viscous fluxes */
    flux(0) = 0;
    flux(1) = tauxx*nvec(0) + tauxy*nvec(1);
    flux(2) = tauxy*nvec(0) + tauyy*nvec(1);
    flux(3) = 0.5e0*(fluxEleft + fluxEright);




#ifndef RELEASE
    T f0 = flux(0);
    T f1 = flux(1);
    T f2 = flux(2);
    T f3 = flux(3);
#endif
}
