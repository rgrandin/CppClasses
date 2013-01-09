/**
 * @file CFD2.h
 * @author Robert Grandin
 * @date 25 March 2012
 * @brief Definition of CFD2 namespace.
 *
 * @section Class Description & Notes
 *
 * This collection of functions is for use in solving hyperbolic partial differential
 * equations.  The member functions of this class are targed towards their use in solving
 * computational fluid dynamics (CFD) problems.
 *
 * These functions are intended for 2-dimensional problems of the form
 * @f$ \frac{\partial Q}{\partial t} + \frac{\partial E}{\partial x} + \frac{\partial F}{\partial y} = 0 @f$,
 * where Q, E, and F are 4-element vectors corresponding to density, x-velocity, y-velocity, and total energy.
 *
 * The physical quantities associated with the Grid2D object are expected to be as follows (unless otherwise noted
 * in the documentation for specific functions):
 *  - 1st scalar quantity: density
 *  - 2nd scalar quantity: total energy
 *  - 3rd scalar quantity: pressure
 *  - 1st vector quantity: velocity (2D velocity, with 3rd component equal to 0)
 *  - 2nd vector quantity: momentum (2D momentum, with 3rd component equal to 0)
 *
 * The density, total energy, and momentum are used in the calculations, with pressure and velocity calculated
 * later.
 *
 * These equations are formulated using a node-centered approach.
 *
 * MPI commands are used within this class.  To enable MPI support, USEMPI must be defined at compile time
 * (-DUSEMPI using the GCC compiler).
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * <b>CURRENT ASSUMPTIONS</b>
 *
 * These functions assume that the "lower" boundary (y = ymin for all x) is a line of symmetry
 * in the problem while the "upper" boundary (y = ymax for all x) is a hard physical boundary
 * similar to a wall (i.e., no normal flow allowed while tangential flow is acceptable).  The
 * first and last column in the grid are assumed to be the applied boundary conditions and are not
 * modified by these routines.
 *
 *
 * @section Revisions
 *
 * @date 25 March 2012
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

#ifndef CFD2_
#define CFD2_

#include <omp.h>

#ifdef USEMPI
#include <mpi.h>
#endif

#include <Array1D.h>
#include <Array2D.h>
#include <Grid2D.h>
#include <MPI_Custom.h>


/**
  @brief Collection of functions for solving 2D CFD problems.  See the documentation for CFD2.h
    for more detailed information regarding this collection of functions.
*/
namespace CFD2 {

/**
  @brief Calculate the time-derivatives of the physical quantities for each cell using a
    first-order spatial discretization.
  @param qtygrid Grid2D object containing the current values for the physical quantities.
  @param dgrid Grid2D object containing the time derivatives of the physical quantities at
    each grid point.  This is the "output" from this function and can be used by numeric
    integrators to integrate the flow field through time.
  @param t Current solution time.
  @return None.
*/
template <class T>
void dQdt_FirstOrder(Grid2D<T> &qtygrid, Grid2D<T> &dgrid, const T t);


/**
  @brief Calculate the time-derivatives of the physical quantities for each cell using a
    second-order spatial discretization.  MPI parallelization is used.
  @param qtygrid Grid2D object containing the current values for the physical quantities.
  @param dgrid Grid2D object containing the time derivatives of the physical quantities at
    each grid point.  This is the "output" from this function and can be used by numeric
    integrators to integrate the flow field through time.
  @param t Current solution time.
  @return None.
*/
template <class T>
void dQdt_SecondOrder(Grid2D<T> &qtygrid, Grid2D<T> &dgrid, const T t);


/**
  @brief Calculate the time-derivatives of the physical quantities for each cell using a
    second-order spatial discretization.  MPI parallelization is used.
  @param qtygrid Grid2D object containing the current values for the physical quantities.
  @param dvals Array3D object containing the time derivatives of the physical quantities at
    each grid point.
  @param t Current solution time.
  @param dependency Array2D containing dependency information for MPI.
  @param nrecv Number of MPI RECV() required by local process when sharing data.
  @param mu Dynamic viscosity.
  @param kthermal Thermal conductivity.
  @return None.
*/
template <class T>
void dQdt_SecondOrder_Viscous_MPI(Grid2D<T> &qtygrid, Array3D<T> &dvals, const T t,
                                  Array2D<int> &dependency, const int nrecv,
                                  const T mu, const T kthermal);


/**
  @brief Calculate the time-derivatives of the physical quantities for each cell using a
    third-order spatial discretization.
  @param qtygrid Grid2D object containing the current values for the physical quantities.
  @param dgrid Grid2D object containing the time derivatives of the physical quantities at
    each grid point.  This is the "output" from this function and can be used by numeric
    integrators to integrate the flow field through time.
  @param t Current solution time.
  @return None.
*/
template <class T>
void dQdt_ThirdOrder(Grid2D<T> &qtygrid, Grid2D<T> &dgrid, const T t);


/**
  @brief Calculate the time-derivatives of the physical quantities for each cell in parallel.
  @param order Order of accuracy of the reconstruction.
  @param i Row-index of cell for which flux is to be calculated.
  @param j Column-index of cell for which flux is to be calculated.
  @param qtygrid Grid2D object containing the current values for the physical quantities.
  @param dgrid Grid2D object containing the time derivatives of the physical quantities at
    each grid point.  This is the "output" from this function and can be used by numeric
    integrators to integrate the flow field through time.
  @param t Current solution time.
  @return None.
  @warning It is left to the user to ensure that accesses to dgrid are appropriate from multiple
    threads.  No thread-safety measures are implemented in this function.
*/
template <class T>
void dQdt_Parallel(const int order, const int i, const int j, Grid2D<T> &qtygrid, Grid2D<T> &dgrid, const T t);


/**
  @brief Calculate the time-derivatives of the physical quantities for each cell in parallel.

  This function
    is for a viscous flow, and calculates the fluxes associated with the inviscid terms.
  @param order Order of accuracy of the reconstruction.
  @param i Row-index of cell for which flux is to be calculated.
  @param j Column-index of cell for which flux is to be calculated.
  @param qtygrid Grid2D object containing the current values for the physical quantities.
  @param dvals Array3D object containing the time derivatives of the physical quantities at
    each grid point.
  @param t Current solution time.
  @param mu Dynamic viscosity.
  @param kthermal Thermal conductivity.
  @return None.
  @warning It is left to the user to ensure that accesses to dgrid are appropriate from multiple
    threads.  No thread-safety measures are implemented in this function.
  @warning This function assumes that qtygrid has 4 scalars: density, x-momentum, y-momentum, and total energy.
*/
template <class T>
void dQdt_Parallel_Viscous(const int order, const int i, const int j, Grid2D<T> &qtygrid,
                           Array3D<T> &dvals, const T t, const T mu, const T kthermal);


/**
  @brief Calculate Rusanov flux.
  @param ql Array1D object containing the left-state vector.
  @param qr Array1D object containing the right-state vector.
  @param nvec Array1D object containing the unit-normal vector to the face (outward normal).
  @param flux Flux associated with each quantity in the ql and qr vectors.  This flux is intended
    to be added to the current cell (i.e., inflow).
  @return None.
  @warning The ql, qr, and flux arrays are expected to contain the quantities: [density, x-momentum, y-momentum, total energy].
*/
template <class T>
void RusanovFlux(Array1D<T> &ql, Array1D<T> &qr, Array1D<T> &nvec, Array1D<T> &flux);


/**
  @brief Calculate diffusive viscous flux terms.

  The interface flux value is the simple average of the flux
    quantities of each cell.  It is assumed that the reconstruction of conserved quantities are 2nd-order
    accurate (piecewise linear), making the shear stresses piecewise constant within the cell.
  @param q1 Conserved quantity values at first point used when calculating gradients.
  @param q2 Conserved quantity values at second point used when calculating gradients.
  @param qfl Conserved quantity values on left side of interface.
  @param qfr Conserved quantity values on right side of interface.
  @param nvec Unit-normal vector for cell face, outward normal.
  @param flux Viscous fluxes associated with each conserved quantity at the interface.
  @param x1 X-coordinate of first point used when calculating gradients.
  @param y1 Y-coordinate of first point used when calculating gradients.
  @param x2 X-coordinate of second point used when calculating gradients.
  @param y2 Y-coordinate of second point used when calculating gradients.
  @param i Row-index of current cell in grid.
  @param j Column-index of current cell in grid.
  @param ir Row-index of adjacent cell in grid.
  @param jr Column-index of adjacent cell in grid.
  @param qtygrid Grid2D object containing the data.
  @param mu Dynamic viscosity.
  @param k Thermal conductivity.
  @return None.
  @warning The ql, qr, and flux arrays are expected to contain the quantities: [density, x-momentum, y-momentum, total energy].
  @warning This routine only handles boundary conditions in the i-index direction.  Like with the other routines in this
    namespace, the first and last column values are left constant as boundary values, so the j-index is expected to vary
    between [1,nx-2].
*/
template <class T>
void ViscousFlux(Array1D<T> &q1, Array1D<T> &q2, Array1D<T> &qfl, Array1D<T> &qfr,
                 Array1D<T> &nvec, Array1D<T> &flux, const T x1, const T y1, const T x2, const T y2,
                 const int i, const int j, const int ir, const int jr, Grid2D<T> &qtygrid,
                 const T mu, const T k);


/**
  @brief Calculate Rusanov flux when the flow is viscous.
  @param ql Array1D object containing the left-state vector.
  @param qr Array1D object containing the right-state vector.
  @param nvec Array1D object containing the unit-normal vector to the face (outward normal).
  @param flux Flux associated with each quantity in the ql and qr vectors.  This flux is intended
    to be added to the current cell (i.e., inflow).
  @return None.
  @warning The ql, qr, and flux arrays are expected to contain the quantities: [density, x-momentum, y-momentum, total energy].
*/
template <class T>
void RusanovFlux_Viscous(Array1D<T> &ql, Array1D<T> &qr, Array1D<T> &nvec, Array1D<T> &flux);


/**
  @brief Perform reconstruction of physical quantities in Grid2D object.

  This is to determine the
    appropriate values for each quantity at the cell interfaces.  Regardless of the chosen order, the
    reconstruction is done only considering one grid-direction at a time (i.e., along rows or columns).
  @param order Order of reconstruction to be used.
    - 1: First-order (piecewise constant values)
    - 2: Second-order (piecewise linear values)
    - 3: Third-order (piecewise parabolic)
  @param qtygrid Grid2D object containing quantity values for node-centered cells.
  @param i Row-index of current cell.
  @param j Column-index of current cell.
  @param di Difference in row-index of interface and current cell.
  @param dj Difference in column-index of interface and current cell.
  @param nvec Array1D containing the (nx,ny) outward normal vector for the interface.
  @param qvec Array1D of the physical quantities at the interface.  This is intended for use in the typical
    CFD equations, making this a 4-element array containing: [density, x-momentum, y-momentum, total energy].
  @warning The scalar and vector quantities discussed in the detailed description of CFD2.h must
    be present with their expected indices (e.g., density is the first scalar quantity).
  @warning It is assumed that the 0th row lies on a symmetry boundary while the upper-most row of the grid
    is adjacent to an impermeable wall.  Further, it is assumed that the outer-most columns remain fixed-value
    as boundary conditions.  For example, consider the computational grid for one half of a nozzle cross-section,
    with fixed inlet and exit conditions (admittedly very limited in applicability, but useful for now).
*/
template <class T>
void Reconstruction(int order, Grid2D<T> &qtygrid, int i, int j, T di, T dj,
                    Array1D<T> &nvec, Array1D<T> &qvec);


/**
  @brief Perform reconstruction of physical quantities in Grid2D object.

  This is to determine the
    appropriate values for each quantity at the cell interfaces.  Regardless of the chosen order, the
    reconstruction is done only considering one grid-direction at a time (i.e., along rows or columns).
  @param order Order of reconstruction to be used.
    - 1: First-order (piecewise constant values)
    - 2: Second-order (piecewise linear values)
    - 3: Third-order (piecewise parabolic)
  @param qtygrid Grid2D object containing quantity values for node-centered cells.
  @param i Row-index of current cell.
  @param j Column-index of current cell.
  @param di Difference in row-index of interface and current cell.
  @param dj Difference in column-index of interface and current cell.
  @param nvec Array1D containing the (nx,ny) outward normal vector for the interface.
  @param qvec Array1D of the physical quantities at the interface.  This is intended for use in the typical
    CFD equations, making this a 4-element array containing: [density, x-momentum, y-momentum, total energy].
  @warning This function assumes that qtygrid has 4 scalars: density, x-momentum, y-momentum, and total energy.
  @warning It is assumed that the 0th row lies on a symmetry boundary while the upper-most row of the grid
    is adjacent to an impermeable wall.  Further, it is assumed that the outer-most columns remain fixed-value
    as boundary conditions.  For example, consider the computational grid for one half of a nozzle cross-section,
    with fixed inlet and exit conditions (admittedly very limited in applicability, but useful for now).
*/
template <class T>
void Reconstruction_v2(int order, Grid2D<T> &qtygrid, int i, int j, T di, T dj,
                    Array1D<T> &nvec, Array1D<T> &qvec);


/**
  @brief Minmod flux limiter.
  @param x First argument.
  @param y Second argument.
  @return Two possibilities:
    - x*y < 0: 0
    - x*y > 0: sign(x)*min(abs(x),abs(y))
*/
template <class T>
T MinMod(T x, T y);


/**
  @brief Calculate the sign of a number.

  Implementation code found at
    http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c/4609795#4609795.
  @param val Value for which its sign is desired.
  @return Integer value identifying if number is negative, zero, or positive (-1, 0, or 1, respectively).
*/
template <class T>
int sgn(T val);


/**
  @brief Calculate the cell vertex coordinates for the specified grid point.

  This assumes that the grid is node-centered
    (i.e., that the grid points are the cell centers) and that the outer edges of the grid are the true boundaries.
    Cell vertices are labeled {A, B, C, D} in a counter-clockwise motion starting with "A" in the upper-right corner.
  @param grid Grid2D object containing the nodal coordinates.
  @param i Row-index of grid point which is the cell-center.
  @param j Column-index of grid point which is the cell-center.
  @param ax X-position of point A.
  @param ay Y-position of point A.
  @param bx X-position of point B.
  @param by Y-position of point B.
  @param cx X-position of point C.
  @param cy Y-position of point C.
  @param dx X-position of point D.
  @param dy Y-position of point D.
  @return None.
*/
template <class T>
void CalculateCellVertices(Grid2D<T> &grid, const int i, const int j, T &ax, T &ay,
                           T &bx, T &by, T &cx, T &cy, T &dx, T &dy);

}

#include "CFD2.cpp"

#endif /* CFD2_ */
