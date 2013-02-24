/**
 * @file Grid2D.h
 * @author Robert Grandin
 * @date 26 November 2010
 * @brief Definition of Grid2D class.
 *
 * @section Class Description & Notes
 *
 * This class represents a 2-dimensional regular grid for use in
 * computational mechanics.  In this class the  ollowing direction-notation is used:
 * 	- Horizontal direction
 * 		- Also known as "x" in spatial coordinates (in the event that the grid
 * 		  is aligned with spatial directions).
 * 		- Uses index "j" in the computational space (column index).
 * 		- Uses greek letter "zi" in the determination of derivatives and the
 * 		  relating of spatial and computational derivatives.
 *
 *  - Vertical direction
 * 		- Also known as "y" in spatial coordinates (in the event that the grid
 * 		  is aligned with spatial directions).
 * 		- Uses index "i" in the computational space (row index).
 * 		- Uses greek letter "eta" in the determination of derivatives and the
 * 		  relating of spatial and computational derivatives.
 *
 * Member arrays containing grid information, such as coordinates of grid points
 * and grid metrics, are public members to simplify their use in solution algorithms.
 * Solution algorithms are deliberately left out of this class and the focus of this
 * class is kept on the storing of grid information and the data associated with each
 * grid point.
 *
 * Bounds-checking during array access is performed if RELEASE is not defined.
 * Defining RELEASE, either manually via 'define RELEASE' or via the compiler
 * using '-DRELEASE' on gcc, will disable bounds-checking.
 *
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 26 November 2010
 *	- Creation date.
 * @date 20 January 2012
 *  - Made updates to improve class in preparation for AerE 547.
 * @date 22 February 2012
 *  - Updated O-grid generation routine to use more-flexible function pointers.
 *  - Added routine for elliptic smoothing of the grid.
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2010-2012, Robert Grandin
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

#ifndef Grid2D_
#define Grid2D_

#include <cmath>
#include <fstream>
#include <Array1D.h>
#include <Array2D.h>
#include <Array3D.h>
#include <PArray1D.h>

#ifdef USEMPI
#include <mpi.h>
#include <MPI_Custom.h>
#endif


/**
 * @brief Two-dimensional grid generation and manipulation for use with
 * 			computational mechanics applications.
 * @warning C++11 features, such as move-constructor and move-assignment, require the symbol
 *  "CXX11" to be defined.
 */
template <class T> class Grid2D{

public:
	/**
     * @brief Create Grid2D object.  No grid parameters defined.  No arrays created
     *     for storing data at grid points.
	 * @pre Sufficient memory exists for object creation.
	 * @post Grid2D object created.
	 * @return None.
	 */
	Grid2D();


    /**
     * @brief Copy constructor.
     * @param a Reference to existing Grid2D object to be copied.
     */
    Grid2D(const Grid2D<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing Grid2D object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    Grid2D(Grid2D<T> &&a);
#endif


	/**
	 * @brief Destructor.
	 * @pre Grid2D object exists.
	 * @post Grid2D object destroyed.
	 * @return None.
	 */
    virtual ~Grid2D();


	/**
     * @brief Create O-Grid.
     *
     * The functions defining coordinates of the inner
     *  and outer surfaces are defined elsewhere and pointers are passed
     *  to this routine.  This allows greater flexibility for this class
     *  as no surfaces are defined in the class itself.
     * @pre Grid2D object exists, and sufficient memory is available for the
     *  grid and its associated properties.
	 * @param nrad Number of radial positions in the grid.
	 * @param ntheta Number of angular positions in the grid.
     * @param fx Pointer to function defining the x-coordinates of the ellipse.
     * @param fy Pointer to function defining the x-coordinates of the ellipse.
     * @param rotouter Rotation of outer ellipse, in radians.
     * @param rotinner Rotation of inner ellipse, in radians.
     * @param rotoffset Rotational-offset, in radians, used when mapping points.
     * @param ao Semi-major axis of outer ellipse.
     * @param bo Semi-minor axis of outer ellipse.
     * @param ai Semi-major axis of inner ellipse.
     * @param bi Semi-minor axis of inner ellipse.
     * @param Pp Adjustable parameter for inner wall.  Recommended value is 0.1.
     * @param Qq adjustable parameter for outer wall.  Recommended value is 0.1.
	 * @post O-Grid created.
	 * @return None.
     * @warning Functions which generate outer grid points expect the following arguments:
     *  angular position, rotation of ellipse, major semi-axis, minor semi-axis.  Note that
     *  semi-major and semi-minor refer to the x and y directions in the unrotated ellipses.
	 */
    void CreateOGridFromFunctions(int nrad, int ntheta, T(*fx)(T, T, T, T),
                                  T(*fy)(T, T, T, T), T rotouter, T rotinner, T rotoffset,
                                  T ao, T bo, T ai, T bi, T Pp, T Qq);


	/**
     * @brief Create H-Grid.
     *
     * The functions defining coordinates of the inner
     *  and outer surfaces are defined elsewhere and pointers are passed
     *  to this routine.  This allows greater flexibility for this class
     *  as no surfaces are defined in the class itself.
     * @pre Grid2D object exists, and sufficient memory is available for the
     *  grid and its associated properties.
	 * @param nx Number of x positions in the grid.
	 * @param ny Number of y positions in the grid.
	 * @param fxin Pointer to function defining the x-coordinates of the inner
	 * 				surface.  This function is paramaterized by a single input
	 * 				value (angular position), and defined elsewhere.
	 * @param fxout Pointer to function defining the x-coordinates of the outer
	 * 				surface.  This function is paramaterized by a single input
	 * 				value (angular position), and defined elsewhere.
	 * @param fyin Pointer to function defining the y-coordinates of the inner
	 * 				surface.  This function is paramaterized by a single input
	 * 				value (angular position), and defined elsewhere.
	 * @param fyout Pointer to function defining the y-coordinates of the outer
	 * 				surface.  This function is paramaterized by a single input
	 * 				value (angular position), and defined elsewhere.
     * @param hermite True if Hermite interpolation is to be used for spacing internal
     *      points.  If true, parameters Pp and Qq will be used.
     * @param Pp Adjustable parameter for inner wall.  Recommended value is 0.1.
     * @param Qq adjustable parameter for outer wall.  Recommended value is 0.1.
     * @param xxmin Minimum x-coordinate of grid.
     * @param xxmax Maximum x-coordinate of grid.
     * @param yymin Minimum y-coordinate of grid.
     * @param yymax Maximum y-coordinate of grid.
     * @param stretch1 Stretch the grid as-described on pages 335-336 of Computational Fluid
     *      Mechanics and Heat Transfer, 2nd Edition.
     * @param alpha Controls stretching location.  Set to 0 for stretching at ymax boundary only.
     *      Set to 1/2 for stretching equally at ymin and ymax boundaries.
     * @param beta Controls scaling of stretching. @f$ \beta = \left( 1 - \frac{\delta}{h} \right)^{-1/2} @f$.
     *      @f$ \delta @f$ is the height of the boundary layer.
	 * @post H-Grid created.
	 * @return None.
	 */
	void CreateHGridFromFunctions(int nx, int ny, T(*fxin)(T),
			T(*fxout)(T), T(*fyin)(T), T(*fyout)(T), const bool hermite, T Pp,
			T Qq, T xxmin, T xxmax, T yymin, T yymax, const bool stretch1, const T alpha, const T beta);


    /**
      @brief Create generic H grid without using boundary functions.
      @pre Grid2D object exist.
      @param nnx Number of points in the x-direction.
      @param nny Number of points in the y-direction.
      @param xxmin Minimum x-extent.
      @param xxmax Maximum x-extent.
      @param yymin Minimum y-extent.
      @param yymax Minimum y-extent.
      @post H-grid created.
      @return None.
    */
    void CreateGenericHGrid(const int nnx, const int nny, const T xxmin, const T xxmax, const T yymin, const T yymax);


    /**
      @brief Smooth grid using elliptic grid generation.
      @pre Grid2D object exists and grid has been initialized.
      @param niterations Number of iterations to be used.  Set to a negative number to
        disable this control.
      @param tol Relative error required to achieve in order to terminate iterations.
        Set to a negative number to disable this control.
      @param iterationsused Number of iterations used in the smoothing process.
      @param tolachieved Relative error achieved in the smoothing process.
      @param fixbound1 Sets if boundary 1 (constant @f$ \eta = 0 @f$ in computational domain) points can be modified.
      @param fixbound2 Sets if boundary 1 (constant @f$ \eta = 1 @f$ in computational domain) points can be modified.
      @param fixbound3 Sets if boundary 1 (constant @f$ \xi = 0 @f$ in computational domain) points can be modified.
      @param fixbound4 Sets if boundary 1 (constant @f$ \xi = 1 @f$ in computational domain) points can be modified.
      @post Grid smoothed.  Original grid points are lost.
      @return None.
      @warning If both niterations and tol are negative (i.e., both methods of control are disabled),
        a tolerance of 1.0e-4 will be used and a maximum-allowable number of iterations of 100,000.
     */
    void EllipticSmoothing(const int niterations, const T tol, int &iterationsused, T &tolachieved,
                           const bool fixbound1, const bool fixbound2, const bool fixbound3, const bool fixbound4);


	/**
     * @brief Output grid to data file for plotting.
     *
     * The resulting file output
	 * 			in VTK Structured Grid format (for flexibility between 2-
	 * 			dimensional grid types).
	 * @pre Grid2D object exists and the grid has been defined.
	 * @param file String containing the filename to contain the data.
	 * @post No object data changed.  File written to disc.
	 * @return None.
	 */
	void WriteGrid(std::string file);


	/**
     * @brief Output solution.  All quantities will be included in the data file.
	 * @pre Grid2D object exists and solution has been found.
	 * @param file String containing the filename to contain the data.
	 * @post No object data changed.  File written to disc.
	 * @return None.
	 */
	void WriteSolution(std::string file);


    /**
	 * @brief Get the number of points in the specified direction.  Number of
	 * 			points is set during the creation/definition of the grid.
	 * @pre Grid2D object exists.
     * @param dim Dimension to be returned.
            - 1:"y"
            - 2:"x"
	 * @post No change to object attributes.
	 * @return None.
	 */
	int GetSize(int dim);


	/**
	 * @brief Calculate the metrics needed to relate the spatial derivatives to
	 * 			the computational domain.
	 * @pre Grid2D object exists and grid has been generated.
	 * @post Metric arrays set.
	 * @return None.
	 */
	void ComputeMetrics();


    /**
      @brief Add an array to hold a scalar quantity at each grid point.
      @pre Grid2D object exists.
      @param name Name of scalar quantity added.
      @post New scalar quantity has been added.
      @return None.
     */
    void AddScalarQuantity(const std::string name);


    /**
      @brief Remove array containing scalar quantity at each grid point.
      @pre Grid2D object exists and scalar quantity exists.
      @param qty 0-based index of quantity to be removed.
      @post Scalar quantity removed.  Remaining quantities remain in their
            same relative order, but their index numbers are updated to
            remain contiguous.
      @return None.
      @warning All data stored in the quantity will be lost.
      @warning Be aware of the index renumbering.  Order remains the same, but
            index numbers are updated to remain contiguous.
     */
    void RemoveScalarQuantity(const int qty);


    /**
      @brief Set the name of a scalar quantity.
      @pre Grid2D object exists and scalar quantity exists.
      @param qty Index of scalar quantity to be named.
      @param name Name to use with scalar quantity.
      @post Scalar quantity name updated.
      @return None.
     */
    void setScalarQuantityName(const int qty, const std::string name);


    /**
      @brief Get the name of a scalar quantity.
      @pre Grid2D object exists and scalar quantity exists.
      @param qty Index of scalar quantity.
      @post No changes to object.
      @return Name of scalar quantity.
     */
    std::string ScalarQuantityName(const int qty);


    /**
      @brief Add an array to hold a vector quantity at each grid point.  It is assumed
            that the vector quantity has 3 components.
      @pre Grid2D object exists.
      @param name Name of vector quantity added.
      @post New vector quantity has been added.
      @return None.
     */
    void AddVectorQuantity(const std::string name);


    /**
      @brief Remove array containing vector quantity at each grid point.
      @pre Grid2D object exists and vector quantity exists.
      @param qty 0-based index of quantity to be removed.
      @post Vector quantity removed.  Remaining quantities remain in their
            same relative order, but their index numbers are updated to
            remain contiguous.
      @return None.
      @warning All data stored in the quantity will be lost.
      @warning Be aware of the index renumbering.  Order remains the same, but
            index numbers are updated to remain contiguous.
     */
    void RemoveVectorQuantity(const int qty);


    /**
      @brief Set the name of a vector quantity.
      @pre Grid2D object exists and vector quantity exists.
      @param qty Index of vector quantity to be named.
      @param name Name to use with vector quantity.
      @post Vector quantity name updated.
      @return None.
     */
    void setVectorQuantityName(const int qty, const std::string name);


    /**
      @brief Get the name of a vector quantity.
      @pre Grid2D object exists and vector quantity exists.
      @param qty Index of vector quantity.
      @post No changes to object.
      @return Name of vector quantity.
     */
    std::string VectorQuantityName(const int qty);


    /**
      @brief Get the number of scalar quantities in the grid.
      @pre Grid2D object exists.
      @post No changes to object.
      @return Number of scalar quantities in the grid.
    */
    int GetNumScalars() const;


    /**
      @brief Get the number of vector quantities in the grid.
      @pre Grid2D object exists.
      @post No changes to object.
      @return Number of vector quantities in the grid.
    */
    int GetNumVectors() const;


    /**
      @brief Set name for grid.  This is used as a descriptor in the VTK output files.
      @pre Grid2D object exists.
      @param name Name to be used for the grid.
      @post Grid name set.
      @return None.
     */
    void setGridName(const std::string name);


    /**
      @brief Get the name for the grid.  This is used as a descriptor in the VTK output files.
      @pre Grid2D object exists.
      @post No changes to object.
      @return Name of the grid.
     */
    std::string GridName() const;


    /**
     * @brief Overload () operator to access 2D quantities.
     * @pre Grid2D object exists.
     * @param ind1 Value of first index.
     * @param ind2 Value of second index.
     * @param qty Quantity to be accessed.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    T& operator()(int ind1, int ind2, int qty);


    /**
     * @brief Overload () operator to access 2D quantities.
     * @pre Grid2D object exists.
     * @param ind1 Value of first index.
     * @param ind2 Value of second index.
     * @param qty Quantity to be accessed.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    const T& operator()(int ind1, int ind2, int qty) const;


    /**
     * @brief Overload () operator to access 3D quantities.
     * @pre Grid2D object exists.
     * @param ind1 Value of first index.
     * @param ind2 Value of second index.
     * @param ind3 Value of third index.
     * @param qty Quantity to be accessed.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    T& operator()(int ind1, int ind2, int ind3, int qty);


    /**
     * @brief Overload () operator to access 3D quantities.
     * @pre Grid2D object exists.
     * @param ind1 Value of first index.
     * @param ind2 Value of second index.
     * @param ind3 Value of third index.
     * @param qty Quantity to be accessed.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    const T& operator()(int ind1, int ind2, int ind3, int qty) const;


    /**
     * @brief Copy-assignment operator.
     * @param a Grid2D object being assigned.
     * @return Reference to instance of Grid2D.
     */
    Grid2D& operator=(Grid2D<T> a);

    /**
      @brief Copy the grid node coordinates from one grid to this grid.
      @pre Grid2D object exists.
      @param sgrid Source grid containing the nodal coordinates which are to be applied
        to this grid.
      @post Grid coordinates of this grid modified.
      @return None.
    */
    void CopyNodalCoordinates(Grid2D<T> &sgrid);


    /**
      @brief Estimate the memory required to store this Grid2D object.
      @pre Grid2D object exists.
      @post No changes to object.
      @return Number of bytes required to store this Grid2D object and member data.
      @warning This routine uses the sizeof() command, which assumes that 1 char variable requires
        1 byte.  If 1 char variable requires a different amount of memory the value returned by this
        function must be scaled accordingly.
    */
    double EstimateMemoryUsage() const;


    /**
      @brief Implement MPI Broadcast for this Grid2D object.

        If MPI is not used, this function has no effect.
        Note that this requires a Grid2D object (with the correct scalar and vector quantities) to exist on each
        MPI process and any values stored in the objects on the non-root process will be lost.
      @pre Grid2D object exists.
      @param src Rank of MPI process which is the source of the data (all other processes receive).
      @post Grid2D data broadcast to all MPI processes.
      @return None.
      @warning This is currently only valid for H-grids.
    */
    void MPI_BCAST(const int src);


    /**
      @brief Implement MPI Broadcast for this Grid2D object.

        If MPI is not used, this function has no effect.
        Note that this requires a Grid2D object (with the correct scalar and vector quantities) to exist on each
        MPI process and any values stored in the objects on the non-root process will be lost.  This function will
        also broadcast the nodal coordinates and metrics.
      @pre Grid2D object exists.
      @param src Rank of MPI process which is the source of the data (all other processes receive).
      @post Grid2D data broadcast to all MPI processes.
      @return None.
      @warning This is currently only valid for H-grids.
    */
    void MPI_BCAST_GRID(const int src);


    /**
      @brief Implement MPI Gather for this Grid2D object.

        If MPI is not used, this function has no effect.
        Note that this requires a Grid2D object (with the correct scalar and vector quantities) to exist on each
        MPI process.  It is assumed that the grid points do not need to be gathered, just the scalar and vector
        quantity data.
      @pre Grid2D object exists.
      @param dest Rank of MPI process which is the destination of the data (all other processes receive).
      @param rank Rank of local MPI process.
      @param np Number of MPI processes.
      @post Grid2D data gathered from all MPI processes.
      @return None.
    */
    void MPI_GATHER(const int dest, const int rank, const int np);


    /** @brief 2D array containing dzi/dx */
    Array2D<T> dzidx;


    /** @brief 2D array containning dzi/dy */
    Array2D<T> dzidy;


    /** @brief 2D array containing deta/dx */
    Array2D<T> detadx;


    /** @brief 2D array containing deta/dy */
    Array2D<T> detady;


    /** @brief 2D array containing the x-coordinates of the nodes of the
     * 			computational grid */
    Array2D<T> xcoords;


    /** @brief 2D array containing the x-coordinates of the nodes of the
     * 			computational grid */
    Array2D<T> ycoords;


    /** @brief 2D array containing the I-blanking parameter for each grid point */
    Array2D<T> iblank;


    /** @brief 1D array containing the residual value as a function of solution
     * 			iteration */
    Array1D<T> converge;

protected:



private:
	/** @brief Number of computational grid points in the "x" direction */
	int nx;

	/** @brief Number of computational grid points in the "y" direction */
	int ny;

	/** @brief Number of iterations used to reach solution */
	int itersused;

	/** @brief Track if grid is O-Grid */
	bool isogrid;

	/** @brief Track if grid is H-Grid */
	bool ishgrid;

	/** @brief Minimum x dimension */
	T xmin;

	/** @brief Maximum x dimension */
	T xmax;

	/** @brief Minimum y dimension */
	T ymin;

	/** @brief Maximum y dimension */
	T ymax;

    /** @brief Value of the mathematical constant PI */
    T pi;

    /** @brief Number of scalar quantities at each point */
    int nscalars;

    /** @brief Number of vector quantities at each point */
    int nvectors;

    /** @brief Array of pointers to scalar quantity arrays */
    PArray1D<Array2D<T>*> pscalars;

    /** @brief Array of pointers to vectory quantity arrays */
    PArray1D<Array3D<T>*> pvectors;

    /** @brief List of names for scalar quantities */
    PArray1D<std::string*> scalar_names;

    /** @brief List of names for vector quantities */
    PArray1D<std::string*> vector_names;

    /** @brief Name of grid */
    std::string gridname;

    /** @brief Size reserved for quantity names */
    int qtysize;

    /** @brief Number of iterations between writing solution "snapshots" to disk */
    int inc_snapshots;

    /** @brief Number of iterations elapsed since last "snapshot" */
    int count_snapshots;

    /** @brief Name to be used for solution "snapshots" */
    std::string name_snapshots;


    /**
     * @brief Grid2DSwap swaps member information between two Grid2D objects.
     * @param first First Grid2D object.
     * @param second Second Grid2D object.
     */
    friend void Grid2DSwap(Grid2D<T> &first, Grid2D<T> &second)
    {
        std::swap(first.dzidx, second.dzidx);
        std::swap(first.detadx, second.detadx);
        std::swap(first.dzidy, second.dzidy);
        std::swap(first.detady, second.detady);
        std::swap(first.xcoords, second.xcoords);
        std::swap(first.ycoords, second.ycoords);
        std::swap(first.iblank, second.iblank);
        std::swap(first.converge, second.converge);
        std::swap(first.nx, second.nx);
        std::swap(first.ny, second.ny);
        std::swap(first.isogrid, second.isogrid);
        std::swap(first.ishgrid, second.ishgrid);
        std::swap(first.itersused, second.itersused);
        std::swap(first.xmin, second.xmin);
        std::swap(first.xmax, second.xmax);
        std::swap(first.ymin, second.ymin);
        std::swap(first.ymax, second.ymax);
        std::swap(first.pi, second.pi);
        std::swap(first.nscalars, second.nscalars);
        std::swap(first.nvectors, second.nvectors);
        std::swap(first.pscalars, second.pscalars);
        std::swap(first.pvectors, second.pvectors);
        std::swap(first.scalar_names, second.scalar_names);
        std::swap(first.vector_names, second.vector_names);
        std::swap(first.gridname, second.gridname);
        std::swap(first.qtysize, second.qtysize);
        std::swap(first.inc_snapshots, second.inc_snapshots);
        std::swap(first.count_snapshots, second.count_snapshots);
        std::swap(first.name_snapshots, second.name_snapshots);
    }


};


#include "Grid2D.cpp"

#endif /* Grid2D_ */
