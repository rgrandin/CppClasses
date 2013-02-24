/**
 * @file Grid1D.h
 * @author Robert Grandin
 * @date 26 November 2010
 * @brief Definition of Grid1D class.
 *
 * @section Class Description & Notes
 *
 * This class represents a 1-dimensional regular grid for use in
 * computational mechanics.
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
 * @date 23 January 2012
 *  - Made updates to improve class in preparation for AerE 547.
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

#ifndef Grid1D_
#define Grid1D_

#include <cmath>
#include <fstream>
#include <Array1D.h>
#include <Array2D.h>
#include <PArray1D.h>


/**
 * @brief One-dimensional grid generation and manipulation for use with
 * 			computational mechanics applications.
 * @warning C++11 features, such as move-constructor and move-assignment, require the symbol
 *  "CXX11" to be defined.
 */
template <class T> class Grid1D{

public:
	/**
     * @brief Create Grid1D object.  No grid parameters defined.  No arrays created
     *     for storing data at grid points.
	 * @pre Sufficient memory exists for object creation.
     * @post Grid1D object created.
	 * @return None.
	 */
    Grid1D();


    /**
     * @brief Copy constructor.
     * @param a Reference to existing Grid1D object to be copied.
     */
    Grid1D(const Grid1D<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing Grid1D object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    Grid1D(Grid1D<T> &&a);
#endif


	/**
	 * @brief Destructor.
     * @pre Grid1D object exists.
     * @post Grid1D object destroyed.
	 * @return None.
	 */
    virtual ~Grid1D();


    /**
      @brief Create uniformly-spaced grid.
      @pre Grid1D object exists.
      @param npts Number of points to be used in the grid.  This includes the
            endpoints, making npts-1 "elements" within the grid.
      @param xmin Minimium spatial coordinate of grid.
      @param xmax Maximum spatial coordinate of grid.
      @post Grid created.
      @return None.
     */
    void CreateUniformGrid(const int npts, const T xmin, const T xmax);


    /**
	 * @brief Output grid to data file for plotting.  The resulting file output
     * 			in VTK Structured Grid format.
     * @pre Grid1D object exists and the grid has been defined.
	 * @param file String containing the filename to contain the data.
	 * @post No object data changed.  File written to disc.
	 * @return None.
	 */
	void WriteGrid(std::string file);


	/**
     * @brief Output solution.  All quantities will be included in the data file.
     * @pre Grid1D object exists and solution has been found.
	 * @param file String containing the filename to contain the data.
	 * @post No object data changed.  File written to disc.
	 * @return None.
	 */
	void WriteSolution(std::string file);


	/**
	 * @brief Output convergence plot to comma-separated-value file for plotting.
     * @pre Grid1D object exists and has been solved.
	 * @param file String containing the filename to contain the data.
	 * @post No object data changed.  File written to disc.
	 * @return None.
	 */
	void WriteConverge(std::string file);


	/**
     * @brief Get the number of points in the grid.
     * @pre Grid1D object exists.
	 * @post No change to object attributes.
	 * @return None.
	 */
    int GetSize();


    /**
      @brief Add an array to hold a scalar quantity at each grid point.
      @pre Grid1D object exists.
      @param name Name of scalar quantity added.
      @post New scalar quantity has been added.
      @return None.
     */
    void AddScalarQuantity(const std::string name);


    /**
      @brief Remove array containing scalar quantity at each grid point.
      @pre Grid1D object exists and scalar quantity exists.
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
      @pre Grid1D object exists and scalar quantity exists.
      @param qty Index of scalar quantity to be named.
      @param name Name to use with scalar quantity.
      @post Scalar quantity name updated.
      @return None.
     */
    void setScalarQuantityName(const int qty, const std::string name);


    /**
      @brief Get the name of a scalar quantity.
      @pre Grid1D object exists and scalar quantity exists.
      @param qty Index of scalar quantity.
      @post No changes to object.
      @return Name of scalar quantity.
     */
    std::string ScalarQuantityName(const int qty);


    /**
      @brief Add an array to hold a vector quantity at each grid point.  The vector
            quantity is assumed to have 3 components.
      @pre Grid1D object exists.
      @param name Name of vector quantity added.
      @post New vector quantity has been added.
      @return None.
     */
    void AddVectorQuantity(const std::string name);


    /**
      @brief Remove array containing vector quantity at each grid point.
      @pre Grid1D object exists and vector quantity exists.
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
      @pre Grid1D object exists and vector quantity exists.
      @param qty Index of vector quantity to be named.
      @param name Name to use with vector quantity.
      @post Vector quantity name updated.
      @return None.
     */
    void setVectorQuantityName(const int qty, const std::string name);


    /**
      @brief Get the name of a vector quantity.
      @pre Grid1D object exists and vector quantity exists.
      @param qty Index of vector quantity.
      @post No changes to object.
      @return Name of vector quantity.
     */
    std::string VectorQuantityName(const int qty);


    /**
      @brief Set name for grid.  This is used as a descriptor in the VTK output files.
      @pre Grid1D object exists.
      @param name Name to be used for the grid.
      @post Grid name set.
      @return None.
     */
    void setGridName(const std::string name);


    /**
      @brief Get the name for the grid.  This is used as a descriptor in the VTK output files.
      @pre Grid1D object exists.
      @post No changes to object.
      @return Name of the grid.
     */
    std::string GridName() const;


    /**
     * @brief Overload () operator to access scalar quantities.
     * @pre Grid1D object exists.
     * @param ind1 Value of first index.
     * @param qty Quantity to be accessed.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    T& operator()(int ind1, int qty);


    /**
     * @brief Overload () operator to access scalar quantities.
     * @pre Grid1D object exists.
     * @param ind1 Value of first index.
     * @param qty Quantity to be accessed.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    const T& operator()(int ind1, int qty) const;


    /**
     * @brief Overload () operator to access vector quantities.
     * @pre Grid1D object exists.
     * @param ind1 Value of first index.
     * @param ind2 Value of second index.
     * @param qty Quantity to be accessed.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    T& operator()(int ind1, int ind2, int qty);


    /**
     * @brief Overload () operator to access vector quantities.
     * @pre Grid1D object exists.
     * @param ind1 Value of first index.
     * @param ind2 Value of second index.
     * @param qty Quantity to be accessed.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    const T& operator()(int ind1, int ind2, int qty) const;


    /**
     * @brief Copy-assignment operator.
     * @param a Grid1D object being assigned.
     * @return Reference to instance of Grid1D.
     */
    Grid1D& operator=(Grid1D<T> a);


    /** @brief 2D array containing the x-coordinates of the nodes of the
     * 			computational grid */
    Array2D<T> xcoords;


    /** @brief 2D array containing the I-blanking parameter for each grid point */
    Array2D<T> iblank;


    /** @brief 1D array containing the residual value as a function of solution
     * 			iteration */
    Array1D<T> converge;

protected:



private:
	/** @brief Number of computational grid points in the "x" direction */
	int nx;

	/** @brief Number of iterations used to reach solution */
	int itersused;

	/** @brief Minimum x dimension */
	T xmin;

	/** @brief Maximum x dimension */
	T xmax;

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
     * @brief Grid1DSwap swaps member information between two Grid1D objects.
     * @param first First Grid1D object.
     * @param second Second Grid1D object.
     */
    friend void Grid1DSwap(Grid1D<T> &first, Grid1D<T> &second)
    {
        std::swap(first.xcoords, second.xcoords);
        std::swap(first.iblank, second.iblank);
        std::swap(first.converge, second.converge);
        std::swap(first.nx, second.nx);
        std::swap(first.itersused, second.itersused);
        std::swap(first.xmin, second.xmin);
        std::swap(first.xmax, second.xmax);
        std::swap(first.pi, second.pi);
        std::swap(first.nscalars, second.nscalars);
        std::swap(first.nvectors, second.nvectors);
        std::swap(first.pscalars, second.pvectors);
        std::swap(first.scalar_names, second.scalar_names);
        std::swap(first.vector_names, second.vector_names);
        std::swap(first.gridname, second.gridname);
        std::swap(first.qtysize, second.qtysize);
        std::swap(first.inc_snapshots, second.inc_snapshots);
        std::swap(first.count_snapshots, second.count_snapshots);
        std::swap(first.name_snapshots, second.name_snapshots);
    }

};


#include "Grid1D.cpp"

#endif /* Grid1D_ */
