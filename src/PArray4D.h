/**
 * @file PArray4D.h
 * @author 	Robert Grandin
 * @date 27 October 2010
 * @brief Definition of PArray4D class.
 *
 *
 * @section Revisions
 *
 * @date 27 October 2010
 * 	- Class created.
 * @date 17 August 2011
 * 	- Updated to be derived from ArrayBase class.
 * @date 5 October 2011
 * 	- Made modifications to Array1D to make this PArray1D
 *
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2010-2011, Robert Grandin
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


#ifndef  PArray4D_
#define  PArray4D_

// REQUIRED INCLUDE FILES
#include <stdlib.h>
#include <cstdlib>
#include <memory>
#include <cmath>

#include <PArrayBase.h>

/**
 * @brief Class definition for storing 4-dimensional array of pointers.
 *
 * This class contains the definition of a 4-dimensional array.  The low-level
 * memory management functions are all handled internally and array elements are
 * accessed via public member functions and overloaded operators.  The object can
 * be queried for its size, both in number of points and memory required to hold
 * the data and its attributes.
 *
 * This array is intended to organize pointers, not numeric data.
 *
 * Bounds-checking during array access is performed if RELEASE is not defined.
 * Defining RELEASE, either manually via 'define RELEASE' or via the compiler
 * using '-DRELEASE' on gcc, will disable bounds-checking.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 * @warning C++11 features, such as move-constructor and move-assignment, require the symbol
 *  "CXX11" to be defined.
 */
template <class T>
class PArray4D : public PArrayBase<T>{
  
  public:
	// CONSTRUCTORS

	/**
	 * @brief Constructor for 3D array using default initial values.
	 * @return None.
	 * @post Array object created and initialized to default values.
	 */
	PArray4D();


	/**
	 * @brief Constructor for 3D array using default initial value.
	 * @pre None.
	 * @param dim1 Size of first dimension of the array.
	 * @param dim2 Size of second dimension of the array.
	 * @param dim3 Size of the third dimension of the array.
	 * @param dim4 Size of the fourth dimension of the array.
	 * @return None.
	 * @post Array object created and initialized to default value.
	 */
    PArray4D(size_t dim1, size_t dim2, size_t dim3, size_t dim4);


    /**
     * @brief Copy constructor.
     * @param a Reference to existing PArray4D object to be copied.
     */
    PArray4D(const PArray4D<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing PArray4D object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    PArray4D(PArray4D<T> &&a);
#endif


	// DECONSTRUCTOR
	/**
	 * @brief Deconstructor for 3D array.
	 * @pre Array object exists.
	 * @return None.
	 * @post Array object destroyed.
	 */
    virtual ~PArray4D();

	
	
	// FUNCTIONS TO ACCESS AND MANIPULATE DATA
	/**
	 * @brief Retrieve the array extent in the specified dimension.
	 * @pre Array object exists.
	 * @param dim Dimension for which extent is desired.
	 * @return Array extent in the specified dimension.
	 * @post Array extents remain unchanged.
	 */
    size_t GetDim(int dim) const;

	
	/**
	 * @brief Reset all array points to a single value.
	 * @pre Array object exists.
	 * @param initvalue Value to be placed at all array points.
	 * @return None.
	 * @post All array points set to 'initvalue'.
	 */
	void ResetVal(const T initvalue);


	/**
	 * @brief Reset array size. All elements set to NULL.
	 * @pre Array object exists.
	 * @warning Any data already existing in the volume will be lost.
	 * @param dim1 New size of the first dimension.
	 * @param dim2 New size of the second dimension.
	 * @param dim3 New size of the third dimension.
	 * @param dim4 New size of the fourth dimension.
	 * @return None.
	 * @post Array size changed to dim1xdim2xdim3 and new points initialized to
	 * 			'initvalue'.
	 */
    void ResetSize(size_t dim1, size_t dim2, size_t dim3, size_t dim4);


	/**
	 * @brief Overload () operator.
	 * @pre Array object exists.
	 * @param dim1 Value of first index.
	 * @param dim2 Value of second index.
	 * @param dim3 Value of third index.
	 * @param dim4 Value of fourth index.
	 * @post No changes to object.
	 * @return Value stored at supplied indices.
	 */
    T& operator()(size_t dim1, size_t dim2, size_t dim3, size_t dim4);


	/**
	 * @brief Overload () operator.
	 * @pre Array object exists.
	 * @param dim1 Value of first index.
	 * @param dim2 Value of second index.
	 * @param dim3 Value of third index.
	 * @param dim4 Value of fourth index.
	 * @post No changes to object.
	 * @return Value stored at supplied indices.
	 */
    const T& operator()(size_t dim1, size_t dim2, size_t dim3, size_t dim4) const;


    /**
     * @brief Copy-assignment operator.
     * @param a Reference to PArray4D object being assigned.
     * @return Reference to instance of PArray4D.
     */
    PArray4D& operator=(PArray4D<T> a);

	
	/**
	 * @brief Get the memory occupied by this object.
	 * @pre Array object exists.
	 * @post No changes made to object.
	 * @return Memory, in bytes, used by this object.
	 * @warning Type 'double' is used to ensure that the size can be accurately
	 * 		rendered.
	 */
    double GetMemoryUsage() const;


	
	
protected:
	// VARIABLES
	/** @brief Number of points along the first dimension. */
    size_t size1;

	/** @brief Number of points along the second dimension. */
    size_t size2;

	/** @brief Number of points along the third dimension. */
    size_t size3;

	/** @brief Number of points along the fourth dimension. */
    size_t size4;



private:

    /**
     * @brief PArray4DSwap swaps member information between two PArray4D objects.
     * @param first First PArray4D object.
     * @param second Second PArray4D object.
     */
    friend void PArray4DSwap(PArray4D<T> &first, PArray4D<T> &second)
    {
        std::swap(first.size1, second.size1);
        std::swap(first.size2, second.size2);
        std::swap(first.size3, second.size3);
        std::swap(first.size4, second.size4);
    }


  
};

#include "PArray4D.cpp"






#endif
