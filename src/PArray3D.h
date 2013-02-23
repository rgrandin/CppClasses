/**
 * @file PArray3D.h
 * @author Robert Grandin
 * @date 15 Oct 2010
 * @brief Definition of PArray3D class.
 *
 * @section Class Description & Notes
 *
 * This class contains the definition of a 3-dimensional array.  The low-level
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
 *
 *
 * @section Revisions
 *
 * @date 15 October 2010
 *	- Creation date.
 * @date 17 August 2011
 * 	- Changed to be derived from ArrayBase class.
 * @date 5 October 2011
 * 	- Made modifications to Array1D to make this PArray1D
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

#ifndef  PArray3D_
#define  PArray3D_

// REQUIRED INCLUDE FILES
#include <stdlib.h>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <PArrayBase.h>

/**
 * @brief Class definition for storing 3-dimensional array of pointers.
 * @warning C++11 features, such as move-constructor and move-assignment, require the symbol
 *  "CXX11" to be defined.
 */
template <class T>
class PArray3D : public PArrayBase<T>{
  
public:
	/**
	 * @brief Constructor for 3D array using default size and initial value.
	 * @return None.
	 * @post Array object created and initialized to default value.
	 */
	PArray3D();


	/**
	 * @brief Constructor for 3D array using default initial value.
	 * @param dim1 Size of first dimension of the array.
	 * @param dim2 Size of the second dimension of the array.
	 * @param dim3 Size of the third dimension of the array.
	 * @return None.
	 * @post Array object created and initialized to default value.
	 */
    PArray3D(size_t dim1, size_t dim2, size_t dim3);


    /**
     * @brief Copy constructor.
     * @param a Reference to existing PArray3D object to be copied.
     */
    PArray3D(PArray3D<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Reference to existing PArray3D object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    PArray3D(PArray3D<T> &&a);
#endif


	/**
	 * @brief Deconstructor for 3D array.
	 * @pre Array object exists.
	 * @return None.
	 * @post Array object destroyed.
	 */
    virtual ~PArray3D();

	
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
	 * @return None.
	 * @post Array size changed to dim1xdim2xdim3 and new points initialized to
	 * 			NULL.
	 */
    void ResetSize(size_t dim1, size_t dim2, size_t dim3);


	/**
	 * @brief Overload () operator.
	 * @pre Array object exists.
	 * @param dim1 Value of first index.
	 * @param dim2 Value of second index.
	 * @param dim3 Value of third index.
	 * @post No changes to object.
	 * @return Value stored at supplied indices.
	 */
    T& operator()(size_t dim1, size_t dim2, size_t dim3);


	/**
	 * @brief Overload () operator.
	 * @pre Array object exists.
	 * @param dim1 Value of first index.
	 * @param dim2 Value of second index.
	 * @param dim3 Value of third index.
	 * @post No changes to object.
	 * @return Value stored at supplied indices.
	 */
    const T& operator()(size_t dim1, size_t dim2, size_t dim3) const;


    /**
     * @brief Copy-assignment operator.
     * @param a Reference to PArray3D object being assigned.
     * @return Reference to instance of PArray3D.
     */
    PArray3D& operator=(PArray3D<T> a);


#ifdef CXX11
    /**
     * @brief Move-assignment operator (C++11).
     * @param a Reference to PArray3D object being assigned.
     * @return Reference to instance of PArray3D.
     * @warning This function requires C++11 compiler support.
     */
    PArray3D& operator=(PArray3D<T> &&a);
#endif


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
	int size1;

	/** @brief Number of points along the second dimension. */
	int size2;

	/** @brief Number of points along the third dimension. */
	int size3;



private:

    /**
     * @brief PArray3DSwap swaps member information between two PArray3D objects.
     * @param first First PArray3D object.
     * @param second Second PArray3D object.
     */
    friend void PArray3DSwap(PArray3D<T> &first, PArray3D<T> &second)
    {
        std::swap(first.npoints, second.npoints);
        std::swap(first.array, second.array);
        std::swap(first.size1, second.size1);
        std::swap(first.size2, second.size2);
        std::swap(first.size3, second.size3);
    }



  
};

#include "PArray3D.cpp"



#endif
