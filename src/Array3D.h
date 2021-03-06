/**
 * @file Array3D.h
 * @author Robert Grandin
 * @date 15 Oct 2010
 * @brief Definition of Array3D class.
 *
 *
 * @section Revisions
 *
 * @date 15 October 2010
 *	- Creation date.
 * @date 17 August 2011
 * 	- Changed to be derived from ArrayBase class.
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

#ifndef  Array3D_
#define  Array3D_

// REQUIRED INCLUDE FILES
#include <stdlib.h>
#include <cstdlib>
#include <memory>
#include <cmath>

#ifdef GENCLASSES_FFTW_TRANSPOSE
#include <fftw3.h>
#endif

#include <ArrayBase.h>



/**
 * @brief Class definition for storing 3-dimensional data.
 *
 * This class contains the definition of a 3-dimensional array.  The low-level
 * memory management functions are all handled internally and array elements are
 * accessed via public member functions and overloaded operators.  The object can
 * be queried for its size, both in number of points and memory required to hold
 * the data and its attributes.
 *
 * Bounds-checking during array access is performed if RELEASE is not defined.
 * Defining RELEASE, either manually via 'define RELEASE' or via the compiler
 * using '-DRELEASE' on gcc, will disable bounds-checking.
 *
 * The transpose function makes use of the FFTW library for single, double,
 * and long-double precisions.  This means that the following libraries must be
 * linked by the compiler: -lfftw3 -lfftw3f -lfftw3l.  To use the FFTW library to
 * perform the transposition, the symbol 'FFTW_TRANSPOSE' must be defined.  Failure
 * to define this symbol will result in a full-copy of the data being made, which
 * requires more memory than the FFTW operations.  The data-copy method is also
 * approximately 4x slower than the FFTW-based method.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 * @warning C++11 features, such as move-constructor and move-assignment, require the symbol
 *  "CXX11" to be defined.
 */
template <class T>
class Array3D : public ArrayBase<T>{
  
  public:
	/**
	 * @brief Constructor for 3D array using default size and initial value.
	 * @return None.
	 * @post Array object created and initialized to default value.
	 */
	Array3D();


	/**
	 * @brief Constructor for 3D array using default initial value.
	 * @param dim1 Size of first dimension of the array.
	 * @param dim2 Size of the second dimension of the array.
	 * @param dim3 Size of the third dimension of the array.
	 * @return None.
	 * @post Array object created and initialized to default value.
	 */
    Array3D(size_t dim1, size_t dim2, size_t dim3);


	/**
	 * @brief Constructor for 3D array using a user-specified initial value.
	 * @param dim1 Size of first dimension of the array.
	 * @param dim2 Size of the second dimension of the array.
	 * @param dim3 Size of the third dimension of the array.
	 * @param initvalue Value to be placed at all array locations.
	 * @return None.
	 * @post Array object created and initialized to default value.
	 */
    Array3D(size_t dim1, size_t dim2, size_t dim3, const T initvalue);


    /**
     * @brief Copy constructor.
     * @param a Reference to existing Array3D object to be copied.
     */
    Array3D(const Array3D<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing Array3D object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    Array3D(Array3D<T> &&a);
#endif

	
	/**
	 * @brief Deconstructor for 3D array.
	 * @pre Array object exists.
	 * @return None.
	 * @post Array object destroyed.
	 */
    virtual ~Array3D();

	
	/**
	 * @brief Retrieve the array extent in the specified dimension.
	 * @pre Array object exists.
	 * @param dim Dimension for which extent is desired.
     *      - dim = 1: number of rows
     *      - dim = 2: number of columns
     *      - dim = 3: number of slices
     *      - other value: error
	 * @return Array extent in the specified dimension.
	 * @post Array extents remain unchanged.
	 */
	size_t GetDim(int dim) const;

	
	/**
	 * @brief Retrieve the value stored in array.
	 * @pre Array object exists.
	 * @param ind1 Index of array point in first dimension.
	 * @param ind2 Index of array point in second dimension.
	 * @param ind3 Index of array point in third dimension.
	 * @return Value stored in array at position (ind1,ind2,ind3).
	 * @post Array remains unchanged.
	 */
	T GetVal(int ind1, int ind2, int ind3) const;


	/**
	 * @brief Set the value stored in array.
	 * @pre Array object exists.
	 * @param ind1 Index of array point in first dimension.
	 * @param ind2 Index of array point in second dimension.
	 * @param ind3 Index of array point in third dimension.
	 * @param value Value to be placed in array at specified indices.
	 * @return None.
	 * @post Value at array point (ind1,ind2,ind3) set to 'value'.
	 */
	void SetVal(int ind1, int ind2, int ind3, const T value);

	
	/**
	 * @brief Reset all array points to a single value.
	 * @pre Array object exists.
	 * @param initvalue Value to be placed at all array points.
	 * @return None.
	 * @post All array points set to 'initvalue'.
	 */
	void ResetVal(const T initvalue);


	/**
     * @brief Reset array size.
     *
     * Previously-existing points retain their previous
     * 	value and newly-created points are initialized to 'initvalue'.
     * 	If input dimensions match existing array dimensions the array is
     * 	set to 'initvalue' at all points (same behavior as ResetVal
     * 	member function).
	 * @pre Array object exists.
	 * @warning Any data already existing in the volume will be lost.
	 * @param dim1 New size of the first dimension.
	 * @param dim2 New size of the second dimension.
	 * @param dim3 New size of the third dimension.
	 * @param initvalue Value to be placed at new array points.
	 * @return None.
	 * @post Array size changed to dim1xdim2xdim3 and new points initialized to
	 * 			'initvalue'.
	 */
    void ResetSize(size_t dim1, size_t dim2, size_t dim3, const T initvalue);


	/**
     * @brief Reset array size.
     *
     * Previously-existing points retain their previous
     * 	value and newly-created points are initialized to 0.0.
     * 	If input dimensions match existing array dimensions the array is
     * 	set to 'initvalue' at all points (same behavior as ResetVal
     * 	member function).
	 * @pre Array object exists.
	 * @warning Any data already existing in the volume will be lost.
	 * @param dim1 New size of the first dimension.
	 * @param dim2 New size of the second dimension.
	 * @param dim3 New size of the third dimension.
	 * @return None.
	 * @post Array size changed to dim1xdim2xdim3 and new points initialized to
	 * 			'initvalue'.
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
     * @param a Array3D object being assigned.
     * @return Reference to instance of Array3D.
     */
    Array3D& operator=(Array3D<T> a);


	/**
	 * @brief Get the memory occupied by this object.
	 * @pre Array object exists.
	 * @post No changes made to object.
	 * @return Memory, in bytes, used by this object.
	 * @warning Type 'double' is used to ensure that the size can be accurately
	 * 		rendered.
	 */
	double GetMemoryUsage() const;


	/**
	 * @brief Find the minimum value in the array, as well as its location.
	 * @pre Array object exists.
	 * @param loc1 Reference to variable which will receive the location value
	 * 		corresponding to the first array dimension.
	 * @param loc2 Reference to the variable which will receive the location value
	 * 		corresponding to the second array dimension.
	 * @param loc3 Reference to the variable which will receive the location value
	 * 		corresponding to the third array dimension.
	 * @post No changes to object.
	 * @return Minimum value in the array.  Location is returned via the referenced
	 * 		input variable.
	 * @warning This is not an absolute-value minimum, so "-3" is considered less
	 * 		than "-2".
	 */
    virtual T MinVal(size_t &loc1, size_t &loc2, size_t &loc3) const;


	/**
	 * @brief Find the maximum value in the array, as well as its location.
	 * @pre Array object exists.
	 * @param loc1 Reference to variable which will receive the location value
	 * 		corresponding to the first array dimension.
	 * @param loc2 Reference to the variable which will receive the location value
	 * 		corresponding to the second array dimension.
	 * @param loc3 Reference to the variable which will receive the location value
	 * 		corresponding to the third array dimension.
	 * @post No changes to object.
	 * @return Maximum value in the array.  Location is returned via the referenced
	 * 		input variable.
	 * @warning This is not an absolute-value maximum, so "-2" is considered greater
	 * 		than "-3".
	 */
    virtual T MaxVal(size_t &loc1, size_t &loc2, size_t &loc3) const;


    /**
     * @brief Transpose specified array dimensions.  Order of specified dimensions does not
     *  matter.
     *
     * Implementation taken from http://agentzlerich.blogspot.com/2010/01/using-fftw-for-in-place-matrix.html.
     * @param dim1 First dimension to be switched.  Must be 0, 1, or 2.
     * @param dim2 Second dimension to be switched.  Must be 0, 1, or 2.
     * @warning This function requires the use of numerical data.  Supported datatypes must fit on of the
     *  following criteria:
     *      - sizeof(datatype) == sizeof(float)
     *          - Integers should work in this case since sizeof(int) == sizeof(float)
     *      - sizeof(datatype) == sizeof(double)
     *      - sizeof(datatype) == sizeof(long double)
     *  The float routines require -lfftw3f, double routines require -lfftw3 (this is the default for FFTW),
     *  and long-double routines require -lfftw3l.  If any one of these libraries is not present on your system,
     *  the corresponding FFTW routines will not be available and unexpected results may occur.
     */
    void Transpose(const int dim1, const int dim2);


    /**
     * @brief Manually set the pointer to an existing array.
     * @param p_data Pointer to start of array.
     * @param npts1 Number of points in first dimension of array.
     * @param npts2 Number of points in second dimension of array.
     * @param npts3 Number of points in third dimension of array.
     * @param useFree If true, 'free()' will be used to free the memory rather than 'delete'.
     *  This is set based on how the array pointed-to by 'p_data' was allocated.
     */
    void SetArrayPointer(T* p_data, size_t npts1, size_t npts2, size_t npts3, bool useFree);


    /**
     * @brief Manually set array size, using existing pointer to serialized array.
     * @param npts1 Number of points in first dimension of array.
     * @param npts2 Number of points in second dimension of array.
     * @param npts3 Number of points in third dimension of array.
     * @param useFree If true, 'free()' will be used to free the memory rather than 'delete'.
     *  This is set based on how the array pointed-to by 'p_data' was allocated.
     */
    void SetArraySize(size_t npts1, size_t npts2, size_t npts3, bool useFree);


	
	
	
protected:
	// VARIABLES
	/** @brief Number of points along the first dimension. */
    size_t size1;

	/** @brief Number of points along the second dimension. */
    size_t size2;

	/** @brief Number of points along the third dimension. */
    size_t size3;



    /**
     * @brief Array3DSwap swaps member information between two ArrayBase objects.
     * @param first First ArrayBase object.
     * @param second Second ArrayBase object.
     */
    friend void Array3DSwap(Array3D<T> &first, Array3D<T> &second)
    {
        std::swap(first.size1, second.size1);
        std::swap(first.size2, second.size2);
        std::swap(first.size3, second.size3);
    }



  private:

  
};

#include "Array3D.cpp"



#endif
