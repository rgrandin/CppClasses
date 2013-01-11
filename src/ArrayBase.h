/**
 * @file ArrayBase.h
 * @author Robert Grandin
 * @date 17 Aug 2011
 * @brief Definition of ArrayBase class.
 *
 * @section Class Description & Notes
 *
 * This class contains the base members and functions for arrays.  A container for
 * the data is provided, along with functions for accessing the data and querying
 * its properties.
 *
 * Bounds-checking during array access is performed if RELEASE is not defined.
 * Defining RELEASE, either manually via 'define RELEASE' or via the compiler
 * using '-DRELEASE' on gcc, will disable bounds-checking.
 *
 * Assert statements are enabled by default, and can be disabled using
 * 'define NDEBUG' or '-DNDEBUG'.  When RELEASE is defined, assert statements
 * are not used.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 *
 * @section Revisions
 *
 * @date 17 August 2011
 * 	- Creation date
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2011, Robert Grandin
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

#ifndef  ArrayBase_
#define  ArrayBase_


// REQUIRED INCLUDE FILES
#include <stdlib.h>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <string>
#include <assert.h>
#include <limits>




// ==================================================================
// ================
// ================    FUNCTION PROTOTYPES AND DESCRIPTIONS
// ================

/**
 * @brief Class definition for storing 1-dimensional data.
 */
template <class T> class ArrayBase{

public:
	/**
	 * @brief Create the array, initialized as a single element initialized to 0.
	 * @return None.
	 * @post Array object created and initialized.
	 */
	ArrayBase();


	/**
	 * @brief Create the array of specified size, initialized to 0.
	 * @param dim1 Number of elements in the array.
	 * @return None.
	 * @post Array object created and initialized.
	 */
	ArrayBase(int dim1);


	/**
	 * @brief Create the array of specified size and initial value.
	 * @param dim1 Size of the array.
	 * @param initvalue Value to be placed at all array locations.
	 * @return None.
	 * @post Array object created and initialized.
	 */
	ArrayBase(int dim1, const T initvalue);


	/**
	 * @brief Destroy array.
	 * @pre Array object exists.
	 * @return None.
	 * @post Array object destroyed.
	 * @warning This is a deep-delete.  All object data is deleted.
	 */
    virtual ~ArrayBase();


	/**
     * @brief Reset array size.
     *
     * If input dimension matches the existing array
     * 	dimension, the array is set to 0 at all points (same
     * 	behavior as ResetVal member function).
	 * @pre Array object exists.
	 * @param dim1 New size of the first dimension.
	 * @return None.
	 * @post Array size changed to dim1 and new points initialized to
	 * 			'initvalue'.
	 * @warning Any data already existing in the array will be lost.
	 * @warning New size must be non-zero and positive.  This is checked with an
	 * 		'assert' statement and as such the check will not be performed when
	 * 		compiled with NDEBUG defined.
	 */
	void ResetSize(long long dim1);


	/**
     * @brief Reset array size.
     *
     * If input dimension matches the existing array
     * 	dimension, the array is set to 'initvalue' at all points (same
     * 	behavior as ResetVal member function).
	 * @pre Array object exists.
	 * @warning Any data already existing in the volume will be lost.
	 * @param dim1 New size of the first dimension.
	 * @param initvalue Value to be placed at new array points.
	 * @return None.
	 * @post Array size changed to dim1 and new points initialized to
	 * 			'initvalue'.
	 */
	void ResetSize(size_t dim1, const T initvalue);


	/**
	 * @brief Reset all array points to a single value.
	 * @pre Array object exists.
	 * @param initvalue Value to be placed at all array points.
	 * @return None.
	 * @post All array points set to 'initvalue'.
	 */
	void ResetVal(const T initvalue);


	/**
	 * @brief Calculate the arithmetic mean of the elements in the array.
	 * @pre Array object exists.
	 * @post No changes to array object.
	 * @return Arithmetic mean of elements in array.
	 * @warning An array of non-floating point integers will have the result
	 * 		typecast for the return value.  For example, an array of integers will
	 * 		require that this function return an integer, so any decimal precision
	 * 		would be lost.  If a non-integer mean is desired when working with
	 * 		integer arrays, see function MeanFloat().
	 * @warning This function has no meaning for non-numeric data (e.g., strings or
	 * 		chars).  Results may be erratic on such data.
	 */
    T Mean() const;


	/**
	 * @brief Calculate the arithmetic mean of the elements in the array.
	 * @pre Array object exists.
	 * @post No changes to array object.
	 * @return Arithmetic mean of elements in array.
	 * @warning This version returns a type-'float' value, which could be useful
	 * 		for integer arrays.
	 * @warning This function has no meaning for non-numeric data (e.g., strings or
	 * 		chars).  Results may be erratic on such data.
	 */
    float MeanFloat() const;


	/**
	 * @brief Calculate the population variance of the elements in the array.
	 * @pre Array object exists.
	 * @post No changes to array object.
	 * @return Variance of elements in array.
	 * @warning An array of non-floating point integers will have the result
	 * 		typecast for the return value.  For example, an array of integers will
	 * 		require that this function return an integer, so any decimal precision
	 * 		would be lost.  If non-integer variance is desired when working with
	 * 		integer arrays, see VarianceFloat().
	 * @warning This function has no meaning for non-numeric data (e.g., strings or
	 * 		chars).  Results may be erratic on such data.
	 */
    T Variance() const;


	/**
	 * @brief Calculate the population variance of the elements in the array.
	 * @pre Array object exists.
	 * @post No changes to array object.
	 * @return Variance of elements in array.
	 * @warning This version returns a type-'float' value, which could be useful
	 * 		for integer arrays.
	 * @warning This function has no meaning for non-numeric data (e.g., strings or
	 * 		chars).  Results may be erratic on such data.
	 */
    float VarianceFloat() const;


	/**
	 * @brief Get the memory occupied by this object.
	 * @pre Array object exists.
	 * @post No changes made to object.
	 * @return Memory, in bytes, used by this object.
	 * @warning Type 'double' is used to ensure that the size can be accurately
	 * 		rendered.
	 */
    virtual double GetMemoryUsage() const;


	/**
	 * @brief Find the minimum value in the array.
	 * @pre Array object exists.
	 * @post No changes to object.
	 * @return Minimum value in the array.
	 * @warning This is not an absolute-value minimum, so "-3" is considered less
	 * 		than "-2".
	 */
    virtual T MinVal() const;


	/**
	 * @brief Find the minimum value in the array, as well as its location.
	 * @pre Array object exists.
	 * @param loc Reference to variable which will receive the location value.
	 * @post No changes to object.
	 * @return Minimum value in the array.  Location is returned via the referenced
	 * 		input variable.
	 * @warning This is not an absolute-value minimum, so "-3" is considered less
	 * 		than "-2".
	 */
    virtual T MinVal(size_t &loc) const;


	/**
	 * @brief Find the maximum value in the array.
	 * @pre Array object exists.
	 * @post No changes to object.
	 * @return Maximum value in the array.
	 * @warning This is not an absolute-value minimum, so "-2" is considered greater
	 * 		than "-3".
	 */
    virtual T MaxVal() const;


	/**
	 * @brief Find the maximum value in the array, as well as its location.
	 * @pre Array object exists.
	 * @param loc Reference to variable which will receive the location value.
	 * @post No changes to object.
	 * @return Maximum value in the array.  Location is returned via the referenced
	 * 		input variable.
	 * @warning If all values in the array are negative, this will not return the
	 * 		least-negative value as the maximum.  Instead, it will return the smallest-
	 * 		expressible non-negative number for the datatype.
	 */
	/* @warning This is not an absolute-value maximum, so "-2" is considered greater
	 * 		than "-3".
	 */
    virtual T MaxVal(size_t &loc) const;


	/**
	 * @brief Find the median value in the array.
	 * @pre Array object exists.
	 * @post No changes to object.
	 * @return Median value in the array.
	 */
    virtual T MedianVal() const;


    /**
     * @brief MemoryRequired calculates the memory required by an instance of this class.
     * @return Memory required, in bytes, of an instance of this class.
     * @warning Units of bytes is contingent upon sizeof(char) = 1.  If type char is of a
     *      different size the value returned by this function will need to be adjusted
     *      accordingly.
     */
    size_t MemoryRequired() const;





protected:

	/** @brief Pointer to actual array.  'Protected' status to allow access by
	 * derived classes. */
	T *array;


	/** @brief Number of points in the array */
    size_t npoints;


	/**
	 * @brief Initialize the array.
	 * @pre Sufficient memory exists for allocation.
	 * @param dim1 Number of elements in the array.
	 * @param initvalue Value to which all elements are initialized.
	 * @return None.
	 * @post Array created with specified size and initial value.
	 */
	void initialize(int dim1, const T initvalue);


	/**
     * @brief Function for comparing two input values (required by qsort).
     *
     * This simply identifies which input value is larger, with no speical
     * 	processing of the numbers.  Must be declared "static" to work
     * 	within this class.
	 * @pre Input data created.
	 * @param a First input value.
	 * @param b Second input value.
	 * @post No data modified.
	 * @return Integer identifying which input value is larger.
	 */
	static int compare(const void *a, const void *b);


private:

};

#include "ArrayBase.cpp"

#endif
