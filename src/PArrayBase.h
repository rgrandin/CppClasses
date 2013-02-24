/**
 * @file PArrayBase.h
 * @author Robert Grandin
 * @date 17 Aug 2011
 * @brief Definition of PArrayBase class.
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

#ifndef  PArrayBase_
#define  PArrayBase_


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
 * @brief Base class for arrays of pointers.
 * @warning C++11 features, such as move-constructor and move-assignment, require the symbol
 *  "CXX11" to be defined.
 */
template <class T> class PArrayBase{

public:
	/**
	 * @brief Create the array, initialized as a single element initialized to NULL.
	 * @return None.
	 * @post Array object created and initialized.
	 */
	PArrayBase();


	/**
	 * @brief Create the array of specified size, initialized to NULL.
	 * @param dim1 Number of elements in the array.
	 * @return None.
	 * @post Array object created and initialized.
	 */
    PArrayBase(size_t dim1);


    /**
     * @brief Copy constructor.
     * @param a Reference to existing PArrayBase object to be copied.
     */
    PArrayBase(const PArrayBase<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing PArrayBase object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    PArrayBase(PArrayBase<T> &&a);
#endif


	/**
	 * @brief Destroy array.
	 * @pre Array object exists.
	 * @return None.
	 * @post Array object destroyed.
	 * @warning This is a deep-delete.  All object data is deleted.
	 */
    virtual ~PArrayBase();


	/**
     * @brief Reset array size.
     *
     * If input dimension matches the existing array
     * 	dimension, the array is set to NULL at all points.
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
    void ResetSize(size_t dim1);


	/**
	 * @brief Reset all array points to a single value.
	 * @pre Array object exists.
	 * @param initvalue Value to be placed at all array points.
	 * @return None.
	 * @post All array points set to 'initvalue'.
	 */
	void ResetVal(const T initvalue);


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
     * @brief Copy-assignment operator.
     * @param a Reference to PArrayBase object being assigned.
     * @return Reference to instance of PArrayBase.
     */
    PArrayBase& operator=(PArrayBase<T> a);


    /**
     * @brief Array subscription operator.
     * @param idx Index to be accessed.
     * @return Value stored at array index.
     */
    T& operator [](const size_t idx);


    /**
     * @brief Array subscription operator.
     * @param idx Index to be accessed.
     * @return Value stored at array index.
     */
    const T& operator [](const size_t idx) const;




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
	 * @return None.
	 * @post Array created with specified size and initial value.
	 */
	void initialize(int dim1);



private:

    /**
     * @brief PArrayBaseSwap swaps member information between two PArrayBase objects.
     * @param first First PArrayBase object.
     * @param second Second PArrayBase object.
     */
    friend void PArrayBaseSwap(PArrayBase<T> &first, PArrayBase<T> &second)
    {
        std::swap(first.npoints, second.npoints);
        std::swap(first.array, second.array);
    }


};

#include "PArrayBase.cpp"

#endif
