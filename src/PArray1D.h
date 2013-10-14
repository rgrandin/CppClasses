/**
 * @file PArray1D.h
 * @author Robert Grandin
 * @date 15 Oct 2010
 * @brief Definition of PArray1D class.
 *
 *
 * @section Revisions
 *
 * @date 15 October 2010
 *	- Creation date.
 * @date 4 August 2011
 * 	- Improving documentation.
 * 	- Adding bounds-checking.
 * @date 17 August 2011
 * 	- Made bounds-checking conditional upon the presence of RELEASE compiler symbol.
 * 	  Defining RELEASE will disable bounds-checking.
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

#ifndef  PArray1D_
#define  PArray1D_


// REQUIRED INCLUDE FILES
#include <stdlib.h>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <string>
#include <assert.h>
#include <PArrayBase.h>



// ==================================================================
// ================
// ================    FUNCTION PROTOTYPES AND DESCRIPTIONS
// ================

/**
 * @brief 1-dimensional array of pointers.
 *
 * This class contains the definition of a 1-dimensional array.  The low-level
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
class PArray1D : public PArrayBase<T>{

public:
	/**
	 * @brief Constructor for 1D array.  Created array contains 1 element
	 * 		initialized to NULL.
	 * @return None.
	 * @post Array object created and initialized to default value.
	 */
	PArray1D();


	/**
	 * @brief Constructor for 1D array.  Array elements initialized to NULL.
	 * @return None.
	 * @post Array object created and initialized to default value.
	 */
    PArray1D(size_t dim1);


    /**
     * @brief Copy constructor.
     * @param a Reference to existing PArray1D object to be copied.
     */
    PArray1D(const PArray1D<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing PArray1D object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    PArray1D(PArray1D<T> &&a);
#endif


	/**
	 * @brief Destructor for 1D array.
     * @pre PArray1D object exists.
	 * @return None.
	 * @post Array object destroyed.  All object data is deleted.
	 * @warning This is a deep-delete.
	 */
    virtual ~PArray1D();


	/**
	 * @brief Retrieve the number of elements in the array.
	 * @pre Array object exists.
	 * @return Number of elements in the array.
	 * @post Array object remains unchanged.
	 */
    size_t GetDim() const;


	/**
	 * @brief Reset all array points to a single value.
	 * @pre Array object exists.
	 * @param initvalue Value to be placed at all array points.
	 * @return None.
	 * @post All array points set to 'initvalue'.
	 */
	void ResetVal(const T initvalue);


	/**
	 * @brief Access the array using indices contained within parentheses: ().
	 * @pre Array object exists.
	 * @param ind1 Value of first index.
	 * @post No changes to object.
	 * @return Value stored at supplied indices.
	 * @warning If ind1 does not fall within valid array bounds, program execution
	 * 		will terminate and output an error message to cerr.  Bounds-checking
	 * 		is disabled if 'NDEBUG' or 'RELEASE' has been defined.
	 */
    T& operator()(size_t ind1);


	/**
	 * @brief Access the array using indices contained within parentheses: ().
	 * @pre Array object exists.
	 * @param ind1 Value of first index.
	 * @post No changes to object.
	 * @return Value stored at supplied indices.
	 * @warning If ind1 does not fall within valid array bounds, program execution
	 * 		will terminate and output an error message to cerr.  Bounds-checking
	 * 		is disabled if 'NDEBUG' or 'RELEASE' has been defined.
	 */
    const T& operator()(size_t ind1) const;


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
     * @brief Copy-assignment operator.
     * @param a Reference to PArray1D object being assigned.
     * @return Reference to instance of PArray1D.
     */
    PArray1D& operator=(PArray1D<T> a);


private:

    /**
     * @brief PArray1DSwap swaps member information between two PArray1D objects.
     * @param first First PArray1D object.
     * @param second Second PArray1D object.
     */
    friend void PArray1DSwap(PArray1D<T> &first, PArray1D<T> &second)
    {
        /* No member functions for this derived class. */
    }


};

#include "PArray1D.cpp"

#endif
