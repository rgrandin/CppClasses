/**
 * @file PArray2D.h
 * @author Robert Grandin
 * @date 15 Oct 2010
 * @brief Definition of PArray2D class.
 *
 * @section Class Description & Notes
 *
 * This class contains the definition of a 2-dimensional array.  The low-level
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
 * 	- Modified to be derived from ArrayBase
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

#ifndef  PArray2D_
#define  PArray2D_


// REQUIRED INCLUDE FILES
#include <stdlib.h>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <assert.h>
#include <PArrayBase.h>



// ==================================================================
// ================
// ================    FUNCTION PROTOTYPES AND DESCRIPTIONS 
// ================

/**
 * @brief Class definition for storing 2-dimensional data.
 */
template <class T>
class PArray2D : public PArrayBase<T>{
  
  public:
	// CONSTRUCTOR
	/**
	 * @brief Default constructor for 2D array.  Default size is 1x1 and array
	 * 		elements are initialized to NULL
	 * @return None.
	 * @post Array object created and initialized to default value.
	 */
	PArray2D();


	/**
	 * @brief Constructor for 2D array initialized to NULL.
	 * @param dim1 Size of first dimension of the array.
	 * @param dim2 Size of the second dimension of the array.
	 * @return None.
	 * @post Array object created and initialized to default value.
	 */
	PArray2D(int dim1, int dim2);
	
	
	// DECONSTRUCTOR
	/**
	 * @brief Deconstructor for 2D array.
	 * @pre Array object exists.
	 * @return None.
	 * @post Array object destroyed.
	 */
    virtual ~PArray2D();

	
	
	// FUNCTIONS TO ACCESS AND MANIPULATE DATA
	/**
	 * @brief Retrieve the array extent in the specified dimension.
	 * @pre Array object exists.
	 * @param dim Dimension for which extent is desired.
	 * @return Array extent in the specified dimension.
	 * @post Array extents remain unchanged.
	 */
    int GetDim(int dim) const;

	
	/**
	 * @brief Reset array size. All points set to NULL.
	 * @pre Array object exists.
	 * @warning Any data already existing in the volume will be lost.
	 * @param dim1 New size of the first dimension.
	 * @param dim2 New size of the second dimension.
	 * @return None.
	 * @post Array size changed to dim1xdim2 and new points initialized to NULL.
	 */
	void ResetSize(int dim1, int dim2);


	/**
	 * @brief Reset all array points to a single value.
	 * @pre Array object exists.
	 * @param initvalue Value to be placed at all array points.
	 * @return None.
	 * @post All array points set to 'initvalue'.
	 */
	void ResetVal(const T initvalue);


	/**
	 * @brief Overload () operator.
	 * @pre Array object exists.
	 * @param ind1 Value of first index.
	 * @param ind2 Value of second index.
	 * @post No changes to object.
	 * @return Value stored at supplied indices.
	 */
	T& operator()(int ind1,int ind2);


	/**
	 * @brief Overload () operator.
	 * @pre Array object exists.
	 * @param ind1 Value of first index.
	 * @param ind2 Value of second index.
	 * @post No changes to object.
	 * @return Value stored at supplied indices.
	 */
	const T& operator()(int ind1,int ind2) const;


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
	/** @brief Number of points along the first dimension */
	int size1;

	/** @brief Number of points along the second dimension */
	int size2;


  
};

#include "PArray2D.cpp"

#endif
