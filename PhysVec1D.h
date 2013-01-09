/**
 * @file PhysVec1D.h
 * @author Robert Grandin
 * @date 25 August 2011
 * @brief Definition of PhysVec1D class.
 *
 * @section Class Description & Notes
 *
 * This class is a specialization of 1D arrays that represents physical vectors.
 * Such vectors occur often in science and engineering, and common vector operations
 * are possible via the member functions of this class and overloaded operators.
 *
 * Note that in this class the term 'vector' refers to these physical vectors
 * and not the C++ container.
 *
 * Bounds-checking and error-checking during execution is performed if RELEASE
 * is not defined.  Defining RELEASE, either manually via 'define RELEASE' or
 * via the compiler using '-DRELEASE' on gcc, will disable bounds-checking.
 *
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 25 August 2011
 *	- Creation date.
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

#ifndef PhysVec1D_
#define PhysVec1D_


/*
 * CLASS DEFINITION HERE
 *
 */
// REQUIRED INCLUDE FILES
#include <stdlib.h>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <string>
#include <assert.h>
#include <Array1D.h>


/**
 * @brief Class to represent physical vectors, such as those that occur in science
 * and engineering.
 */
template <class T>
class PhysVec1D : public Array1D<T> {

public:
	/**
	 * @brief Default constructor creates a single-element vector (i.e., a scalar).
	 * @pre Sufficient memory exists.
	 * @post Single-element vector created.
	 * @return None.
	 */
	PhysVec1D();


	/**
	 * @brief Copy constructor.  Created vector matches the input vector.
	 * @pre Input vector has been created and sufficeint memory exists for the
	 * 		new object.
	 * @param vec Vector from which values will be copied.
	 * @post New object created with the same parameters and contents as the
	 * 		input vector.
	 * @return None.
	 * @warning This is a true copy.  After this object is created, changes to the
	 * 		input vector will not be reflected in this vector.
	 */
	PhysVec1D(PhysVec1D<T> &vec);


	/**
	 * @brief Generalized constructor.
	 * @pre Sufficient memory exists.
	 * @param n Number of elements in the vector.
	 * @param initval Value to which all elements are initialized.
	 * @post Vector created as-specified.
	 * @return None.
	 */
	PhysVec1D(int n, const T initval);


	/**
	 * @brief Deconstructor.
	 * @pre PhysVec1D object exists.
	 * @post PhysVec1D object destroyed.
	 * @return None.
	 * @warning This is a deep-delete.  All data is lost and any references to this
	 * 		object will be broken.
	 */
    virtual ~PhysVec1D();


	/**
	 * @brief Assignment operator.
	 * @pre PhysVec1D object exists.
	 * @param vec Reference to vector from which values are assigned.
	 * @post Values assigned.
	 * @return Reference to self, which has been updated to match vectors.
	 */
	PhysVec1D<T>& operator=(const PhysVec1D<T>& vec);


	/**
	 * @brief Calculate the norm of the vector.
	 * @pre PhysVec1D object exists.
	 * @param p Norm to be calculated.  Use any negative number for infinity norm.
	 * @post No changes to object.
	 * @return Specified norm of the vector.
	 */
	const T Norm(const T p) const;


	/**
     * @brief Calculate the dot-product of this vector with another.
     *
     * Specifically, the operation is (this) dot (other_vector).
	 * @pre PhysVec1D objects exist.
	 * @param vec Vector with which this vector is dotted.
	 * @post No changes to objects.
	 * @return Dot product result.
	 */
	const T Dot(const PhysVec1D<T> &vec) const;


	/**
     * @brief Calculate the cross-product of this vector with another.
     *
     * Specifically, the operation is (this) cross (other_vector).
	 * @pre PhysVec1D objects exits.
	 * @param ivec Vector with which this is crossed.
	 * @param ovec Vector which contains result.
	 * @post No changes to this vector or input vector.
	 * @return None.
	 * @warning If input vectors are 2-element vectors, the result contains 3-elements.
	 * 		Use of 3-element input vectors produce a 3-element output vector.
	 */
	void Cross(const PhysVec1D<T> &ivec, PhysVec1D<T> &ovec);


	/**
     * @brief Rotate the vector in three dimensions.
     *
     * Vector magnitude is unchanged.  Rotations
     * 	are applied about the first axis, then the second axis, and finally the third axis.
	 * @pre PhysVec1D object exists.
	 * @param alpha Rotation, in radians, about the first (commonly 'x') axis.
	 * @param beta Rotation, in radians, about the second (commonly 'y') axis.
	 * @param gamma Rotation, in radians, about the third (commonly 'z') axis.
	 * @post Vector components updated to reflect the new orientation.
	 * @warning Original component values are lost.
	 * @warning Positive rotations defined using the right-hand rule.
	 * @return None.
	 */
	void Rotate(const T alpha, const T beta, const T gamma);



};


#include "PhysVec1D.cpp"

#endif /* PhysVec1D_ */
