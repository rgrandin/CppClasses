/**
 * @file Array2D.h
 * @author Robert Grandin
 * @date 15 Oct 2010
 * @brief Definition of Array2D class.
 *
 * @section Revisions
 *
 * This list may not be complete.
 *
 * @date 15 October 2010
 *	- Creation date.
 * @date 17 August 2011
 * 	- Modified to be derived from ArrayBase
 * @date February 2013
 *  - Addition of additional constructors and assignment operators.
 * @date 9 March 2013
 *  - Addition of transpose operation.
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2010-2013, Robert Grandin
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

#ifndef  Array2D_
#define  Array2D_


// REQUIRED INCLUDE FILES
#include <stdlib.h>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>

#ifdef GENCLASSES_FFTW_TRANSPOSE
#include <fftw3.h>
#endif

#include <ArrayBase.h>
#include <PArray1D.h>




// ==================================================================
// ================
// ================    FUNCTION PROTOTYPES AND DESCRIPTIONS 
// ================

/**
 * @brief Class definition for storing 2-dimensional data.
 *
 * This class contains the definition of a 2-dimensional array.  The low-level
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
 * perform the transposition, the symbol 'GENCLASSES_FFTW_TRANSPOSE' must be defined.  Failure
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
class Array2D : public ArrayBase<T>{

public:
    // CONSTRUCTOR
    /**
     * @brief Constructor for 2D array using default initial size and value.
     * @return None.
     * @post Array object created and initialized to default value.
     */
    Array2D();


    /**
     * @brief Constructor for 2D array using default initial value.
     * @param dim1 Size of first dimension of the array.
     * @param dim2 Size of the second dimension of the array.
     * @param initvalue Value to be placed at all array locations.
     * @return None.
     * @post Array object created and initialized to default value.
     */
    Array2D(size_t dim1, size_t dim2, const T initvalue);


    /**
     * @brief Copy constructor.
     * @param a Reference to existing Array2D object to be copied.
     */
    Array2D(const Array2D<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing Array2D object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    Array2D(Array2D<T> &&a);
#endif


    // DECONSTRUCTOR
    /**
     * @brief Deconstructor for 2D array.
     * @pre Array object exists.
     * @return None.
     * @post Array object destroyed.
     */
    virtual ~Array2D();



    // FUNCTIONS TO ACCESS AND MANIPULATE DATA
    /**
     * @brief Retrieve the array extent in the specified dimension.
     * @pre Array object exists.
     * @param dim Dimension for which extent is desired.
     *      - dim = 1: number of rows
     *      - dim = 2: number of columns
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
     * @return Value stored in array at position (ind1,ind2).
     * @post Array remains unchanged.
     */
    T GetValue(size_t ind1, size_t ind2) const;


    /**
     * @brief Set the value stored in array.
     * @pre Array object exists.
     * @param ind1 Index of array point in first dimension.
     * @param ind2 Index of array point in second dimension.
     * @param value Value to be placed in array at specified indices.
     * @return None.
     * @post Value at array point (ind1,ind2) set to 'value'.
     */
    void SetValue(size_t ind1, size_t ind2, const T value);


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
     * @return None.
     * @post Array size changed to dim1xdim2 and new points initialized to
     * 			'initvalue'.
     */
    void ResetSize(size_t dim1, size_t dim2);


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
     * @param initvalue Value to be placed at new array points.
     * @return None.
     * @post Array size changed to dim1xdim2 and new points initialized to
     * 			'initvalue'.
     */
    void ResetSize(size_t dim1, size_t dim2, const T initvalue);


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
    T& operator()(size_t ind1, size_t ind2);


    /**
     * @brief Overload () operator.
     * @pre Array object exists.
     * @param ind1 Value of first index.
     * @param ind2 Value of second index.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    const T& operator()(size_t ind1, size_t ind2) const;


    /**
     * @brief Overload () operator.
     * @pre Array object exists.
     * @param ind1 Value of first index.
     * @param ind2 Value of second index.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    const T& operator()(int ind1, int ind2) const;


    /**
     * @brief Overload () operator.
     * @pre Array object exists.
     * @param ind1 Value of first index.
     * @param ind2 Value of second index.
     * @post No changes to object.
     * @return Value stored at supplied indices.
     */
    T& operator()(int ind1, int ind2);


#ifdef CXX11
    /**
     * @brief Copy-assignment operator.
     * @param a Array2D object being assigned.
     * @return Reference to instance of Array2D.
     */
    Array2D& operator=(Array2D<T> a);
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


    /**
     * @brief Find the minimum value in the array, as well as its location.
     * @pre Array object exists.
     * @param loc1 Reference to variable which will receive the location value
     * 		corresponding to the first array dimension.
     * @param loc2 Reference to the variable which will receive the location value
     * 		corresponding to the second array dimension.
     * @post No changes to object.
     * @return Minimum value in the array.  Location is returned via the referenced
     * 		input variable.
     * @warning This is not an absolute-value minimum, so "-3" is considered less
     * 		than "-2".
     */
    virtual T MinVal(size_t &loc1, size_t &loc2) const;


    /**
     * @brief Find the maximum value in the array, as well as its location.
     * @pre Array object exists.
     * @param loc1 Reference to variable which will receive the location value
     * 		corresponding to the first array dimension.
     * @param loc2 Reference to the variable which will receive the location value
     * 		corresponding to the second array dimension.
     * @post No changes to object.
     * @return Maximum value in the array.  Location is returned via the referenced
     * 		input variable.
     * @warning This is not an absolute-value maximum, so "-2" is considered greater
     * 		than "-3".
     */
    virtual T MaxVal(size_t &loc1, size_t &loc2) const;


    /**
     * @brief Read comma-separated values (CSV) file into array.
     * @pre Array2D object exists.
     * @param filename Name of file to be read in.
     * @param nheader Number of header rows to be skipped.
     * @param mincols Minimum number of columns to be used.
     * @param defval Default value to be placed in unused locations.
     * @param linewidth Number of characters on a line.  This <B>must</B> be greater
     * 		than the longest line.  Default value is 4096 characters, which should
     * 		be large enough for most applications.
     * @post Data read-in to array.
     * @return None.
     * @warning Array resized to match the size of the input data.  Minimum number
     * 		of columns is the greater of user-specified value and number of columns
     * 		contained in the data file.  Any data currently stored in the array
     * 		will be lost.
     */
    void ReadCSVFile(const std::string filename, const int nheader, const int mincols,
                     const T defval, const int linewidth);


    /**
      @brief Write array values to disk as a comma-separated values (CSV) file.
      @pre Array2D object exists.
      @param filename Name of file to be written.  This is the absolute filename (i.e., includes
        the full path if the file is to be written outside the current directory).
      @param labels Reference to PArray1D object containing pointers to strings containing the labels for
        each column.  The ith string pointed-to by this array corresponds to the ith column of the array.
      @param ndec Number of decimals to use with scientific notation in output.  If set to 0 or a negative number,
        this parameter is ignored and the output format will be determined by the standard library.  If set to
        a positive number, the output numbers will be expressed in scientific notation with 'ndec' decimals used,
        regardless of if all the decimals are required (e.g., '0.5' with 4 decimal points would be '0.5000').
      @post No changes to object.
      @return None.
    */
    void WriteCSVFile(const std::string filename, const PArray1D<std::string *> &labels, const int ndec) const;


    /**
     * @brief Transpose transpose the array.
     *
     *  Implementation taken from http://agentzlerich.blogspot.com/2010/01/using-fftw-for-in-place-matrix.html.
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
    void Transpose();


    /**
     * @brief Manually set the pointer to an existing array.
     * @param p_data Pointer to start of array.
     * @param npts1 Number of points in first dimension of array.
     * @param npts2 Number of points in second dimension of array.
     * @param useFree If true, 'free()' will be used to free the memory rather than 'delete'.
     *  This is set based on how the array pointed-to by 'p_data' was allocated.
     */
    void SetArrayPointer(T* p_data, size_t npts1, size_t npts2, bool useFree);


protected:
    // VARIABLES
    /** @brief Number of points along the first dimension */
    size_t size1;

    /** @brief Number of points along the second dimension */
    size_t size2;




    /**
     * @brief Array2DSwap swaps member information between two ArrayBase objects.
     * @param first First ArrayBase object.
     * @param second Second ArrayBase object.
     */
    friend void Array2DSwap(Array2D<T> &first, Array2D<T> &second)
    {
        std::swap(first.size1, second.size1);
        std::swap(first.size2, second.size2);
    }


private:


};

#include "Array2D.cpp"

#endif
