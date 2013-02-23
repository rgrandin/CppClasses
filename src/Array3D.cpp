/**
 * @file Array3D.cpp
 * @author Robert Grandin
 * @brief Implementation of Array3D class.
 */


#include <Array3D.h>


// ==================================================================
// ================
// ================    PRIVATE FUNCTIONS
// ================

/** @brief Macro to map 3D indices to serialized index. */
#define IND1_IND2_IND3 ind3*size1*size2+ind1*size2+ind2


// ==================================================================
// ================
// ================    PUBLIC FUNCTIONS
// ================

// CONSTRUCTORS AND DECONSTRUCTORS
template <class T>
Array3D<T>::Array3D()
{
    size1 = 1;
    size2 = 1;
    size3 = 1;
    npoints = size1*size2*size3;
}


template <class T>
Array3D<T>::Array3D(int dim1, int dim2, int dim3)
{
    size1 = dim1;
    size2 = dim2;
    size3 = dim3;
    npoints = size1*size2*size3;
    ArrayBase<T>::ResetSize(npoints,(T)0.0e0);
}


template <class T>
Array3D<T>::Array3D(int dim1, int dim2, int dim3, const T initvalue)
{
    size1 = dim1;
    size2 = dim2;
    size3 = dim3;
    npoints = size1*size2*size3;
    ArrayBase<T>::ResetSize(npoints,initvalue);
}


template <class T>
Array3D<T>::Array3D(Array3D<T> &a) : Array3D()
{
    Array3DSwap(*this, a);
}


#ifdef CXX11
template <class T>
Array3D<T>::Array3D(Array3D<T> &&a) : Array3D()
{
    Array3DSwap(*this, a);
}
#endif


template <class T>
Array3D<T>::~Array3D()
{
}


// () OPERATOR
template < class T > inline
T& Array3D<T>::operator()(size_t ind1, size_t ind2, size_t ind3)
{
#ifndef RELEASE
    /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
    assert(ind1 < size1);
    assert(ind2 < size2);
    assert(ind3 < size3);
    return ArrayBase<T>::array[IND1_IND2_IND3];
#else
    /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
    return ArrayBase<T>::array[IND1_IND2_IND3];
#endif
}

template < class T > inline
const T& Array3D<T>::operator()(size_t ind1, size_t ind2, size_t ind3) const
{
#ifndef RELEASE
    /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
    assert(ind1 >= 0 && ind1 < size1);
    assert(ind2 >= 0 && ind2 < size2);
    assert(ind3 >= 0 && ind3 < size3);
    return ArrayBase<T>::array[IND1_IND2_IND3];
#else
    /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
    return ArrayBase<T>::array[IND1_IND2_IND3];
#endif
}


template <class T>
Array3D<T>& Array3D<T>::operator=(Array3D<T> a)
{
    Array3DSwap(*this, a);
    return *this;
}


#ifdef CXX11
template <class T>
Array3D<T>& Array3D<T>::operator=(Array3D<T> &&a)
{
    Array3DSwap(*this, a);
    return *this;
}
#endif



// DATA ACCESS AND MODIFICATION FUNCTIONS
template <class T>
size_t Array3D<T>::GetDim(int dim) const
{
	#ifndef RELEASE
		assert(dim > 0 && dim < 4);
        size_t retval = 0;
		if(dim == 1){retval = size1;}
		if(dim == 2){retval = size2;}
		if(dim == 3){retval = size3;}
		return retval;
	#else
        size_t retval = 0;
		if(dim == 1){retval = size1;}
		if(dim == 2){retval = size2;}
		if(dim == 3){retval = size3;}
		return retval;
	#endif
}

template <class T> inline
T Array3D<T>::GetVal(int ind1, int ind2, int ind3) const
{
#ifndef RELEASE
    /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
    assert(ind1 >= 0 && ind1 < size1);
    assert(ind2 >= 0 && ind2 < size2);
    assert(ind3 >= 0 && ind3 < size3);
    return ArrayBase<T>::array[IND1_IND2_IND3];
#else
    /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
    return ArrayBase<T>::array[IND1_IND2_IND3];
#endif
}

template <class T>
void Array3D<T>::ResetVal(const T initval)
{
    for(size_t i=0; i<npoints; i++){
        ArrayBase<T>::array[i] = initval;
    }
}

template <class T> inline
void Array3D<T>::SetVal(int ind1, int ind2, int ind3,
                        T value)
{
#ifndef RELEASE
    /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
    assert(ind1 >= 0 && ind1 < size1);
    assert(ind2 >= 0 && ind2 < size2);
    assert(ind3 >= 0 && ind3 < size3);
    ArrayBase<T>::array[IND1_IND2_IND3] = value;
#else
    /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
    ArrayBase<T>::array[IND1_IND2_IND3] = value;
#endif
}

template <class T>
void Array3D<T>::ResetSize(size_t dim1, size_t dim2, size_t dim3)
{
    // CALL ResetSize WITHIN INITIAL VALUE SET TO "0".
    ResetSize(dim1,dim2,dim3,0);
}

template <class T>
void Array3D<T>::ResetSize(size_t dim1, size_t dim2, size_t dim3,
                           const T initvalue)
{
    // CHECK THAT INPUT BOUNDS ARE INDEED DIFFERENT THAN CURRENT BOUNDS BEFORE
    // ATTEMPTING TO RESIZE THE ARRAY.
    if(dim1 != size1 || dim2 != size2 ||  dim3 != size3){
        size1 = dim1;
        size2 = dim2;
        size3 = dim3;
        npoints = size1*size2*size3;

        ArrayBase<T>::ResetSize(npoints,initvalue);
    } else {
        // IF INPUT BOUNDS MATCH EXISTING BOUNDS, RESET ALL ARRAY POINTS TO
        // 'initvalue'
        for(size_t i=0; i<npoints; i++){
            ArrayBase<T>::array[i] = initvalue;
        }
    }
}


template <class T>
double Array3D<T>::GetMemoryUsage() const
{
    double retval = 0.0e0;

    // MEMORY REQUIREMENT FOR DATA AND npoints (MEMBERS FROM SOURCE CLASS)
    retval += ArrayBase<T>::GetMemoryUsage();

    // SPACE REQUIRED FOR PARAMETERS NOT INCLUDED IN SOURCE CLASS
    retval += (double)sizeof(size1);
    retval += (double)sizeof(size2);
    retval += (double)sizeof(size3);

    return retval;
}


template <class T>
T Array3D<T>::MinVal(size_t &loc1, size_t &loc2, size_t &loc3) const
{
    size_t loc = 0;
    T min = ArrayBase<T>::MinVal(loc);

    // CONVERT SERIALIZED LOCATION TO 3-INDEX LOCATION
    loc3 = (size_t)(loc/(size1*size2));
    loc1 = (size_t)((loc - loc3*size1*size2)/size2);
    loc2 = loc - loc3*size1*size2 - loc1*size2;

    return min;
}


template <class T>
T Array3D<T>::MaxVal(size_t &loc1, size_t &loc2, size_t &loc3) const
{
    size_t loc = 0;
    T max = ArrayBase<T>::MaxVal(loc);

    // CONVERT SERIALIZED LOCATION TO 3-INDEX LOCATION
    loc3 = (size_t)(loc/(size1*size2));
    loc1 = (size_t)((loc - loc3*size1*size2)/size2);
    loc2 = loc - loc3*size1*size2 - loc1*size2;

    return max;
}



