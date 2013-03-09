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
    ArrayBase<T>::npoints = size1*size2*size3;
}


template <class T>
Array3D<T>::Array3D(int dim1, int dim2, int dim3)
{
    size1 = dim1;
    size2 = dim2;
    size3 = dim3;
    ArrayBase<T>::npoints = size1*size2*size3;
    ArrayBase<T>::ResetSize(ArrayBase<T>::npoints,(T)0.0e0);
}


template <class T>
Array3D<T>::Array3D(int dim1, int dim2, int dim3, const T initvalue)
{
    size1 = dim1;
    size2 = dim2;
    size3 = dim3;
    ArrayBase<T>::npoints = size1*size2*size3;
    ArrayBase<T>::ResetSize(ArrayBase<T>::npoints,initvalue);
}


template <class T>
Array3D<T>::Array3D(const Array3D<T> &a) : ArrayBase<T>(a),
    size1(a.size1),
    size2(a.size2),
    size3(a.size3)
{
}


#ifdef CXX11
template <class T>
Array3D<T>::Array3D(Array3D<T> &&a) : ArrayBase<T>(std::move(a))
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
    return ArrayBase<T>::p_array[IND1_IND2_IND3];
#else
    /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
    return ArrayBase<T>::p_array[IND1_IND2_IND3];
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
    return ArrayBase<T>::p_array[IND1_IND2_IND3];
#else
    /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
    return ArrayBase<T>::p_array[IND1_IND2_IND3];
#endif
}


template <class T>
Array3D<T>& Array3D<T>::operator=(Array3D<T> a)
{
    ArrayBase<T>::operator=(static_cast<ArrayBase<T>>(a));
    Array3DSwap(*this, a);
    return *this;
}


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
    return ArrayBase<T>::p_array[IND1_IND2_IND3];
#else
    /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
    return ArrayBase<T>::p_array[IND1_IND2_IND3];
#endif
}

template <class T>
void Array3D<T>::ResetVal(const T initval)
{
    for(size_t i=0; i<ArrayBase<T>::npoints; i++){
        ArrayBase<T>::p_array[i] = initval;
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
    ArrayBase<T>::p_array[IND1_IND2_IND3] = value;
#else
    /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
    ArrayBase<T>::p_array[IND1_IND2_IND3] = value;
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
        size_t npts = size1*size2*size3;

        ArrayBase<T>::ResetSize(npts, initvalue);
    } else {
        // IF INPUT BOUNDS MATCH EXISTING BOUNDS, RESET ALL ARRAY POINTS TO
        // 'initvalue'
        for(size_t i=0; i<ArrayBase<T>::npoints; i++){
            ArrayBase<T>::p_array[i] = initvalue;
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


template <class T>
void Array3D<T>::Transpose(const int dim1, const int dim2)
{
    /* Require valid values for dimensions to be transposed. */
    assert(dim1 >= 0 && dim1 < 3);
    assert(dim2 >= 0 && dim2 < 3);

    T tmpval;
    size_t tmpsize;
    size_t idx1, idx2;


    /* Transpose 1st and 2nd dimensions. */
    if((dim1 == 0 && dim2 == 1) || (dim1 == 1 && dim2 == 0)){

        for(size_t k=0; k<size3; k++){
            for(size_t i=0; i<size1; i++){
                for(size_t j=0; j<size2; j++){

                    idx1 = k*size1*size2 + i*size2 + j;
                    idx2 = k*size1*size2 + j*size1 + i;
                    tmpval = ArrayBase<T>::p_array[idx1];
                    ArrayBase<T>::p_array[idx1] = ArrayBase<T>::p_array[idx2];
                    ArrayBase<T>::p_array[idx2] = tmpval;

                }
            }
        }

        tmpsize = size1;
        size1 = size2;
        size2 = tmpsize;
    }


    /* Transpose 1st and 3rd dimensions. */
    if((dim1 == 0 && dim2 == 2) || (dim1 == 2 && dim2 == 0)){

        for(size_t k=0; k<size3; k++){
            for(size_t i=0; i<size1; i++){
                for(size_t j=0; j<size2; j++){

                    idx1 = k*size1*size2 + i*size2 + j;
                    idx2 = i*size3*size2 + k*size2 + j;
                    tmpval = ArrayBase<T>::p_array[idx1];
                    ArrayBase<T>::p_array[idx1] = ArrayBase<T>::p_array[idx2];
                    ArrayBase<T>::p_array[idx2] = tmpval;

                }
            }
        }

        tmpsize = size1;
        size1 = size3;
        size3 = tmpsize;
    }


    /* Transpose 2nd and 3rd dimensions. */
    if((dim1 == 1 && dim2 == 2) || (dim1 == 2 && dim2 == 1)){

        for(size_t k=0; k<size3; k++){
            for(size_t i=0; i<size1; i++){
                for(size_t j=0; j<size2; j++){

                    idx1 = k*size1*size2 + i*size2 + j;
                    idx2 = j*size1*size3 + i*size3 + k;
                    tmpval = ArrayBase<T>::p_array[idx1];
                    ArrayBase<T>::p_array[idx1] = ArrayBase<T>::p_array[idx2];
                    ArrayBase<T>::p_array[idx2] = tmpval;

                }
            }
        }

        tmpsize = size2;
        size2 = size3;
        size3 = tmpsize;
    }
}
