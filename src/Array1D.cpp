/**
 * @file Array1D.cpp
 * @author Robert Grandin
 * @brief Implementation of Array1D class.
 */


#include <Array1D.h>



// ==================================================================
// ================
// ================    PRIVATE FUNCTIONS
// ================




// ==================================================================
// ================
// ================    PUBLIC FUNCTIONS
// ================

// CONSTRUCTORS AND DECONSTRUCTORS
template <class T>
Array1D<T>::Array1D()
{
}


template <class T>
Array1D<T>::Array1D(size_t dim1)
{
	ArrayBase<T>::ResetSize(dim1,(T)0.0e0);
}


template <class T>
Array1D<T>::Array1D(size_t dim1, const T initvalue)
{
	ArrayBase<T>::ResetSize(dim1,initvalue);
}


template <class T>
Array1D<T>::Array1D(const Array1D<T> &a) : ArrayBase<T>(a)
{
    /* No actions required here.  Calling the copy constructor of ArrayBase
     * takes care of the rest of the work. */
}


#ifdef CXX11
template <class T>
Array1D<T>::Array1D(Array1D<T> &&a) : ArrayBase<T>(std::move(a))
{
    Array1DSwap(*this, a);
}
#endif


template <class T>
Array1D<T>::~Array1D()
{
}


// () OPERATOR
template <class T> inline
T& Array1D<T>::operator()(size_t ind1)
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
        assert(ind1 < ArrayBase<T>::npoints);
        return ArrayBase<T>::p_array[ind1];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
        return ArrayBase<T>::p_array[ind1];
	#endif
}


template <class T> inline
const T& Array1D<T>::operator()(size_t ind1) const
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < ArrayBase<T>::npoints);
        return ArrayBase<T>::p_array[ind1];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
        return ArrayBase<T>::p_array[ind1];
	#endif
}


template <class T>
Array1D<T>& Array1D<T>::operator=(Array1D<T> a)
{
#ifdef CXX11
    ArrayBase<T>::operator=(static_cast<ArrayBase<T>>(a));
    Array1DSwap(*this, a);
    return *this;
#else
    std::cerr << "Array1D<T>::operator= ERROR: Not implemented without C++11 support!" << std::endl;
#endif
}


template <class T> inline
size_t Array1D<T>::GetDim() const
{
	return ArrayBase<T>::npoints;
}

template <class T> inline
T Array1D<T>::GetVal(size_t ind1) const
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < ArrayBase<T>::npoints);
        return ArrayBase<T>::p_array[ind1];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
        return ArrayBase<T>::p_array[ind1];
	#endif
}

template <class T> inline
void Array1D<T>::SetVal(size_t ind1, T value)
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < ArrayBase<T>::npoints);
        ArrayBase<T>::p_array[ind1] = value;
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
        ArrayBase<T>::p_array[ind1] = value;
	#endif
}

template <class T>
void Array1D<T>::ResetVal(const T initval)
{
	ArrayBase<T>::ResetVal(initval);
}

template <class T>
void Array1D<T>::ResetSize(size_t dim1)
{
	ArrayBase<T>::ResetSize(dim1,0);
}

template <class T>
void Array1D<T>::ResetSize(size_t dim1, const T initvalue)
{
	ArrayBase<T>::ResetSize(dim1,initvalue);
}


template <class T>
double Array1D<T>::GetMemoryUsage() const
{
	double retval = 0.0e0;

	// NUMBER OF VALUES STORED IN THE ARRAY
	retval += ArrayBase<T>::GetMemoryUsage();

	return retval;
}


template <class T>
void Array1D<T>::SetArrayPointer(T *p_data, size_t npts, bool useFree)
{
    ArrayBase<T>::SetArrayPointer(p_data, npts, useFree);
}
