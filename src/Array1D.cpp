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
Array1D<T>::Array1D(Array1D<T> &ab) : Array1D()
{
    Array1DSwap(*this, ab);
}


template <class T>
Array1D<T>::Array1D(Array1D<T> &&ab) : Array1D()
{
    Array1DSwap(*this, ab);
}

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
		return ArrayBase<T>::array[ind1];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return ArrayBase<T>::array[ind1];
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
		return ArrayBase<T>::array[ind1];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return ArrayBase<T>::array[ind1];
	#endif
}


template <class T>
Array1D<T>& Array1D<T>::operator=(const Array1D<T> &ab)
{
    Array1DSwap(*this, ab);
    return *this;
}


template <class T>
Array1D<T>& Array1D<T>::operator=(const Array1D<T> &&ab)
{
    Array1DSwap(*this, ab);
    return *this;
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
		return ArrayBase<T>::array[ind1];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return ArrayBase<T>::array[ind1];
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
		ArrayBase<T>::array[ind1] = value;
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		ArrayBase<T>::array[ind1] = value;
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
