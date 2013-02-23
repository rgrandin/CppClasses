/**
 * @file PArray1D.cpp
 * @author Robert Grandin
 * @brief Implementation of PArray1D class.
 */


#include <PArray1D.h>



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
PArray1D<T>::PArray1D()
{
	;
}

template <class T>
PArray1D<T>::PArray1D(size_t dim1)
{
	PArrayBase<T>::ResetSize(dim1);
}


template <class T>
PArray1D<T>::PArray1D(PArray1D<T> &a) : PArray1D()
{
    PArray1DSwap(*this, a);
}


template <class T>
PArray1D<T>::PArray1D(PArray1D<T> &&a) : PArray1D()
{
    PArray1DSwap(*this, a);
}


template <class T>
PArray1D<T>::~PArray1D()
{
	;
}


// () OPERATOR
template <class T> inline
T& PArray1D<T>::operator()(size_t ind1)
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
        assert(ind1 < PArrayBase<T>::npoints);
		return PArrayBase<T>::array[ind1];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return PArrayBase<T>::array[ind1];
	#endif
}


template <class T> inline
const T& PArray1D<T>::operator()(size_t ind1) const
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < PArrayBase<T>::npoints);
		return PArrayBase<T>::array[ind1];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return PArrayBase<T>::array[ind1];
	#endif
}


template <class T>
PArray1D<T>& PArray1D<T>::operator=(const PArray1D<T> &a)
{
    PArray1DSwap(*this, a);
    return *this;
}


template <class T>
PArray1D<T>& PArray1D<T>::operator=(const PArray1D<T> &&a)
{
    PArray1DSwap(*this, a);
    return *this;
}


template <class T> inline
size_t PArray1D<T>::GetDim() const
{
	return PArrayBase<T>::npoints;
}


template <class T>
void PArray1D<T>::ResetVal(const T initval)
{
	PArrayBase<T>::ResetVal(initval);
}

template <class T>
void PArray1D<T>::ResetSize(size_t dim1)
{
	PArrayBase<T>::ResetSize(dim1);
}


template <class T>
double PArray1D<T>::GetMemoryUsage() const
{
	double retval = 0.0e0;

	// NUMBER OF VALUES STORED IN THE ARRAY
	retval += PArrayBase<T>::GetMemoryUsage();

	return retval;
}
