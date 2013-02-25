/**
 * @file PArray2D.cpp
 * @author Robert Grandin
 * @brief Implementation of PArray2D class.
 */


#include <PArray2D.h>


// ==================================================================
// ================
// ================    PRIVATE FUNCTIONS
// ================

/** @brief Macro to map 2D indices to serialized index. */
#define IND1_IND2 ind1*size2+ind2



// ==================================================================
// ================
// ================    PUBLIC FUNCTIONS
// ================

// CONSTRUCTORS AND DESTRUCTORS
template <class T>
PArray2D<T>::PArray2D()
{
	size1 = 1;
	size2 = 1;
	PArrayBase<T>::npoints = size1*size2;
}

template <class T>
PArray2D<T>::PArray2D(size_t dim1, size_t dim2)
{
	size1 = dim1;
	size2 = dim2;
	PArrayBase<T>::ResetSize(size1*size2);
}


template <class T>
PArray2D<T>::PArray2D(const PArray2D<T> &a) : PArrayBase<T>(a),
    size1(a.size1), size2(a.size2), npoints(a.npoints)
{
}


#ifdef CXX11
template <class T>
PArray2D<T>::PArray2D(PArray2D<T> &&a) : PArrayBase<T>(std::move(a))
{
    PArray2DSwap(*this, a);
}
#endif


template <class T>
PArray2D<T>::~PArray2D()
{
	;
}


// () OPERATOR
template < class T > inline
T& PArray2D<T>::operator()(size_t ind1, size_t ind2)
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < size1);
		assert(ind2 >= 0 && ind2 < size2);
        return PArrayBase<T>::p_array[IND1_IND2];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
        return PArrayBase<T>::p_array[IND1_IND2];
	#endif
}

template < class T > inline
const T& PArray2D<T>::operator()(size_t ind1, size_t ind2) const
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < size1);
		assert(ind2 >= 0 && ind2 < size2);
        return PArrayBase<T>::p_array[IND1_IND2];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
        return PArrayBase<T>::p_array[IND1_IND2];
	#endif
}


template <class T>
PArray2D<T>& PArray2D<T>::operator=(PArray2D<T> a)
{
    PArrayBase<T>::operator=(static_cast<PArrayBase<T>>(a));
    PArray2DSwap(*this, a);
    return *this;
}



// DATA ACCESS AND MODIFICATION FUNCTIONS
template <class T>
size_t PArray2D<T>::GetDim(int dim) const
{
	#ifndef RELEASE
		assert(dim > 0 && dim < 3);
        size_t retval = 0;
		if(dim == 1){retval = size1;}
		if(dim == 2){retval = size2;}
		return retval;
	#else
        size_t retval = 0;
		if(dim == 1){retval = size1;}
		if(dim == 2){retval = size2;}
		return retval;
	#endif
}


template <class T>
void PArray2D<T>::ResetVal(const T initval)
{
	for(int i=0; i<PArrayBase<T>::npoints; i++){
        PArrayBase<T>::p_array[i] = initval;
	}
}


template <class T>
void PArray2D<T>::ResetSize(size_t dim1, size_t dim2)
{
	// CHECK THAT INPUT BOUNDS ARE INDEED DIFFERENT THAN CURRENT BOUNDS BEFORE
	// ATTEMPTING TO RESIZE THE ARRAY.
	if(dim1 != size1 || dim2 != size2){
		size1 = dim1;
		size2 = dim2;

		PArrayBase<T>::ResetSize(size1*size2);
	} else {
		// IF INPUT BOUNDS MATCH EXISTING BOUNDS, RESET ALL ARRAY POINTS TO
		// NULL
		for(int i=0; i<PArrayBase<T>::npoints; i++){
            PArrayBase<T>::p_array[i] = NULL;
		}
	}
}


template <class T>
double PArray2D<T>::GetMemoryUsage() const
{
	double retval = 0.0e0;

	// MEMORY REQUIREMENT FOR DATA AND npoints (MEMBERS FROM SOURCE CLASS)
	retval += PArrayBase<T>::GetMemoryUsage();

	// SPACE REQUIRED FOR PARAMETERS NOT INCLUDED IN SOURCE CLASS
	retval += (double)sizeof(size1);
	retval += (double)sizeof(size2);

	return retval;
}

