/**
 * @file PArray3D.cpp
 * @author Robert Grandin
 * @brief Implementation of PArray3D class.
 */


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
PArray3D<T>::PArray3D()
{
	size1 = 1;
	size2 = 1;
	size3 = 1;
	PArrayBase<T>::npoints = size1*size2*size3;
}


template <class T>
PArray3D<T>::PArray3D(int dim1, int dim2, int dim3)
{
	size1 = dim1;
	size2 = dim2;
	size3 = dim3;
	PArrayBase<T>::ResetSize(size1*size2*size3);
}


template <class T>
PArray3D<T>::~PArray3D()
{
	;
}


// () OPERATOR
template < class T > inline
T& PArray3D<T>::operator()(int ind1, int ind2, int ind3)
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < size1);
		assert(ind2 >= 0 && ind2 < size2);
		assert(ind3 >= 0 && ind3 < size3);
		return PArrayBase<T>::array[IND1_IND2_IND3];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return PArrayBase<T>::array[IND1_IND2_IND3];
	#endif
}

template < class T > inline
const T& PArray3D<T>::operator()(int ind1, int ind2, int ind3) const
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < size1);
		assert(ind2 >= 0 && ind2 < size2);
		assert(ind3 >= 0 && ind3 < size3);
		return PArrayBase<T>::array[IND1_IND2_IND3];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return PArrayBase<T>::array[IND1_IND2_IND3];
	#endif
	}



// DATA ACCESS AND MODIFICATION FUNCTIONS
template <class T>
int PArray3D<T>::GetDim(int dim) const
{
	#ifndef RELEASE
		assert(dim > 0 && dim < 4);
		int retval = 0;
		if(dim == 1){retval = size1;}
		if(dim == 2){retval = size2;}
		if(dim == 3){retval = size3;}
		return retval;
	#else
		int retval = 0;
		if(dim == 1){retval = size1;}
		if(dim == 2){retval = size2;}
		if(dim == 3){retval = size3;}
		return retval;
	#endif
}


template <class T>
void PArray3D<T>::ResetVal(const T initval)
{
	for(int i=0; i<PArrayBase<T>::npoints; i++){
		PArrayBase<T>::array[i] = initval;
	}
}


template <class T>
void PArray3D<T>::ResetSize(int dim1, int dim2, int dim3)
{
	// CHECK THAT INPUT BOUNDS ARE INDEED DIFFERENT THAN CURRENT BOUNDS BEFORE
	// ATTEMPTING TO RESIZE THE ARRAY.
	if(dim1 != size1 || dim2 != size2 ||  dim3 != size3){
		size1 = dim1;
		size2 = dim2;
		size3 = dim3;

		PArrayBase<T>::ResetSize(size1*size2*size3);
	} else {
		// IF INPUT BOUNDS MATCH EXISTING BOUNDS, RESET ALL ARRAY POINTS TO
		// 'initvalue'
		for(int i=0; i<PArrayBase<T>::npoints; i++){
			PArrayBase<T>::array[i] = NULL;
		}
	}
}


template <class T>
double PArray3D<T>::GetMemoryUsage() const
{
	double retval = 0.0e0;

	// MEMORY REQUIREMENT FOR DATA AND npoints (MEMBERS FROM SOURCE CLASS)
	retval += PArrayBase<T>::GetMemoryUsage();

	// SPACE REQUIRED FOR PARAMETERS NOT INCLUDED IN SOURCE CLASS
	retval += (double)sizeof(size1);
	retval += (double)sizeof(size2);
	retval += (double)sizeof(size3);

	return retval;
}
