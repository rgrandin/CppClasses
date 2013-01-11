/**
 * @file PArray4D.cpp
 * @author Robert Grandin
 * @brief Implementation of PArray4D class.
 */




// ==================================================================
// ================
// ================    PRIVATE FUNCTIONS
// ================

/** @brief Macro to map 4D indices to serialized index. */
#define IND1_IND2_IND3_IND4 ind4*size1*size2*size3+ind3*size1*size2+ind1*size2+ind2



// ==================================================================
// ================
// ================    PUBLIC FUNCTIONS
// ================

// CONSTRUCTORS AND DECONSTRUCTORS
template <class T>
PArray4D<T>::PArray4D()
{
	size1 = 1;
	size2 = 1;
	size3 = 1;
	size4 = 1;
	PArrayBase<T>::npoints = size1*size2*size3*size4;
}


template <class T>
PArray4D<T>::PArray4D(int dim1, int dim2, int dim3, int dim4)
{
	size1 = dim1;
	size2 = dim2;
	size3 = dim3;
	size4 = dim4;
	PArrayBase<T>::ResetSize(size1*size2*size3*size4);
}


template <class T>
PArray4D<T>::~PArray4D()
{
	;
}


// () OPERATOR
template < class T > inline
T& PArray4D<T>::operator()(int ind1, int ind2, int ind3, int ind4)
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < size1);
		assert(ind2 >= 0 && ind2 < size2);
		assert(ind3 >= 0 && ind3 < size3);
		assert(ind4 >= 0 && ind4 < size4);
		return PArrayBase<T>::array[IND1_IND2_IND3_IND4];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return PArrayBase<T>::array[IND1_IND2_IND3_IND4];
	#endif
}

template < class T > inline
const T& PArray4D<T>::operator()(int ind1, int ind2, int ind3,
		int ind4) const
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < size1);
		assert(ind2 >= 0 && ind2 < size2);
		assert(ind3 >= 0 && ind3 < size3);
		assert(ind4 >= 0 && ind4 < size4);
		return PArrayBase<T>::array[IND1_IND2_IND3_IND4];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return PArrayBase<T>::array[IND1_IND2_IND3_IND4];
	#endif
}



// DATA ACCESS AND MODIFICATION FUNCTIONS
template <class T>
int PArray4D<T>::GetDim(int dim) const

{
	#ifndef RELEASE
		assert(dim > 0 && dim < 5);
		int retval = 0;
		if(dim == 1){retval = size1;}
		if(dim == 2){retval = size2;}
		if(dim == 3){retval = size3;}
		if(dim == 4){retval = size4;}
		return retval;
	#else
		int retval = 0;
		if(dim == 1){retval = size1;}
		if(dim == 2){retval = size2;}
		if(dim == 3){retval = size3;}
		if(dim == 4){retval = size4;}
		return retval;
	#endif
}


template <class T>
void PArray4D<T>::ResetVal(const T initval)
{
	for(int i=0; i<PArrayBase<T>::npoints; i++){
		PArrayBase<T>::array[i] = initval;
	}
}


template <class T>
void PArray4D<T>::ResetSize(int dim1, int dim2, int dim3, int dim4)
{
	// CHECK THAT INPUT BOUNDS ARE INDEED DIFFERENT THAN CURRENT BOUNDS BEFORE
	// ATTEMPTING TO RESIZE THE ARRAY.
	if(dim1 != size1 || dim2 != size2 ||  dim3 != size3 || dim4 != size4){
		size1 = dim1;
		size2 = dim2;
		size3 = dim3;
		size4 = dim4;

		PArrayBase<T>::ResetSize(size1*size2*size3*size4);
	} else {
		// IF INPUT BOUNDS MATCH EXISTING BOUNDS, RESET ALL ARRAY POINTS TO
		// NULL
		for(int i=0; i<PArrayBase<T>::npoints; i++){
			PArrayBase<T>::array[i] = NULL;
		}
	}
}


template <class T>
double PArray4D<T>::GetMemoryUsage() const
{
	double retval = 0.0e0;

	// MEMORY REQUIREMENT FOR DATA AND npoints (MEMBERS FROM SOURCE CLASS)
	retval += PArrayBase<T>::GetMemoryUsage();

	// SPACE REQUIRED FOR PARAMETERS NOT INCLUDED IN SOURCE CLASS
	retval += (double)sizeof(size1);
	retval += (double)sizeof(size2);
	retval += (double)sizeof(size3);
	retval += (double)sizeof(size4);

	return retval;
}
