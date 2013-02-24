/**
 * @file PArrayBase.cpp
 * @author Robert Grandin
 * @brief Implementation of PArrayBase class.
 */


#include <PArrayBase.h>



// ==================================================================
// ================
// ================    PRIVATE FUNCTIONS
// ================

// ARRAY INITIALIZATION
template <class T>
void PArrayBase<T>::initialize(int dim1)
{
	npoints = dim1;
	if(!npoints){							// CHECK npoints VALUE SET
		array = (T*)NULL;					// SET POINTER TO NULL
	} else {								// ALLOCATE MEMORY FOR array
		array = new T[npoints];

        for(size_t i=0; i<npoints; i++){
			array[i] = NULL;				// INITIALIZE TO NULL
		}
	}
}



// ==================================================================
// ================
// ================    PUBLIC FUNCTIONS
// ================

// CONSTRUCTORS AND DESTRUCTORS
template <class T>
PArrayBase<T>::PArrayBase()
{
	PArrayBase<T>::initialize(1);
}

template <class T>
PArrayBase<T>::PArrayBase(size_t dim1)
{
	initialize(dim1);
}


template <class T>
PArrayBase<T>::PArrayBase(const PArrayBase<T> &a) : npoints(a.npoints),
    array(npoints ? new T[npoints] : 0)
{
    std::copy(a.array, a.array + npoints, array);
}


#ifdef CXX11
template <class T>
PArrayBase<T>::PArrayBase(PArrayBase<T> &&a) : PArrayBase<T>()
{
    PArrayBaseSwap(*this, a);
}
#endif


template <class T>
PArrayBase<T>::~PArrayBase()
{
    for(size_t i=0; i<npoints; i++){
		delete array[i];
	}
	delete [] array;
}


template <class T>
void PArrayBase<T>::ResetSize(size_t dim1)
{
	// REQUIRE THE NUMBER OF POINTS TO BE POSITIVE AND NON-ZERO
	assert(dim1 > 0);

	// CHECK FOR NEW DIMENSIONS MATCHING EXISTING DIMENSIONS.  IF NOT, RESET
	// THE ARRAY SIZE AS REQUIRED.
	if(dim1 != npoints){
        for(size_t i=0; i<npoints; i++){
			if(array[i] != NULL){
				delete array[i];
			}
		}

		npoints = dim1;
		delete [] array;
		array = new T[npoints];

        for(size_t i=0; i<npoints; i++){
			array[i] = NULL;
		}
	} else {
		// IF dim1 IS THE CURRENT ARRAY SIZE, SET VALUE AT ALL POINTS TO
		// NULL.
        for(size_t i=0; i<npoints; i++){
			array[i] = NULL;
		}
	}
}


template <class T>
void PArrayBase<T>::ResetVal(const T initvalue)
{
	for(int i=0; i<npoints; i++){
		array[i] = initvalue;
	}
}


template <class T>
double PArrayBase<T>::GetMemoryUsage() const
{
	double retval = 0.0e0;

	// NUMBER OF VALUES STORED IN THE ARRAY
	retval += (double)npoints*(double)sizeof(T);

	// SPACE REQUIRED FOR ARRAY PARAMETERS
	retval += (double)sizeof(npoints);

	return retval;
}


template <class T>
PArrayBase<T>& PArrayBase<T>::operator=(PArrayBase<T> a)
{
    PArrayBaseSwap(*this, a);
    return *this;
}


template <class T>
T& PArrayBase<T>::operator [](const size_t idx)
{
    return array[idx];
}


template <class T>
const T& PArrayBase<T>::operator [](const size_t idx) const
{
    return array[idx];
}
