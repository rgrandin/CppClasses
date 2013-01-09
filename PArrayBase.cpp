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
PArrayBase<T>::PArrayBase(long long dim1)
{
	initialize(dim1);
}


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
