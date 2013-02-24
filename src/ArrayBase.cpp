/**
 * @file ArrayBase.cpp
 * @author Robert Grandin
 * @brief Implementation of ArrayBase class.
 */


#include "ArrayBase.h"

// ==================================================================
// ================
// ================    PRIVATE FUNCTIONS
// ================

// ARRAY INITIALIZATION
template <class T>
void ArrayBase<T>::initialize(size_t dim1, T initvalue)
{
	npoints = dim1;
	if(!npoints){							// CHECK npoints VALUE SET
        array = (T*)NULL;					// SET POINTER TO NULL
	} else {								// ALLOCATE MEMORY FOR array
        array = new T[npoints];

        for(size_t i=0; i<npoints; i++){
            array[i] = initvalue;           // INITIALIZE TO SPECIFIED VALUE
		}

	}
}


template <class T>
int ArrayBase<T>::compare(const void *a, const void *b)
{
	return (int)( *(T*)a - *(T*)b );
}



// ==================================================================
// ================
// ================    PUBLIC FUNCTIONS
// ================

// CONSTRUCTORS AND DESTRUCTORS
template <class T>
ArrayBase<T>::ArrayBase()
{
	ArrayBase<T>::initialize(1,0);
}


template <class T>
ArrayBase<T>::ArrayBase(size_t dim1)
{
	initialize(dim1,(T)0.0e0);
}


template <class T>
ArrayBase<T>::ArrayBase(size_t dim1, const T initvalue)
{
	initialize(dim1,initvalue);
}


template <class T>
ArrayBase<T>::ArrayBase(const ArrayBase<T> &ab) : npoints(ab.npoints), array(npoints ? new T[npoints] : 0)
{
    std::copy(ab.array, ab.array + npoints, array);
}


#ifdef CXX11
template <class T>
ArrayBase<T>::ArrayBase(ArrayBase<T> &&ab) : ArrayBase()
{
    ArrayBaseSwap(*this, ab);
}
#endif


template <class T>
ArrayBase<T>::~ArrayBase()
{
    if(array != NULL){
        delete [] array;
    }
    array = NULL;
}


template <class T>
ArrayBase<T>& ArrayBase<T>::operator=(ArrayBase<T> ab)
{
    ArrayBaseSwap(*this, ab);
    return *this;
}


#ifdef CXX11
template <class T>
ArrayBase<T>& ArrayBase<T>::operator=(ArrayBase<T> &&ab)
{
    ArrayBaseSwap(*this, ab);
    return *this;
}
#endif


template <class T>
T& ArrayBase<T>::operator [](const size_t idx)
{
    return array[idx];
}


template <class T>
const T& ArrayBase<T>::operator [](const size_t idx) const
{
    return array[idx];
}


template <class T>
void ArrayBase<T>::ResetSize(size_t dim1)
{
	if(typeid(T) == typeid(std::string) ||
		typeid(T) == typeid(char)){
		// CHECK FOR NEW DIMENSIONS MATCHING EXISTING DIMENSIONS.  IF NOT, RESET
		// THE ARRAY SIZE AS REQUIRED.
		if(dim1 != npoints){
            delete [] array;
            array = new T[dim1];
			npoints = dim1;
		}
	} else {
		// CALL ResetSize WITH INITIAL VALUE SPECIFIED TO DEFAULT VALUE OF "0".
		ResetSize(dim1,(T)0.0e0);
	}
}

template <class T>
void ArrayBase<T>::ResetSize(size_t dim1, const T initvalue)
{
	//std::cout << "Resetting array size to " << dim1 << " points and initialized to " <<
	//		initvalue << std::endl;

	// REQUIRE THE NUMBER OF POINTS TO BE POSITIVE AND NON-ZERO
	assert(dim1 > 0);

	// CHECK FOR NEW DIMENSIONS MATCHING EXISTING DIMENSIONS.  IF NOT, RESET
	// THE ARRAY SIZE AS REQUIRED.
	if(dim1 != npoints){
		npoints = dim1;
        delete [] array;
        array = new T[npoints];

        for(size_t i=0; i<npoints; i++){
            array[i] = initvalue;
		}
	} else {
		// IF dim1 IS THE CURRENT ARRAY SIZE, SET VALUE AT ALL POINTS TO
		// initvalue.
        for(size_t i=0; i<npoints; i++){
            array[i] = initvalue;
		}
	}
}


template <class T>
void ArrayBase<T>::ResetVal(const T initval)
{
    assert(array);
    for(size_t i=0; i<npoints; i++){
        ArrayBase<T>::array[i] = initval;
	}
}


template <class T>
T ArrayBase<T>::Mean() const
{
    T retval = (T)0.0e0;

    T sum = (T)0.0e0;
    for(size_t i=0; i<npoints; i++){
        sum += array[i];
    }
    retval = (T)sum/(T)npoints;

    return retval;
}


template <class T>
float ArrayBase<T>::MeanFloat() const
{
	float retval = (T)0.0e0;
	float sum = 0.0e0;

	for(int i=0; i<npoints; i++){
        sum += (float)array[i];
	}

	retval = sum/(float)npoints;

	return retval;
}


template <class T>
T ArrayBase<T>::Variance() const
{
    T retval = (T)0.0e0;
    T mean = ArrayBase<T>::Mean();

    T diffsq = (T)0.0e0;
    for(size_t i=0; i<npoints; i++){
        diffsq += (array[i] - mean)*(array[i] - mean);
    }
    retval = diffsq/(T)npoints;

    return retval;
}


template <class T>
float ArrayBase<T>::VarianceFloat() const
{
	float retval = (T)0.0e0;
	float mean = ArrayBase<T>::Mean();
	float diffsq = 0.0e0;
	for(int i=0; i<npoints; i++){
        diffsq += ((float)array[i] - mean)*((float)array[i] - mean);
	}
	retval = diffsq/(float)npoints;

	return retval;
}


template <class T>
double ArrayBase<T>::GetMemoryUsage() const
{
	double retval = 0.0e0;

	// NUMBER OF VALUES STORED IN THE ARRAY
	retval += (double)npoints*(double)sizeof(T);

	// SPACE REQUIRED FOR ARRAY PARAMETERS
	retval += (double)sizeof(npoints);

	return retval;
}


template <class T>
T ArrayBase<T>::MinVal() const
{
    size_t loc = 0;
	return ArrayBase<T>::MinVal(loc);
}


template <class T>
T ArrayBase<T>::MinVal(size_t &loc) const
{
	// INITIALIZE TO MAXIMUM VALUE REPRESENTED BY VARIABLE TYPE, T
    T min = std::numeric_limits<T>::max();

	// LOOP THROUGH ARRAY AND CHECK IF ARRAY VALUE IS LESS THAN CURRENT MINIMUM
    for(size_t i=0; i<npoints; i++){
        if(array[i] < min){
            min = array[i];
			loc = i;
		}
	}

	return min;
}


template <class T>
T ArrayBase<T>::MaxVal() const
{
    size_t loc = 0;
	return ArrayBase<T>::MaxVal(loc);
}


template <class T>
T ArrayBase<T>::MaxVal(size_t &loc) const
{
	// INITIALIZE TO MAXIMUM VALUE REPRESENTED BY VARIABLE TYPE, T
    T max = -std::numeric_limits<T>::max();

	// LOOP THROUGH ARRAY AND CHECK IF ARRAY VALUE IS LESS THAN CURRENT MINIMUM
    for(size_t i=0; i<npoints; i++){
        if(array[i] > max){
            max = array[i];
			loc = i;
		}
	}

	return max;
}


template <class T>
T ArrayBase<T>::MedianVal() const
{
    // COPY DATA TO TEMPORARY ARRAY TO AVOID REORDERING DATA IN array
    T *ptmp = new T[npoints];
    for(size_t i=0; i<npoints; i++){
        ptmp[i] = array[i];
    }

    // CALL qsort.  ptmp WILL CONTAIN THE ASCENDINGLY-SORTED VALUES
    qsort(ptmp,npoints,sizeof(T),ArrayBase<T>::compare);

    /*
     * DETERMINE MEDIAN VALUE FROM SORTED VECTOR
     */
    T medval = (T)0.0e0;
    size_t offset = npoints/2;
    if(npoints%2 == 0){
        // MUST AVERAGE MIDDLE TWO VALUES
        medval = (T)(0.5e0*((T)ptmp[offset] + (T)ptmp[offset+1]));
    } else {
        // TAKE MIDDLE VALUE
        medval = (T)ptmp[offset];
    }

    delete [] ptmp;

    return medval;
}


template <class T>
size_t ArrayBase<T>::MemoryRequired() const
{
    return ArrayBase<T>::GetMemoryUsage();
}


template <class T>
T ArrayBase<T>::RMS() const
{
    T xrms = (T)0.0e0;

    for(size_t i=0; i<npoints; i++){
        xrms += array[i]*array[i];
    }

    return (T)sqrt((double)xrms/(double)npoints);
}


template <class T>
void ArrayBase<T>::Test(std::string &result)
{
    result = "SUCCESS";
    bool errorfound = false;

    /* Create initial array. */
    size_t npts = 15;
    ArrayBase<T> array1(npts);
    for(size_t i=0; i<npts; i++){
        array1[i] = (T)i;
    }


    /* Test copy constructor. */
    ArrayBase<T> array2(array1);
    bool compare = true;
    T eps = 1.0e-10;
    for(size_t i=0; i<npts; i++){
        if(fabs(array1[i] - array2[i]) > eps){
            compare = false;
        }
    }
    if(!compare){
        if(!errorfound){
            errorfound = true;
            result = "";
        }
        result += "  - Copy constructor failed \n";
    }


    /* Test copy assignment. */
    ArrayBase<T> array3;
    array3 = array1;
    compare = true;
    for(size_t i=0; i<npts; i++){
        if(fabs(array1[i] - array3[i]) > eps){
            compare = false;
        }
    }
    if(!compare){
        if(!errorfound){
            errorfound = true;
            result = "";
        }
        result += "  - Copy assignment failed \n";
    }



}
