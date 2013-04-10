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
    use_free = false;   /* Use 'new' and 'delete' by default. */

    array_name = "NoName";

	npoints = dim1;
	if(!npoints){							// CHECK npoints VALUE SET
        p_array = (T*)NULL;					// SET POINTER TO NULL
	} else {								// ALLOCATE MEMORY FOR array
        p_array = new T[npoints];

        for(size_t i=0; i<npoints; i++){
            p_array[i] = initvalue;           // INITIALIZE TO SPECIFIED VALUE
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
ArrayBase<T>::ArrayBase(const ArrayBase<T> &ab) : npoints(ab.npoints),
    p_array(npoints ? new T[npoints] : 0)
{
    //std::copy(ab.p_array, ab.p_array + npoints, p_array);
    for(size_t i=0; i<npoints; i++){
        p_array[i] = ab.p_array[i];
    }
}


#ifdef CXX11
template <class T>
ArrayBase<T>::ArrayBase(ArrayBase<T> &&ab)
{
    ArrayBaseSwap(*this, ab);
}
#endif


template <class T>
ArrayBase<T>::~ArrayBase()
{
    FreeMemory();
    p_array = NULL;
}


template <class T>
void ArrayBase<T>::FreeMemory()
{
    if(p_array != NULL){
        if(use_free){
            free(p_array);
        } else {
            delete [] p_array;
        }

        p_array = NULL;
    }
}


template <class T>
ArrayBase<T>& ArrayBase<T>::operator=(ArrayBase<T> ab)
{
    ArrayBaseSwap(*this, ab);
    return *this;
}


template <class T>
T& ArrayBase<T>::operator [](const size_t idx)
{
    return p_array[idx];
}


template <class T>
const T& ArrayBase<T>::operator [](const size_t idx) const
{
    return p_array[idx];
}


template <class T>
void ArrayBase<T>::ResetSize(size_t dim1)
{
	if(typeid(T) == typeid(std::string) ||
		typeid(T) == typeid(char)){
		// CHECK FOR NEW DIMENSIONS MATCHING EXISTING DIMENSIONS.  IF NOT, RESET
		// THE ARRAY SIZE AS REQUIRED.
		if(dim1 != npoints){
            FreeMemory();
            p_array = new T[dim1];
			npoints = dim1;
            use_free = false;
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
//	if(dim1 != npoints){
		npoints = dim1;
        FreeMemory();
        p_array = new T[npoints];
        use_free = false;

        for(size_t i=0; i<npoints; i++){
            p_array[i] = initvalue;
		}
//	} else {
//		// IF dim1 IS THE CURRENT ARRAY SIZE, SET VALUE AT ALL POINTS TO
//		// initvalue.
//        for(size_t i=0; i<npoints; i++){
//            p_array[i] = initvalue;
//		}
//	}
}


template <class T>
void ArrayBase<T>::ResetVal(const T initval)
{
    assert(p_array);
    for(size_t i=0; i<npoints; i++){
        p_array[i] = initval;
	}
}


template <class T>
T ArrayBase<T>::Mean() const
{
    T retval = (T)0.0e0;

    T sum = (T)0.0e0;
    for(size_t i=0; i<npoints; i++){
        sum += p_array[i];
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
        sum += (float)p_array[i];
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
        diffsq += (p_array[i] - mean)*(p_array[i] - mean);
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
        diffsq += ((float)p_array[i] - mean)*((float)p_array[i] - mean);
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
        if(p_array[i] < min){
            min = p_array[i];
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
    T max = (T)0.0e0;
    if(std::numeric_limits<T>::is_signed){
        max = -std::numeric_limits<T>::max();
    }

	// LOOP THROUGH ARRAY AND CHECK IF ARRAY VALUE IS LESS THAN CURRENT MINIMUM
    for(size_t i=0; i<npoints; i++){
        if(p_array[i] > max){
            max = p_array[i];
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
        ptmp[i] = p_array[i];
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
        xrms += p_array[i]*p_array[i];
    }

    return (T)sqrt((double)xrms/(double)npoints);
}


template <class T>
size_t ArrayBase<T>::NPts()
{
    return npoints;
}


template <class T>
void ArrayBase<T>::SetArrayPointer(T *p_data, size_t npts, bool useFree)
{
    /* Reset size to delete any existing data and free memory. */
    FreeMemory();

    p_array = p_data;
    npoints = npts;
    use_free = useFree;
}


template <class T>
void ArrayBase<T>::SetArrayPointerNULL()
{
    npoints = 0;
    p_array = NULL;
    use_free = false;
}


template <class T>
void ArrayBase<T>::setName(const std::string name)
{
    array_name = name;
}


template <class T>
std::string ArrayBase<T>::Name() const
{
    return array_name;
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
        array1[i] = (T)i + (T)1.0e0;
    }


    /* Test copy constructor. */
    ArrayBase<T> array2(array1);
    bool compare = true;
    T eps = (T)1.0e-10;
    for(size_t i=0; i<npts; i++){
        if(fabs((float)(array1[i] - array2[i])) > eps){
            compare = false;
        }
    }
    if(!compare){
        if(!errorfound){
            errorfound = true;
            result = "FAILED \n";
        }
        result += "  - Copy constructor failed \n";
    }


    /* Test copy assignment. */
    ArrayBase<T> array3;
    array3 = array1;
    compare = true;
    for(size_t i=0; i<npts; i++){
        if(fabs((float)(array1[i] - array3[i])) > eps){
            compare = false;
        }
    }
    if(!compare){
        if(!errorfound){
            errorfound = true;
            result = "FAILED \n";
        }
        result += "  - Copy assignment failed \n";
    }


    /* Test array statistics calculations. */
    T truesum = (T)npts*((T)npts + (T)1.0)*(T)0.5e0;
    T truemean = truesum/(T)array1.NPts();
    T testmean = array1.Mean();
    if(fabs((float)(truemean - testmean)) > eps){
        if(!errorfound){
            errorfound = true;
            result = "FAILED \n";
        }
        result += "  - Mean calculation failed \n";
    }

    T truemedian = (T)0.0e0;
    if(npts % 2 == 1){
        truemedian = array1[npts/2];
    } else {
        truemedian = (T)0.5e0*(array1[(npts-1)/2] + array1[npts/2]);
    }
    T testmedian = array1.MedianVal();
    if(fabs((float)(truemedian - testmedian)) > eps){
        if(!errorfound){
            errorfound = true;
            result = "FAILED \n";
        }
        result += "  - Median calculation failed \n";
    }

    T truemin = (T)1.0e0;
    T truemax = (T)npts;
    T testmin = array1.MinVal();
    T testmax = array1.MaxVal();
    if(fabs((float)(truemin - testmin)) > eps){
        if(!errorfound){
            errorfound = true;
            result = "FAILED \n";
        }
        result += "  - Min-value calculation failed \n";
    }
    if(fabs((float)(truemax - testmax)) > eps){
        if(!errorfound){
            errorfound = true;
            result = "FAILED \n";
        }
        result += "  - Max-value calculation failed \n";
    }
}
