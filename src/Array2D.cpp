/**
 * @file Array2D.cpp
 * @author Robert Grandin
 * @brief Implementation of Array2D class.
 */

#include <Array2D.h>



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
Array2D<T>::Array2D()
{
	size1 = 1;
	size2 = 1;
	npoints = size1*size2;
}

template <class T>
Array2D<T>::Array2D(size_t dim1, size_t dim2, const T initvalue)
{
	size1 = dim1;
	size2 = dim2;
	npoints = size1*size2;
	ArrayBase<T>::ResetSize(npoints,initvalue);
}

template <class T>
Array2D<T>::~Array2D()
{
	;
}


/*
// ASSIGNMENT OPERATORS
template <class T>
Array2D<T>& Array2D<T>::operator=(const Array2D<T>& a)
{
  if(this != &a){
	int temp = a.rows*a.cols;
	if(temp <= npoints){
	  for(int i=0; i<tmp; i++){
		array[i] = a.array[i];
	  }
	} else {
	  delete [] array;
	  initialize(a.size1,a.size2,
	}
  }

  return *this;
}
*/

// () OPERATOR
template < class T > inline
T& Array2D<T>::operator()(size_t ind1, size_t ind2)
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
        assert(ind1 < size1);
        assert(ind2 < size2);
		return ArrayBase<T>::array[IND1_IND2];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return ArrayBase<T>::array[IND1_IND2];
	#endif
}

template < class T > inline
const T& Array2D<T>::operator()(size_t ind1, size_t ind2) const
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < size1);
		assert(ind2 >= 0 && ind2 < size2);
		return ArrayBase<T>::array[IND1_IND2];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return ArrayBase<T>::array[IND1_IND2];
	#endif
}


template < class T > inline
T& Array2D<T>::operator()(int ind1i, int ind2i)
{
    size_t ind1 = (size_t)ind1i;
    size_t ind2 = (size_t)ind2i;
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 < size1);
        assert(ind2 < size2);
        return ArrayBase<T>::array[IND1_IND2];
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return ArrayBase<T>::array[IND1_IND2];
    #endif
}

template < class T > inline
const T& Array2D<T>::operator()(int ind1i, int ind2i) const
{
    size_t ind1 = (size_t)ind1i;
    size_t ind2 = (size_t)ind2i;
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 >= 0 && ind1 < size1);
        assert(ind2 >= 0 && ind2 < size2);
        return ArrayBase<T>::array[IND1_IND2];
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return ArrayBase<T>::array[IND1_IND2];
    #endif
}



// DATA ACCESS AND MODIFICATION FUNCTIONS
template <class T>
size_t Array2D<T>::GetDim(int dim) const
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

template <class T> inline
T Array2D<T>::GetValue(size_t ind1, size_t ind2) const
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < size1);
		assert(ind2 >= 0 && ind2 < size2);
		return ArrayBase<T>::array[IND1_IND2];
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		return ArrayBase<T>::array[IND1_IND2];
	#endif
}

template <class T> inline
void Array2D<T>::SetValue(size_t ind1, size_t ind2, T value)
{
	#ifndef RELEASE
		/*
		 * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
		 */
		assert(ind1 >= 0 && ind1 < size1);
		assert(ind2 >= 0 && ind2 < size2);
		ArrayBase<T>::array[IND1_IND2] = value;
	#else
		/*
		 * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
		 */
		ArrayBase<T>::array[IND1_IND2] = value;
	#endif
}

template <class T>
void Array2D<T>::ResetVal(const T initval)
{
	for(int i=0; i<npoints; i++){
		ArrayBase<T>::array[i] = initval;
	}
}

template <class T>
void Array2D<T>::ResetSize(size_t dim1, size_t dim2)
{
	// CALL ResetSize WITHIN INITIAL VALUE SET TO "0".
	ResetSize(dim1,dim2,0);
}



template <class T>
void Array2D<T>::ResetSize(size_t dim1, size_t dim2, const T initvalue)
{
	// CHECK THAT INPUT BOUNDS ARE INDEED DIFFERENT THAN CURRENT BOUNDS BEFORE
	// ATTEMPTING TO RESIZE THE ARRAY.
	if(dim1 != size1 || dim2 != size2){
		size1 = dim1;
		size2 = dim2;
        npoints = size1*size2;

		ArrayBase<T>::ResetSize(npoints,initvalue);
	} else {
		// IF INPUT BOUNDS MATCH EXISTING BOUNDS, RESET ALL ARRAY POINTS TO
		// 'initvalue'
        for(size_t i=0; i<npoints; i++){
			ArrayBase<T>::array[i] = initvalue;
		}
	}
}


template <class T>
double Array2D<T>::GetMemoryUsage() const
{
	double retval = 0.0e0;

	// MEMORY REQUIREMENT FOR DATA AND npoints (MEMBERS FROM SOURCE CLASS)
	retval += ArrayBase<T>::GetMemoryUsage();

	// SPACE REQUIRED FOR PARAMETERS NOT INCLUDED IN SOURCE CLASS
	retval += (double)sizeof(size1);
	retval += (double)sizeof(size2);

	return retval;
}


template <class T>
T Array2D<T>::MinVal(size_t &loc1, size_t &loc2) const
{
    size_t loc = 0;
	T min = ArrayBase<T>::MinVal(loc);

	// CONVERT SERIALIZED LOCATION TO 2-INDEX LOCATION
    loc1 = (size_t)(loc/size2);
	loc2 = loc - loc1*size2;

	return min;
}


template <class T>
T Array2D<T>::MaxVal(size_t &loc1, size_t &loc2) const
{
    size_t loc = 0;
	T max = ArrayBase<T>::MaxVal(loc);

	// CONVERT SERIALIZED LOCATION TO 2-INDEX LOCATION
    loc1 = (size_t)(loc/size2);
	loc2 = loc - loc1*size2;

	return max;
}


template <class T>
void Array2D<T>::ReadCSVFile(const std::string filename, const int nheader, const int mincols,
		const T defval, const int linewidth=4096)
{
	std::fstream file;
	std::stringstream sstmp1;
	std::stringstream sstmp2;
	char vals1[linewidth];
	char vals2[linewidth];
	long filestartlocation;
	T qty_val = 0.0e0;

	// Open file and save starting location of the file.
    file.open(filename.c_str(),std::ios::in);
	filestartlocation = file.tellg();

	// Determine number of lines in the file and maximum number of columns.
	int ncols = 0;
	int tmpncols = 0;
	int linecount = 0;
	while(!file.eof()){
		tmpncols = 0;
		file.getline(vals1,linewidth);
		sstmp1 << vals1;
		while(sstmp1.getline(vals2,linewidth,',')){
			tmpncols++;
		}
		if(tmpncols > ncols){
			ncols = tmpncols;
		}
		linecount++;
	}
	file.clear();
	file.seekg(filestartlocation);

	// Subtract 1 due to extra row being added prior to encountering EOF in
	// above loop.  'linecount' is now the total number of lines in the file.
	linecount -= 1;

	// Subtract the user-specified number of header rows.
	linecount -= nheader;

	// If the number of columns in the file is less than the user-specified
	// minimum number of columns, use the user-specified value.  Otherwise, use
	// the number of columns in the text file.
	if(ncols < mincols){
		ncols = mincols;
	}

	// Resize array to hold values.  Array initialized to default value.
	Array2D<T>::ResetSize(linecount,ncols,defval);

	// Read header row and discard.  Then loop through remaining lines.
	//Array2D<T> tmparray(linecount,ncols,defvalue);
	int qty_count = 0;
	int ind1 = 0;
	int ind2 = 0;
	file.getline(vals1,linewidth);
	for(int i=0; i<linecount; i++){
		file.getline(vals1,linewidth);
		sstmp1.clear(); sstmp1.str("");
		sstmp1 << vals1;
		qty_count = 0;
		while(sstmp1.getline(vals2,linewidth,',')){
			sstmp2.clear(); sstmp2.str("");
			sstmp2 << vals2;
			sstmp2 >> qty_val; sstmp2.str("");

			// Set dummy variables to make use of macro defined above
			ind1 = i;
			ind2 = qty_count;
			ArrayBase<T>::array[IND1_IND2] = qty_val;

			qty_count++;
		}
	}
}


template <class T>
void Array2D<T>::WriteCSVFile(const std::string filename, const PArray1D<std::string*> &labels, const int ndec) const
{
    int ncols = GetDim(2);

    if(labels.GetDim() != ncols){
        std::cerr << "ERROR: Incorrect number of labels supplied" << std::endl;
        std::cerr << "       (Array2D<T>::WriteCSVFile)" << std::endl;
        assert(labels.GetDim() == ncols);
    }

    std::fstream file;
    file.open(filename.c_str(),std::ios::out);

    if(ndec > 0){
        file.setf ( std::ios::scientific, std::ios::floatfield );
        file.precision(ndec);
    }

    /* Write labels to file */
    for(int i=0; i<ncols-1; i++){
        file << labels(i)->substr() << ", ";
    }
    file << labels(ncols-1)->substr() << std::endl;

    /* Write array data */
    int ind1 = 0;
    int ind2 = 0;
    for(ind1=0; ind1<GetDim(1); ind1++){
        for(ind2=0; ind2<ncols-1; ind2++){
            file << ArrayBase<T>::array[IND1_IND2] << ", ";
        }
        ind2 = ncols-1;
        file << ArrayBase<T>::array[IND1_IND2] << std::endl;
    }

    file.close();

}
