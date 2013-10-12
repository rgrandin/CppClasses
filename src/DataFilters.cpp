/**
 * @file DataFilters.cpp
 * @author Robert Grandin
 * @brief Implementation of DataFilters class.
 */



#include "DataFilters.h"




// ==================================================================
// ================
// ================    PRIVATE FUNCTIONS
// ================
template <class T>
int DataFilters<T>::compare(const void *a, const void *b)
{
	return (int)( *(T*)a - *(T*)b );
}


template <class T>
void DataFilters<T>::Sinc1D_GetFFT(size_t nfft, T rolloff, T expatten)
{
	/*
	 * CALCULATE THE FFT OF THE SINC FILTER
	 */

    T normfactor = expatten*0.5e0*exp(1.0e0);

	if(n_sinc_fft != nfft){

		delete [] sinc_fft;

		sinc_fft = new fftw_complex[nfft];
		double *tmp = new double[nfft];
		n_sinc_fft = nfft;

		T wint = 2.0*PI/(T)nfft;
		T w, rn1, rn2, rd;

		// CALCULATE FILTER BASED ON THESE FREQUENCIES
        for(size_t i=0; i<nfft; i++){
			w = (T)i*wint - PI;
            if(w > PI){ w -= (T)2.0*PI; }
			rn1 = abs(2.0/rolloff*sin(rolloff*w/2.0));
			rn2 = sin(rolloff*w/2.0);
			rd = rolloff*w/2.0;
            tmp[i] = 2.0e0*PI*rn1*(rn2/rd)*(rn2/rd)*normfactor;
		}

		// SHIFT FILTER TO RANGE FROM 0 TO +2*PI
        size_t j = nfft/2;
        for(size_t i=0; i<nfft; i++){
			sinc_fft[j][0] = tmp[i];
			sinc_fft[j][1] = 0.0;
			j++;
			if(j == nfft){ j = 0; }
		}

		delete [] tmp;

		return;
	} else {
		// NO ACTION REQUIRED
	}
}


template <class T>
void DataFilters<T>::Sinc2D_SetPlan(int plan, fftw_complex *input, fftw_complex *output)
{
    /** @brief Define custom filter flags for FFTW. */
    #define PLANFLAGS FFTW_MEASURE|FFTW_UNALIGNED
    //#define PLANFLAGS FFTW_ESTIMATE|FFTW_UNALIGNED
	//#define PLANFLAGS FFTW_PATIENT|FFTW_UNALIGNED
	//#define PLANFLAGS FFTW_EXHAUSTIVE|FFTW_UNALIGNED

	// CHECK plan VALUE
	#ifndef RELEASE
		assert(plan == 1 || plan == -1);
    #endif

	// DEFINE FORWARD PLAN
	if(plan == 1 && sinc2d_forward_plan_set == false){
        sinc2d_forward_plan = fftw_plan_dft_2d((int)n_sinc_fft1,(int)n_sinc_fft2,input,output,
				FFTW_FORWARD, PLANFLAGS);
		sinc2d_forward_plan_set = true;
	}

	// DEFINE BACKWARD PLAN
	if(plan == -1 && sinc2d_backward_plan_set == false){
        sinc2d_backward_plan = fftw_plan_dft_2d((int)n_sinc_fft1,(int)n_sinc_fft2,input,output,
				FFTW_BACKWARD, PLANFLAGS);
		sinc2d_backward_plan_set = true;
	}

}


template <class T>
T DataFilters<T>::CustomFilter2D_GetAmplitudeReal(const T w1, const T w2) const
{
	// Require that filter-definition file has been read
	#ifndef RELEASE
		assert(custom2d_fileread == true);
	#endif

	//T eps = 1.0e-15;		// Tolerance for determining DC offset.
	fftw_complex retval;
	retval[0] = 1.0e0;		// Default return value does not affect amplitude
	retval[1] = 0.0e0;		// Default return value does not affect phase

	/*
	 * Check that supplied frequencies are not DC offset (i.e., w1 = 0 & w2 = 0).
	 * Both frequencies must be within (-eps,eps) in order for them to meet the
	 * requirement of being DC offset.  Such a window must be used since exact
	 * comparison of floating point numbers is not allowed.  This window should
	 * be small.
	 *
	 * If frequencies are DC offset, the filter returns {1,0} so that the offset
	 * is not affected either in amplitude or phase.
	 *
	 * Note: this w1,w2 ~= 0 conditional is commented-out based on testing which
	 * showed that its inclusion (even with eps = 1.0e-15) introduced artifacts
	 * and its exclusion produced identical results to the hard-coded CT filter
	 * when a ramp filter was defined via this custom mechanism.
	 *
	 * If frequencies are not DC offset, the appropriate filter value is determined
	 * and returned.
	 *
	 * If the first user-defined frequency in either direction has a magnitude
	 * greater than 2*PI, a radially-symmetric filter is used with
	 * w = sqrt(w1*w1 + w2*w2).
	 *
	 * If the first user-defined frequency in both direction has a magnitude
	 * less than 2*PI, the filter is defined by finding the interpolated value
	 * in each direction and then using the product of the two values.
	 */
	//if((abs(w1) > eps) && (abs(w2) > eps)){
		if(abs(custom2d_points(0,0)) > 2.0e0*PI || abs(custom2d_points(0,2)) > 2.0e0*PI){
			retval[0] = CustomFilter2D_GetAmplitudeRadialReal(w1,w2);
		} else {

			T w1abs = abs(w1);
			T w2abs = abs(w2);

			T a1 = CustomFilter_InterpolateMagnitude(0,w1abs);
			T a2 = CustomFilter_InterpolateMagnitude(0,w2abs);

			retval[0] = a1*a2;

		} // End of conditional that first user-defined frequency < 2*PI

	//} else {
	//	retval[0] = 1.0e0;
	//} // End of conditional that w1,w2 =/= 0

	return retval[0];
}


template <class T>
T DataFilters<T>::CustomFilter2D_GetAmplitudeRadialReal(const T w1, const T w2) const
{
	int interpdir = -1;
	if(abs(custom2d_points(0,0)) > 2.0e0*PI){
		interpdir = 2;
	} else {
		interpdir = 0;
	}

	T w = sqrt(w1*w1 + w2*w2);

	return CustomFilter_InterpolateMagnitude(interpdir,w);
}


template <class T>
T DataFilters<T>::CustomFilter_InterpolateMagnitude(const int dir, const T w,
		const int method=0) const
{
	T retval = 0.0e0;

	/*
	 * Piecewise-linear interpolation
	 */
	if(method == 0){
		int limitindex1_lower = 0;
		int limitindex1_upper = 1;
		int nterms = custom2d_points.GetDim(1);

		// Loop through terms read-in from filter-definition file to find bounding
		// points on either side of |w1| and |w2|
		for(int i=0; i<nterms-1; i++){
			if(w >= custom2d_points(limitindex1_upper,dir)){
				// Increment indices and check again
				limitindex1_lower++;
				limitindex1_upper++;
			}
		}

		// Apply piecewise-linear model using the bounding points.
		fftw_complex c1;	// Interpolated coefficient

		// Determine interpolated value.  If the desired frequency is a quarter of
		// the way between the user-specified points, a quarter of the magnitude
		// change in that interval is added to the magnitude at the start of the
		// interval.

		// Only interpolate if desired frequency is bounded on both sides by
		// user-specified values.  If the desired frequency is greater than the
		// maximum user-specified value, the magnitude corresponding to this
		// maximum-defined frequency is used.
		if(limitindex1_upper < custom2d_points.GetDim(1)){
			// Determine change in magnitude across interval
			T dc1 = custom2d_points(limitindex1_upper,dir+1) - custom2d_points(limitindex1_lower,dir+1);

			// Determine location within user-specified frequency range and normalize
			// to [0.0, 1.0]
			T dw1 = custom2d_points(limitindex1_upper,dir) - custom2d_points(limitindex1_lower,dir);
			T ndw1 = (w - custom2d_points(limitindex1_lower,0))/dw1;

			// Perform the actual interpolation
			c1[0] = custom2d_points(limitindex1_lower,dir+1) + dc1*ndw1;
			c1[1] = 0.0e0;

			retval = c1[0];

			//printf("Frequency: %f\n",w);
			//printf("  Lower-bound: %f -- %f \n",custom2d_points(limitindex1_lower,dir),
			//		custom2d_points(limitindex1_lower,dir+1));
			//printf("  Upper-bound: %f -- %f \n",custom2d_points(limitindex1_upper,dir),
			//		custom2d_points(limitindex1_upper,dir+1));
			//printf("  Interpolated value: %f \n",retval);

		} else {
			retval = custom2d_points(custom2d_points.GetDim(1)-1,dir+1);
		}

	}

	return retval;
}

// ==================================================================
// ================
// ================    PUBLIC FUNCTIONS
// ================

template <class T>
DataFilters<T>::DataFilters()
{
	// DEFINE INITIAL VALUES
	sinc_fft = new fftw_complex[1];
	n_sinc_fft = 1;
	n_sinc_fft1 = 1;
	n_sinc_fft2 = 1;
	sinc2d_forward_plan_set = false;
	sinc2d_backward_plan_set = false;
	custom2d_forward_plan_set = false;
	custom2d_backward_plan_set = false;
	custom2d_fileread = false;
	custom2d_nfft1 = 1;
	custom2d_nfft2 = 1;
	custom2d_points.ResetSize(1,1);
	custom2d_fft = new fftw_complex[1];

	// CREATE FFTW PLANS TO ALLOW FOR VALID DESTROY ON OBJECT DESTRUCTION.  THIS IS
	// TO AVOID ERRORS WHEN fftw_destroy_plan() IS CALLED.
	fftw_complex *tmp = new fftw_complex[4];
	fftw_complex *dtmp = new fftw_complex[4];
	sinc2d_forward_plan = fftw_plan_dft_2d(2,2,dtmp,tmp,FFTW_FORWARD, FFTW_ESTIMATE);
	sinc2d_backward_plan = fftw_plan_dft_2d(2,2,tmp,dtmp,FFTW_BACKWARD, FFTW_ESTIMATE);
	custom2d_forward_plan = fftw_plan_dft_2d(2,2,dtmp,tmp,FFTW_FORWARD, FFTW_ESTIMATE);
	custom2d_backward_plan = fftw_plan_dft_2d(2,2,tmp,dtmp,FFTW_BACKWARD, FFTW_ESTIMATE);
	delete [] tmp;
	delete [] dtmp;
}


template <class T>
DataFilters<T>::DataFilters(const DataFilters<T> &a) :
    n_sinc_fft(a.n_sinc_fft), n_sinc_fft1(a.n_sinc_fft1), n_sinc_fft2(a.n_sinc_fft2),
    sinc2d_backward_plan(a.sinc2d_backward_plan), sinc2d_backward_plan_set(a.sinc2d_backward_plan_set),
    sinc2d_forward_plan(a.sinc2d_forward_plan), sinc2d_forward_plan_set(a.sinc2d_forward_plan_set),
    custom2d_backward_plan(a.custom2d_backward_plan),
    custom2d_backward_plan_set(a.custom2d_backward_plan_set),
    custom2d_fileread(a.custom2d_fileread),
    custom2d_forward_plan(a.custom2d_forward_plan),
    custom2d_forward_plan_set(a.custom2d_forward_plan_set),
    custom2d_nfft1(a.custom2d_nfft1), custom2d_nfft2(a.custom2d_nfft2),
    custom2d_points(a.custom2d_points),
    sinc_fft(n_sinc_fft ? new fftw_complex[n_sinc_fft] : 0),
    custom2d_fft(custom2d_nfft1*custom2d_nfft2 ? new fftw_complex[custom2d_nfft1*custom2d_nfft2] : 0)
{
    //std::copy(a.custom2d_fft, a.custom2d_fft + custom2d_nfft1*custom2d_nfft2, custom2d_fft);
    //std::copy(a.sinc_fft, a.sinc_fft + n_sinc_fft, sinc_fft);
    size_t npts = custom2d_nfft1*custom2d_nfft2;
    if(custom2d_fft){
        for(size_t i=0; i<npts; i++){
            custom2d_fft[i][0] = a.custom2d_fft[i][0];
            custom2d_fft[i][1] = a.custom2d_fft[i][1];
        }
    }

    if(sinc_fft){
        for(size_t i=0; i<n_sinc_fft; i++){
            sinc_fft[i][0] = a.sinc_fft[i][0];
            sinc_fft[i][1] = a.sinc_fft[i][1];
        }
    }
}


#ifdef CXX11
template <class T>
DataFilters<T>::DataFilters(DataFilters<T> &&a) : DataFilters<T>(std::move(a))
{
    DataFiltersSwap(*this, a);
}
#endif


template <class T>
DataFilters<T>::~DataFilters()
{
	// FREE MEMORY
    if(sinc_fft){
        delete [] sinc_fft;
    }
    if(custom2d_fft){
        delete [] custom2d_fft;
    }
    if(sinc2d_forward_plan_set){
        fftw_destroy_plan(sinc2d_forward_plan);
    }
    if(sinc2d_backward_plan_set){
        fftw_destroy_plan(sinc2d_backward_plan);
    }
    if(custom2d_forward_plan_set){
        fftw_destroy_plan(custom2d_forward_plan);
    }
    if(custom2d_backward_plan_set){
        fftw_destroy_plan(custom2d_backward_plan);
    }
}


template <class T>
DataFilters<T>& DataFilters<T>::operator=(DataFilters<T> a)
{
    DataFiltersSwap(*this, a);
    return *this;
}



template <class T>
void DataFilters<T>::Median1D(Array1D<T> &vec, int kernalsize,
		int endcond, bool overwrite)
{
	// GET SIZE OF INPUT ARRAY
    int n = (int)vec.GetDim();

	// CREATE ARRAY TO HOLD FILTERED VALUES
	Array1D<T> filtered(2,0.0);
	if(overwrite == false){ filtered.ResetSize(n);}

	/*
	 * CHECK IF KERNAL SIZE IS ODD OR EVEN AND DETERMINE OFFSET SIZE ACCORDINGLY
	 */
	int offset;
	bool isodd;
	if(kernalsize%2 == 0){  // EVEN
		offset = kernalsize/2;
		isodd = false;
	} else { // ODD
		offset = (kernalsize - 1)/2;
		isodd = true;
	}

	Array1D<T> tmp(kernalsize,0.0);
    for(int i=0; i<n; i++){
        int j = i-offset;
        for(int k=0; k<kernalsize; k++){
            int l = j + k;
			/*
			 * THIS SECTION OF CODE USES A PERIODIC BOUNDARY CONDITION ON THE
			 * TEMPORARY 1D ARRAY.
			 */
			if(endcond == 1){
				if(l < 0){tmp(k) = vec(n+l);}
				if(l > n-1){tmp(k) = vec(l-n);}
			}

			/*
			 * THIS SECTION OF CODE EXTENDS THE END POINTS AS NECESSARY WHEN
			 * GENERATING THE TEMPORARY 1D ARRAY
			 */
			if(endcond == 2){
				if(l < 0){tmp(k) = vec(0);}
				if(l > n-1){tmp(k) = vec(n-1);}
			}

			/*
			 * THIS LINE IS UNAFFECTED BY THE CHOICE OF BOUNDARY CONDITIONS
			 */
			if(l >= 0 && l <= n-1){tmp(k) = vec(l);}
		}

		// MUST USE POINTER METHOD FOR qsort
		T *ptmp = new T[kernalsize];
		for(int j=0; j<kernalsize; j++){ptmp[j] = tmp(j);}

		// CALL qsort.  ptmp WILL CONTAIN THE ASCENDINGLY-SORTED VALUES
		qsort(ptmp,kernalsize,sizeof(T),DataFilters<T>::compare);

		/*
		 * DETERMINE MEDIAN VALUE FROM SORTED VECTOR
		 */
		T medval;
		if(isodd == false){	// MUST AVERAGE MIDDLE TWO VALUES
			medval = 0.5*((T)ptmp[offset] + (T)ptmp[offset+1]);
		} else {	// TAKE MIDDLE VALUE
			medval = (T)ptmp[offset];
		}

		// DELETE THE POINTER-BASED ARRAY USED FOR qsort
		delete [] ptmp;

		if(overwrite == true){
			// PLACE MEDIAN VALUE INTO INPUT VECTOR
			vec(i) = medval;
		} else {
			// PLACE MEDIAN VALUE INTO STORAGE VECTOR.  INPUT VECTOR IS
			// *NOT* OVERWRITTEN
			filtered(i) = medval;

		}

	} // END OF PASS THROUGH INPUT VECTOR

	/*
	 * IF OVERWRITING WAS NOT ALLOWED, TRANSFER THE STORED VALUES FROM filtered
	 * ARRAY TO THE INPUT ARRAY.  IF OVERWRITING WAS ALLOWED, NO FURTHER ACTION
	 * NEEDS TO BE TAKEN.
	 */
	if(overwrite == false){
        for(int i=0; i<n; i++){
			vec(i) = (T)filtered(i);
		}
	}
}



template <class T>
void DataFilters<T>::Sinc1D_Real(Array1D<T> &data, T rolloff, T expatten, size_t nfft)
{
	/*
	 * CONVOLVE REAL 1D DATA WITH A SINC FILTER.
	 *
	 * THIS FUNCTION BASED ON SAMPLE CODE FOUND AT
	 * http://people.sc.fsu.edu/~jburkardt/c_src/fftw3/fftw3.html
	 * ON 2 APRIL 2011.
	 */

    size_t n = data.GetDim();	// NUMBER OF DATA POINTS
	double *data_sp;
	double *data_filtered_sp;
	fftw_complex *data_fft;
	fftw_complex *data_filtered_fft;
	fftw_plan data_plan_forward;
	fftw_plan data_plan_backward;



	/*
	 * CREATE ARRAY TO HOLD INPUT DATA AND COPY VALUES INTO IT.
	 */
	data_sp = new double[nfft];

    for (size_t i=0; i<nfft; i++)
	{
		if(i < n) {
			data_sp[i] = (double)data(i);
		} else {
			data_sp[i] = 0.0;
		}
    }

	/*
	 * CREATE ARRAY TO HOLD THE HAMMING FILTER IN THE SPATIAL DOMAIN AND
	 * DEFINE ITS VALUES.  ALSO, CALCULATE ITS FFT.  IF A HAMMING FILTER WITH
	 * THE DESIRED VALUES (n,m, & nfft) HAS ALREADY BEEN CALCULATED, IT IS NOT
	 * RE-CALCULATED.  IF ANY OF THE DESIRED VALUES HAVE CHANGED, THE FILTER
	 * IS RE-CALCULATED ACCORDINGLY AND THE PREVIOUS FILTER IS LOST.
	 */
    DataFilters<T>::Sinc1D_GetFFT(nfft,rolloff,expatten);

	/*
	 * SETUP AN ARRAY TO HOLD THE TRANSFORMED DATA, CREATE THE FFTW "PLAN", AND
	 * EXECUTE THE PLAN TO OBTAIN THE FFT TRANSFORM OF THE DATA.
	 *
	 * ALSO, EXECUTE THE SAME PLAN ON THE HAMMING FILTER TO GET ITS
	 * REPRESENTATION IN THE FREQUENCY-DOMAIN.
	 */

	data_fft = new fftw_complex[nfft];
	data_filtered_fft = new fftw_complex[nfft];
    int infft = (int)n_sinc_fft;
    data_plan_forward = fftw_plan_dft_r2c_1d(infft,data_sp,data_fft,
			FFTW_ESTIMATE);
	fftw_execute (data_plan_forward);

	double scale = 1.0/(double)n_sinc_fft;

    for(size_t i=0; i<nfft; i++){
		data_filtered_fft[i][0] = (data_fft[i][0]*sinc_fft[i][0] -
				data_fft[i][1]*sinc_fft[i][1])*scale;
		data_filtered_fft[i][1] = (data_fft[i][0]*sinc_fft[i][1] +
				data_fft[i][1]*sinc_fft[i][0])*scale;
	}



	/*
	 * CREATE ARRAY TO HOLD THE INVERSE-TRANSFORMED FILTERED DATA AND PERFORM
	 * THE INVERSION.
	 */
	data_filtered_sp = new double[nfft];

    data_plan_backward = fftw_plan_dft_c2r_1d(infft,data_filtered_fft,
			data_filtered_sp,FFTW_ESTIMATE );

	fftw_execute(data_plan_backward);

	// REPLACE INPUT DATA WITH FILTERED RESULT
    for(size_t i=0; i<n; i++){
		data(i) = (T)data_filtered_sp[i];
	}


	/*
	 * FREE MEMORY ALLOCATED BY FFTW ROUTINES
	 */
	fftw_destroy_plan(data_plan_forward);
	fftw_destroy_plan(data_plan_backward);

	delete [] data_sp;
	delete [] data_fft;
	delete [] data_filtered_sp;
	delete [] data_filtered_fft;

	return;
}


template <class T>
void DataFilters<T>::Sinc2D_GetFFT(const size_t nfft1, const size_t nfft2,
                                   const T expatten1, const T expatten2)
{
    bool debugoutputtodisk = false;

	/*
	 * CALCULATE THE FFT OF THE SINC FILTER
	 */

	if(n_sinc_fft1 != nfft1 || n_sinc_fft2 != nfft2){

        /* Since FFT size has changed, force recalculation of FFTW plans. */
		sinc2d_forward_plan_set = false;
		sinc2d_backward_plan_set = false;

		n_sinc_fft1 = nfft1;
		n_sinc_fft2 = nfft2;
		n_sinc_fft = nfft1*nfft2;

		delete [] sinc_fft;
		sinc_fft = new fftw_complex[n_sinc_fft];

		T wint1 = 2.0e0*PI/(T)nfft1;
		T wint2 = 2.0e0*PI/(T)nfft2;
        T w1, w2;
//        T expof1 = exp(1.0e0);

        /* Calculate filter. */
//        double scale = 1.0e0/((double)n_sinc_fft);
        for(size_t i=0; i<nfft2; i++){
            w2 = (T)i*wint2;
            if(w2 > PI){ w2 -= 2.0e0*PI; }

            for(size_t j=0; j<nfft1; j++){
                w1 = (T)j*wint1;
                if(w1 > PI){ w1 -= 2.0e0*PI; }


				if(i == 0 && j == 0){
                    sinc_fft[0][0] = scale;
                    sinc_fft[0][1] = 0.0e0;
				} else {
                    /* Calculate frequency magnitude. */
                    T r = sqrt(w1*w1 + w2*w2);

                    /* Weighted sum of exponential attenuation values. */
                    T expval = w1*w1/r/r*expatten1 + w2*w2/r/r*expatten2;

//                    /* Calculate normalization factors so maximum filter amplitude is 1.0. */
//                    T normfactor = expval*0.5e0*expof1;     /* Only valid if 1/expval > PI */

//                    if(expval < 1.0e0/PI){
//                        normfactor = 0.5e0*expf(PI*r)/PI;   /* For case where 1/expval < PI */
//                    }

//                    normfactor = (T)1.0e0;      /* Reset to disable scaling effect. */


                    /* Set filter magnitude. */
                    sinc_fft[i*nfft1+j][0] = 2.0e0*r*exp(-r*expval);//*scale*normfactor;
                    sinc_fft[i*nfft1+j][1] = 0.0e0;

                } /* Conditional on i == 0 && j == 0 */
            } /* Loop over j */
        } /* Loop over i */



        /* Print frequency domain (debugging) */
        if(debugoutputtodisk == true){
            std::fstream ssfile;
            std::string fname("filter.txt");
            ssfile.open(fname.c_str(),std::ios::out);
            for(size_t i=0; i<nfft1; i++){
                for(size_t j=0; j<nfft2; j++){
                    ssfile << sinc_fft[i*nfft2+j][0] << " " << std::endl;
                }
                ssfile << std::endl;
            }
            ssfile.close();
        }


	} else {
        /* No action required.  The frequency domain of the filter has already been calculated. */
	}
}


template <class T>
void DataFilters<T>::Sinc2D_Real(Array2D<T> &data)
{
	/*
	 * CONVOLVE REAL 2D DATA WITH A SINC FILTER.
	 *
	 * THIS FUNCTION BASED ON SAMPLE CODE FOUND AT
	 * http://people.sc.fsu.edu/~jburkardt/c_src/fftw3/fftw3.html
	 * ON 2 APRIL 2011.
	 */

    size_t n1 = data.GetDim(1);	// NUMBER OF DATA POINTS
    size_t n2 = data.GetDim(2);	// NUMBER OF DATA POINTS


	// REQUIRE THAT FFTW PLANS HAVE BEEN SET
	assert(sinc2d_forward_plan_set == true && sinc2d_backward_plan_set == true);


	/*
	 * CREATE ARRAYS TO HOLD DATA FOR FFTW ROUTINES
	 */
	fftw_complex *data_sp = new fftw_complex[n_sinc_fft];
	fftw_complex *data_fft = new fftw_complex[n_sinc_fft];
	fftw_complex *data_filtered_fft = new fftw_complex[n_sinc_fft];
	fftw_complex *data_filtered_sp = new fftw_complex[n_sinc_fft];

    for(size_t i=0; i<n_sinc_fft; i++){
		data_sp[i][0] = 0.0e0;
		data_sp[i][1] = 0.0e0;
		data_filtered_sp[i][0] = 0.0e0;
		data_filtered_sp[i][1] = 0.0e0;
		data_fft[i][0] = 0.0e0;
		data_fft[i][1] = 0.0e0;
		data_filtered_fft[i][0] = 0.0e0;
		data_filtered_fft[i][1] = 0.0e0;
	}




    /* "Smooth padding"
     *
     * Interpolation scheme used:
     *
     *                      nfftx
     *      |------------|                              |
     *      |------------|                              |
     *      |------------|                              |
     *   n  |-- DATA ----|               I              |
     *   f  |------------|                              |
     *   f  |------------|                              |
     *   t  |------------|                              |
     *   y  |------------|..............................|
     *      |            |                              |
     *      |            |                              |
     *      |            |                              |
     *      |            |                              |
     *      |            |                              |
     *      |     II     |             III              |
     *      |            |                              |
     *      |            |                              |
     *      |            |                              |
     *      |            |                              |
     *
     *  DATA: This region corresponds to the input data to be filtered (array 'rdata').
     *  I:    Padding along the first direction ('nfftx' direction).  This is linearly-
     *        interpolated along the rows so that the last element (right-most) matches the
     *        first element of the row (left-most in the DATA region).
     *  II:   Padding along the second direction ('nffty' direction).  This is linearly-
     *        interpolated along the columns so that the last element (bottom-most) matches
     *        the first element of the column (top-most in the DATA region).
     *  III:  This region is padded using the average of the upper and left boundary values.
     *        The upper-boundary value is taken from the linear interpolation in the bottom
     *        row of region I.  The left-boundary value is taken from the linear interpolation
     *        in the right-most column of region II.
     */

    for(size_t i=0; i<n_sinc_fft2; i++){
        for(size_t j=0; j<n_sinc_fft1; j++){
            if(i < n1 && j < n2){
                data_sp[i*n_sinc_fft1+j][0] = (double)data(i,j);
            } else {
                if(i < n1){
                    /* (i,j) is in region I. */
                    double dval = (double)(data(i,(size_t)0) - data(i,n2-(size_t)1))/(double)(n_sinc_fft1 - n2);
                    data_sp[i*n_sinc_fft1+j][0] = (double)data(i,n2-(size_t)1) + dval*(double)(j - n2);
                }
                if(j < n2){
                    /* (i,j) is in region II. */
                    double dval = (double)(data((size_t)0,j) - data(n1-(size_t)1,j))/(double)(n_sinc_fft2 - n1);
                    data_sp[i*n_sinc_fft1+j][0] = (double)data(n1-(size_t)1,j) + dval*(double)(i - n1);
                }
                if(i >= n1 && j >= n2){
                    /* (i,j) is in region III. */
                    double dval = (double)(data(n1-(size_t)1,(size_t)0) - data(n1-(size_t)1,n2-(size_t)1))/
                            (double)(n_sinc_fft1 - n2);
                    double data1 = (double)data(n1-(size_t)1,n2-(size_t)1) + dval*(double)(j - n2);

                    dval = (double)(data((size_t)0,n2-(size_t)1) - data(n1-(size_t)1,n2-(size_t)1))/
                            (double)(n_sinc_fft2 - n1);
                    double data2 = (double)data(n1-(size_t)1,n2-(size_t)1) + dval*(double)(i - n1);

                    data_sp[i*n_sinc_fft1+j][0] = 0.5e0*(data1 + data2);
                }
            }
        }
    }




	/*
	 * EXECUTE THE FORWARD PLAN
	 */

	//fftw_execute(sinc2d_forward_plan);
	fftw_execute_dft(sinc2d_forward_plan,data_sp,data_fft);


    /*
     * POINT-WISE MULTIPLICATION (SCALING FACTOR IN sinc_fft CALCULATION)
	 */
    for(size_t i=0; i<n_sinc_fft; i++){
		data_filtered_fft[i][0] = (data_fft[i][0]*sinc_fft[i][0] -
                data_fft[i][1]*sinc_fft[i][1]);
		data_filtered_fft[i][1] = (data_fft[i][0]*sinc_fft[i][1] +
                data_fft[i][1]*sinc_fft[i][0]);
	}


	/*
	 * EXECUTE THE BACKWARD PLAN
	 */

	//fftw_execute(sinc2d_backward_plan);
	fftw_execute_dft(sinc2d_backward_plan,data_filtered_fft,data_filtered_sp);

	/*
	 * PLACE DATA INTO INPUT ARRAY, OVERWRITING EXISTING VALUES
	 */
    for(size_t i=0; i<n1; i++){
        for(size_t j=0; j<n2; j++){
            data(i,j) = (T)data_filtered_sp[i*n_sinc_fft1+j][0];
		}
	}

	delete [] data_fft;
	delete [] data_sp;
	delete [] data_filtered_fft;
	delete [] data_filtered_sp;

	return;
}


template <class T>
void DataFilters<T>::FFTW2D_PreparePlans(size_t nfft1, size_t nfft2)
{
	// SET NUMBER OF FFT FREQUENCIES
	n_sinc_fft1 = nfft1;
	n_sinc_fft2 = nfft2;
	n_sinc_fft = nfft1*nfft2;

	// ALLOCATE MEMORY
	fftw_complex *data_sp;
	fftw_complex *data_fft;
	fftw_complex *data_filtered_fft;
	fftw_complex *data_filtered_sp;

	data_sp = new fftw_complex[n_sinc_fft];
	data_filtered_sp = new fftw_complex[n_sinc_fft];

	data_fft = new fftw_complex[n_sinc_fft];
	data_filtered_fft = new fftw_complex[n_sinc_fft];

    for(size_t i=0; i<n_sinc_fft; i++){
		data_sp[i][0] = 0.0e0;
		data_sp[i][1] = 0.0e0;
		data_filtered_sp[i][0] = 0.0e0;
		data_filtered_sp[i][1] = 0.0e0;
		data_fft[i][0] = 0.0e0;
		data_fft[i][1] = 0.0e0;
		data_filtered_fft[i][0] = 0.0e0;
		data_filtered_fft[i][1] = 0.0e0;
	}

	// CREATE FFTW PLANS
	DataFilters<T>::Sinc2D_SetPlan(1,data_sp,data_fft);
	DataFilters<T>::Sinc2D_SetPlan(-1,data_filtered_fft,data_filtered_sp);

	delete [] data_sp;
	delete [] data_fft;
	delete [] data_filtered_sp;
	delete [] data_filtered_fft;
}


template <class T>
void DataFilters<T>::CustomFilter2D_Read(std::string filename)
{
	custom2d_points.ReadCSVFile(filename,1,4,-10.0e0);
	custom2d_fileread = true;
}


template <class T>
void DataFilters<T>::CustomFilter2D_Prepare(size_t nfft1, size_t nfft2, std::string dir,
		std::string filename)
{
	// Create FFTW plans
	custom2d_nfft1 = nfft1;
	custom2d_nfft2 = nfft2;

	fftw_complex *input = new fftw_complex[custom2d_nfft1*custom2d_nfft2];
	fftw_complex *output = new fftw_complex[custom2d_nfft1*custom2d_nfft2];

	delete [] custom2d_fft;
	custom2d_fft = new fftw_complex[custom2d_nfft1*custom2d_nfft2];

    /** @brief Define custom filter flags for FFTW. */
    #define CUSTOM2DFLAGS FFTW_MEASURE|FFTW_UNALIGNED
	//#define CUSTOM2DFLAGS FFTW_ESTIMATE|FFTW_UNALIGNED
	//#define CUSTOM2DFLAGS FFTW_PATIENT|FFTW_UNALIGNED
	//#define CUSTOM2DFLAGS FFTW_EXHAUSTIVE|FFTW_UNALIGNED

	// Define forward plan
	if(custom2d_forward_plan_set == false){
		custom2d_forward_plan = fftw_plan_dft_2d(custom2d_nfft1,custom2d_nfft2,input,output,
				FFTW_FORWARD, CUSTOM2DFLAGS);
		custom2d_forward_plan_set = true;
	}

	// Define backward plan
	if(custom2d_backward_plan_set == false){
		custom2d_backward_plan = fftw_plan_dft_2d(custom2d_nfft1,custom2d_nfft2,input,output,
				FFTW_BACKWARD, CUSTOM2DFLAGS);
		custom2d_backward_plan_set = true;
	}

	delete [] input;
	delete [] output;


	// Define filter in the frequency domain
	T wint1 = 2.0e0*PI/(T)custom2d_nfft1;
	T wint2 = 2.0e0*PI/(T)custom2d_nfft2;
	T w1 = 0.0e0;
	T w2 = 0.0e0;
	fftw_complex filter_amplitude;
	filter_amplitude[0] = 1.0e0;
	filter_amplitude[1] = 0.0e0;
	for(int i=0; i<custom2d_nfft1; i++){
		w1 = (T)i*wint1;
		if(w1 > PI){
			w1 -= 2.0e0*PI;
		}

		for(int j=0; j<custom2d_nfft2; j++){
			w2 = (T)j*wint2;
			if(w2 > PI){
				w2 -= 2.0e0*PI;
			}

			filter_amplitude[0] = CustomFilter2D_GetAmplitudeReal(w1,w2);

			custom2d_fft[i*custom2d_nfft2+j][0] = filter_amplitude[0];
			custom2d_fft[i*custom2d_nfft2+j][1] = filter_amplitude[1];
		}
	}

    // Copy filter data to UniformVolumeCore object and write VTK file to disk.
	// Note that this assumes the filter is real.  A complex filter will need an
    // alternate approach since UniformVolumeCore doesn't support direct access
	// to vector data.
    UniformVolume<T> filteroutput(custom2d_nfft1,custom2d_nfft2,1,0.0e0);
	filteroutput.SetDir(dir);
	filteroutput.SetFileStem(filename);
	filteroutput.SetImageOutput();
	filteroutput.SetIntensityOutput();
	filteroutput.SetXmin(-PI);
	filteroutput.SetXmax(PI);
	filteroutput.SetYmin(-PI);
	filteroutput.SetYmax(PI);
	filteroutput.SetZmin(0.0e0);
	filteroutput.SetZmax(0.0e0);

	// Shift filter so it spans [-PI,PI] rather than [0,2*PI] and the zero-fequency
	// position is centered in the plane.
	int ii = custom2d_nfft1/2;
	int jj = custom2d_nfft2/2;
	for(int i=0; i<custom2d_nfft1; i++){
		jj = custom2d_nfft2/2;
		for(int j=0; j<custom2d_nfft2; j++){
			filteroutput(i,j,0) = custom2d_fft[ii*custom2d_nfft2+jj][0];
			jj++;
			if(jj == custom2d_nfft2){ jj = 0; }
		}
		ii++;
		if(ii == custom2d_nfft1){ ii = 0; }
	}
	filteroutput.VTKWrite();
}


template <class T>
void DataFilters<T>::CustomFilter2D_Real(Array2D<T> &data)
{
	/*
	 * CONVOLVE REAL 2D DATA WITH A CUSTOM FILTER.
	 */

	int n1 = data.GetDim(1);	// NUMBER OF DATA POINTS
	int n2 = data.GetDim(2);	// NUMBER OF DATA POINTS


	// REQUIRE THAT FFTW PLANS HAVE BEEN SET
	assert(custom2d_forward_plan_set == true && custom2d_backward_plan_set == true);


	/*
	 * CREATE ARRAYS TO HOLD DATA FOR FFTW ROUTINES
	 */
	fftw_complex *data_sp = new fftw_complex[custom2d_nfft1*custom2d_nfft2];
	fftw_complex *data_fft = new fftw_complex[custom2d_nfft1*custom2d_nfft2];
	fftw_complex *data_filtered_fft = new fftw_complex[custom2d_nfft1*custom2d_nfft2];
	fftw_complex *data_filtered_sp = new fftw_complex[custom2d_nfft1*custom2d_nfft2];

	for(int i=0; i<custom2d_nfft1*custom2d_nfft2; i++){
		data_sp[i][0] = 0.0e0;
		data_sp[i][1] = 0.0e0;
		data_filtered_sp[i][0] = 0.0e0;
		data_filtered_sp[i][1] = 0.0e0;
		data_fft[i][0] = 0.0e0;
		data_fft[i][1] = 0.0e0;
		data_filtered_fft[i][0] = 0.0e0;
		data_filtered_fft[i][1] = 0.0e0;
	}

	for(int i=0; i<custom2d_nfft1; i++){
		for(int j=0; j<custom2d_nfft2; j++){
			if(i < n1 && j < n2){
				data_sp[i*custom2d_nfft2+j][0] = (double)data(i,j);
			}
		}
	}


	/*
	 * EXECUTE THE FORWARD PLAN
	 */
	fftw_execute_dft(custom2d_forward_plan,data_sp,data_fft);

	double scale = 1.0e0/(double)(custom2d_nfft1*custom2d_nfft2);

	/*
	 * POINT-WISE MULTIPLICATION
	 */
	for(int i=0; i<custom2d_nfft1*custom2d_nfft2; i++){
		data_filtered_fft[i][0] = (data_fft[i][0]*custom2d_fft[i][0] -
				data_fft[i][1]*custom2d_fft[i][1])*scale;
		data_filtered_fft[i][1] = (data_fft[i][0]*custom2d_fft[i][1] +
				data_fft[i][1]*custom2d_fft[i][0])*scale;
	}


	/*
	 * EXECUTE THE BACKWARD PLAN
	 */
	fftw_execute_dft(custom2d_backward_plan,data_filtered_fft,data_filtered_sp);

	/*
	 * PLACE DATA INTO INPUT ARRAY, OVERWRITING EXISTING VALUES
	 */
	for(int i=0; i<n1; i++){
		for(int j=0; j<n2; j++){
			data(i,j) = (T)data_filtered_sp[i*custom2d_nfft2+j][0];
		}
	}

	delete [] data_fft;
	delete [] data_sp;
	delete [] data_filtered_fft;
	delete [] data_filtered_sp;

	return;
}


template <class T>
void DataFilters<T>::Test(std::string &result)
{
    result = "SUCCESS";
    bool errorfound = false;

    /* Create test object. */
    DataFilters<T> df1;
    size_t nfft1 = 8;
    size_t nfft2 = 16;
    T exp1 = (T)1.0e0;
    T exp2 = (T)2.0e0;
    df1.Sinc2D_GetFFT(nfft1, nfft2, exp1, exp2);



    /* Test copy constructor. */
    DataFilters<T> df2(df1);


    /* Test copy assignment. */
    DataFilters<T> df3;
    df3 = df1;


    if(!errorfound){
        errorfound = true;
    }

}
