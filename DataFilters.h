/**
 * @file DataFilters.h
 * @author 	Robert Grandin
 * @date 18 November 2010
 * @brief Definition of DataFilters class.
 *
 * @section Class Description & Notes
 *
 * This class contains routines necessary for applying filters to datasets.
 *
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 18 November 2010
 *      - Creation date.
 * @date 5 October 2012
 *      - Removed unused filter functions as part of CT reconstruction program
 *        rewrite.
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2010-2012, Robert Grandin
 * All rights reserved.
 *
 * Redistribution and use of this file is permitted provided that the following
 * conditions are met:
 * 	-# 	Redistributions must produce the above copyright notice, this list of
 * 		conditions, and the following disclaimer in the documentation and/or
 * 		other materials provided with the distribution.
 * 	-#	Neither the name of the organization nor the names of its contributors
 * 		may be used to endorse or promote products derived from this software
 * 		without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
 * EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING BUT NOT
 * LIMITING TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 */


#ifndef PI_
#define PI_

/** @brief Define mathematical constant @f$ \pi @f$ */
#define PI (T)3.14159265358979e0

#endif


#ifndef DataFilters_
/** @brief Define symbol for DataFilters class. */
#define DataFilters_



#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <complex>
#include <Array1D.h>
#include <Array2D.h>
#include <UniformVolume3D.h>


// ==================================================================
// ================
// ================    FUNCTION PROTOTYPES AND DESCRIPTIONS
// ================

/**
 * @brief Class definition for filtering routines.
 *
 * Functions referring to a sinc function are a ramp filter in the frequency domain.
 * A sinc function in the spatial/temporal domain corresponds to a ramp function in
 * the frequency domain.
 */
template <class T>
class DataFilters {

public:

	/**
	 * @brief Constructor.
	 * @pre Memory exists for temporary data storage.
	 * @post DataFilters object exists.
	 * @return None.
	 */
	DataFilters();


	/**
	 * @brief Deconstructor.
	 * @pre DataFilters object exists.
	 * @post DataFilters object destroyed.
	 * @return None.
	 */
    virtual ~DataFilters();


	/**
	 * @brief One-dimensional median filter.
	 * @pre Input data created.
	 * @param vec Reference to 1-D array to be filtered.
	 * @param kernalsize Size of median kernel.
	 * @param endcond End condition.  1:Periodic, 2:Repeated-end-value
	 * @param overwrite Boolean controlling if filtered and unfiltered values
	 * 					can be mixed during the filtering process.
	 * @post No data external to this routine modified.
	 * @return None.  Input array is overwritten by filtered data.
	 */
    void Median1D(Array1D<T> &vec, int kernalsize, int endcond,
			bool overwrite);


	/**
     * @brief 1-Dimensional sinc filter as described online at
     *  http://www.clear.rice.edu/elec431/projects96/DSP/sincfilt.html
     *  on 3 April 2011.
     * @pre Input data created.
	 * @param data Reference to the 1D input data.
     * @param rolloff Number between 0 and 1 that controls the high-frequency
     *  roll-off.
     * @param expatten Exponential attenuation factor.
	 * @param nfft Number of FFT freqencies to be used.
	 * @post Input data is replaced by filter data.
	 * @return None.  Input data is overwritten by the filter results.
     * @image html sinc_rolloff_a0.png "Rolloff value of 0"
     * @image html sinc_rolloff_a1.png "Rolloff value of 1"
     * @image html sinc_rolloff_a2.png "Rolloff value of 2"
     * @image html sinc_rolloff_a3.png "Rolloff value of 3"
     * @image latex sinc_rolloff_a0.png "Rolloff value of 0" width=\textwidth
     * @image latex sinc_rolloff_a1.png "Rolloff value of 1" width=\textwidth
     * @image latex sinc_rolloff_a2.png "Rolloff value of 2" width=\textwidth
     * @image latex sinc_rolloff_a3.png "Rolloff value of 3" width=\textwidth
	 */
    void Sinc1D_Real(Array1D<T> &data, T rolloff, T expatten, size_t nfft);


	/**
     * @brief 2-Dimensional sinc filter as described online at
     *  http://www.clear.rice.edu/elec431/projects96/DSP/sincfilt.html
     *  on 3 April 2011.
	 * @pre Input data created.
	 * @param data Reference to the 2D input data.
	 * @post Input data is replaced by filter data.
	 * @return None.  Input data is overwritten by the filter results.
	 */
	void Sinc2D_Real(Array2D<T> &data);


	/**
     * @brief Determine FFT of 2D sinc filter.
     *
     * The filter is a multiplication of a
     *  linear ramp filter and an exponential decay (both defined in the frequency
     *  domain as a function of the magnitude of the frequency, |w1,w2|).
	 * @pre Input data created and sinc filter defined.
     * @param nfft1 Number of FFT points to be used in the first dimension (across
     *  rows).
     * @param nfft2 Number of FFT points to be used in the second dimension
     *  (across columns).
     * @param expatten1 Exponential attenuation associated with the first dimension.
     * @param expatten2 Exponential attenuation associated with the second dimension.
	 * @post FFT performed and defined.
	 * @return None.
	 */
    void Sinc2D_GetFFT(const size_t nfft1, const size_t nfft2, const T expatten1, const T expatten2);


	/**
	 * @brief Perform preparations for 2D FFT.
	 * @pre DataFilters object exists.
	 * @param nfft1 Number of FFT frequencies in the first dimension.
	 * @param nfft2 Number of FFT frequencies in the second dimension.
	 * @post FFTW plans created for forward and backward transforms.
	 */
	void FFTW2D_PreparePlans(size_t nfft1, size_t nfft2);


	/**
	 * @brief Read CSV file containing definition of custom 2D filter.
	 * @pre DataFilters object exists.
	 * @param filename Name of file to be read.
	 * @post Filter definition read.
	 * @return None.
	 * @warning First row in data file must be column headers.  Subsequent rows are
	 * 		expected to be in ascending order by frequency.  Negative frequency values
	 * 		will cause that frequency/amplitude pair to be ignored.
	 */
	void CustomFilter2D_Read(std::string filename);


	/**
     * @brief Prepare custom 2D filter.
     *
     * Number of FFT frequencies is defined and
     *  the filter is calculated in the frequency-domain.  Frequency-domain
     *  representation is output as a VTK file.
	 * @pre DataFilters object exists.
	 * @param nfft1 Number of FFT frequencies in the first dimension.
	 * @param nfft2 Number of FFT frequencies in the second dimension.
	 * @param dir Directory to output filter VTK file.
	 * @param filename Name to use for filter VTK file.
	 * @post Filter defined and ready for application.
	 * @return None.
	 */
    void CustomFilter2D_Prepare(size_t nfft1, size_t nfft2, std::string dir,
			std::string filename);


	/**
	 * @brief Apply custom 2D filter.
	 * @pre DataFilters object exists.
	 * @param data Reference to 2D array of data which is to be filtered.
	 * @post Filtered result written to input array.
	 * @warning Input data is lost.
	 * @return None.
	 */
	void CustomFilter2D_Real(Array2D<T> &data);




private:

	/**
     * @brief Function for comparing two input values (required by qsort).
     *
     * This
     *  simply identifies which input value is larger, with no speical
     *  processing of the numbers.  Must be declared "static" to work
     *  within this class.
	 * @pre Input data created.
	 * @param a First input value.
	 * @param b Second input value.
	 * @post No data modified.
	 * @return Integer identifying which input value is larger.
	 */
	static int compare(const void *a, const void *b);


	/**
	 * @brief Determine FFT of 1D sinc filter.
     * @pre Input data created and sinc filter defined.
     * @param nfft Number of FFT points to be used.
     * @param rolloff Number between 0 and 1 that controls the high-frequency
     *  roll-off.
     * @param expatten Exponential attenuation factor.
	 * @post FFT performed and defined.
     * @image html sinc_rolloff_a0.png "Rolloff value of 0"
     * @image html sinc_rolloff_a1.png "Rolloff value of 1"
     * @image html sinc_rolloff_a2.png "Rolloff value of 2"
     * @image html sinc_rolloff_a3.png "Rolloff value of 3"
     * @image latex sinc_rolloff_a0.png "Rolloff value of 0" width=\textwidth
     * @image latex sinc_rolloff_a1.png "Rolloff value of 1" width=\textwidth
     * @image latex sinc_rolloff_a2.png "Rolloff value of 2" width=\textwidth
     * @image latex sinc_rolloff_a3.png "Rolloff value of 3" width=\textwidth
	 * @return None.
	 */
    void Sinc1D_GetFFT(size_t nfft, T rolloff, T expatten);



	/**
	 * @brief Determine FFTW plan for the 2D sinc filter.
	 * @pre DataFilters object exists.
	 * @param plan Identify which plan to set:
	 * 		- 1:Foward plan
	 * 		- -1:Backward plan
	 * @param input Data input to transform.
	 * @param output Transform result.
	 * @post FFTW plan defined.
	 * @return None.
	 */
	void Sinc2D_SetPlan(int plan, fftw_complex *input, fftw_complex *output);


	/**
	 * @brief Determine custom 2D filter amplitude for given frequencies.
	 * @pre DataFilters object exists.
	 * @param w1 Angular frequency in the first dimension.
	 * @param w2 Angular frequency in the second dimension.
	 * @post No object attributes changed.
	 * @return Filter amplitude for given frequency components.  This is the REAL
	 * 		portion of the filter.
	 * @warning Input frequencies are expected to be in the range [-PI,PI] and
	 * 		be in unis of [radians/second] (i.e., angular frequency).
	 */
	T CustomFilter2D_GetAmplitudeReal(const T w1, const T w2) const;


	/**
	 * @brief Called by CustomFilter2D_GetAmplitudeReal when a radially-symmetric
	 * 		filter is required.
	 * @pre DataFilters object exists.
	 * @param w1 Angular frequency in the first dimension.
	 * @param w2 Angular frequency in the second dimension.
	 * @post No object attributes changed.
	 * @return Filter amplitude for given frequency components.  This is the REAL
	 * 		portion of the filter.
	 * @warning Input frequencies are expected to be in the range [-PI,PI] and
	 * 		be in unis of [radians/second] (i.e., angular frequency).
	 */
	T CustomFilter2D_GetAmplitudeRadialReal(const T w1, const T w2) const;


	/**
	 * @brief Determine interpolated filter amplitude.  This is for a single dimension.
	 * @pre DataFilters object exists.
	 * @param dir Direction with which interpolation is to be done.
	 * 		- 0: First user-defined dimention
	 * 		- 2: Second user-defined dimension
	 * 		- 4: Third user-defined dimension
	 * 		- ...and so on
	 * @param w Frequency at which filter magnitude is to be interpolated.  This is
	 * 		expected to be in the range [0,PI].
	 * @param method Method of interpolation used.  Default is piecewise-linear.
	 * 		- 0: Piecewise-linear
	 * @post No object attributes changed.
	 * @return Filter amplitude for given frequency.  This is the REAL portion of
	 * 		the filter.
	 * @warning Input frequencies are expected to be in the range [0,PI] and
	 * 		be in unis of [radians/second] (i.e., angular frequency).
	 */
	T CustomFilter_InterpolateMagnitude(const int dir, const T w, const int method) const;



    /** @brief Number of elements in the frequency domain sinc filter */
    size_t n_sinc_fft;

	/** @brief Number of elements in the frequency domain sinc filter, along the
	 * first dimension */
    size_t n_sinc_fft1;

	/** @brief Number of elements in the frequency domain sinc filter, along the
	 * second dimension */
    size_t n_sinc_fft2;

	/** @brief Sinc filter in frequency domain */
	fftw_complex *sinc_fft;

	/** @brief Sinc 2D filter FFTW forward-plan */
	fftw_plan sinc2d_forward_plan;

	/** @brief Sinc 2D filter FFTW backward-plan */
	fftw_plan sinc2d_backward_plan;

	/** @brief Track if sinc 2D filter FFTW forward-plan has been set */
	bool sinc2d_forward_plan_set;

	/** @brief Track if sinc 2D filter FFTW backward-plan has been set */
	bool sinc2d_backward_plan_set;

	/** @brief Track if custom 2D filter FFTW forward-plan has been set */
	bool custom2d_forward_plan_set;

	/** @brief Track if custom 2D filter FFTW backward-plan has been set */
	bool custom2d_backward_plan_set;

	/** @brief Number of FFT frequencies in first dimension for custom 2D filter */
    size_t custom2d_nfft1;

	/** @brief Number of FFT frequencies in second dimension for custom 2D filter */
    size_t custom2d_nfft2;

	/** @brief 2D array storing (frequency,amplitude) pairs for 2D filter */
	Array2D<T> custom2d_points;

	/** @brief Frequency-domain definition of the custom 2D filter */
	fftw_complex *custom2d_fft;

	/** @brief FFTW forward plan for custom 2D filter */
	fftw_plan custom2d_forward_plan;

	/** @brief FFTW backward plan for custom 2D filter */
	fftw_plan custom2d_backward_plan;

	/** @brief Track if custom filter definition file has been read */
	bool custom2d_fileread;


};

#include "DataFilters.cpp"



#endif /* DataFilters_ */
