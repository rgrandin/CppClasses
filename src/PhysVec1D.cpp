/**
 * @file PhysVec1D.cpp
 * @author Robert Grandin
 * @brief Implementation of PhysVec1D class.
 */


#include <PhysVec1D.h>



template <class T>
PhysVec1D<T>::PhysVec1D()
{

}


template <class T>
PhysVec1D<T>::PhysVec1D(PhysVec1D<T> &a) : PhysVec1D()
{
    PhysVec1DSwap(*this, a);
}


#ifdef CXX11
template <class T>
PhysVec1D<T>::PhysVec1D(PhysVec1D<T> &&a) : PhysVec1D()
{
    PhysVec1DSwap(*this, a);
}
#endif


template <class T>
PhysVec1D<T>::PhysVec1D(size_t n, const T initval)
{
    ArrayBase<T>::initialize(n,(T)initval);
}


template <class T>
PhysVec1D<T>::~PhysVec1D()
{

}


template <class T>
PhysVec1D<T>& PhysVec1D<T>::operator=(PhysVec1D<T> a)
{
    PhysVec1DSwap(*this, a);
    return *this;
}


#ifdef CXX11
template <class T>
PhysVec1D<T>& PhysVec1D<T>::operator=(PhysVec1D<T> &&a)
{
    PhysVec1DSwap(*this, a);
    return *this;
}
#endif


template <class T>
const T PhysVec1D<T>::Norm(const T p) const
{
	T retval = 0.0e0;
	if(p > 0){
		T sum = (T)0.0e0;

		for(int i=0; i<ArrayBase<T>::npoints; i++){
			sum += pow(ArrayBase<T>::array[i],p);
		}

		retval = pow(sum,1.0e0/p);
	}

	if(p < 0){
		retval = ArrayBase<T>::MaxVal();
	}

	return retval;
}


template <class T>
const T PhysVec1D<T>::Dot(const PhysVec1D<T> &vec) const
{
	T sum = (T)0.0e0;

	for(int i=0; i<ArrayBase<T>::npoints; i++){
		sum += ArrayBase<T>::array[i]*vec(i);
	}

	return sum;
}


template <class T>
PhysVec1D& PhysVec1D<T>::Cross(const PhysVec1D<T> &ivec)
{
	// CHECK VECTOR DIMENSIONS
	assert(ArrayBase<T>::npoints == 2 || ArrayBase<T>::npoints == 3);
	assert(ivec.GetDim() == 2 || ivec.GetDim() == 3);

    PhysVec1D result(ivec);

	T a1 = (T)0.0e0;
	T a2 = (T)0.0e0;
	T a3 = (T)0.0e0;
	T b1 = (T)0.0e0;
	T b2 = (T)0.0e0;
	T b3 = (T)0.0e0;
	a1 = ArrayBase<T>::array[0];
	a2 = ArrayBase<T>::array[1];
	if(ArrayBase<T>::npoints == 3){ a3 = ArrayBase<T>::array[2]; }
	b1 = ivec(0);
	b2 = ivec(1);
	if(ivec.GetDim() == 3){ b3 = ivec(2); }

    result(0) = a2*b3 - a3*b2;
    result(1) = a3*b1 - a1*b3;
    result(2) = a1*b2 - a2*b1;

    return &result;
}


template <class T>
void PhysVec1D<T>::Rotate(const T alpha, const T beta, const T gamma)
{
	T px0 = ArrayBase<T>::array[0];
	T py0 = ArrayBase<T>::array[1];
	T pz0 = ArrayBase<T>::array[2];



	T cosalpha = cos(alpha);
	T sinalpha = sin(alpha);
	T cosbeta = cos(beta);
	T sinbeta = sin(beta);
	T cosgamma = cos(gamma);
	T singamma = sin(gamma);

	ArrayBase<T>::array[0] = px0*cosbeta*cosgamma - py0*cosbeta*singamma + pz0*sinbeta;
	ArrayBase<T>::array[1] = px0*(cosalpha*singamma + cosgamma*sinalpha*sinbeta) +
			py0*(cosalpha*cosgamma - sinalpha*sinbeta*singamma) -
			pz0*cosbeta*sinalpha;
	ArrayBase<T>::array[2] = px0*(sinalpha*singamma - cosalpha*cosgamma*sinbeta) +
			py0*(cosgamma*sinalpha + cosalpha*sinbeta*singamma) +
			pz0*cosalpha*cosbeta;

	/*
	COPIED FROM PREVIOUS WORK IN THE GENERAL-MOTION RECONSTRUCTION ROUTINE

	alpha <-> gamma
	beta <-> psi
	gamma <-> theta

	px = px0*cospsi*costheta - py0*cospsi*sintheta + pz0*sinpsi;
	py = px0*(cosgamma*sintheta + costheta*singamma*sinpsi) +
			py0*(cosgamma*costheta - singamma*sinpsi*sintheta) -
			pz0*cospsi*singamma;
	pz = px0*(singamma*sintheta - cosgamma*costheta*sinpsi) +
			py0*(costheta*singamma + cosgamma*sinpsi*sintheta) +
			pz0*cosgamma*cospsi;
	*/
}
