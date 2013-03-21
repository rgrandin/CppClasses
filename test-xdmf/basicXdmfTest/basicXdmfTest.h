#ifndef basicXdmfTest_H
#define basicXdmfTest_H

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include <hdf5.h>
#include <Xdmf.h>
#include <XdmfH5Driver.h>

#include <Array3D.h>
#include <StringManip.h>
#include <UniformVolume.h>

#include <omp.h>


namespace basicXdmfTest{


void writeXDMFFile();


void convertVolume();


}


#endif // basicXdmfTest_H
