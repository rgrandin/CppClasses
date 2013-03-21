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


namespace basicXdmfTest{


void writeXDMFFile(size_t size1, size_t size2, size_t size3);

}


#endif // basicXdmfTest_H
