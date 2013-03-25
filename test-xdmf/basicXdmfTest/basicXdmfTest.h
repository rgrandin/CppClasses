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

#include <Array1D.h>
#include <Array2D.h>
#include <Array3D.h>
#include <Array4D.h>
#include <StringManip.h>
#include <UniformVolume.h>
#include <XdmfIO.h>

#include <omp.h>


/** @brief Namespace demonstrating basic XDMF usage. */
namespace basicXdmfTest{


/**
 * @brief Write a single array to disk using XDMF/HDF5 format.
 *
 *  User will be asked for array dimensions and file name during function execution.
 */
void writeXDMFFile();


/**
 * @brief Write a multiple arrays to disk using XDMF/HDF5 format.
 *
 *  User will be asked for array quantity, array dimensions and file name during
 *  function execution.
 */
void writeXDMFFile_MultiArray();


/**
 * @brief Read an XDMF/HDF5 file and write data information to std::cout.
 *
 *  User will be asked for filename during function execution.
 */
void readXDMFFile();


/**
 * @brief Convert a volume file (VTK or CNDE VOL format) to XDMF/HDF5 format.
 *
 *  User will be asked for filenames during function execution.
 */
void convertVolume();


/**
 * @brief Write XDMF/HDF5 file with multiple arrays of multiple dimensionality, size, and datatype.
 */
void writeXDMFFile_MultiTest();



}


#endif // basicXdmfTest_H
