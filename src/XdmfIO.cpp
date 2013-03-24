/**
 * @file XdmfIO.cpp
 * @author Robert Grandin
 * @brief Implementation of XdmfIO routines for array writing/reading using XDMF/HDF5 format.
 */


#include "XdmfIO.h"

void XdmfIO::writeUniformArrays(std::string filename,
                                XdmfIO::data_info<char> *char_data,
                                XdmfIO::data_info<unsigned char> *uchar_data,
                                XdmfIO::data_info<short> *short_data,
                                XdmfIO::data_info<unsigned short> *ushort_data,
                                XdmfIO::data_info<int> *int_data,
                                XdmfIO::data_info<unsigned int> *uint_data,
                                XdmfIO::data_info<long> *long_data,
                                XdmfIO::data_info<float> *float_data,
                                XdmfIO::data_info<double> *double_data,
                                int compression)

