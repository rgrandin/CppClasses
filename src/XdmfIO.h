/**
 * @file XdmfIO.h
 * @author 	Robert Grandin
 * @date 23 March 2013
 * @brief Definition of namespace to facilitate data IO with XDMF/HDF5 format.
 *
 * @section Description & Notes
 *
 * This class implements functions to read/write data contained within the Array1D,
 * Array2D, Array3D, and Array4D data structures.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 23 March 2013
 *	- Creation date.
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2013, Robert Grandin
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



#ifndef XDMFIO_H
#define XDMFIO_H


#include <stdlib.h>
#include <iostream>
#include <stdio.h>

#include <Xdmf.h>

#include <PArray1D.h>
#include <Array1D.h>


/** @brief Namespace containing functions for reading and writing Array data using the XDMF file format. */
namespace XdmfIO
{


template <typename T>
struct data_info {
    /** @brief Array of pointers to data arrays. */
    PArray1D<T*> *data;

    /** @brief Array of pointers to Array1D objects containing actual data dimensions in row-major (IJK) order. */
    PArray1D<Array1D<size_t>*> *dims;

    /** @brief Array of pointers to Array1D objects containing origin values for each data dimension. */
    PArray1D<Array1D<float>*> *origin;

    /** @brief Array of pointers to Array1D objects containing spacing values for each data dimension. */
    PArray1D<Array1D<float>*> *spacing;

    /** @brief Array of names for data arrays. */
    PArray1D<std::string*> *data_name;
};



/**
 * @brief Define uniform data array and insert into tree grid.
 * @param data_struct Pointer to XdmfIO::data_info structure containing data to be written.
 * @param tree_grid Pointer to XdmfGrid which is to contain the data.
 * @param xarray Pointer to array of pointers for XdmfArray arrays to be defined and written.
 * @param grid Pointer to array of pointers for XdmfGrid objects for each XdmfArray.
 * @param topo Pointer to array of pointers for XdmfTopology objects for each XdmfGrid object.
 * @param geo Pointer to array of pointers for XdmfGeometry objects for each XdmfGrid object.
 * @param array_count Number of arrays already created.  This is used to determine the appropriate offsets to
 *  use when accessing the arrays of pointers so as to not overwrite previously-defined data.
 * @param filename Name of output file, without extension.
 * @param compression Level of compression to be used.  Value must be in the range [0,9], with 0 being no compression and
 *  9 being maximum compression.
 * @warning As each array is added, the pointer-to-data stored in data_struct is set to NULL.
 */
template <class T>
void addUniformArrays(XdmfIO::data_info<T> *data_struct, XdmfGrid *tree_grid, XdmfArray **xarray,
                      XdmfGrid **grid, XdmfTopology **topo, XdmfGeometry **geo, int &array_count,
                      std::string filename, int compression)
{
    int narrays = (int)data_struct->data->GetDim();


    for(int i=0; i<narrays; i++){

        int array_number = i + array_count;

        /* Create an XdmfArray object to represent the data.  No values are actually moved out of the
         * pre-existing array (in this case the Array3D object populated above).  Deletion of the XdmfArray
         * will not affect the actual data, either.  The XdmfArray object only facilitates the interface
         * between the actual data and the data-output. */

        XdmfPointer data_pointer = data_struct->data->operator()(i);
        std::string arrayname = data_struct->data_name->operator()(i)->substr();


        /* Check number of quantities defined for data dimensions, origin, and spacing.  If only one set of
         * geometry information is specified, use it for all data arrays.  If multiple sets of geometry info
         * are specified, use the set corresponding to the data array. */
        Array1D<size_t> *data_dims;
        Array1D<float> *data_origin;
        Array1D<float> *data_spacing;
        if(data_struct->dims->GetDim() == 1 && data_struct->origin->GetDim() == 1 && data_struct->spacing->GetDim() == 1){
            data_dims = data_struct->dims->operator()(0);
            data_origin = data_struct->origin->operator()(0);
            data_spacing = data_struct->spacing->operator()(0);
        } else {
            data_dims = data_struct->dims->operator()(i);
            data_origin = data_struct->origin->operator()(i);
            data_spacing = data_struct->spacing->operator()(i);
        }




        size_t rank = data_dims->GetDim();
        XdmfInt64 *shape = new XdmfInt64[rank];
        for(size_t r=0; r<rank; r++){
            shape[r] = (XdmfInt64)data_dims->operator()(rank-r-1);
        }

        XdmfInt32 number_type = XDMF_UNKNOWN_TYPE;
        if(typeid(T) == typeid(char)){
            number_type = XDMF_INT8_TYPE;
        }
        if(typeid(T) == typeid(unsigned char)){
            number_type = XDMF_UINT8_TYPE;
        }
        if(typeid(T) == typeid(short)){
            number_type = XDMF_INT16_TYPE;
        }
        if(typeid(T) == typeid(unsigned short)){
            number_type = XDMF_UINT16_TYPE;
        }
        if(typeid(T) == typeid(int)){
            number_type = XDMF_INT32_TYPE;
        }
        if(typeid(T) == typeid(unsigned int)){
            number_type = XDMF_UINT32_TYPE;
        }
        if(typeid(T) == typeid(float)){
            number_type = XDMF_FLOAT32_TYPE;
        }
        if(typeid(T) == typeid(double)){
            number_type = XDMF_FLOAT64_TYPE;
        }


        xarray[array_number] = new XdmfArray;
        xarray[array_number]->SetAllowAllocate(false);                 /* Prevent the XdmfArray from allocating any memory. */
        xarray[array_number]->SetShape((XdmfInt32)rank, shape);        /* Set rank and dimension sizes of array. */
        xarray[array_number]->SetDataPointer(data_pointer);            /* Set pointer to previously-existing data. */
        xarray[array_number]->SetNumberType(number_type);              /* Identify the data type. */


        grid[array_number] = new XdmfGrid;                  /* Create a grid for the data. */
        grid[array_number]->SetName("Structured Grid");     /* Name the grid.  This is user-choice. */
        grid[array_number]->SetGridType(XDMF_GRID_UNIFORM); /* Set grid type.  This example is for image data, which is
                                                             * on a uniform grid.  Grid choice can be other types based
                                                             * on needs. */

        topo[array_number] = new XdmfTopology;                      /* Create a topology object to help define grid. */
        topo[array_number] = grid[array_number]->GetTopology();     /* Set to be grid topology. */

        topo[array_number]->SetTopologyType(XDMF_3DCORECTMESH);     /* Set grid type, again this is for image data and other
                                                                     * options are possible for other applications. */
        topo[array_number]->GetShapeDesc()->SetShape((XdmfInt32)rank, shape);  /* Set shape to match that of the data array. */


        geo[array_number] = new XdmfGeometry;                   /* Create a geometry object to complete grid
                                                                 * structure definition. */
        geo[array_number] = grid[array_number]->GetGeometry();  /* Set to be the grid geometry. */
        if(data_origin && data_spacing){
            /* If both data_origin and data_spacing is defined, use those values when defining the geometry. */

            XdmfFloat64 *origin = new XdmfFloat64[rank];
            for(size_t r=0; r<rank; r++){
                origin[r] = data_origin->operator()(rank-r-1);
            }

            XdmfFloat64 *spacing = new XdmfFloat64[rank];
            for(size_t r=0; r<rank; r++){
                spacing[r] = data_spacing->operator()(rank-r-1);
            }

            geo[array_number]->SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ);   /* Set for image data.  Again, more-complex options exist. */
            geo[array_number]->SetOrigin(origin);                              /* Set origin. */
            geo[array_number]->SetDxDyDz(spacing);                             /* Set point spacing. */
        } else {
            /* If origin/spacing information is missing, assume values.  Reading the produced file into ParaView
             * fails if origin/spacing information is not provided. */

            XdmfFloat64 origin[3] = { 0.0, 0.0, 0.0 };
            XdmfFloat64 spacing[3] = { 1.0, 1.0, 1.0 };

            geo[array_number]->SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ);   /* Set for image data.  Again, more-complex options exist. */
            geo[array_number]->SetOrigin(origin);                              /* Set origin. */
            geo[array_number]->SetDxDyDz(spacing);                             /* Set point spacing. */
        }


        /* Insert this data grid into the super-grid. */
        tree_grid->Insert(grid[array_number]);

        XdmfAttribute *attrib = new XdmfAttribute;              /* Create attribute to describe data. */
        attrib->SetName(arrayname.c_str());                     /* Set descriptive name for data. */
        attrib->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE); /* Define data to exist at nodes (rather than within
                                                                 * cells or on cell faces). */
        attrib->SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);   /* Define data to be scalar. */
        attrib->SetValues(xarray[array_number]);                              /* Define data values to be those stored in the
                                                                 * XdmfArray object created above (which in-turn takes
                                                                 * us to the c-array represented by the XdmfArray object. */

        /* Build name of file/path to be used with HDF5 heavy-data file.  This is set so that the
         * HDF5 filename matches that of the XML file.  The array name is used to identify the data within the
         * HDF5 file.  All arrays will be placed in the same HDF5 file. */
        std::string heavydataname(filename);
        heavydataname = heavydataname + ".h5:/" + arrayname;

        xarray[array_number]->SetHeavyDataSetName(heavydataname.c_str());  /* Set name used for heavy-data.  This will cause the
                                                              * actual data values to be written to the specified
                                                              * HDF5 file/array.  If the array is small-enough,
                                                              * values will be directly-inserted into the XML.  The
                                                              * default value for this threshold is 100, but can be
                                                              * modified by the user.  If data values are located
                                                              * within the XML file, this function call has no
                                                              * effect. */

        xarray[array_number]->SetCompression(compression);                 /* Set compression level to be used by HDF5 file.  This
                                                              * is also ignored if data is embedded in the XML file.
                                                              * Compression method will be noted in the HDF5 file and
                                                              * decompression will happen automatically when the
                                                              * HDF5 file is opened in a suitable viewer or these
                                                              * Xdmf____ classes. */

        /* Insert attribute into grid.  Note that this grid is a child of the "super-grid". */
        grid[array_number]->Insert(attrib);

        /* Set pointer to data to NULL.  Actual data is not affected.  This prevents attempted double-deletion of
         * the data when both the PArray1D object (pointers to data) and Array object (the actual data) are destroyed. */
        data_struct->data->operator()(i) = NULL;

    } /* Loop through arrays to be added. */

    /* Increment array offset. */
    array_count += narrays;

} /* XdmfIO::addUniformArrays() */





/**
 * @brief Define uniform data array and insert into tree grid.
 * @param data_struct Pointer to XdmfIO::data_info structure containing data to be written.
 * @param grid Pointer to XdmfGrid which is to contain the data.
 * @param xarray Pointer to array of pointers for XdmfArray arrays to be defined and written.
 * @param array_count Number of arrays already created.  This is used to determine the appropriate offsets to
 *  use when accessing the arrays of pointers so as to not overwrite previously-defined data.
 * @param filename Name of output file, without extension.
 * @param compression Level of compression to be used.  Value must be in the range [0,9], with 0 being no compression and
 *  9 being maximum compression.
 * @warning As each array is added, the pointer-to-data stored in data_struct is set to NULL.
 */
template <class T>
void addUniformArraysExistingGrid(XdmfIO::data_info<T> *data_struct, XdmfGrid *grid, XdmfArray **xarray,
                      int &array_count,
                      std::string filename, int compression)
{
    int narrays = (int)data_struct->data->GetDim();


    for(int i=0; i<narrays; i++){

        int array_number = i + array_count;

        /* Create an XdmfArray object to represent the data.  No values are actually moved out of the
         * pre-existing array (in this case the Array3D object populated above).  Deletion of the XdmfArray
         * will not affect the actual data, either.  The XdmfArray object only facilitates the interface
         * between the actual data and the data-output. */

        XdmfPointer data_pointer = data_struct->data->operator()(i);
        std::string arrayname = data_struct->data_name->operator()(i)->substr();


        XdmfInt32 number_type = XDMF_UNKNOWN_TYPE;
        if(typeid(T) == typeid(char)){
            number_type = XDMF_INT8_TYPE;
        }
        if(typeid(T) == typeid(unsigned char)){
            number_type = XDMF_UINT8_TYPE;
        }
        if(typeid(T) == typeid(short)){
            number_type = XDMF_INT16_TYPE;
        }
        if(typeid(T) == typeid(unsigned short)){
            number_type = XDMF_UINT16_TYPE;
        }
        if(typeid(T) == typeid(int)){
            number_type = XDMF_INT32_TYPE;
        }
        if(typeid(T) == typeid(unsigned int)){
            number_type = XDMF_UINT32_TYPE;
        }
        if(typeid(T) == typeid(float)){
            number_type = XDMF_FLOAT32_TYPE;
        }
        if(typeid(T) == typeid(double)){
            number_type = XDMF_FLOAT64_TYPE;
        }


        size_t rank = data_struct->dims->operator()(0)->GetDim();
        XdmfInt64 *shape = new XdmfInt64[rank];
        for(size_t r=0; r<rank; r++){
            shape[r] = (XdmfInt64)data_struct->dims->operator()(0)->operator()(rank-r-1);
        }

        xarray[array_number] = new XdmfArray;
        xarray[array_number]->SetAllowAllocate(false);                 /* Prevent the XdmfArray from allocating any memory. */
        xarray[array_number]->SetShape((XdmfInt32)rank, shape);        /* Set rank and dimension sizes of array. */
        xarray[array_number]->SetDataPointer(data_pointer);            /* Set pointer to previously-existing data. */
        xarray[array_number]->SetNumberType(number_type);              /* Identify the data type. */




        XdmfAttribute *attrib = new XdmfAttribute;              /* Create attribute to describe data. */
        attrib->SetName(arrayname.c_str());                     /* Set descriptive name for data. */
        attrib->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE); /* Define data to exist at nodes (rather than within
                                                                 * cells or on cell faces). */
        if(rank <= 3){
            attrib->SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);   /* Define data to be scalar. */
        } else {
            attrib->SetAttributeType(XDMF_ATTRIBUTE_TYPE_VECTOR);   /* Define data to be vector. */
        }
        attrib->SetValues(xarray[array_number]);                /* Define data values to be those stored in the
                                                                 * XdmfArray object created above (which in-turn takes
                                                                 * us to the c-array represented by the XdmfArray object. */

        /* Build name of file/path to be used with HDF5 heavy-data file.  This is set so that the
         * HDF5 filename matches that of the XML file.  The array name is used to identify the data within the
         * HDF5 file.  All arrays will be placed in the same HDF5 file. */
        std::string heavydataname(filename);
        heavydataname = heavydataname + ".h5:/" + arrayname;

        xarray[array_number]->SetHeavyDataSetName(heavydataname.c_str());  /* Set name used for heavy-data.  This will cause the
                                                              * actual data values to be written to the specified
                                                              * HDF5 file/array.  If the array is small-enough,
                                                              * values will be directly-inserted into the XML.  The
                                                              * default value for this threshold is 100, but can be
                                                              * modified by the user.  If data values are located
                                                              * within the XML file, this function call has no
                                                              * effect. */

        xarray[array_number]->SetCompression(compression);                 /* Set compression level to be used by HDF5 file.  This
                                                              * is also ignored if data is embedded in the XML file.
                                                              * Compression method will be noted in the HDF5 file and
                                                              * decompression will happen automatically when the
                                                              * HDF5 file is opened in a suitable viewer or these
                                                              * Xdmf____ classes. */

        /* Insert attribute into grid.  Note that this grid is a child of the "super-grid". */
        grid->Insert(attrib);

        /* Set pointer to data to NULL.  Actual data is not affected.  This prevents attempted double-deletion of
         * the data when both the PArray1D object (pointers to data) and Array object (the actual data) are destroyed. */
        data_struct->data->operator()(i) = NULL;

    } /* Loop through arrays to be added. */

    /* Increment array offset. */
    array_count += narrays;

} /* XdmfIO::addUniformArraysExistingGrid() */







/**
 * @brief Write arrays to disk as a single uniform grid.
 *
 *  This function is written to support several datatypes.  For each datatype, data is expected to be contained within
 *  an Array1D object.  Pointers to these arrays are stored in PArray1D objects to support the use of multiple arrays.
 *  For multi-dimensional arrays, data is expected to be serialized in row-major format (i.e., "C-style").
 *  The actual dimensions for each array are specified in another Array1D object.  Indices are expected to be in "IJK"
 *  order, meaning that the fastest-varying index comes first.  The origin and spacing of the array are also specified
 *  using Array1D objects.
 *
 *  If not data of a given type is to be written, pass a NULL pointer for the XdmfIO::data_info struct.
 *
 * @param filename Name if file to be written, without extension.  This name will be used for both the XDMF file as well as the
 *  HDF5 file used to store heavy data.
 * @param char_data XdmfIO::data_info struct for char data.  This is considered synonomous with XDMF 8-bit integers.
 * @param uchar_data XdmfIO::data_info struct for char data.  This is considered synonomous with XDMF 8-bit
 *  unsigned integers.
 * @param short_data XdmfIO::data_info struct for short data.  This is considered synonomous with XDMF 16-bit integers.
 * @param ushort_data XdmfIO::data_info struct for unsigned short data.  This is considered synonomous with XDMF 16-bit
 *  unsigned integers.
 * @param int_data XdmfIO::data_info struct for integer data.  This is considered synonomous with XDMF 32-bit integers.
 * @param uint_data XdmfIO::data_info struct for unsigned integer data.  This is considered synonomous with XDMF 32-bit
 *  unsigned integers.
 * @param long_data XdmfIO::data_info struct for long integer data.  This is considered synonomous with XDMF 64-bit integers.
 * @param float_data XdmfIO::data_info struct for float data.  This is considered synonomous with XDMF 32-bit floats.
 * @param double_data XdmfIO::data_info struct for double data.  This is considered synonomous with XDMF 64-bit floats.
 * @param compression Level of compression to be used.  Value must be in the range [0,9], with 0 being no compression and
 *  9 being maximum compression.
 * @warning It is expected that only a single grid is to be written.  Multiple data arrays can correspond to the grid.
 */
inline
void writeSingleUniformGrid(std::string filename,
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
{
    /* Determine the total number of arrays to be written.  This only refers to the data arrays and
     * does not include the origin/spacing data for each array.
     *
     * Also, get the information needed to define the grid.  If a multiple datatypes are specified, and
     * each datatype defines the grid, the last grid definition will be used.  Only a single grid
     * definition is allowed, and it is left to the programmer to ensure that a valid grid definition
     * is passed-in.  If conflicting grid definitions are passed-in, the last-defined grid will be use.  There
     * are no checks for duplicate/differing grids. */
    int narrays = 0;
    Array1D<size_t> *data_dims = NULL;
    Array1D<float> *data_origin = NULL;
    Array1D<float> *data_spacing = NULL;

    if(char_data){
        narrays += (int)char_data->data->GetDim();

        if(char_data->dims){
            data_dims = char_data->dims->operator()(0);
        }
        if(char_data->origin){
            data_origin = char_data->origin->operator()(0);
        }
        if(char_data->spacing){
            data_spacing = char_data->spacing->operator()(0);
        }
    }
    if(uchar_data){
        narrays += (int)uchar_data->data->GetDim();

        if(char_data->dims){
            data_dims = char_data->dims->operator()(0);
        }
        if(char_data->origin){
            data_origin = char_data->origin->operator()(0);
        }
        if(char_data->spacing){
            data_spacing = char_data->spacing->operator()(0);
        }
    }
    if(short_data){
        narrays += (int)short_data->data->GetDim();

        if(short_data->dims){
            data_dims = short_data->dims->operator()(0);
        }
        if(short_data->origin){
            data_origin = short_data->origin->operator()(0);
        }
        if(short_data->spacing){
            data_spacing = short_data->spacing->operator()(0);
        }
    }
    if(ushort_data){
        narrays += (int)ushort_data->data->GetDim();

        if(ushort_data->dims){
            data_dims = ushort_data->dims->operator()(0);
        }
        if(ushort_data->origin){
            data_origin = ushort_data->origin->operator()(0);
        }
        if(ushort_data->spacing){
            data_spacing = ushort_data->spacing->operator()(0);
        }
    }
    if(int_data){
        narrays += (int)int_data->data->GetDim();

        if(int_data->dims){
            data_dims = int_data->dims->operator()(0);
        }
        if(int_data->origin){
            data_origin = int_data->origin->operator()(0);
        }
        if(int_data->spacing){
            data_spacing = int_data->spacing->operator()(0);
        }
    }
    if(uint_data){
        narrays += (int)uint_data->data->GetDim();

        if(uint_data->dims){
            data_dims = uint_data->dims->operator()(0);
        }
        if(uint_data->origin){
            data_origin = uint_data->origin->operator()(0);
        }
        if(uint_data->spacing){
            data_spacing = uint_data->spacing->operator()(0);
        }
    }
    if(long_data){
        narrays += (int)long_data->data->GetDim();

        if(long_data->dims){
            data_dims = long_data->dims->operator()(0);
        }
        if(long_data->origin){
            data_origin = long_data->origin->operator()(0);
        }
        if(long_data->spacing){
            data_spacing = long_data->spacing->operator()(0);
        }
    }
    if(float_data){
        narrays += (int)float_data->data->GetDim();

        if(float_data->dims){
            data_dims = float_data->dims->operator()(0);
        }
        if(float_data->origin){
            data_origin = float_data->origin->operator()(0);
        }
        if(float_data->spacing){
            data_spacing = float_data->spacing->operator()(0);
        }
    }
    if(double_data){
        narrays += (int)double_data->data->GetDim();

        if(double_data->dims){
            data_dims = double_data->dims->operator()(0);
        }
        if(double_data->origin){
            data_origin = double_data->origin->operator()(0);
        }
        if(double_data->spacing){
            data_spacing = double_data->spacing->operator()(0);
        }
    }


    /* Check that grid information has been defined.  If origin/spacing are not specified,
     * default values are assumed below.  Dimensionality is required, though. */
    if(!data_dims){
        std::cerr << "XdmfIO::writeSingleUniformGrid()  ERROR: Incomplete grid specified!" << std::endl;
        std::cerr << "                                         Aborting data-writing." << std::endl;
        return;
    }




    size_t rank = data_dims->GetDim();
    XdmfInt64 *shape = new XdmfInt64[rank];
    for(size_t r=0; r<rank; r++){
        shape[r] = (XdmfInt64)data_dims->operator()(rank-r-1);
    }


    /* Create grid. */
    XdmfGrid *grid;
    grid = new XdmfGrid;                  /* Create a grid for the data. */
    grid->SetName("Structured Grid");     /* Name the grid.  This is user-choice. */
    grid->SetGridType(XDMF_GRID_UNIFORM); /* Set grid type.  This example is for image data, which is
                                           * on a uniform grid.  Grid choice can be other types based
                                           * on needs. */

    XdmfTopology *topo;
    topo = new XdmfTopology;                      /* Create a topology object to help define grid. */
    topo = grid->GetTopology();                   /* Set to be grid topology. */

    topo->SetTopologyType(XDMF_3DCORECTMESH);     /* Set grid type, again this is for image data and other
                                                   * options are possible for other applications. */
    topo->GetShapeDesc()->SetShape((XdmfInt32)rank, shape);  /* Set shape to match that of the data array. */


    XdmfGeometry *geo;
    geo = new XdmfGeometry;                   /* Create a geometry object to complete grid
                                               * structure definition. */
    geo = grid->GetGeometry();                /* Set to be the grid geometry. */
    if(data_origin && data_spacing){
        /* If both data_origin and data_spacing is defined, use those values when defining the geometry. */

        XdmfFloat64 *origin = new XdmfFloat64[rank];
        for(size_t r=0; r<rank; r++){
            origin[r] = data_origin->operator()(rank-r-1);
        }

        XdmfFloat64 *spacing = new XdmfFloat64[rank];
        for(size_t r=0; r<rank; r++){
            spacing[r] = data_spacing->operator()(rank-r-1);
        }

        geo->SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ);   /* Set for image data.  Again, more-complex options exist. */
        geo->SetOrigin(origin);                              /* Set origin. */
        geo->SetDxDyDz(spacing);                             /* Set point spacing. */
    } else {
        /* If origin/spacing information is missing, assume values.  Reading the produced file into ParaView
         * fails if origin/spacing information is not provided. */

        XdmfFloat64 origin[3] = { 0.0, 0.0, 0.0 };
        XdmfFloat64 spacing[3] = { 1.0, 1.0, 1.0 };

        geo->SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ);   /* Set for image data.  Again, more-complex options exist. */
        geo->SetOrigin(origin);                              /* Set origin. */
        geo->SetDxDyDz(spacing);                             /* Set point spacing. */
    }


    XdmfDOM *d = new XdmfDOM;
    XdmfRoot *root = new XdmfRoot;
    XdmfDomain *domain = new XdmfDomain;
    XdmfInformation *info = new XdmfInformation;

    XdmfArray **xarray = new XdmfArray*[narrays];


    /* Perform "boiler plate" operations for the output file. */
    root->SetDOM(d);
    root->SetVersion(2.2);
    root->Build();

    info->SetName("InfoName");
    info->SetValue("1.23");
    root->Insert(info);

    root->Insert(domain);

    domain->Insert(grid);


    /* For each input data structure, write array data to disk.  If passed-in pointer is NULL, no action is taken
     * for that datatype.  Domain is passed to the 'add' function rather than a tree grid so that the single
     * grid to be created is added directly to the domain. */
    int array_count = 0;
    if(char_data){
        addUniformArraysExistingGrid(char_data, grid, xarray, array_count, filename, compression);
    }
    if(uchar_data){
        addUniformArraysExistingGrid(uchar_data, grid, xarray, array_count, filename, compression);
    }
    if(short_data){
        addUniformArraysExistingGrid(short_data, grid, xarray, array_count, filename, compression);
    }
    if(ushort_data){
        addUniformArraysExistingGrid(ushort_data, grid, xarray, array_count, filename, compression);
    }
    if(int_data){
        addUniformArraysExistingGrid(int_data, grid, xarray, array_count, filename, compression);
    }
    if(uint_data){
        addUniformArraysExistingGrid(uint_data, grid, xarray, array_count, filename, compression);
    }
    if(long_data){
        addUniformArraysExistingGrid(long_data, grid, xarray, array_count, filename, compression);
    }
    if(float_data){
        addUniformArraysExistingGrid(float_data, grid, xarray, array_count, filename, compression);
    }
    if(double_data){
        addUniformArraysExistingGrid(double_data, grid, xarray, array_count, filename, compression);
    }


    /* Now that all data has been inserted, build the XML document. */
    root->Build();

    /* Add file extension to user-specified file stem and write data to the disk.  This will cause the XML file
     * to be written as well as the HDF5 file, if necessary.  If multiple HDF5 files are specified above, each will
     * be written at this point. */
    filename = filename + ".xmf";
    d->Write(filename.c_str());


    /* Cleanup pointers. */
    for(int i=0; i<narrays; i++){
        delete xarray[i];
    }
    delete [] xarray;  /* Deletion of XmdfDOM object does not delete array, so 'xarray'
                        * must be manually deleted here. */
    delete d;       /* Deleting XmdfDOM object deletes its children, so other
                     * Xmdf____ objects do not need to be manually deleted here.
                     * This behavior was verified using Valgrind to locate memory
                     * leaks, and this implementation produced no errors.  Valgrind shows that there
                     * may still be lost memory, but this appears to be within the Xdmf____ classes
                     * and thus is not under the programmer's control (without digging deep into source
                     * code). */

} /* XdmfIO::writeSingleUniformGrid() */




/**
 * @brief Write arrays to disk as uniform grids.
 *
 *  This function is written to support several datatypes.  For each datatype, data is expected to be contained within
 *  an Array1D object.  Pointers to these arrays are stored in PArray1D objects to support the use of multiple arrays.
 *  For multi-dimensional arrays, data is expected to be serialized in row-major format (i.e., "C-style").
 *  The actual dimensions for each array are specified in another Array1D object.  Indices are expected to be in "IJK"
 *  order, meaning that the fastest-varying index comes first.  The origin and spacing of the array are also specified
 *  using Array1D objects.
 *
 *  If not data of a given type is to be written, pass a NULL pointer for the XdmfIO::data_info struct.
 *
 * @param filename Name if file to be written, without extension.  This name will be used for both the XDMF file as well as the
 *  HDF5 file used to store heavy data.
 * @param char_data XdmfIO::data_info struct for char data.  This is considered synonomous with XDMF 8-bit integers.
 * @param uchar_data XdmfIO::data_info struct for char data.  This is considered synonomous with XDMF 8-bit
 *  unsigned integers.
 * @param short_data XdmfIO::data_info struct for short data.  This is considered synonomous with XDMF 16-bit integers.
 * @param ushort_data XdmfIO::data_info struct for unsigned short data.  This is considered synonomous with XDMF 16-bit
 *  unsigned integers.
 * @param int_data XdmfIO::data_info struct for integer data.  This is considered synonomous with XDMF 32-bit integers.
 * @param uint_data XdmfIO::data_info struct for unsigned integer data.  This is considered synonomous with XDMF 32-bit
 *  unsigned integers.
 * @param long_data XdmfIO::data_info struct for long integer data.  This is considered synonomous with XDMF 64-bit integers.
 * @param float_data XdmfIO::data_info struct for float data.  This is considered synonomous with XDMF 32-bit floats.
 * @param double_data XdmfIO::data_info struct for double data.  This is considered synonomous with XDMF 64-bit floats.
 * @param compression Level of compression to be used.  Value must be in the range [0,9], with 0 being no compression and
 *  9 being maximum compression.
 */
inline
void writeUniformGrid(std::string filename,
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
{
    /* Determine the total number of arrays to be written.  This only refers to the data arrays and
     * does not include the origin/spacing data for each array. */
    int narrays = 0;
    if(char_data){
        narrays += (int)char_data->data->GetDim();
    }
    if(uchar_data){
        narrays += (int)uchar_data->data->GetDim();
    }
    if(short_data){
        narrays += (int)short_data->data->GetDim();
    }
    if(ushort_data){
        narrays += (int)ushort_data->data->GetDim();
    }
    if(int_data){
        narrays += (int)int_data->data->GetDim();
    }
    if(uint_data){
        narrays += (int)uint_data->data->GetDim();
    }
    if(long_data){
        narrays += (int)long_data->data->GetDim();
    }
    if(float_data){
        narrays += (int)float_data->data->GetDim();
    }
    if(double_data){
        narrays += (int)double_data->data->GetDim();
    }

    int ngrids = 0;
    if(char_data){
        ngrids += (int)char_data->dims->GetDim();
    }
    if(uchar_data){
        ngrids += (int)uchar_data->dims->GetDim();
    }
    if(short_data){
        ngrids += (int)short_data->dims->GetDim();
    }
    if(ushort_data){
        ngrids += (int)ushort_data->dims->GetDim();
    }
    if(int_data){
        ngrids += (int)int_data->dims->GetDim();
    }
    if(uint_data){
        ngrids += (int)uint_data->dims->GetDim();
    }
    if(long_data){
        ngrids += (int)long_data->dims->GetDim();
    }
    if(float_data){
        ngrids += (int)float_data->dims->GetDim();
    }
    if(double_data){
        ngrids += (int)double_data->dims->GetDim();
    }

    /* If only a single grid was defined, call the function to write this single grid.  This results in an XDMF file
     * which does not contain a hierarchal tree structure. */
    if(ngrids == 1){
        writeSingleUniformGrid(filename, char_data, uchar_data, short_data, ushort_data, int_data,
                               uint_data, long_data, float_data, double_data, compression);
        return;
    }



    XdmfDOM *d = new XdmfDOM;
    XdmfRoot *root = new XdmfRoot;
    XdmfDomain *domain = new XdmfDomain;
    XdmfInformation *info = new XdmfInformation;

    XdmfArray **xarray = new XdmfArray*[narrays];
    XdmfGrid **grid = new XdmfGrid*[narrays];
    XdmfTopology **topo = new XdmfTopology*[narrays];
    XdmfGeometry **geo = new XdmfGeometry*[narrays];


    /* Perform "boiler plate" operations for the output file. */
    root->SetDOM(d);
    root->SetVersion(2.2);
    root->Build();

    info->SetName("InfoName");
    info->SetValue("1.23");
    root->Insert(info);

    root->Insert(domain);



    /* Define a grid which is the sole member of the domain.  When multiple arrays are present, their grids
     * are children of this "super-grid".  If only one array is to be used, the "super-grid" can be omitted,
     * and then the data grid is inserted into the domain.  See writeXDMFFile() for a single-arrayexample. */
    XdmfGrid *tree_grid = new XdmfGrid;
    tree_grid->SetName("ParentGrid");
    tree_grid->SetGridType(XDMF_GRID_TREE);
    domain->Insert(tree_grid);



    /* For each input data structure, write array data to disk.  If passed-in pointer is NULL, no action is taken
     * for that datatype. */
    int array_count = 0;
    if(char_data){
        addUniformArrays(char_data, tree_grid, xarray, grid, topo, geo, array_count, filename, compression);
    }
    if(uchar_data){
        addUniformArrays(uchar_data, tree_grid, xarray, grid, topo, geo, array_count, filename, compression);
    }
    if(short_data){
        addUniformArrays(short_data, tree_grid, xarray, grid, topo, geo, array_count, filename, compression);
    }
    if(ushort_data){
        addUniformArrays(ushort_data, tree_grid, xarray, grid, topo, geo, array_count, filename, compression);
    }
    if(int_data){
        addUniformArrays(int_data, tree_grid, xarray, grid, topo, geo, array_count, filename, compression);
    }
    if(uint_data){
        addUniformArrays(uint_data, tree_grid, xarray, grid, topo, geo, array_count, filename, compression);
    }
    if(long_data){
        addUniformArrays(long_data, tree_grid, xarray, grid, topo, geo, array_count, filename, compression);
    }
    if(float_data){
        addUniformArrays(float_data, tree_grid, xarray, grid, topo, geo, array_count, filename, compression);
    }
    if(double_data){
        addUniformArrays(double_data, tree_grid, xarray, grid, topo, geo, array_count, filename, compression);
    }


    /* Now that all data has been inserted, build the XML document. */
    root->Build();

    /* Add file extension to user-specified file stem and write data to the disk.  This will cause the XML file
     * to be written as well as the HDF5 file, if necessary.  If multiple HDF5 files are specified above, each will
     * be written at this point. */
    filename = filename + ".xmf";
    d->Write(filename.c_str());


    /* Cleanup pointers. */
    for(int i=0; i<narrays; i++){
        delete xarray[i];
    }
    delete [] xarray;  /* Deletion of XmdfDOM object does not delete array, so 'xarray'
                        * must be manually deleted here. */
    delete d;       /* Deleting XmdfDOM object deletes its children, so other
                     * Xmdf____ objects do not need to be manually deleted here.
                     * This behavior was verified using Valgrind to locate memory
                     * leaks, and this implementation produced no errors.  Valgrind shows that there
                     * may still be lost memory, but this appears to be within the Xdmf____ classes
                     * and thus is not under the programmer's control (without digging deep into source
                     * code). */

} /* XdmfIO::writeUniformGrid() */



} /* namespace XdmfIO */



//#include "XdmfIO.cpp"

#endif // XDMFIO_H
