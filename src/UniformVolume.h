/**
 * @file UniformVolume.h
 * @author Robert Grandin
 * @date 15 Oct 2010
 * @brief Definition of UniformVolume class.
 *
 * @section Class Description & Notes
 *
 * This class represents a 3-dimensional volume discretized into uniformly-sized
 * voxels (volume elements).  Grid spacing is uniform along each direction, but
 * does not need to be uniform between dimensions (i.e., x-spacing and y-spacing
 * can be different, but along the x-direction the spacing is uniform).
 *
 * An arbitrary number of scalar and vector quantities can be represented at each
 * voxel.  The number of quantities can be specified when an instance of of this
 * is declared, as well as anytime after that.  Quantities can also be deleted.
 *
 * Spatial parameters for the minimum and maximum extents can be specified,
 * making it possible to relate voxel position and size with a physical region in
 * space.
 *
 * The Qt Framework is used to provide functionality with graphical user interfaces
 * (GUIs).  Member functions which may require significant time to execute will emit
 * signals containing their progress so that a GUI can be updated accordingly.
 * If this class is to be used in a non-graphical application, two compilation
 * options exist.  In either case, no changes need to be made to this class.
 *  1 If the program is compiled with the Qt libraries, no special action needs
 *    to be take with respect to this class.  The normal Qt compilation proceedure
 *    should work.
 *  1 If Qt libraries are not available, or are not to be used, define the symbol
 *    NOQT with the compiler to disable Qt functionality for the entire application
 *    (and by extension, this class).
 *
 * Qt-specific functionality is handled by the parent class of this object, QtIntermediary.
 * QtIntermediary handles the details related to Qt vs. non-Qt compilation.  This class
 * does not require modification when either enabling or disabling Qt.
 *
 * By default, bounds-checking and error-checking measures are enabled for this class.
 * These can be disabled by defining RELEASE at compile-time.  Disabling them is
 * recommended for performance reasons once the program has been debugged.
 *
 * Some functionality changes are required when compiling this on Windows with the MSVC
 * compiler and linking against the VTK libraries.  Any such changes are controlled by
 * the definition of the symbol 'WIN_MSVC'.
 *
 * Reading and writing volume data is handled using the VTK library.  The VTK shared library
 * files (e.g., *.so and *.dll) must be present when linking an application using this class.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 *
 * @section Revisions
 *
 * @date 15 October 2010
 *	- Creation date.
 *
 * @date 30 November 2010
 * 	- Changed from Volume3D to UniformVolume3D.  No change to the functionality,
 * 	  just a name-change to better-reflect the purpose of this class.
 * @date 7 April 2011
 * 	- Separated into Core and GUI versions.
 * @date 27 January 2012
 *  - Code cleanup and switch to using Qt
 * @date October 2012
 *  - Added support for binary files in legacy VTK format.
 * @date 19 February 2013
 *  - Rename from "UniformVolume3D" to "UniformVolume".  Inclusion of "3D" seemed redundant
 *    since "volume" implies 3 dimensions.
 * @date March 2013
 *  - Implement VTK libraries for reading and writing files.
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2010-2013, Robert Grandin
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


#ifndef  UniformVolume_
#define  UniformVolume_


/* Non-Qt includes, standard C++ */
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string.h>

/* Non-Qt includes, my classes */
#include <Array1D.h>
#include <Array3D.h>
#include <Array4D.h>
#include <PArray1D.h>
#include <StringManip.h>

/* VTK classes for file IO */
#include <vtkVersion.h>
#include <vtkSmartPointer.h>

#include <vtkStructuredPoints.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPImageDataWriter.h>
#include <vtkImageImport.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPImageDataReader.h>
#include <vtkImageExport.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkXdmfWriter.h>
#include <vtkXdmfReader.h>


/* Include the base class which serves as an intermediary between this class
  and QObject.  It requires NOQT or USEQT to be defined for the entire application.
 */
#include <QtIntermediaryBase.h>





/**
 * @brief Representation of 3-dimensional volumetric data.
 * @warning C++11 features, such as move-constructor and move-assignment, require the symbol
 *  "CXX11" to be defined.
 */
template <class T>
class UniformVolume : public QtIntermediaryBase
{

public:
    /**
      @brief Default constructor.
      @pre Sufficient memory exists for object.
      @post UniformVolume object created.  Default values are
            - Spatial extents: -1.0 to 1.0 in each direction
            - 1 scalar quantity, 0 vector quantities
      @return None.
     */
    UniformVolume();


    /**
    @brief Copy constructor.
    @pre Sufficient memory exists for object.
    @param uc3d Reference to UniformVolume object which is to be copied.
    @post Current object updated to reflect all data and parameters of
        referenced object.
    @return None.
    */
    UniformVolume(const UniformVolume<T> &uc3d);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing UniformVolume object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    UniformVolume(UniformVolume<T> &&a);
#endif


    /**
      @brief Destructor.
      @pre UniformVolume object exists.
      @post Object destroyed.  All object data deleted.
      @return None.
     */
    virtual ~UniformVolume();


    /**
     * @brief Copy-assignment operator.
     * @param a UniformVolume object being assigned.
     * @return Reference to instance of UniformVolume.
     */
    UniformVolume& operator=(UniformVolume<T> a);


    /**
    @brief Retrieve the number of scalar quantities present in the object.
    @pre UniformVolume object exists.
    @post No changes to object.
    @return Number of scalar quantities in the object.
    */
    size_t NumScalarQuantities() const;


    /**
    @brief Retrieve the number of vector quantities present in the object.
    @pre UniformVolume object exists.
    @post No changes to object.
    @return Number of vector quantities in the object.
    */
    size_t NumVectorQuantities() const;


    /**
    @brief Add an array to hold a scalar quantity at each grid point.
    @pre UniformVolume object exists.
    @param name Name of scalar quantity added.
    @post New scalar quantity has been added.
    @return None.
    */
    void AddScalarQuantity(const std::string name);


    /**
    @brief Add an array to hold a scalar quantity at each grid point.
    @pre UniformVolume object exists.
    @param name Name of scalar quantity added.
    @param data Pointer to data.
    @post New scalar quantity has been added.
    @return None.
    */
    void AddScalarQuantity(const std::string name, Array3D<T> *data);


    /**
    @brief Remove array containing scalar quantity at each grid point.

    This removes the reference to the data array identified by quantity index 'qty'.
    @pre UniformVolume object exists and scalar quantity exists.
    @param qty 0-based index of quantity to be removed.
    @post Scalar quantity removed.  Remaining quantities remain in their
          same relative order, but their index numbers are updated to
          remain contiguous.
    @return None.
    @warning All data stored in the quantity will be lost.
    @warning Be aware of the index renumbering.  Order remains the same, but
          index numbers are updated to remain contiguous.
    @warning No data is actually deleted.  The pointer to the data array is simply
        forgotten and the original data array is untouched.
     */
    void RemoveScalarQuantityRef(const size_t qty);


    /**
    @brief Remove array containing scalar quantity at each grid point.
    @pre UniformVolume object exists and scalar quantity exists.
    @param qty 0-based index of quantity to be removed.
    @post Scalar quantity removed.  Remaining quantities remain in their
          same relative order, but their index numbers are updated to
          remain contiguous.
    @return None.
    @warning All data stored in the quantity will be lost.
    @warning Be aware of the index renumbering.  Order remains the same, but
          index numbers are updated to remain contiguous.
    */
    void RemoveScalarQuantity(const size_t qty);


    /**
     * @brief PointerToScalarData retrieves the pointer to the scalar data identified
     *  by index 'qty'.
     * @param qty 0-based index of scalar quantity to which a pointer is desired.
     * @return Pointer to Array3D object containing data.
     * @warning The data in this object is not forgotten by this object when this function
     *  is called.  It is left to the programmer to ensure that the memory occupied by the
     *  Array3D object containing data is only freed once.
     * @see RemoveScalarQuantityRef() for information about having this UniformVolume object
     *  "forget" about the data array.
     */
    Array3D<T>* PointerToScalarData(const size_t qty);


    /**
    @brief Set the name of a scalar quantity.
    @pre UniformVolume object exists and scalar quantity exists.
    @param qty Index of scalar quantity to be named.
    @param name Name to use with scalar quantity.
    @post Scalar quantity name updated.
    @return None.
    */
    void setScalarQuantityName(const size_t qty, const std::string name);


    /**
    @brief Get the name of a scalar quantity.
    @pre UniformVolume object exists and scalar quantity exists.
    @param qty Index of scalar quantity.
    @post No changes to object.
    @return Name of scalar quantity.
    */
    std::string ScalarQuantityName(const size_t qty);


    /**
    @brief Add an array to hold a vector quantity at each grid point.
    @pre UniformVolume object exists.
    @param name Name of vector quantity added.
    @param ncomp Number of components in the vector.
    @post New vector quantity has been added.
    @return None.
    */
    void AddVectorQuantity(const std::string name, const size_t ncomp);



    /**
    @brief Remove array containing vector quantity at each grid point.
    @pre UniformVolume object exists and vector quantity exists.
    @param qty 0-based index of quantity to be removed.
    @post Vector quantity removed.  Remaining quantities remain in their
          same relative order, but their index numbers are updated to
          remain contiguous.
    @return None.
    @warning All data stored in the quantity will be lost.
    @warning Be aware of the index renumbering.  Order remains the same, but
          index numbers are updated to remain contiguous.
    */
    void RemoveVectorQuantity(const size_t qty);


    /**
    @brief Set the name of a vector quantity.
    @pre UniformVolume object exists and vector quantity exists.
    @param qty Index of vector quantity to be named.
    @param name Name to use with vector quantity.
    @post Vector quantity name updated.
    @return None.
    */
    void setVectorQuantityName(const size_t qty, const std::string name);


    /**
    @brief Get the name of a vector quantity.
    @pre UniformVolume object exists and vector quantity exists.
    @param qty Index of vector quantity.
    @post No changes to object.
    @return Name of vector quantity.
    */
    std::string VectorQuantityName(const size_t qty);


    /**
     * @brief RemoveAllData deletes all data arrays (scalar and vector).
     */
    void RemoveAllData();


    /**
    * @brief Overload () operator to access scalar quantities.
    * @pre UniformVolume object exists.
    * @param ind1 Value of first index.
    * @param ind2 Value of second index.
    * @param ind3 Value of third index.
    * @param qty Quantity to be accessed.
    * @post No changes to object.
    * @return Value stored at supplied indices.
    */
    T& operator()(const size_t ind1, const size_t ind2, const size_t ind3, const size_t qty);


    /**
    * @brief Overload () operator to access scalar quantities.
    * @pre UniformVolume object exists.
    * @param ind1 Value of first index.
    * @param ind2 Value of second index.
    * @param ind3 Value of third index.
    * @param qty Quantity to be accessed.
    * @post No changes to object.
    * @return Value stored at supplied indices.
    */
    const T& operator()(const size_t ind1, const size_t ind2, const size_t ind3, const size_t qty) const;


    /**
    * @brief Overload () operator to access vector quantities.
    * @pre UniformVolume object exists.
    * @param ind1 Value of first index.
    * @param ind2 Value of second index.
    * @param ind3 Value of third index.
    * @param comp Component of vector.
    * @param qty Quantity to be accessed.
    * @post No changes to object.
    * @return Value stored at supplied indices.
    */
    T& operator()(const size_t ind1, const size_t ind2, const size_t ind3,
                  const size_t comp, const size_t qty);


    /**
    * @brief Overload () operator to access vector quantities.
    * @pre UniformVolume object exists.
    * @param ind1 Value of first index.
    * @param ind2 Value of second index.
    * @param ind3 Value of third index.
    * @param comp Component of vector.
    * @param qty Quantity to be accessed.
    * @post No changes to object.
    * @return Value stored at supplied indices.
    */
    const T& operator()(const size_t ind1, const size_t ind2, const size_t ind3,
                        const size_t comp, const size_t qty) const;


    /**
    @brief Estmiate the memory used by the volume.

    This only
        accounts for the data arrays (scalar and vector quantities)
        and not individual member variables.
    @pre UniformVolume object exists.
    @post No changes to object.
    @return Number of bytes used by this object.  Type 'double' is used
        to avoid possible overflow issues.
    @warning The number of bytes is calculated using the 'sizeof()'
        function and it is assumed that a single 'char' variable
        requires 1 byte.  If 'char' requires a different amount, the
        value returned by this function will need to be scaled
        accordingly.
    */
    double EstimateMemoryUsage() const;


    /**
    * @brief Reset values for a scalar data quantity.
    * @pre Volume object exists.
    * @param qty Index of quantity to be reset.
    * @param initvalue Value to be placed at index locations.
    * @return None.
    * @post Scalar values reset.
    */
    void ResetScalarVals(const int qty, const T initvalue);


    /**
    * @brief Reset values for a vector data quantity.  All components of
    *    vector quantity are reset.
    * @pre Volume object exists.
    * @param qty Index of quantity to be reset.
    * @param initvalue Value to be placed at index locations.
    * @return None.
    * @post Vector values reset.
    */
    void ResetVectorVals(const int qty, const T initvalue);


    /**
    @brief Provide the number components in the specified vector quantity.
    @pre UniformVolume object exists.
    @param qty Index of quantity to be queried.
    @post No changes to object.
    @return Number of components in vector.
    */
    int VectorQuantityComponents(const int qty) const;


    /**
    @brief Set the specified spatial extent for the volume.
    @pre UniformVolume object exists.
    @param dir Direction for which extent is being set.
        - 1: x-direction
        - 0: y-direction
        - 2: z-direction
    @param maxmin Identify if maximum or minimum value is being set.
        - positive values: maximum value being set
        - negative values: minimum value being set
    @param val Value of spatial extent to be set.
    @post Spatial extent updated.
    @return None.
    */
    void setSpatialExtent(const int dir, const int maxmin, const T val);


    /**
    @brief Retrieved the specified spatial extent for the volume.
    @pre UniformVolume object exists.
    @param dir Direction for which extent is being retrieved.
        - 1: x-direction
        - 0: y-direction
        - 2: z-direction
    @param maxmin Identify if maximum or minimum value is being retrieved.
        - positive values: maximum value being set
        - negative values: minimum value being set
    @post No change to object.
    @return Spatial extent specified.
    */
    T SpatialExtent(const int dir, const int maxmin) const;


    /**
    @brief Get the spacing between points in the specified direction.
    @pre UniformVolume object exists.
    @param dir Direction for which spacing is to be retrieved.
        - 0: x-direction
        - 1: y-direction
        - 2: z-direction
    @post No change to object.
    @return Spacing in the desired direction.
    */
    T PointSpacing(const int dir) const;


    /**
    @brief Set if data should be written to disk using the VTK Image Data
        format.

    This has no effect on VOL-format output.  This does not
        write the data to the disk.  This has no effect on which quantities
        will be written to the disk.
    @pre UniformVolume object exists.
    @param state Desired setting for Image Data output.  If set to true, all
        other VTK output formats are disabled.
    @post Setting updated.
    @return None.
    */
    void setVTKImageDataOutput(const bool state);


    /**
    @brief Set if data should be written to disk using the VTK Rectilinear Grid
        format.

    This has no effect on VOL-format output.  This does not
        write the data to the disk.  This has no effect on which quantities
        will be written to the disk.
    @pre UniformVolume object exists.
    @param state Desired setting for Rectilinear Grid output.  If set to true, all
        other VTK output formats are disabled.
    @post Setting updated.
    @return None.
    */
    void setVTKRectDataOutput(const bool state);


    /**
    @brief Return if VTK Image Data output is enabled/disabled.
    @pre UniformVolume object exists.
    @post No changes to object.
    @return Boolean representing if VTK Image Data output is enabled (true)
        or disabled (false).
    */
    bool VTKImageDataOutput() const;


    /**
    @brief Return if VTK Rectilinear Grid output is enabled/disabled.
    @pre UniformVolume object exists.
    @post No changes to object.
    @return Boolean representing if VTK Rectilinear Grid output is enabled (true)
        or disabled (false).
    */
    bool VTKRectDataOutput() const;


    /**
     * @brief Set if XDMF/HDF5 output format is to be used.
     * @param state True: Use XDMF/HDF5 format.  False: do not use XDMF/HDF5 format.
     */
    void setXDMFOutput(const bool state);


    /**
     * @brief Return if XDMF/HDF5 output format is to be used.
     * @return True: Use XDMF/HDF5 format.  False: do not use XDMF/HDF5 format.
     */
    bool XDMFOutput() const;


    /**
    @brief Set the stem of the filename to use when volume data is written to
        the disk.

    This should not include any directory information.  The
        appropriate file extension will be added by the file-writint functions.
    @pre UniformVolume object exists.
    @param stem String to be used for the filename stem.
    @post Stem value set.
    @return None.
    */
    void setFilenameStem(const std::string stem);


    /**
    @brief Retrieve the stem to be used for output filename.
    @pre UniformVolume object exists.
    @post No changes to object.
    @return String to be used for filename stem.
    */
    std::string FilenameStem() const;


    /**
    @brief Set the directory to which files are to be written.
    @pre UniformVolume object exists.
    @param dir Directory to be used for file output.
    @post Directory value set.
    @return None.
    */
    void setOutputDirectory(const std::string dir);


    /**
    @brief Retrieve the stem to be used for output files.
    @pre UniformVolume object exists.
    @post No changes to object.
    @return Directory to which output files will be written.
    */
    std::string OutputDirectory() const;


    /**
    @brief Reset the resolution of the volume.

    This will reset the
        resolution for all data quantities.  Vector quantities will maintian
        their current number of components.
    @pre UniformVolume object exists.
    @param nx Number of points in the x-direction.
    @param ny Number of points in the y-direction.
    @param nz Number of points in the z-direction.
    @param initval Value to which volume data is to be initialized.
    @post All data quantities resized and all values reset.
    @return None.
    @warning All data stored in this object will be lost.
    */
    void ResetResolution(const size_t ny, const size_t nx, const size_t nz, const T initval);


    /**
    @brief Get the volume resolution in the specified direction.
    @pre UniformVolume object exists.
    @param dir Direction in which resolution is desired.
        - 1: x-direction
        - 0: y-direction
        - 2: z-direction
    @post No changes to object.
    @return Number of points in specified direction.
    */
    size_t GetResolution(const int dir);


    /**
    @brief Write the volume data to the disk as a VTK file.

    File type
        is set by VTKImageDataOutput() and VTKRectDataOutput() functions.
        All data quantities in this volume will be written to the disk.
    @pre UniformVolume object exists and sufficient disk space exists
        for the output file.
    @post No changes to object.  Data written to disk.
    @return None.
    */
    void VTKWrite();


    /**
     * @brief Write volume data to disk as VTK image data file.  Output file will use
     *  VTK's XML format and file writing is performed by the VTK libraries.
     */
    void VTKWriteImageData();


    /**
    @brief Write the volume data to the disk as a binary VTK file.

    File type is set by VTKImageDataOutput() and VTKRectDataOutput() functions.
        All data quantities in this volume will be written to the disk.
    @pre UniformVolume object exists and sufficient disk space exists
        for the output file.
    @post No changes to object.  Data written to disk.
    @return None.
    @ warning This routine flips the byte-order (endian-ness) during output.  This
        is done to produce a big-endian file on a little-endian machine.  ParaView
        assumes that legacy VTK format binary files are big-endian.
    */
    void VTKWriteBinaryBitFlip();


    /**
    @brief Write the volume data to the disk as a binary VTK file.

    File type is set by VTKImageDataOutput() and VTKRectDataOutput() functions.
        All data quantities in this volume will be written to the disk.
    @pre UniformVolume object exists and sufficient disk space exists
        for the output file.
    @post No changes to object.  Data written to disk.
    @return None.
    @ warning This routine does not reverse the byte-order (endian-ness).
    */
    void VTKWriteBinary();


    /**
     * @brief VTKWriteBigEndian writes the volume data to a VTK binary file using
     *      Big Endian byte ordering (required ordering for use with ParaView).
     *
     * The endian-ness of the system is checked at run-time when this function
     *      is called and the appropriate binary-output file routine (VTKWriteBinary()
     *      or VTKWriteBinaryBitFlip()) is called.
     */
    void VTKWriteBinaryBigEndian();


    /**
    @brief Write volume data to the disk as a CNDE VOL file.

     Only scalar quantities can be written in this fashion.

     The coordinate system shown below is assumed to be the basis for this volume object.  When
     writing the data as a CNDE VOL file, the following axis-switching is done.  Capital letters
     denote the axes shown below while lower-case letters denote their corresponding axis in the
     CNDE VOL output file (as-labeled in the CNDE 3D visualization program).
        - X becomes -y
        - Y becomes -z
        - Z becomes +x
     @image html UV3D_coordsys.png "Coordinate system used in QGenRecon"
     @image latex UV3D_coordsys.png "Coordinate system used in QGenRecon"
    @pre UniformVolume object exists and sufficient disk space exists
        for the output file.
    @param qty Identifier for scalar quantity to be output.
    @param allowScaling Specify if data may be scaled to span a minimum range of values
        if it does not do so already.  This is to combat a bug in the CNDE 3D Visualization
        application where data is cast-to-int when adjusting the lookup table.  A small
        range of data values will get mapped to too few grayscales on the display, causing
        the visualization to be ineffective.
    @param minRange Minimum range data is to span.  If allowScaling is false, this parameter
        has no effect.
    @post No changes to object.  Data written to disk.
    @return None.
    @warning Although IEEE 754 requires elementary arithmetic operations to be within
        0.5 <a href="http://en.wikipedia.org/wiki/Unit_in_the_last_place>ULP</a> of the
        mathematically-exact result, the scaling and subsequent un-scaling may result in
        a loss of precision.  Check results accordingly.
    */
    void VOLWrite(const int qty, const bool allowScaling, const T minRange);


    /**
    @brief Read data from disk into UniformVolume object.
    @pre UniformVolume object exists and sufficient memory exists to hold
        the read-in data.
    @param filename Name of file to be read in.  If file is located in a different
        directory than the working directory of the application, a path must
        be provided.  Providing an absolute path is recommended.
    @param isBigEndian Identifies byte-ordering of binary file.  Ignored if
           datafile is ASCII.
           - True: data file to be read has Big Endian byte-ordering.
           - False: data file to be read has Little Endian byte-ordering.
    @post File read-in from file.
    @return None.
    @warning All data within this object will be lost.  Data within this object
        will match that in the input file after read-in.
    */
    void ReadFile(const std::string filename, const bool isBigEndian);


    /**
     * @brief Specify an existing array into which the next scalar quantity is to be
     *  read.
     * @param array Array into which data is to be read.
     * @param array_size Number of points in array.
     */
    void setScalarDestinationArray(T* array, const size_t array_size);


    /**
     * @brief Get the number of points read during the last scalar quantity.
     * @return points_read Number of points read.
     */
    size_t NumberScalarPointsRead();


    /**
    @brief Calculate the arithmetic mean of the specified data quantity.
    @pre UniformVolume object exists.
    @param scalarqty Identifier for scalar quantity to be processed.  Ignored
        if a negative number.
    @param vectorqty Identifier for vector quantity to be processed.  Ignored
        if a negative number.
    @param vectorcomp Component of vector quantity to be processed.  If negative,
        the magnitude of the vector is used.
    @post No changes to object.
    @return Mean value of specified quantity.
    */
    T Mean(const int scalarqty, const int vectorqty, const int vectorcomp) const;


    /**
    @brief Calculate the median of the specified data quantity.
    @pre UniformVolume object exists.
    @param scalarqty Identifier for scalar quantity to be processed.  Ignored
        if a negative number.
    @param vectorqty Identifier for vector quantity to be processed.  Ignored
        if a negative number.
    @param vectorcomp Component of vector quantity to be processed.  If negative,
        the magnitude of the vector is used.
    @post No changes to object.
    @return Median value of specified quantity.
    */
    T Median(const int scalarqty, const int vectorqty, const int vectorcomp) const;


    /**
    @brief Calculate the variance of the specified data quantity.
    @pre UniformVolume object exists.
    @param scalarqty Identifier for scalar quantity to be processed.  Ignored
        if a negative number.
    @param vectorqty Identifier for vector quantity to be processed.  Ignored
        if a negative number.
    @param vectorcomp Component of vector quantity to be processed.  If negative,
        the magnitude of the vector is used.
    @post No changes to object.
    @return Variance value of specified quantity.
    */
    T Variance(const int scalarqty, const int vectorqty, const int vectorcomp) const;


    /**
    @brief Locate the minimum of the specified data quantity.
    @pre UniformVolume object exists.
    @param scalarqty Identifier for scalar quantity to be processed.  Ignored
        if a negative number.
    @param vectorqty Identifier for vector quantity to be processed.  Ignored
        if a negative number.
    @param vectorcomp Component of vector quantity to be processed.  If negative,
        the magnitude of the vector is used.
    @param loc_x X-direction index of minimum-valued point.  Set by this function.
    @param loc_y Y-direction index of minimum-valued point.  Set by this function.
    @param loc_z Z-direction index of minimum-valued point.  Set by this function.
    @post No changes to object.
    @return Minimum value of specified quantity.
    */
    T MinVal(const int scalarqty, const int vectorqty, const int vectorcomp,
             int &loc_x, int &loc_y, int &loc_z) const;


    /**
    @brief Locate the maximum of the specified data quantity.
    @pre UniformVolume object exists.
    @param scalarqty Identifier for scalar quantity to be processed.  Ignored
        if a negative number.
    @param vectorqty Identifier for vector quantity to be processed.  Ignored
        if a negative number.
    @param vectorcomp Component of vector quantity to be processed.  If negative,
        the magnitude of the vector is used.
    @param loc_x X-direction index of maximum-valued point.  Set by this function.
    @param loc_y Y-direction index of maximum-valued point.  Set by this function.
    @param loc_z Z-direction index of maximum-valued point.  Set by this function.
    @post No changes to object.
    @return Maximum value of specified quantity.
    */
    T MaxVal(const int scalarqty, const int vectorqty, const int vectorcomp,
             int &loc_x, int &loc_y, int &loc_z) const;


    /**
     * @brief MemoryRequired calculates the memory required by an instance of this class.
     * @return Memory required, in bytes, of an instance of this class.
     * @warning Units of bytes is contingent upon sizeof(char) = 1.  If type char is of a
     *      different size the value returned by this function will need to be adjusted
     *      accordingly.
     */
    size_t MemoryRequired() const;


    /** @brief QtIntermediaryBase object to provide signals for CUDA function wrappers. */
    QtIntermediaryBase *qtsignals;

    /** @brief Format-conversion options that must be set prior to calling PerformCalculations. */
    struct conv_opts{
        /** @brief File to be read-in. */
        std::string input_file;

        /** @brief Indicate if input file has Big Endian byte ordering. */
        bool isBigEndian;

        /** @brief Indicate if output should be CNDE VOL format. */
        bool output_CNDEVOL;

        /** @brief Indicate if output should be VTK Rectilinear Grid format. */
        bool output_VTKRectGrid;

        /** @brief Indicate if output should be VTK Image Data format. */
        bool output_VTKImageData;

        /** @brief Indicate if data is allowed to be scaled with CNDE VOL output. */
        bool output_CNDEVOL_allowScaling;

        /** @brief Indicate if output should be XDMF/HDF5 format. */
        bool output_xdmf;

        /** @brief Minimum data range to be used if data-scaling is allowed for CNDE VOL output. */
        T output_CNDEVOL_scalingRange;
    };

    /** @brief Instance of UniformVolume::conv_opts. */
    conv_opts conversion_options;

    /**
     * @brief PointCoordinates returns the (x,y,z) coordinates of the volume point (row,col,slice).
     * @param row 0-based row index of point.
     * @param col 0-based column index of point.
     * @param slice 0-based slice index of point.
     * @param y Y-coordinate of point spatial position.
     * @param x X-coordinate of point spatial position.
     * @param z Z-coordinate of point spatial position.
     */
    void PointCoordinates(const size_t row, const size_t col, const size_t slice, T &y, T &x, T &z) const;


    /**
     * @brief Write data to disk in XDMF/HDF5 format.
     * @param compression Level of compression to be used.  Value must be in the range 0 - 9, 0 is no
     *  compression.
     */
    void WriteXdmf(const int compression);


    /**
     * @brief Free a memory array which was originally allocated by the VTK libraries.
     * @param data_ptr Pointer to array.
     * @param npts1 Number of points in first dimension of array.
     * @param npts2 Number of points in second dimension of array.
     * @param npts3 Number of points in third dimension of array.
     */
    void FreeMemoryVTK(T* data_ptr, size_t npts1, size_t npts2, size_t npts3);


    /**
     * @brief Identify if data contained in this object was read-in using the VTK libraries.
     * @return True if data resides in memory allocated by VTK libraries, False if data resides in memory
     *  allocated by this object.
     */
    bool DataFromVTK() const;








#ifdef USEQT
public slots:
#endif
    /**
     * @brief PerformCalculations is a overload of QtIntermediary::PerformCalculations.
     *
     * In UniformVolume
     *      this function handles the reading-in of a volume file and the writing of that data in a
     *      user-specified format.  This is done to allow volume-conversion operations to happen in its
     *      own thread.
     */
    void PerformCalculations();


protected:
    /** @brief Minimum x-coordinate of volume */
    T xmin;

    /** @brief Maximum x-coordinate of volume */
    T xmax;

    /** @brief Minimum y-coordinate of volume */
    T ymin;

    /** @brief Maximum y-coordinate of volume */
    T ymax;

    /** @brief Minimum z-coordinate of volume */
    T zmin;

    /** @brief Maximum z-coordinate of volume */
    T zmax;

    /** @brief Spatial resolution in the x-direction */
    T xspacing;

    /** @brief Spatial resolution in the y-direction */
    T yspacing;

    /** @brief Spatial resolution in the z-direction */
    T zspacing;

    /** @brief Tracks if xmin has been set */
    bool xminset;

    /** @brief Tracks if xmax has been set */
    bool xmaxset;

    /** @brief Tracks if ymin has been set */
    bool yminset;

    /** @brief Tracks if ymax has been set */
    bool ymaxset;

    /** @brief Tracks if zmin has been set */
    bool zminset;

    /** @brief Tracks if zmax has been set */
    bool zmaxset;

    /** @brief Stem of filename to be used for output */
    std::string filenamestem;

    /** @brief Output directory */
    std::string outputdir;

    /** @brief Tracks if VTK output is in Image Data format */
    bool imageoutput;

    /** @brief Tracks if VTK output is in Rectilinear Grid format */
    bool rectoutput;

    /** @brief Tracks if XDMF/HDF5 output is to be used. */
    bool xdmfoutput;

    /** @brief Name of the volume.  This is the stem of any output files which
    * 			use this volume */
    std::string volname;

    /** @brief Tracks if volume name has been set. */
    bool volnameset;

    /** @brief Number of rows in the reconstruction volume */
    size_t vrows;

    /** @brief Number of columns in the reconstruction volume */
    size_t vcols;

    /** @brief Number of slices in the reconstruction volume */
    size_t vslices;

    /** @brief Units used for length quantities */
    std::string unitslength;

    /** @brief Tracks if the gradient magnitude is to be subtracted from the intensity */
    bool intensitysubtractgradient;

    /** @brief Tracks if the laplacian magnitude is to be subtracted from the intensity */
    bool intensitysubtractlaplacian;

    /** @brief Number of scalar quantities at each point */
    size_t nscalars;

    /** @brief Number of vector quantities at each point */
    size_t nvectors;

    /** @brief Array of pointers to scalar quantity arrays */
    PArray1D<Array3D<T>*> pscalars;

    /** @brief Array of pointers to vectory quantity arrays */
    PArray1D<Array4D<T>*> pvectors;

    /** @brief List of names for scalar quantities */
    PArray1D<std::string*> scalar_names;

    /** @brief List of names for vector quantities */
    PArray1D<std::string*> vector_names;

    /** @brief Size reserved for quantity names */
    int qtysize;

    /** @brief Human-readable form of datatype */
    std::string dtypename;

    /** @brief Identifies if the current binary file-writing is Big-Endian. */
    bool writebigendian;

    /** @brief Array to be used for reading-in scalar data. */
    T* scalar_data;

    /** @brief Number of points allowed for reading scalar data into pre-existing array. */
    size_t scalar_data_size;

    /** @brief Number of points actually read on most-recent scalar data read. */
    size_t scalar_data_points_read;

    /** @brief Flags if data was read-in using VTK library.  */
    bool data_from_vtk;



private:
    /**
    @brief Initialization of object.
    @pre UniformVolume object exists.
    @param nx Number of points in x-direction.
    @param ny Number of points in y-direction.
    @param nz Number of points in z-direction.
    @param minx Minimum x-coordinate.
    @param maxx Maximum x-coordinate.
    @param miny Minimum y-coordinate.
    @param maxy Maximum y-coordinate.
    @param minz Minimum z-coordinate.
    @param maxz Maximum z-coordinate.
    @param initval Initial value for all points.  Used for all scalar quantities and all
        components of all vector quantities.
    @param n_scalars Number of scalar quantities to create.
    @param n_vectors Number of vector quantities to create.
    @param qty_label_size Size to be reserved for scalar and vector labels.
    @post Object initialized.
    @return None.
    */
    void Initialize(const int nx, const int ny, const int nz,
                    const T minx, const T maxx, const T miny, const T maxy,
                    const T minz, const T maxz, const T initval,
                    const int n_scalars, const int n_vectors, const int qty_label_size);


    /**
    * @brief Read VOL file.
    * @pre UniformVolume object exists.
    * @param filename Name of file to be read.
    * @post VOL file contents read into object.
    * @return None.
    */
    void ReadVOLFile(std::string filename);


    /**
    * @brief Read VTK file.
    * @pre UniformVolume object exists.
    * @param filename Name of file to be read.
    * @param isBigEndian Identifies byte-ordering of binary file.  Ignored if
    *       datafile is ASCII.
    *       - True: data file to be read has Big Endian byte-ordering.
    *       - False: data file to be read has Little Endian byte-ordering.
    * @warning It is assumed that the scalar intensity data is the first quantity
    * 		represented in the VTK file.
    * @post VTK file contents read into object.
    * @return None.
    */
    void ReadVTKFile(std::string filename, const bool isBigEndian);


    /**
     * @brief Read XDMF file and place data into this object.
     * @param filename Name of file to be read.
     */
    void ReadXDMFFile(std::string filename);


    /**
     * @brief Read a legacy-format VTK file using the VTK library functions.
     * @param filename Name of file to be read.
     * @warning Only structured points data (i.e. image data) is supported.
     * @warning Only single-quantity (i.e., one scalar) files are supported.
     */
    void ReadLegacyVTKFile(std::string filename);


    /**
     * @brief Read a XML-format VTK file using the VTK library functions.
     * @param filename Name of file to be read.
     * @warning Only structured points data (i.e. image data) is supported.
     * @warning Only single-quantity (i.e., one scalar) files are supported.
     */
    void ReadXMLVTKFile(std::string filename);


    /**
    * @brief Function that writes the actual VOL file.
    *
    * This is not called directly by the user, but rather by VOLWrite().
    * @pre UniformVolumeCore object exists.
    * @param fileName Name of the output file, including extension.
    * @param xx Number of elements in the 2nd dimension.  Set by this function.
    * @param yy Number of elements in the 1st dimension.  Set by this function.
    * @param zz Number of elements in the 3rd dimension.  Set by this function.
    * @param voxelIntensity Pointer to serialized volume array.
    * @post Volume contents written to disk in the form of a VOL file.
    * @return Status code.
    * 		- 0: Success
    * 		- Any other value: Failure
    */
    int WriteVOLFile(std::string fileName, size_t xx, size_t yy, size_t zz,
                     T *voxelIntensity);


    /**
     * @brief ByteSwap reverses the byte-order of the input value.
     *
     * Code from http://stackoverflow.com/a/3824338.
     * @param value Pointer to value to have its bytes swapped.
     * @param numbytes Size of value to be swapped, in bytes.
     * @warning Value will be modified!
     */
    void ByteSwap(void *value, size_t numbytes);


    /**
     * @brief IsBigEndian checks if the system endian-ness is Big Endian.
     * @return 'true' if system is Big Endian.
     */
    bool IsBigEndian();


    /**
     * @brief VTKReadLegacy reads legacy VTK formatted files.
     * @param file std::fstream object of opened file.
     * @param isBigEndian Indiciates endian-ness of file.
     *      - True: File to be read is Big Endian.
     *      - False: File to be read is Little Endian.
     * @warning File is not closed by this function.
     * @warning It is expected that the first line read by this function contains the header keyword
     *      "DATASET", and the second line "DIMENSIONS".  Remainder of header is handled internal to
     *      this function.  This requires the first three lines to be read in the calling function (or
     *      farther up the hierarchy).
     */
    void VTKReadLegacyBinary(std::fstream &file, const bool isBigEndian);


    /**
     * @brief VTKReadLegacy reads legacy VTK formatted files.
     * @param file std::fstream object of opened file.
     * @warning It is expected that the first line read by this function contains the header keyword
     *      "DATASET", and the second line "DIMENSIONS".  Remainder of header is handled internal to
     *      this function.  This requires the first three lines to be read in the calling function (or
     *      farther up the hierarchy).
     */
    void VTKReadLegacyASCII(std::fstream &file);


    /**
    @brief Write a portion of the volume data to the disk as a binary VTK file.

    File type is set by VTKImageDataOutput() and VTKRectDataOutput() functions.
        All data quantities in this volume will be written to the disk.
    @param first_slice 0-based index of first slice to be written.
    @param num_slices Number of slices to be written.
    @param slice_min Z-coordinate value of first slice.
    @warning This routine does not reverse the byte-order (endian-ness).
    */
    void VTKWriteBinaryPartial(const size_t first_slice, const size_t num_slices, const float slice_min);


    /**
     * @brief VTKWriteBigEndian writes a portion of the volume data to a VTK binary file using
     *      Big Endian byte ordering (required ordering for use with ParaView).
     *
     * The endian-ness of the system is checked at run-time when this function
     *      is called and the appropriate binary-output file routine (VTKWriteBinary()
     *      or VTKWriteBinaryBitFlip()) is called.
     * @param first_slice 0-based index of first slice to be written.
     * @param num_slices Number of slices to be written.
     * @param slice_min Z-coordinate value of first slice.
     */
    void VTKWriteBinaryBitFlipPartial(const size_t first_slice, const size_t num_slices, const float slice_min);


    /**
     * @brief UniformVolumeSwap swaps member information between two UniformVolume objects.
     * @param first First UniformVolume object.
     * @param second Second UniformVolume object.
     */
    friend void UniformVolumeSwap(UniformVolume<T> &first, UniformVolume<T> &second)
    {
        std::swap(first.conversion_options, second.conversion_options);
        std::swap(first.xmin, second.xmin);
        std::swap(first.xmax, second.xmax);
        std::swap(first.ymin, second.ymin);
        std::swap(first.ymax, second.ymax);
        std::swap(first.zmin, second.zmin);
        std::swap(first.zmax, second.zmax);
        std::swap(first.xspacing, second.xspacing);
        std::swap(first.yspacing, second.yspacing);
        std::swap(first.zspacing, second.zspacing);
        std::swap(first.xminset, second.xminset);
        std::swap(first.xmaxset, second.xmaxset);
        std::swap(first.yminset, second.yminset);
        std::swap(first.ymaxset, second.ymaxset);
        std::swap(first.zminset, second.zminset);
        std::swap(first.zmaxset, second.zmaxset);
        std::swap(first.filenamestem, second.filenamestem);
        std::swap(first.outputdir, second.outputdir);
        std::swap(first.imageoutput, second.imageoutput);
        std::swap(first.rectoutput, second.rectoutput);
        std::swap(first.volname, second.volname);
        std::swap(first.volnameset, second.volnameset);
        std::swap(first.vrows, second.vrows);
        std::swap(first.vcols, second.vrows);
        std::swap(first.vslices, second.vslices);
        std::swap(first.unitslength, second.unitslength);
        std::swap(first.intensitysubtractgradient, second.intensitysubtractgradient);
        std::swap(first.intensitysubtractlaplacian, second.intensitysubtractlaplacian);
        std::swap(first.nscalars, second.nscalars);
        std::swap(first.nvectors, second.nvectors);
        std::swap(first.pscalars, second.pscalars);
        std::swap(first.pvectors, second.pvectors);
        std::swap(first.scalar_names, second.scalar_names);
        std::swap(first.vector_names, second.vector_names);
        std::swap(first.qtysize, second.qtysize);
        std::swap(first.dtypename, second.dtypename);
        std::swap(first.writebigendian, second.writebigendian);
        std::swap(first.scalar_data, second.scalar_data);
        std::swap(first.scalar_data_size, second.scalar_data_size);
        std::swap(first.scalar_data_points_read, second.scalar_data_points_read);
    }


    /**
     * @brief Load data contained within VTK Image Data dataset into this object's data structures.
     * @param dataset Pointer to vtkImageData object containing data to be loaded.
     * @param reader Pointer to vtkDataReader used to read the input file.
     */
    void LoadVTKDataset(vtkImageData *dataset, vtkAlgorithm *reader);

};

#include "UniformVolume.cpp"
#endif /* UniformVolume_ */
