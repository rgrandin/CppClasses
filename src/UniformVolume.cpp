/**
 * @file UniformVolume.cpp
 * @author Robert Grandin
 * @brief Implementation of UniformVolume class.
 */


#include "UniformVolume.h"    /* Included for syntax-hilighting */



/* =======================================================================================
 *
 *
 *                  PRIVATE FUNCTIONS
 *
 *
 *
 */
template <class T>
void UniformVolume<T>::Initialize(const int nx, const int ny, const int nz,
                                    const T minx, const T maxx, const T miny, const T maxy,
                                    const T minz, const T maxz, const T initval, const int n_scalars,
                                    const int n_vectors, const int qty_label_size)
{
    vrows = ny;
    vcols = nx;
    vslices = nz;
    xmin = minx;
    xmax = maxx;
    ymin = miny;
    ymax = maxy;
    zmin = minz;
    zmax = maxz;
    xminset = false;
    xmaxset = false;
    yminset = false;
    ymaxset = false;
    zminset = false;
    zmaxset = false;
    if(vcols > 1){
        xspacing = (xmax - xmin)/((T)vcols - 1.0e0+1);
    } else {
        xspacing = 0.0e0;
    }
    if(vrows > 1){
        yspacing = (ymax - ymin)/((T)vrows - 1.0e0+1);
    } else {
        yspacing = 0.0e0;
    }
    if(vslices > 1){
        zspacing = (zmax - zmin)/((T)vslices - 1.0e0+1);
    } else {
        zspacing = 0.0e0;
    }
    nscalars = 0;
    nvectors = 0;
    qtysize = qty_label_size;

    /* Remove all initial data. */
    UniformVolume<T>::RemoveAllData();

    /* Add specified quantities. */
    std::string scalarlabel("Scalar");
    std::stringstream ssnum;
    for(int i=0; i<n_scalars; i++){
        ssnum.str(""); ssnum << i;
        std::string tmpstr(scalarlabel + ssnum.str());
        UniformVolume<T>::AddScalarQuantity(tmpstr);
        pscalars(i)->ResetVal(initval);
    }

    std::string vectorlabel("Vector");
    for(int i=0; i<n_vectors; i++){
        ssnum.str(""); ssnum << i;
        std::string tmpstr(vectorlabel + ssnum.str());
        UniformVolume<T>::AddVectorQuantity(tmpstr,3);
        pvectors(i)->ResetVal(initval);
    }


    /* Set type name */
    dtypename = "unknown_type";
    if(typeid(T) == typeid(char)){
        dtypename = "char";
    }
    if(typeid(T) == typeid(short)){
        dtypename = "short";
    }
    if(typeid(T) == typeid(int)){
        dtypename = "int";
    }
    if(typeid(T) == typeid(long)){
        dtypename = "long";
    }
    if(typeid(T) == typeid(long long)){
        dtypename = "long long";
    }
    if(typeid(T) == typeid(float)){
        dtypename = "float";
    }
    if(typeid(T) == typeid(double)){
        dtypename = "double";
    }
    if(typeid(T) == typeid(long double)){
        dtypename = "long double";
    }

    qtsignals = new QtIntermediaryBase;

    writebigendian = false;

    xdmfoutput = false;

    conversion_options.input_file = "Not_Set";
    conversion_options.isBigEndian = false;
    conversion_options.output_CNDEVOL = false;
    conversion_options.output_VTKImageData = false;
    conversion_options.output_VTKRectGrid = false;
    conversion_options.output_CNDEVOL_allowScaling = false;
    conversion_options.output_CNDEVOL_scalingRange = (T)1024.0e0;
    conversion_options.output_xdmf = false;

    scalar_data = NULL;
    scalar_data_size = 0;
    scalar_data_points_read = 0;
}



/*
    NOTE: This function is taken straight from the code as-provided by Jia-Dong
        and Joe Gray.  Its calling routine, VOLWrite(), has been modified for
        use with this class.
 */
template <class T>
int UniformVolume<T>::WriteVOLFile(
    std::string	fileName,
    size_t			xx,
    size_t			yy,
    size_t			zz,
    T		* voxelIntensity)
{

    /*
     * Conversion between (X,Y,Z) axes used in my code and their corresponding axes in HRCT, PSCT, and
     * 3D Visualization:
     *  (my coords) --> (corresponding VOL coords)
     *            X --> y
     *            Y --> z
     *            Z --> x
     */


    size_t			len;
    size_t			headerLength;

    /* Get integer representation of resolution.  Convert between coordinate systems here. */
    int x = (int)zz;
    int y = (int)xx;
    int z = (int)yy;

    char		endSliceStr[100], countStr[30],startSliceStr[100];
    char		extraInfoStr[10];
    char		stateValidFlg;
    std::ofstream	ffp;

    T		xRotateAngle=0.0;
    T		yRotateAngle=0.0;

    T		xzoom=1.0;
    T		yzoom=1.0;
    T		zzoom=1.0;
    int			corex=1, corey=1;

    int			radius=0;
    char		* voxelState;
    int			pixelsPerVolume;

    /* Dummy code to remove unused variable warnings. */
    T tmp = voxelIntensity[0];
    tmp += (T)1.0e0;

    pixelsPerVolume = z*x*y;

    stateValidFlg='Y';
    sprintf(countStr,"%dx%dx%d",x, y, z);
    sprintf(endSliceStr, "Volume Size: %s\n", countStr);
    sprintf(extraInfoStr, " ");
    sprintf(startSliceStr, "Start Slice:\ne:\\output\\%s001.asc\n",filenamestem.c_str());

    headerLength = strlen(endSliceStr)+(size_t)1+strlen(startSliceStr)+(size_t)1+strlen(extraInfoStr)+(size_t)1;
    headerLength +=1;									//state part

    len = strlen(startSliceStr);
    len = strlen(endSliceStr);
    len = strlen(extraInfoStr);
    ffp.open(fileName.c_str(),std::ios::out|std::ios::binary);
    if(!ffp){
        std::string errorstring;
        errorstring = "Cannot open " + fileName + "\n";
        std::cerr << errorstring << std::endl;
        return 1;
    }

    int iheaderLength = (int)headerLength;
    ffp.write((char *)&iheaderLength,sizeof(int));
    ffp.write(startSliceStr,sizeof(char)* (strlen(startSliceStr)+(size_t)1));
    ffp.write(endSliceStr,sizeof(char)* (strlen(endSliceStr)+(size_t)1));
    ffp.write(extraInfoStr,sizeof(char)* (strlen(extraInfoStr)+(size_t)1));
    ffp.write(&stateValidFlg,sizeof(char));

    ffp.write((char *)&x, sizeof(int)*1);
    ffp.write((char *)&y, sizeof(int)*1);
    ffp.write((char *)&z, sizeof(int)*1);

    ffp.write((char *)&xRotateAngle,sizeof(T));
    ffp.write((char *)&yRotateAngle,sizeof(T));

    ffp.write((char *)&xzoom,sizeof(T));
    ffp.write((char *)&yzoom,sizeof(T));
    ffp.write((char *)&zzoom,sizeof(T));

    ffp.write((char*) &corex, sizeof(int));
    ffp.write((char*) &corey,sizeof(int));
    ffp.write((char*)&radius, sizeof(int));



    /* Typecast values individually into a 'float' variable.  Then, write to place in the
     * data file.  This is to avoid problems writing a 'float' VOL file
     * from a 'double' UniformVolume object. */
    std::string desc;
    desc = "Writing CNDE VOL File";
    qtsignals->EmitFunctionDesc(desc);
    qtsignals->EmitFunctionProgress(0.0e0);

    float ftmp = (float)0.0e0;
    for(size_t i=0; i<yy; i++){
        for(size_t j=0; j<xx; j++){
            for(size_t k=0; k<zz; k++){
                ftmp = (float)pscalars(0)->operator ()(i,j,k);
                ffp.write(reinterpret_cast<char*>(&ftmp), sizeof(float));
            }
        }
        qtsignals->EmitFunctionProgress((float)i/(float)yy);
    }

    desc = "";
    qtsignals->EmitFunctionDesc(desc);
    qtsignals->EmitFunctionProgress(0.0e0);




    voxelState = new char[pixelsPerVolume];
    for(int iii=0; iii<pixelsPerVolume; iii++)
        voxelState[iii] = 0;
    if(stateValidFlg=='Y')
        ffp.write((char *)voxelState,sizeof(char)*pixelsPerVolume);
    if (voxelState != NULL) {delete [] voxelState; voxelState=NULL;}
    ffp.close();

    /* Dummy code to remove unused variable warning. */
    len++;

    return(0);
}


template <class T>
void UniformVolume<T>::ReadVOLFile(std::string filename)
{
    /*
     * Conversion between (X,Y,Z) axes used in my code and their corresponding axes in HRCT, PSCT, and
     * 3D Visualization:
     *  (my coords) --> (corresponding VOL coords)
     *            X --> y
     *            Y --> z
     *            Z --> x
     */


    /* Delete all data and create a single scalar quantity. */
    UniformVolume<T>::RemoveAllData();
    if(scalar_data == NULL){
        UniformVolume<T>::AddScalarQuantity("cnde_vol_data");
    } else {
        scalar_names(0) = new std::string;
        scalar_names(0)->reserve(qtysize);
        scalar_names(0)->assign("cnde_vol_data");
        nscalars = 1;
    }

    /*
     * NOTE: This routine adapted from code provided by Jia-Dong
     */

    char * voxelState = NULL;
    float xRotateAngle;
    float yRotateAngle;
    unsigned long pixelsPerVolume;
    float xzoom;
    float yzoom;
    float zzoom;
    int corex;
    int corey;
    int radius;
    int x;
    int y;
    int z;

    std::ifstream	fp;
    char		*headerInfo = NULL;
    char		*startSlice;
    char		*endSlice;
    char		*extraInfo;
    char		stateflg;

    int headerLength;

    char *extraInfoStr;
    char *endSliceStr;
    char *startSliceStr;

    fp.open(filename.c_str(),std::ios::in|std::ios::binary);

    fp.read((char *)&headerLength,sizeof(int));

    headerInfo = (char *)new char[headerLength];
    fp.read(headerInfo,sizeof(char)*headerLength);


    startSlice = headerInfo;
    endSlice = startSlice+strlen(startSlice)+1;
    extraInfo = endSlice+strlen(endSlice)+1;

    int len;
    len = (int)strlen(startSlice);
    len = (int)strlen(endSlice);
    len = (int)strlen(extraInfo);
    startSliceStr = startSlice;
    endSliceStr = endSlice;
    extraInfoStr = extraInfo;
    stateflg = *(headerInfo+headerLength-1);


    fp.read((char *)&x,sizeof(int)*1);
    fp.read((char *)&y,sizeof(int)*1);
    fp.read((char *)&z,sizeof(int)*1);

    /* Get resolution values in my (X,Y,Z) coordinate system. */
    size_t XX = (size_t)y;
    size_t YY = (size_t)z;
    size_t ZZ = (size_t)x;

    pixelsPerVolume = (unsigned long)(XX*YY*ZZ);

    /* Set volume resolution */
    UniformVolume<T>::ResetResolution(YY,XX,ZZ,(T)0.0e0);

    fp.read((char *)&xRotateAngle,sizeof(float));
    fp.read((char *)&yRotateAngle,sizeof(float));

    fp.read((char *)&xzoom,sizeof(float));
    fp.read((char *)&yzoom,sizeof(float));
    fp.read((char *)&zzoom,sizeof(float));

    fp.read((char *)&corex,sizeof(int));
    fp.read((char *)&corey,sizeof(int));
    fp.read((char *)&radius,sizeof(int));

    if(voxelState!=NULL)	{
        delete voxelState;
       }
    voxelState = new char[pixelsPerVolume];

    if(voxelState==NULL){
        fp.close();
    }



    /* Dummy code to avoid unused variable warnings. */
    len++;
    extraInfoStr++;
    endSliceStr++;
    startSliceStr++;




    /* Read values individually into a 'float' variable.  Then, typecast to place in the
     * data array for this object.  This is to avoid problems reading a 'float' VOL file
     * into a 'double' UniformVolume object. */
    std::string desc;
    desc = "Reading CNDE VOL File";
    qtsignals->EmitFunctionDesc(desc);
    qtsignals->EmitFunctionProgress(0.0e0);

    if(scalar_data == NULL){
        float ftmp = (float)0.0e0;
        for(size_t i=0; i<YY; i++){
            for(size_t j=0; j<XX; j++){
                for(size_t k=0; k<ZZ; k++){
                    fp.read(reinterpret_cast<char*>(&ftmp),sizeof(float));
                    pscalars(0)->operator ()(i,j,k) = (T)ftmp;
                }
            }
            qtsignals->EmitFunctionProgress((float)i/(float)YY);
        }
    } else {

        scalar_data_points_read = 0;

        float ftmp = (float)0.0e0;
        for(size_t i=0; i<YY; i++){
            for(size_t j=0; j<XX; j++){
                for(size_t k=0; k<ZZ; k++){
                    fp.read(reinterpret_cast<char*>(&ftmp),sizeof(float));
                    scalar_data[scalar_data_points_read] = (T)ftmp;

                    /* Update number of points read.  If this value is greater than the size of the
                     * array, inform the user of the error and return. */
                    scalar_data_points_read++;
                    if(scalar_data_points_read > scalar_data_size){
                        std::cerr << "UniformVolume::ReadVOLFile() ERROR: Insufficient array size provided." << std::endl;
                        return;
                    }
                }
            }
            qtsignals->EmitFunctionProgress((float)i/(float)YY);
        }
    }

    desc = "";
    qtsignals->EmitFunctionDesc(desc);
    qtsignals->EmitFunctionProgress(0.0e0);

    if(stateflg=='Y')
    fp.read((char *)voxelState,sizeof(char)*pixelsPerVolume);

    fp.close();
    delete [] voxelState;



    /* Set other volume parameters */
    UniformVolume<T>::setSpatialExtent(1,-1,0.0e0);
    UniformVolume<T>::setSpatialExtent(1,1,(T)XX);
    UniformVolume<T>::setSpatialExtent(0,-1,0.0e0);
    UniformVolume<T>::setSpatialExtent(0,1,(T)YY);
    UniformVolume<T>::setSpatialExtent(2,-1,0.0e0);
    UniformVolume<T>::setSpatialExtent(2,1,(T)ZZ);
}


template <class T>
void UniformVolume<T>::ReadVTKFile(std::string filename, const bool isBigEndian)
{
    /*
      Delete all scalar and vector quantities
     */
    UniformVolume<T>::RemoveAllData();

    /* Define temporary variables */
    std::fstream vtkfile;
    char discard[256];
    std::string ascii_vs_binary("neither");
    std::stringstream sstmp;
    std::string format("unknown");


    vtkfile.open(filename.c_str(), std::ios::in);

    if(vtkfile.is_open() == true){

        vtkfile.getline(discard,256);           /* VTK version line */

        vtkfile.close();

        /* Check first character of first line.  If '#', file is legacy VTK format.  If '<'
         * file is XML format.  Call appropriate reader. */
        if(strncmp(&discard[0],"#",1) == 0){
            /* File is legacy format. */
//            format = "legacy";
            UniformVolume<T>::ReadLegacyVTKFile(filename);
        }
        if(strncmp(&discard[0],"<",1) == 0){
            /* File is XML format. */
//            format = "XML";
            UniformVolume<T>::ReadXMLVTKFile(filename);
        }


        /* Calls to code written by me.  Deprecated in favor of using VTK libraries. */
        if(false){
            if(format == "legacy"){
                vtkfile.getline(discard,256);           /* Description line */
                vtkfile.getline(discard,256);           /* ASCII vs. Binary */

                sstmp.str(discard);
                ascii_vs_binary = sstmp.str();

                if(ascii_vs_binary == "ASCII"){

                    /* Call ASCII reader. */
                    UniformVolume<T>::VTKReadLegacyASCII(vtkfile);
                    vtkfile.close();

                }
                if(ascii_vs_binary == "BINARY"){

                    /* Close file opened in ASCII mode. */
                    vtkfile.close();

                    /* Re-open file in binary mode. */
                    vtkfile.open(filename.c_str(), std::ios::in | std::ios::binary);

                    if(vtkfile.is_open() == true){

                        /* Re-read initial header lines. */
                        vtkfile.getline(discard,256);           /* VTK version line */
                        vtkfile.getline(discard,256);           /* Description line */
                        vtkfile.getline(discard,256);           /* ASCII vs. Binary */


                        /* Call binary reader for legacy VTK format. */
                        UniformVolume<T>::VTKReadLegacyBinary(vtkfile,isBigEndian);
                        vtkfile.close();

                    } /* End of conditional on successful file open. */
                } /* End of conditional on binary file. */
            } /* End of check for legacy format. */

            if(format == "XML"){
                std::cerr << "UniformVolume::ReadVTKFile() - XML format not supported" << std::endl;

                vtkfile.close();
            }


            /* Close file if it's still open. */
            if(vtkfile.is_open()){
                vtkfile.close();
            }

        }



    } else {
        std::cout << "UniformVolume::ReadVTKFile() - Could not open file " << filename << " for reading" << std::endl;
    }
}


template <class T>
void UniformVolume<T>::ReadLegacyVTKFile(std::string filename)
{
    UniformVolume<T>::RemoveAllData();


    /* imageExport is used to move read-in data to a C-style array controlled by this object.  If the data cannot
     * be moved to an array controlled by this object, it will be copied.  If copied, we want the VTK library to
     * free the memory it allocates.  This means that the dataset which is returned by the file reading function is
     * not actually assigned to the export object until we know that a move (not a copy) is to be performed. */
    vtkSmartPointer<vtkImageExport> imageExport = vtkSmartPointer<vtkImageExport>::New();


    vtkSmartPointer<vtkStructuredPointsReader> reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    reader->SetFileName(filename.c_str());
    reader->UpdateInformation();
    reader->UpdateWholeExtent();

    vtkSmartPointer<vtkImageData> dataset = vtkSmartPointer<vtkImageData>::New();
    dataset = vtkImageData::SafeDownCast(reader->GetOutput());


    int dims[3] = {0, 0, 0};
    double *origin;
    double *spacing;

    dataset->GetDimensions(dims);
    origin = dataset->GetOrigin();
    spacing = dataset->GetSpacing();
    char *name = dataset->GetPointData()->GetScalars()->GetName();
    std::string str_name(name);


    int type_this;
    if(typeid(T) == typeid(char)){
        type_this = 2;
    }

    if(typeid(T) == typeid(signed char)){
        type_this = 15;
    }

    if(typeid(T) == typeid(unsigned char)){
        type_this = 3;
    }

    if(typeid(T) == typeid(short)){
        type_this = 4;
    }

    if(typeid(T) == typeid(unsigned short)){
        type_this = 5;
    }

    if(typeid(T) == typeid(int)){
        type_this = 6;
    }

    if(typeid(T) == typeid(unsigned int)){
        type_this = 7;
    }

    if(typeid(T) == typeid(long)){
        type_this = 8;
    }

    if(typeid(T) == typeid(unsigned long)){
        type_this = 9;
    }

    if(typeid(T) == typeid(float)){
        type_this = 10;
    }

    if(typeid(T) == typeid(double)){
        type_this = 11;
    }




    int type_file = dataset->GetScalarType();

    /* If either data is to be read into a user-supplied array (pointer to which is 'scalar_data'), or the datatypes
     * do not match, the read-in data will have to be copied to its destination array. */
    if(scalar_data || type_this != type_file){

        size_t npts = dims[0]*dims[1]*dims[2];

        if(scalar_data){

            /* If data is to be read into pre-existing array, copy values into that array. */

            for(size_t ii=0; ii<scalar_data_size; ii++){

                scalar_data[ii] = (T)dataset->GetPointData()->GetScalars()->GetComponent(ii, 0);
                scalar_data_points_read++;

                if(scalar_data_points_read == npts){
                    std::cerr << "UniformVolume<T>::ReadLegacyVTKFile()  ERROR: Insufficient destination array size." << std::endl;
                    std::cerr << "                                              Aborting." << std::endl;
                    return;
                }
            }

        } else {

            /* If no pre-existing array is to be used, then we must place data into this object, and typecast the
             * values to match the template type. */

            std::string str_name("nothing");

            UniformVolume<T>::ResetResolution((size_t)dims[1], (size_t)dims[0], (size_t)dims[2], (T)0.0e0);
            UniformVolume<T>::AddScalarQuantity(str_name);

            for(size_t i=0; i<npts; i++){
                pscalars(0)->operator [](i) = (T)dataset->GetPointData()->GetScalars()->GetComponent(i, 0);
            }

        }

    } else {

        /* Data destination is this object, and datatypes match, so the data pointer can simply be assigned
         * and no duplication of data is required. */

        reader->Register(reader->GetOutput());          /* Prevents loss of data when 'reader' goes out of scope.  This
                                                         * is desired since we want this object to control the memory
                                                         * rather than the VTK library functions. */
        imageExport->SetInputData(dataset);             /* Assign data to export object for movement into separate
                                                         * C-style array. */

//        UniformVolume<T>::ResetResolution(1, 1, 1, (T)0.0e0);
        UniformVolume<T>::AddScalarQuantity(str_name);

        vrows = dims[1];
        vcols = dims[0];
        vslices = dims[2];

        /* Setting the pointer like this appears to "move" the data into my array structure.  Valgrind
         * does not show any missing deallocation as a result of doing this, so I think this is OK. */
//        T* dataptr = (T*)imageExport->GetPointerToData();
//        imageExport->Export(&pscalars(0)->operator [](0));
//        T* dataptr;
//        imageExport->Export(dataptr);
//        pscalars(0)->SetArraySize(vrows, vcols, vslices, true);

        imageExport->Update();
        pscalars(0)->SetArrayPointer((T*)imageExport->GetPointerToData(), vrows, vcols, vslices, true);
    }


    xmin = origin[0];
    ymin = origin[1];
    zmin = origin[2];

    xspacing = spacing[0];
    yspacing = spacing[1];
    zspacing = spacing[2];

    xmax = xmin + xspacing*(float)vcols;
    ymax = ymin + yspacing*(float)vrows;
    zmax = zmin + zspacing*(float)vslices;

}


template <class T>
void UniformVolume<T>::ReadXMLVTKFile(std::string filename)
{
    UniformVolume<T>::RemoveAllData();


    vtkSmartPointer<vtkImageExport> imageExport = vtkSmartPointer<vtkImageExport>::New();
    vtkImageData *dataset;

    std::string extension;
    extension = StringManip::DetermFileExt(filename);

    if(extension == "pvti" || extension == "PVTI"){
        vtkSmartPointer<vtkXMLPImageDataReader> reader = vtkSmartPointer<vtkXMLPImageDataReader>::New();
        reader->SetFileName(filename.c_str());
        reader->UpdateInformation();
        reader->UpdateWholeExtent();
        reader->GetOutput()->Register(reader);
        dataset = vtkImageData::SafeDownCast(reader->GetOutput());
    } else {
        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        reader->SetFileName(filename.c_str());
        reader->UpdateInformation();
        reader->UpdateWholeExtent();
        reader->GetOutput()->Register(reader);
        dataset = vtkImageData::SafeDownCast(reader->GetOutput());
    }

    int dims[3] = {0, 0, 0};
    double *origin;
    double *spacing;

    dataset->GetDimensions(dims);
    origin = dataset->GetOrigin();
    spacing = dataset->GetSpacing();
    char *name = dataset->GetPointData()->GetScalars()->GetName();
    std::string str_name(name);

    int type_this;
    if(typeid(T) == typeid(char)){
        type_this = 2;
    }

    if(typeid(T) == typeid(signed char)){
        type_this = 15;
    }

    if(typeid(T) == typeid(unsigned char)){
        type_this = 3;
    }

    if(typeid(T) == typeid(short)){
        type_this = 4;
    }

    if(typeid(T) == typeid(unsigned short)){
        type_this = 5;
    }

    if(typeid(T) == typeid(int)){
        type_this = 6;
    }

    if(typeid(T) == typeid(unsigned int)){
        type_this = 7;
    }

    if(typeid(T) == typeid(long)){
        type_this = 8;
    }

    if(typeid(T) == typeid(unsigned long)){
        type_this = 9;
    }

    if(typeid(T) == typeid(float)){
        type_this = 10;
    }

    if(typeid(T) == typeid(double)){
        type_this = 11;
    }




    int type_file = imageExport->GetDataScalarType();

    /* If either data is to be read into a user-supplied array (pointer to which is 'scalar_data'), or the datatypes
     * do not match, the read-in data will have to be copied to its destination array. */
    if(scalar_data || type_this != type_file){

        size_t npts = dims[0]*dims[1]*dims[2];

        vtkSmartPointer<vtkDoubleArray> data = vtkDoubleArray::SafeDownCast(dataset->GetPointData()->GetScalars());

        if(scalar_data){

            /* If data is to be read into pre-existing array, copy values into that array. */

            for(size_t ii=0; ii<scalar_data_size; ii++){

                scalar_data[ii] = (T)data->GetValue(ii);
                scalar_data_points_read++;

                if(scalar_data_points_read == npts){
                    std::cerr << "UniformVolume<T>::ReadLegacyVTKFile()  ERROR: Insufficient destination array size." << std::endl;
                    std::cerr << "                                              Aborting." << std::endl;
                    return;
                }
            }

        } else {

            /* If no pre-existing array is to be used, then we must place data into this object, and typecast the
             * values to match the template type. */

            UniformVolume<T>::ResetResolution((size_t)dims[1], (size_t)dims[0], (size_t)dims[2], (T)0.0e0);
            UniformVolume<T>::AddScalarQuantity(str_name);

            for(size_t i=0; i<npts; i++){
                pscalars(0)->operator [](i) = (T)data->GetValue(i);
            }

        }

    } else {

        /* Data destination is this object, and datatypes match, so the data pointer can simply be assigned
         * and no duplication of data is required. */

        imageExport->SetInputData(dataset);             /* Assign data to export object for movement into separate
                                                         * C-style array. */

        UniformVolume<T>::AddScalarQuantity(str_name);
        vcols = dims[0];
        vrows = dims[1];
        vslices = dims[2];

        /* Setting the pointer like this appears to "move" the data into my array structure.  Valgrind
         * does not show any missing deallocation as a result of doing this, so I think this is OK. */
        pscalars(0)->SetArrayPointer((T*)imageExport->GetPointerToData(), vrows, vcols, vslices, true);

    }


    xmin = origin[0];
    ymin = origin[1];
    zmin = origin[2];

    xspacing = spacing[0];
    yspacing = spacing[1];
    zspacing = spacing[2];

    xmax = xmin + xspacing*(float)vcols;
    ymax = ymin + yspacing*(float)vrows;
    zmax = zmin + zspacing*(float)vslices;

} /* UniformVolume<T>::ReadXMLVTKFile() */











/* =======================================================================================
 *
 *
 *                  PUBLIC FUNCTIONS
 *
 *
 *
 */
template <class T>
UniformVolume<T>::UniformVolume()
{
    UniformVolume<T>::Initialize(1,1,1, -1.0e0, 1.0e0, -1.0e0, 1.0e0, -1.0e0, 1.0e0, 0.0e0, 1, 0, 256);
}


template <class T>
UniformVolume<T>::UniformVolume(const UniformVolume<T> &uc3d) : QtIntermediaryBase()
{
    /* Set resolution and copy all data */
    nscalars = 0;
    nvectors = 0;
    qtysize = 256;

    int nx = 0; int ny = 0; int nz = 0;
    nx = uc3d.GetResolution(1);
    ny = uc3d.GetResolution(0);
    nz = uc3d.GetResolution(2);
    UniformVolume<T>::ResetResolution(ny,nx,nz,0.0e0);

    int n_scalars = uc3d.NumScalarQuantities();
    int n_vectors = uc3d.NumVectorQuantities();
    std::string qtyname("name");

    for(int q=0; q<n_scalars; q++){
        qtyname = uc3d.ScalarQuantityName(q);
        UniformVolume<T>::AddScalarQuantity(qtyname);
        for(int k=0; k<vslices; k++){
            for(int i=0; i<vrows; i++){
                for(int j=0; j<vcols; j++){
                    pscalars(nscalars-1)->operator ()(i,j,k) = uc3d(i,j,k,q);
                }
            }
        }
    }

    int ncomp = 0;
    for(int q=0; q<n_vectors; q++){
        qtyname = uc3d.VectorQuantityName(q);
        ncomp = uc3d.VectorQuantityComponents(q);
        UniformVolume<T>::AddVectorQuantity(qtyname,ncomp);
        for(int k=0; k<vslices; k++){
            for(int i=0; i<vrows; i++){
                for(int j=0; j<vcols; j++){
                    for(int l=0; l<ncomp; l++){
                        pvectors(nvectors-1)->operator ()(i,j,k,l) = uc3d(i,j,k,l,q);
                    }
                }
            }
        }
    }

    /* Set spatial parameters */
    xmin = uc3d.SpatialExtent(1,-1); xmax = uc3d.SpatialExtent(1,1);
    ymin = uc3d.SpatialExtent(0,-1); ymax = uc3d.SpatialExtent(0,1);
    zmin = uc3d.SpatialExtent(2,-1); zmax = uc3d.SpatialExtent(2,1);
    if(vcols > 1){
        xspacing = (xmax - xmin)/((T)vcols - 1.0e0);
    } else {
        xspacing = 0.0e0;
    }
    if(vrows > 1){
        yspacing = (ymax - ymin)/((T)vrows - 1.0e0);
    } else {
        yspacing = 0.0e0;
    }
    if(vslices > 1){
        zspacing = (zmax - zmin)/((T)vslices - 1.0e0);
    } else {
        zspacing = 0.0e0;
    }

    std::cout << "Spatial limits: x " << xmin << " , " << xmax << "  " << xspacing << std::endl;
    std::cout << "                y " << ymin << " , " << ymax << "  " << yspacing << std::endl;
    std::cout << "                z " << zmin << " , " << zmax << "  " << zspacing << std::endl;

    /* Set output parameters */
    outputdir = uc3d.OutputDirectory();
    filenamestem = uc3d.FilenameStem();
    imageoutput = uc3d.VTKImageDataOutput();
    rectoutput = uc3d.VTKRectDataOutput();
    volname = filenamestem;

    /* Set type name */
    dtypename = "unknown_type";
    if(typeid(T) == typeid(char)){
        dtypename = "char";
    }
    if(typeid(T) == typeid(short)){
        dtypename = "short";
    }
    if(typeid(T) == typeid(int)){
        dtypename = "int";
    }
    if(typeid(T) == typeid(long)){
        dtypename = "long";
    }
    if(typeid(T) == typeid(long long)){
        dtypename = "long long";
    }
    if(typeid(T) == typeid(float)){
        dtypename = "float";
    }
    if(typeid(T) == typeid(double)){
        dtypename = "double";
    }
    if(typeid(T) == typeid(long double)){
        dtypename = "long double";
    }

    writebigendian = false;
}


#ifdef CXX11
template <class T>
UniformVolume<T>::UniformVolume(UniformVolume<T> &&a) : UniformVolume<T>()
{
    UniformVolumeSwap(*this, a);
}
#endif


template <class T>
UniformVolume<T>::~UniformVolume()
{
    if(qtsignals){
        delete qtsignals;
    }
}


template <class T>
UniformVolume<T>& UniformVolume<T>::operator=(UniformVolume<T> a)
{
    UniformVolumeSwap(*this, a);
    return *this;
}


template <class T>
size_t UniformVolume<T>::NumScalarQuantities() const
{
    return nscalars;
}


template <class T>
size_t UniformVolume<T>::NumVectorQuantities() const
{
    return nvectors;
}


template <class T>
void UniformVolume<T>::AddScalarQuantity(const std::string name)
{
    // Increase count of number of scalar quantities
    nscalars++;

    /*
      Increase the size of the PArray1D object containing pointers to
      Array3D objects containing data.  Pointers to existing arrays are
      copied to a temporary PArray1D so that the resize operation does not
      cause existing data to be lost.

      No extra memory is required for this approach since the data is moved
      between PArray1D objects rather than copied.

      After the move is complete, add a new Array3D object to the final
      entry in the resized PArray1D object, set its size equal to the grid
      size, and all points initialized to 0.0.
     */

    /*
      If the first element of pscalars is NULL, no quantity exists to be copied,
      so this first part can be skipped.
     */
    if(pscalars(0) != NULL){
        PArray1D<Array3D<T>*> tmparray(nscalars);
        for(size_t i=0; i<nscalars-1; i++){
            /* Move pointer to data to temporary array.  Set pointer within pscalars to NULL so that
             * the resize operation does not cause data to be lost. */
            tmparray(i) = pscalars(i);
            pscalars(i) = NULL;
        }
        pscalars.ResetSize(nscalars);
        for(size_t i=0; i<nscalars-1; i++){
            /* Move pointer from temporary array back to pscalars. */
            pscalars(i) = tmparray(i);
            tmparray(i) = NULL;
        }

        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nscalars);
        for(size_t i=0; i<nscalars-1; i++){
            tmparray2(i) = scalar_names(i);
            scalar_names(i) = NULL;
        }
        scalar_names.ResetSize(nscalars);
        for(size_t i=0; i<nscalars-1; i++){
            scalar_names(i) = tmparray2(i);
            tmparray2(i) = NULL;
        }
    }

    /* By this point, the last element of both 'pscalars' and 'scalar_names' has not been set.  Create
     * the appropriate objects for each array and initialize values.  */
//    std::stringstream procdesc;
//    procdesc << "Allocating memory for scalar quantity '" << name << "'";
//    qtsignals->EmitFunctionDesc2(procdesc.str());

    pscalars(nscalars-1) = new Array3D<T>;
    pscalars(nscalars-1)->ResetSize(vrows,vcols,vslices);
    pscalars(nscalars-1)->ResetVal((T)0.0e0);

    scalar_names(nscalars-1) = new std::string;
    scalar_names(nscalars-1)->reserve(qtysize);
    scalar_names(nscalars-1)->assign(name);

//    procdesc.str("");
//    qtsignals->EmitFunctionDesc2(procdesc.str());
}


template <class T>
void UniformVolume<T>::AddScalarQuantity(const std::string name, Array3D<T> *data)
{
    // Increase count of number of scalar quantities
    nscalars++;

    /*
      Increase the size of the PArray1D object containing pointers to
      Array3D objects containing data.  Pointers to existing arrays are
      copied to a temporary PArray1D so that the resize operation does not
      cause existing data to be lost.

      No extra memory is required for this approach since the data is moved
      between PArray1D objects rather than copied.

      After the move is complete, add a new Array3D object to the final
      entry in the resized PArray1D object, set its size equal to the grid
      size, and all points initialized to 0.0.
     */

    /*
      If the first element of pscalars is NULL, no quantity exists to be copied,
      so this first part can be skipped.
     */
    if(pscalars(0) != NULL){
        PArray1D<Array3D<T>*> tmparray(nscalars);
        for(size_t i=0; i<nscalars-1; i++){
            /* Move pointer to data to temporary array.  Set pointer within pscalars to NULL so that
             * the resize operation does not cause data to be lost. */
            tmparray(i) = pscalars(i);
            pscalars(i) = NULL;
        }
        pscalars.ResetSize(nscalars);
        for(size_t i=0; i<nscalars-1; i++){
            /* Move pointer from temporary array back to pscalars. */
            pscalars(i) = tmparray(i);
            tmparray(i) = NULL;
        }

        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nscalars);
        for(size_t i=0; i<nscalars-1; i++){
            tmparray2(i) = scalar_names(i);
            scalar_names(i) = NULL;
        }
        scalar_names.ResetSize(nscalars);
        for(size_t i=0; i<nscalars-1; i++){
            scalar_names(i) = tmparray2(i);
            tmparray2(i) = NULL;
        }
    }

    /* By this point, the last element of both 'pscalars' and 'scalar_names' has not been set.  Set the
     * pointer to the data to the passed-in value and set the appropriate name.  */
    pscalars(nscalars-1) = data;

    scalar_names(nscalars-1) = new std::string;
    scalar_names(nscalars-1)->reserve(qtysize);
    scalar_names(nscalars-1)->assign(name);
}


template <class T>
void UniformVolume<T>::RemoveScalarQuantityRef(const size_t qty)
{
    /* Set pointer to NULL to prevent deletion of data. */
    pscalars(qty) = NULL;

    /* Remove quantity from array of pointers. */
    UniformVolume<T>::RemoveScalarQuantity(qty);

    std::cerr << "Removed scalar quantity reference " << qty << std::endl;
}


template <class T>
void UniformVolume<T>::RemoveAllData()
{
    /* Set member variables tracking data quantities to 0.  Reset arrays of pointers to data and data-names
     * to contain a single NULL element (cannot create 0-length array, so create array with no real data). */
    nscalars = 0;
    nvectors = 0;

    pscalars.ResetSize(1);
    pscalars(0) = NULL;

    scalar_names.ResetSize(1);
    scalar_names(0) = NULL;

    pvectors.ResetSize(1);
    pvectors(0) = NULL;

    vector_names.ResetSize(1);
    vector_names(0) = NULL;
}


template <class T>
Array3D<T>* UniformVolume<T>::PointerToScalarData(const size_t qty)
{
    return pscalars(qty);
}


template <class T>
void UniformVolume<T>::RemoveScalarQuantity(const size_t qty)
{
    // Decrease count of number of scalar quantities
    nscalars--;

    /*
      Decrease the size of the PArray1D object containing pointers to
      Array3D objects containing data.  As with the AddScalarQuantity function,
      pointers to the data are moved between temporary arrays to avoid the need
      for extra memory to copy data.
     */

    /*
      If the first element of pscalars is NULL, no scalar quantities exist, thus there
      are none to be removed.
     */
    if(pscalars(0) != NULL && nscalars > 0){
        PArray1D<Array3D<T>*> tmparray(nscalars+1);
        for(size_t i=0; i<nscalars+1; i++){

            /* Save pointer to data only for data not to be removed.  Missing element is set to
             * NULL in the initialization of 'tmparray'. */
            if(i != qty){
                tmparray(i) = pscalars(i);
                pscalars(i) = NULL;
            }
        }
        pscalars.ResetSize(nscalars);
        size_t j = 0;
        for(size_t i=0; i<nscalars+1; i++){

            /* Move pointer to data only if loop iteration is for data to be retained.  Variable
             * 'j' is only incremented when a valid pointer is moved.  When i > qty, j = i-1. */
            if(i != qty){
                pscalars(j) = tmparray(i);
                tmparray(i) = NULL;
                j++;
            }
        }


        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nscalars+1);
        for(size_t i=0; i<nscalars+1; i++){
            tmparray2(i) = scalar_names(i);
            scalar_names(i) = NULL;
        }
        scalar_names.ResetSize(nscalars);
        j = 0;
        for(size_t i=0; i<nscalars+1; i++){
            if(i != qty){
                if(i != qty){
                    scalar_names(j) = tmparray2(i);
                    tmparray2(i) = NULL;
                    j++;
                }
            }
        }
    }

    /*
      If updated value of nscalars is 0, delete pointer to first object and set
      to NULL.
     */
    if(nscalars == 0){
        if(pscalars(0) != NULL){
            delete pscalars(0);
        }
        pscalars(0) = NULL;

        if(scalar_names(0) != NULL){
            delete scalar_names(0);
        }
        scalar_names(0) = NULL;
    }
}


template <class T>
void UniformVolume<T>::AddVectorQuantity(const std::string name, const size_t ncomp)
{
    // Increase count of number of scalar quantities
    nvectors++;

    /*
      Increase the size of the PArray1D object containing pointers to
      Array4D objects containing data.  Pointers to existing arrays are
      copied to a temporary PArray1D so that the resize operation does not
      cause existing data to be lost.

      No extra memory is required for this approach since the data is moved
      between PArray1D objects rather than copied.

      After the move is complete, add a new Array3D object to the final
      entry in the resized PArray1D object, set its size equal to the grid
      size, and all points initialized to 0.0.
     */

    /*
      If the first element of pscalars is NULL, no quantity exists to be copied,
      so this first part can be skipped.
     */
    if(pvectors(0) != NULL){
        PArray1D<Array4D<T>*> tmparray(nvectors);
        for(size_t i=0; i<nvectors-1; i++){
            /* Move pointer to data to temporary array.  Set pointer within pscalars to NULL so that
             * the resize operation does not cause data to be lost. */
            tmparray(i) = pvectors(i);
            pvectors(i) = NULL;
        }
        pvectors.ResetSize(nvectors);
        for(size_t i=0; i<nvectors-1; i++){
            /* Move pointer from temporary array back to pscalars. */
            pvectors(i) = tmparray(i);
            tmparray(i) = NULL;
        }

        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nvectors);
        for(size_t i=0; i<nvectors-1; i++){
            tmparray2(i) = vector_names(i);
            vector_names(i) = NULL;
        }
        vector_names.ResetSize(nvectors);
        for(size_t i=0; i<nvectors-1; i++){
            vector_names(i) = tmparray2(i);
            tmparray2(i) = NULL;
        }
    }

    /* By this point, the last element of both 'pscalars' and 'scalar_names' has not been set.  Create
     * the appropriate objects for each array and initialize values.  */
//    std::stringstream procdesc;
//    procdesc << "Allocating memory for vector quantity '" << name << "'";
//    qtsignals->EmitFunctionDesc2(procdesc.str());

    pvectors(nvectors-1) = new Array4D<T>;
    pvectors(nvectors-1)->ResetSize(vrows,vcols,vslices,ncomp);
    pvectors(nvectors-1)->ResetVal((T)0.0e0);

    vector_names(nvectors-1) = new std::string;
    vector_names(nvectors-1)->reserve(qtysize);
    vector_names(nvectors-1)->assign(name);

//    procdesc.str("");
//    qtsignals->EmitFunctionDesc2(procdesc.str());
}


template <class T>
void UniformVolume<T>::RemoveVectorQuantity(const size_t qty)
{
    // Decrease count of number of vector quantities
    nvectors--;

    /*
      Decrease the size of the PArray1D object containing pointers to
      Array4D objects containing data.  As with the AddVectorQuantity function,
      pointers to the data are moved between temporary arrays to avoid the need
      for extra memory to copy data.
     */

    /*
      If the first element of pvectors is NULL, no vector quantities exist, thus there
      are none to be removed.
     */
    if(pvectors(0) != NULL && nvectors > 0){
        PArray1D<Array4D<T>*> tmparray(nvectors+1);
        for(size_t i=0; i<nvectors+1; i++){

            /* Save pointer to data only for data not to be removed.  Missing element is set to
             * NULL in the initialization of 'tmparray'. */
            if(i != qty){
                tmparray(i) = pvectors(i);
                pvectors(i) = NULL;
            }
        }
        pvectors.ResetSize(nvectors);
        size_t j = 0;
        for(size_t i=0; i<nvectors+1; i++){

            /* Move pointer to data only if loop iteration is for data to be retained.  Variable
             * 'j' is only incremented when a valid pointer is moved.  When i > qty, j = i-1. */
            if(i != qty){
                pvectors(j) = tmparray(i);
                tmparray(i) = NULL;
                j++;
            }
        }


        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nvectors+1);
        for(size_t i=0; i<nvectors+1; i++){
            tmparray2(i) = vector_names(i);
            vector_names(i) = NULL;
        }
        vector_names.ResetSize(nvectors);
        j = 0;
        for(size_t i=0; i<nvectors+1; i++){
            if(i != qty){
                if(i != qty){
                    vector_names(j) = tmparray2(i);
                    tmparray2(i) = NULL;
                    j++;
                }
            }
        }
    }

    /*
      If updated value of nvectors is 0, delete pointer to first object and set
      to NULL.
     */
    if(nvectors == 0){
        if(pvectors(0) != NULL){
            delete pvectors(0);
        }
        pvectors(0) = NULL;

        if(vector_names(0) != NULL){
            delete vector_names(0);
        }
        vector_names(0) = NULL;
    }
}


template <class T>
int UniformVolume<T>::VectorQuantityComponents(const int qty) const
{
    return pvectors(qty)->GetDim(4);
}



template <class T>
void UniformVolume<T>::setScalarQuantityName(const size_t qty, const std::string name)
{
    scalar_names(qty)->assign(name);
}


template <class T>
std::string UniformVolume<T>::ScalarQuantityName(const size_t qty)
{
    return scalar_names(qty)->substr();
}


template <class T>
void UniformVolume<T>::setVectorQuantityName(const size_t qty, const std::string name)
{
    vector_names(qty)->assign(name);
}


template <class T>
std::string UniformVolume<T>::VectorQuantityName(const size_t qty)
{
    return vector_names(qty)->substr();
}



template < class T > inline
T& UniformVolume<T>::operator()(const size_t ind1, const size_t ind2,
                                  const size_t ind3, const size_t qty)
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 < vrows);
        assert(ind2 < vcols);
        assert(ind3 < vslices);
        return pscalars(qty)->operator ()(ind1,ind2,ind3);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pscalars(qty)->operator ()(ind1,ind2,ind3);
    #endif
}


template < class T > inline
const T& UniformVolume<T>::operator()(const size_t ind1, const size_t ind2,
                                        const size_t ind3, const size_t qty) const
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 < vrows);
        assert(ind2 < vcols);
        assert(ind3 < vslices);
        return pscalars(qty)->operator ()(ind1,ind2,ind3);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pscalars(qty)->operator ()(ind1,ind2,ind3);
    #endif
}


template < class T > inline
T& UniformVolume<T>::operator()(const size_t ind1, const size_t ind2,
                                  const size_t ind3, const size_t comp, const size_t qty)
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 < vrows);
        assert(ind2 < vcols);
        assert(ind3 < vslices);
        assert(comp < pvectors(qty)->GetDim(4));
        return pvectors(qty)->operator ()(ind1,ind2,ind3,comp);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pvectors(qty)->operator ()(ind1,ind2,ind3,comp);
    #endif
}


template < class T > inline
const T& UniformVolume<T>::operator()(const size_t ind1, const size_t ind2,
                                        const size_t ind3, const size_t comp, const size_t qty) const
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 < vrows);
        assert(ind2 < vcols);
        assert(ind3 < vslices);
        assert(comp < pvectors(qty)->GetDim(4));
        return pvectors(qty)->operator ()(ind1,ind2,ind3,comp);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pvectors(qty)->operator ()(ind1,ind2,ind3,comp);
    #endif
}


template <class T>
double UniformVolume<T>::EstimateMemoryUsage() const
{
    double retval = 0.0e0;
    for(int i=0; i<nscalars; i++){
        retval += pscalars(i)->GetMemoryUsage();
    }
    for(int i=0; i<nvectors; i++){
        retval += pvectors(i)->GetMemoryUsage();
    }
    return retval;
}


template <class T>
void UniformVolume<T>::ResetScalarVals(const int qty, const T initvalue)
{
    pscalars(qty)->ResetVal(initvalue);
}


template <class T>
void UniformVolume<T>::ResetVectorVals(const int qty, const T initvalue)
{
    pvectors(qty)->ResetVal(initvalue);
}


template <class T>
void UniformVolume<T>::setSpatialExtent(const int dir, const int maxmin, const T val)
{
    if(dir == 1){
        if(maxmin < 0){
            xmin = val;
            xminset = true;
        }
        if(maxmin > 0){
            xmax = val;
            xmaxset = true;
        }
        if(xminset == true && xmaxset == true){
            if(vcols > 1){
                xspacing = (xmax - xmin)/((T)vcols - 1.0e0);
            } else {
                xspacing = 0.0e0;
            }
        }
    }
    if(dir == 0){
        if(maxmin < 0){
            ymin = val;
            yminset = true;
        }
        if(maxmin > 0){
            ymax = val;
            ymaxset = true;
        }
        if(yminset == true && ymaxset == true){
            if(vrows > 1){
                yspacing = (ymax - ymin)/((T)vrows - 1.0e0);
            } else {
                yspacing = 0.0e0;
            }
        }
    }
    if(dir == 2){
        if(maxmin < 0){
            zmin = val;
            zminset = true;
        }
        if(maxmin > 0){
            zmax = val;
            zmaxset = true;
        }
        if(zminset == true && zmaxset == true){
            if(vslices > 1){
                zspacing = (zmax - zmin)/((T)vslices - 1.0e0);
            } else {
                zspacing = 0.0e0;
            }
        }
    }
}


template <class T>
T UniformVolume<T>::SpatialExtent(const int dir, const int maxmin) const
{
    T retval = 0.0e0;

    if(dir == 1){
        if(maxmin < 0){
            retval = xmin;
        }
        if(maxmin > 0){
            retval = xmax;
        }
    }
    if(dir == 0){
        if(maxmin < 0){
            retval = ymin;
        }
        if(maxmin > 0){
            retval = ymax;
        }
    }
    if(dir == 2){
        if(maxmin < 0){
            retval = zmin;
        }
        if(maxmin > 0){
            retval = zmax;
        }
    }

    return retval;
}


template <class T>
T UniformVolume<T>::PointSpacing(const int dir) const
{
    T retval = 0.0e0;
    if(dir == 0){
        retval = xspacing;
    }
    if(dir == 1){
        retval = yspacing;
    }
    if(dir == 2){
        retval = zspacing;
    }
    return retval;
}


template <class T>
void UniformVolume<T>::setVTKImageDataOutput(const bool state)
{
    imageoutput = state;
    if(imageoutput == true){
        rectoutput = false;
    }
}


template <class T>
void UniformVolume<T>::setVTKRectDataOutput(const bool state)
{
    rectoutput = state;
    if(rectoutput == true){
        imageoutput = false;
    }
}


template <class T>
bool UniformVolume<T>::VTKImageDataOutput() const
{
    return imageoutput;
}


template <class T>
bool UniformVolume<T>::VTKRectDataOutput() const
{
    return rectoutput;
}


template <class T>
void UniformVolume<T>::setXDMFOutput(const bool state)
{
    xdmfoutput = state;
}


template <class T>
bool UniformVolume<T>::XDMFOutput() const
{
    return xdmfoutput;
}


template <class T>
void UniformVolume<T>::setFilenameStem(const std::string stem)
{
    filenamestem = stem;
}


template <class T>
std::string UniformVolume<T>::FilenameStem() const
{
    return filenamestem;
}


template <class T>
void UniformVolume<T>::setOutputDirectory(const std::string dir)
{
    outputdir = dir;
}


template <class T>
std::string UniformVolume<T>::OutputDirectory() const
{
    return outputdir;
}


template <class T>
void UniformVolume<T>::ResetResolution(const size_t ny, const size_t nx, const size_t nz, const T initval)
{
    for(size_t i=0; i<nscalars; i++){
        if(pscalars(i) != NULL){
            pscalars(i)->ResetSize(ny,nx,nz,initval);
        }
    }
    for(size_t i=0; i<nvectors; i++){
        if(pvectors(i) != NULL){
            size_t n = pvectors(i)->GetDim(3);
            pvectors(i)->ResetSize(ny,nx,nz,n,initval);
        }
    }
    vrows = ny;
    vcols = nx;
    vslices = nz;

    if(vcols > 1){
        xspacing = (xmax - xmin)/((T)vcols - 1.0e0);
    } else {
        xspacing = 0.0e0;
    }

    if(vrows > 1){
        yspacing = (ymax - ymin)/((T)vrows - 1.0e0);
    } else {
        yspacing = 0.0e0;
    }

    if(vslices > 1){
        zspacing = (zmax - zmin)/((T)vslices - 1.0e0);
    } else {
        zspacing = 0.0e0;
    }
}


template <class T>
size_t UniformVolume<T>::GetResolution(const int dir)
{
    size_t retval = 0;
    if(dir == 1){
        retval = vcols;
    }
    if(dir == 0){
        retval = vrows;
    }
    if(dir == 2){
        retval = vslices;
    }
    return retval;
}


template <class T>
void UniformVolume<T>::VTKWrite()
{

    if(false){
        UniformVolume<T>::VTKWriteImageData();      /* Currently some sort of bug in here. */
    } else {
        std::string filename;
        std::fstream ssfile;
        filename = outputdir + "/" + filenamestem + ".vtk";
        ssfile.open(filename.c_str(),std::ios::out);

        /*
          Variables needed for tracking progress
         */
        float Tchunks = (float)vslices*((float)nscalars + (float)nvectors);
        float chunkscomplete = 0.0e0;
        float completefrac = chunkscomplete/Tchunks;
        std::string descriptor("Writing ASCII VTK File");

        /*
          Emit signals for the beginning of the file-writing process
         */
        qtsignals->EmitFunctionDesc(descriptor);
        qtsignals->EmitFunctionProgress(completefrac);


        /* Temporary stream to reduce I/O calls. */
        std::stringstream tmpss;


        if(ssfile.is_open() == true){
            ssfile << "# vtk DataFile Version 2.0" << std::endl;
            ssfile << filenamestem << std::endl;
            ssfile << "ASCII" << std::endl;

            if(imageoutput == true){
                ssfile << "DATASET STRUCTURED_POINTS" << std::endl;
                ssfile << "DIMENSIONS " << vcols << " " << vrows << " " << vslices << std::endl;
                ssfile << "ORIGIN " << xmin << " " << ymin << " " << zmin << std::endl;
                ssfile << "SPACING " << xspacing << " " << yspacing << " " << zspacing << std::endl;
            }

            if(rectoutput == true){
                ssfile << "DATASET RECTILINEAR_GRID" << std::endl;
                ssfile << "DIMENSIONS " << vcols << " " << vrows << " " << vslices << std::endl;
                ssfile << "X_COORDINATES " << vcols << " float" << std::endl;
                for(int i=0; i<vcols; i++){
                    float x = xmin + float(i)*xspacing;
                    ssfile << x << std::endl;
                }
                ssfile << "Y_COORDINATES " << vrows << " float" << std::endl;
                for(int i=0; i<vrows; i++){
                    float y = ymin + float(i)*yspacing;
                    ssfile << y << std::endl;
                }
                ssfile << "Z_COORDINATES " << vslices << " float" << std::endl;
                for(int i=0; i<vslices; i++){
                    float z = zmin + float(i)*zspacing;
                    ssfile << z << std::endl;
                }
            }

            tmpss.clear();
            tmpss.str("");
            ssfile << "POINT_DATA " << vrows*vcols*vslices << std::endl;
            for(int s=0; s<nscalars; s++){

                descriptor = "Writing ASCII VTK File - scalars '" + ScalarQuantityName(s) + "'";
                qtsignals->EmitFunctionDesc(descriptor);

                chunkscomplete = 0.0e0;
                std::string name(scalar_names(s)->substr());
                ssfile << "SCALARS " << name << " " << dtypename << " 1" << std::endl;
                ssfile << "LOOKUP_TABLE default" << std::endl;
                for(int k=0; k<vslices; k++){
                    for(int i=0; i<vrows; i++){
                        for(int j=0; j<vcols; j++){
                            tmpss << pscalars(s)->operator ()(i,j,k) << std::endl;
                        }
                    }
                    ssfile << tmpss.str();
                    tmpss.clear();
                    tmpss.str("");
                    chunkscomplete += 1.0e0;
                    completefrac = chunkscomplete/Tchunks;
                    qtsignals->EmitFunctionProgress(completefrac, descriptor);
                    qtsignals->EmitFunctionProgress(completefrac);
                }
            }

            tmpss.clear();
            tmpss.str("");
            for(int v=0; v<nvectors; v++){

                descriptor = "Writing ASCII VTK File - vectors '" + VectorQuantityName(v) + "'";
                qtsignals->EmitFunctionDesc(descriptor);

                chunkscomplete = 0.0e0;
                std::string name(vector_names(v)->substr());
                int n = pvectors(v)->GetDim(4);
                if(n > 3){
                    std::cerr << "WARNING: VTK file format requires 3-component vectors." << std::endl;
                    std::cerr << "         Vector quantity " << name << " has " << n << " components." << std::endl;
                    std::cerr << "         File will be written with all components, which may cause" << std::endl;
                    std::cerr << "         problems with some file readers." << std::endl;
                }
                ssfile << "VECTORS " << name << " " << dtypename << std::endl;
                for(int k=0; k<vslices; k++){
                    for(int i=0; i<vrows; i++){
                        for(int j=0; j<vcols; j++){
                            for(int l=0; l<n; l++){
                                tmpss << pvectors(v)->operator ()(i,j,k,l) << " ";
                            }
                            tmpss << std::endl;
                        }
                    }
                    ssfile << tmpss.str();
                    tmpss.clear();
                    tmpss.str("");
                    chunkscomplete += 1.0e0;
                    completefrac = chunkscomplete/Tchunks;
                    qtsignals->EmitFunctionProgress(completefrac, descriptor);
                    qtsignals->EmitFunctionProgress(completefrac);
                }
            }
            ssfile.close();
        } else {
            std::cout << "ERROR!  Cannot open file: " << filename << " for data-output!" << std::endl;
        }
    }
}


template <class T>
void UniformVolume<T>::VTKWriteBinary()
{
    /* Calculate the number of points in the data volume, and compare to the maximum value expressible
     * using type 'int'.  If too many data points are in the volume, automatically write the data as
     * multiple VTK files. */
    std::string filestem_save(filenamestem);
    size_t first_slice = 0;
    size_t num_slices = vslices;
    size_t remainder = 0;
    size_t ntasks = 1;
    float slice_min = (float)zmin;
    float dz = 0.0e0;

    size_t max_points_int = 2147483647;     /* 2^31 - 1 */

    size_t max_slices = max_points_int/(vrows*vcols);

    /* If too many points exist, re-calculate the number of sub-volume-files required and their size.  The
     * last sub-volume-file may be smaller than the others. */
    if(vslices > max_slices){
        size_t loop_limit = vslices;
        size_t loop_counter = 0;
        while(num_slices > max_slices && loop_counter != loop_limit){
            /* Calculate the number of sub-volume-files required. */
            ntasks = vslices/max_slices;    /* INTEGER MATH! */

            /* Add additional task if remainder is non-zero. */
            if(vslices % max_slices != 0){
                ntasks++;
            }

            /* Calculate number of slices per task. */
            num_slices = vslices/ntasks;    /* INTEGER MATH! */

            /* Calculate remainder. */
            remainder = vslices % ntasks;

            loop_counter++;
        }
    }

    /* Calculate z-coordinate increment between each sub-volume-file. */
    dz = (float)num_slices*zspacing;

    /* Loop through the number of tasks and write the sub-volume-file. */
    for(size_t i=0; i<ntasks; i++){

        /* If multiple tasks are to be used, append "_partXX" to the current filename stem. */
        if(ntasks > 1){
            std::stringstream tmpss;
            tmpss << filenamestem << "_part" << StringManip::FormattedNumber<size_t>(i,2,0,'0',false);
            filenamestem = tmpss.str();
        }

        slice_min = (float)zmin + dz*(float)i;
        first_slice = num_slices*i;

        /* If this is the last iteration and the remainder is non-zero, change the number of slices to
         * be equal to the remainder since the final iteration may be smaller than the others. */
        if(i == ntasks-1 && remainder != 0){
            num_slices = remainder;
        }

        /* Build a string to update the user-interface to indicate which file is being written, emit the
         * signal to update the interface, and call the function to write the specified portion of the data. */
        std::stringstream tmpss;
        tmpss << "Writing volume file " << i+1 << " of " << ntasks;
        qtsignals->EmitFunctionDesc2(tmpss.str());

        UniformVolume<T>::VTKWriteBinaryPartial(first_slice, num_slices, slice_min);
    }

    /* Reset filename stem to original value. */
    filenamestem = filestem_save;
}


template <class T>
void UniformVolume<T>::VTKWriteBinaryBitFlip()
{
    /* Calculate the number of points in the data volume, and compare to the maximum value expressible
     * using type 'int'.  If too many data points are in the volume, automatically write the data as
     * multiple VTK files. */
    std::string filestem_save(filenamestem);
    size_t first_slice = 0;
    size_t num_slices = vslices;
    size_t remainder = 0;
    size_t ntasks = 1;
    float slice_min = (float)zmin;
    float dz = 0.0e0;

    size_t max_points_int = 2147483647;     /* 2^31 - 1 -- Maximum value for signed integer. */

    size_t max_slices = max_points_int/(vrows*vcols);

    /* If too many points exist, re-calculate the number of sub-volume-files required and their size.  The
     * last sub-volume-file may be smaller than the others. */
    if(vslices > max_slices){
        size_t loop_limit = vslices;
        size_t loop_counter = 0;
        while(num_slices > max_slices && loop_counter != loop_limit){
            /* Calculate the number of sub-volume-files required. */
            ntasks = vslices/max_slices;    /* INTEGER MATH! */

            /* Add additional task if remainder is non-zero. */
            if(vslices % max_slices != 0){
                ntasks++;
            }

            /* Calculate number of slices per task. */
            num_slices = vslices/ntasks;    /* INTEGER MATH! */

            /* Calculate remainder. */
            remainder = vslices % ntasks;

            loop_counter++;
        }
    }

    /* Calculate z-coordinate increment between each sub-volume-file. */
    dz = (float)num_slices*zspacing;

    /* Loop through the number of tasks and write the sub-volume-file. */
    for(size_t i=0; i<ntasks; i++){

        /* If multiple tasks are to be used, append "_partXX" to the current filename stem. */
        if(ntasks > 1){
            std::stringstream tmpss;
            tmpss << filestem_save << "_part" << StringManip::FormattedNumber<size_t>(i,2,0,'0',false);
            filenamestem = tmpss.str();
        }

        slice_min = (float)zmin + dz*(float)i;
        first_slice = num_slices*i;

        /* If this is the last iteration and the remainder is non-zero, change the number of slices to
         * be equal to the remainder since the final iteration may be smaller than the others. */
        if(i == ntasks-1 && remainder != 0){
            num_slices = remainder;
        }

        /* Build a string to update the user-interface to indicate which file is being written, emit the
         * signal to update the interface, and call the function to write the specified portion of the data. */
        std::stringstream tmpss;
        tmpss << "Writing volume file " << i+1 << " of " << ntasks;
        qtsignals->EmitFunctionDesc2(tmpss.str());

        UniformVolume<T>::VTKWriteBinaryBitFlipPartial(first_slice, num_slices, slice_min);
    }

    /* Reset filename stem to original value. */
    filenamestem = filestem_save;
}


template <class T>
void UniformVolume<T>::VTKWriteBinaryPartial(const size_t first_slice, const size_t num_slices, const float slice_min)
{
    std::string filename;
    std::fstream ssfile;
    if(!UniformVolume<T>::IsBigEndian()){
        filename = outputdir + "/" + filenamestem + ".vtk";
    } else {
        filename = outputdir + "/" + filenamestem + "_BigEndian.vtk";
    }
    ssfile.open(filename.c_str(),std::ios::out|std::ios::binary);

    /*
      Variables needed for tracking progress
     */
    float Tchunks = (float)num_slices*((float)nscalars + (float)nvectors);
    float chunkscomplete = 0.0e0;
    float completefrac = chunkscomplete/Tchunks;
    std::string descriptor("Writing Binary VTK File");

    /*
      Emit signals for the beginning of the file-writing process
     */
    qtsignals->EmitFunctionDesc(descriptor);
    qtsignals->EmitFunctionProgress(completefrac);


    /* Buffer for text values to be written. */
    char buffer[2048];
    int size;

    T tmpfloat;

    /* Use std::stringstream to allow cross-compiler handling of size_t array sizes without
     * requiring separate format specifiers in 'sprintf()'. */
    std::stringstream tmpss;


    if(ssfile.is_open() == true){
        size = sprintf(buffer,"# vtk DataFile Version 2.0\n");
        ssfile.write(buffer,size);
        size = sprintf(buffer,"%s\n",filenamestem.c_str());
        ssfile.write(buffer,size);
        size = sprintf(buffer,"BINARY\n");
        ssfile.write(buffer,size);

        if(imageoutput == true){
            size = sprintf(buffer,"DATASET STRUCTURED_POINTS\n");
            ssfile.write(buffer,size);
            tmpss.clear();  tmpss.str("");
            tmpss << "DIMENSIONS " << vcols << " " << vrows << " " << num_slices << std::endl;
            size = sprintf(buffer,"%s",tmpss.str().c_str());
            ssfile.write(buffer,size);
            size = sprintf(buffer,"ORIGIN %g %g %g\n",xmin,ymin,slice_min);
            ssfile.write(buffer,size);
            size = sprintf(buffer,"SPACING %g %g %g\n",xspacing,yspacing,zspacing);
            ssfile.write(buffer,size);
        }

        if(rectoutput == true){
            size = sprintf(buffer,"DATASET RECTILINEAR_GRID\n");
            ssfile.write(buffer,size);

            tmpss.clear();  tmpss.str("");
            tmpss << "DIMENSIONS " << vcols << " " << vrows << " " << num_slices << std::endl;
            size = sprintf(buffer,"%s",tmpss.str().c_str());
            ssfile.write(buffer,size);

            tmpss.clear();  tmpss.str("");
            tmpss << "X_COORDINATES " << vcols << " float" << std::endl;
            size = sprintf(buffer,"%s",tmpss.str().c_str());
            ssfile.write(buffer,size);
            for(size_t i=0; i<vcols; i++){
                tmpfloat = (float)i*(float)xspacing;
                ssfile.write((char*)&tmpfloat,sizeof(float));
            }

            tmpss.clear();  tmpss.str("");
            tmpss << "Y_COORDINATES " << vrows << " float" << std::endl;
            size = sprintf(buffer,"%s",tmpss.str().c_str());
            ssfile.write(buffer,size);
            for(size_t i=0; i<vrows; i++){
                tmpfloat = (float)i*(float)yspacing;
                ssfile.write((char*)&tmpfloat,sizeof(float));
            }

            tmpss.clear();  tmpss.str("");
            tmpss << "Z_COORDINATES " << num_slices << " float" << std::endl;
            size = sprintf(buffer,"%s",tmpss.str().c_str());
            ssfile.write(buffer,size);
            for(size_t i=0; i<num_slices; i++){
                tmpfloat = (float)i*(float)zspacing + slice_min;
                ssfile.write((char*)&tmpfloat,sizeof(float));
            }
        }

        size_t nchunks = num_slices*(nscalars + nvectors);
        size_t chunks_processed = 0;

        tmpss.clear();  tmpss.str("");
        tmpss << "POINT_DATA " << vcols*vrows*num_slices << std::endl;
        size = sprintf(buffer,"%s",tmpss.str().c_str());
        ssfile.write(buffer,size);

        for(size_t s=0; s<nscalars; s++){

            descriptor = "Writing Binary VTK File - scalars '" + ScalarQuantityName(s) + "'";
            qtsignals->EmitFunctionDesc(descriptor);

            size = sprintf(buffer,"SCALARS %s %s 1\n",scalar_names(s)->substr().c_str(),dtypename.c_str());
            ssfile.write(buffer,size);
            size = sprintf(buffer,"LOOKUP_TABLE default\n");
            ssfile.write(buffer,size);
            for(size_t k=0; k<num_slices; k++){
                size_t kk = k + first_slice;
                for(size_t i=0; i<vrows; i++){
                    for(size_t j=0; j<vcols; j++){
                        tmpfloat = pscalars(s)->operator ()(i,j,kk);
                        ssfile.write((char*)&tmpfloat,sizeof(T));
                    }
                }
                chunks_processed++;
                qtsignals->EmitFunctionProgress((float)chunks_processed/(float)nchunks);
            }
        }

        for(size_t v=0; v<nvectors; v++){

            descriptor = "Writing Binary VTK File - vectors '" + VectorQuantityName(v) + "'";
            qtsignals->EmitFunctionDesc(descriptor);

            std::string name(vector_names(v)->substr());
            size_t n = pvectors(v)->GetDim(4);
            if(n > 3){
                std::cerr << "WARNING: VTK file format requires 3-component vectors." << std::endl;
                std::cerr << "         Vector quantity " << name << " has " << n << " components." << std::endl;
                std::cerr << "         File will be written with all components, which may cause" << std::endl;
                std::cerr << "         problems with some file readers." << std::endl;
            }
            size = sprintf(buffer,"VECTORS %s %s\n",vector_names(v)->substr().c_str(),dtypename.c_str());
            ssfile.write(buffer,size);
            for(size_t k=0; k<num_slices; k++){
                size_t kk = k + first_slice;
                for(size_t i=0; i<vrows; i++){
                    for(size_t j=0; j<vcols; j++){
                        for(size_t l=0; l<n; l++){ /* Write components of vectors. */
                            tmpfloat = pvectors(v)->operator ()(i,j,kk,l);
                            ssfile.write((char*)&tmpfloat,sizeof(T));
                        }
                        for(size_t l=n; l<3; l++){ /* Add 0 for extra dimension(s) if necessary. */
                            tmpfloat = (T)0.0e0;
                            ssfile.write((char*)&tmpfloat,sizeof(T));
                        }
                    }
                }
                chunks_processed++;
                qtsignals->EmitFunctionProgress((float)chunks_processed/(float)nchunks);
            }
        }
        ssfile.close();
    } else {
        std::cout << "ERROR!  Cannot open file: " << filename << " for data-output!" << std::endl;
    }

    descriptor = "";
    qtsignals->EmitFunctionDesc(descriptor);
    qtsignals->EmitFunctionProgress(0.0e0);
}


template <class T>
void UniformVolume<T>::VTKWriteBinaryBitFlipPartial(const size_t first_slice, const size_t num_slices, const float slice_min)
{
    std::string filename;
    std::fstream ssfile;
    if(UniformVolume<T>::IsBigEndian()){
        filename = outputdir + "/" + filenamestem + ".vtk";     /* Flipping bits will make output Little Endian. */
    } else {
        filename = outputdir + "/" + filenamestem + "_BigEndian.vtk";
    }
    ssfile.open(filename.c_str(),std::ios::out|std::ios::binary);

    /*
      Variables needed for tracking progress
     */
    float Tchunks = (float)num_slices*((float)nscalars + (float)nvectors);
    float chunkscomplete = 0.0e0;
    float completefrac = chunkscomplete/Tchunks;
    std::string descriptor("Writing Binary VTK File");

    /*
      Emit signals for the beginning of the file-writing process
     */
    qtsignals->EmitFunctionDesc(descriptor);
    qtsignals->EmitFunctionProgress(completefrac);


    /* Buffer for text values to be written. */
    char buffer[2048];
    int size;

    T flippedfloat;

    /* Use std::stringstream to allow cross-compiler handling of size_t array sizes without
     * requiring separate format specifiers in 'sprintf()'. */
    std::stringstream tmpss;


    if(ssfile.is_open() == true){
        size = sprintf(buffer,"# vtk DataFile Version 2.0\n");
        ssfile.write(buffer,size);
        size = sprintf(buffer,"%s\n",filenamestem.c_str());
        ssfile.write(buffer,size);
        size = sprintf(buffer,"BINARY\n");
        ssfile.write(buffer,size);

        if(imageoutput == true){
            size = sprintf(buffer,"DATASET STRUCTURED_POINTS\n");
            ssfile.write(buffer,size);
            tmpss.clear();  tmpss.str("");
            tmpss << "DIMENSIONS " << vcols << " " << vrows << " " << num_slices << std::endl;
            size = sprintf(buffer,"%s",tmpss.str().c_str());
            ssfile.write(buffer,size);
            size = sprintf(buffer,"ORIGIN %g %g %g\n",xmin,ymin,slice_min);
            ssfile.write(buffer,size);
            size = sprintf(buffer,"SPACING %g %g %g\n",xspacing,yspacing,zspacing);
            ssfile.write(buffer,size);
        }

        if(rectoutput == true){
            size = sprintf(buffer,"DATASET RECTILINEAR_GRID\n");
            ssfile.write(buffer,size);

            tmpss.clear();  tmpss.str("");
            tmpss << "DIMENSIONS " << vcols << " " << vrows << " " << num_slices << std::endl;
            size = sprintf(buffer,"%s",tmpss.str().c_str());
            ssfile.write(buffer,size);

            tmpss.clear();  tmpss.str("");
            tmpss << "X_COORDINATES " << vcols << " float" << std::endl;
            size = sprintf(buffer,"%s",tmpss.str().c_str());
            ssfile.write(buffer,size);
            for(size_t i=0; i<vcols; i++){
                flippedfloat = (float)i*(float)xspacing;
                UniformVolume<T>::ByteSwap(&flippedfloat, sizeof(float));
                ssfile.write(reinterpret_cast<char*>(&flippedfloat),sizeof(float));
            }

            tmpss.clear();  tmpss.str("");
            tmpss << "Y_COORDINATES " << vrows << " float" << std::endl;
            size = sprintf(buffer,"%s",tmpss.str().c_str());
            ssfile.write(buffer,size);
            for(size_t i=0; i<vrows; i++){
                flippedfloat = (float)i*(float)yspacing;
                UniformVolume<T>::ByteSwap(&flippedfloat, sizeof(float));
                ssfile.write(reinterpret_cast<char*>(&flippedfloat),sizeof(float));
            }

            tmpss.clear();  tmpss.str("");
            tmpss << "Z_COORDINATES " << num_slices << " float" << std::endl;
            size = sprintf(buffer,"%s",tmpss.str().c_str());
            ssfile.write(buffer,size);
            for(size_t i=0; i<num_slices; i++){
                flippedfloat = (float)i*(float)zspacing + slice_min;
                UniformVolume<T>::ByteSwap(&flippedfloat, sizeof(float));
                ssfile.write(reinterpret_cast<char*>(&flippedfloat),sizeof(float));
            }
        }

        size_t nchunks = num_slices*(nscalars + nvectors);
        size_t chunks_processed = 0;

        tmpss.clear();  tmpss.str("");
        tmpss << "POINT_DATA " << vcols*vrows*num_slices<< std::endl;
        size = sprintf(buffer,"%s",tmpss.str().c_str());
        ssfile.write(buffer,size);
        for(size_t s=0; s<nscalars; s++){

            descriptor = "Writing Binary VTK File - scalars '" + ScalarQuantityName(s) + "'";
            qtsignals->EmitFunctionDesc(descriptor);

            size = sprintf(buffer,"SCALARS %s %s 1\n",scalar_names(s)->substr().c_str(),dtypename.c_str());
            ssfile.write(buffer,size);
            size = sprintf(buffer,"LOOKUP_TABLE default\n");
            ssfile.write(buffer,size);
            for(size_t k=0; k<num_slices; k++){
                size_t kk = k + first_slice;
                for(size_t i=0; i<vrows; i++){
                    for(size_t j=0; j<vcols; j++){
                        flippedfloat = pscalars(s)->operator ()(i,j,kk);
                        UniformVolume<T>::ByteSwap(&flippedfloat, sizeof(T));
                        ssfile.write(reinterpret_cast<char*>(&flippedfloat),sizeof(T));
                    }
                }
                chunks_processed++;
                qtsignals->EmitFunctionProgress((float)chunks_processed/(float)nchunks);
            }
        }

        for(size_t v=0; v<nvectors; v++){
            std::string name(vector_names(v)->substr());
            size_t n = pvectors(v)->GetDim(4);
            if(n > 3){
                std::cerr << "WARNING: VTK file format requires 3-component vectors." << std::endl;
                std::cerr << "         Vector quantity " << name << " has " << n << " components." << std::endl;
                std::cerr << "         File will be written with all components, which may cause" << std::endl;
                std::cerr << "         problems with some file readers." << std::endl;
            }

            descriptor = "Writing Binary VTK File - vectors '" + VectorQuantityName(v) + "'";
            qtsignals->EmitFunctionDesc(descriptor);

            size = sprintf(buffer,"VECTORS %s %s\n",vector_names(v)->substr().c_str(),dtypename.c_str());
            ssfile.write(buffer,size);
            for(size_t k=0; k<num_slices; k++){
                size_t kk = k + first_slice;
                for(size_t i=0; i<vrows; i++){
                    for(size_t j=0; j<vcols; j++){
                        for(size_t l=0; l<n; l++){ /* Write components of vectors. */
                            flippedfloat = pvectors(v)->operator ()(i,j,kk,l);
                            UniformVolume<T>::ByteSwap(&flippedfloat, sizeof(T));
                            ssfile.write(reinterpret_cast<char*>(&flippedfloat),sizeof(T));
                        }
                        for(size_t l=n; l<3; l++){ /* Add 0 for extra dimension(s) if necessary. */
                            flippedfloat = (T)0.0e0;
                            UniformVolume<T>::ByteSwap(&flippedfloat, sizeof(T));
                            ssfile.write(reinterpret_cast<char*>(&flippedfloat),sizeof(T));
                        }
                    }
                }
                chunks_processed++;
                qtsignals->EmitFunctionProgress((float)chunks_processed/(float)nchunks);
            }
        }
        ssfile.close();
    } else {
        std::cout << "ERROR!  Cannot open file: " << filename << " for data-output!" << std::endl;
    }

    descriptor = "";
    qtsignals->EmitFunctionDesc(descriptor);
    qtsignals->EmitFunctionProgress(0.0e0);
}


template <class T>
void UniformVolume<T>::VTKWriteBinaryBigEndian()
{

    if(false){
        UniformVolume<T>::VTKWriteImageData();      /* Currently some sort of bug in here. */
    } else {

        writebigendian = true;

        if(UniformVolume<T>::IsBigEndian()){
            UniformVolume<T>::VTKWriteBinary();
        } else {
            UniformVolume<T>::VTKWriteBinaryBitFlip();
        }

        writebigendian = false;
    }
}


template <class T>
void UniformVolume<T>::ByteSwap(void *value, size_t numbytes)
{
    unsigned char *memp = reinterpret_cast<unsigned char*>(value);
    std::reverse(memp, memp + numbytes);
}


template <class T>
bool UniformVolume<T>::IsBigEndian()
{
    /* From http://http://www.codeproject.com/Articles/4804/Basic-concepts-on-Endianness */
    short word = 0x4321;
    if((*(char *)& word) != 0x21 ){
        return true;
    } else  {
        return false;
    }
}


template <class T>
void UniformVolume<T>::VOLWrite(const int qty, const bool allowScaling, const T minRange)
{
    int retval = -1;

    std::string filename;
    filename = outputdir + "/" + filenamestem + ".vol";

    /* Multiply all values by a scaling factor to ensure that the data spans multiple integers.
     *
     * This scaling is done to compensate for a bug in the behavior of the CNDE 3D Visualization
     * program.  When adjusting the lookup table in 3D Viz, it is thought that a cast-to-int is
     * applied to the data values when determining the appropriate grayscale to be used in the
     * display.  If data values span a range less than an integer, no remapping is possible because
     * all values get cast to the same integer.  Similarly, if very few integers are spanned, the
     * mapping result will appear very poor due to a small number of integers which are cast to
     * grayscales.  By applying a linear scaling to the data, the numeric values can be made to
     * span a sufficiently-large range of values to prevent the cast-to-int operation from
     * adversely affecting the color-mapping.
     *
     * Note that this addresses a graphics/visualization problem in 3D Viz.  The application does
     * accurately read-in the floating-point data. */
    std::string desc;
    size_t loc1, loc2, loc3;
    T maxval = pscalars(qty)->MaxVal(loc1,loc2,loc3);
    T minval = pscalars(qty)->MinVal(loc1,loc2,loc3);
    T range = maxval - minval;
    bool scalingPerformed = false;

    if(range < minRange && allowScaling){
        /* Apply scaling to data, only if data range must be increased AND scaling is allowed. */
        scalingPerformed = true;

        desc = "Scaling Data For VOL File";
        qtsignals->EmitFunctionDesc(desc);
        qtsignals->EmitFunctionProgress(0.0e0);

        for(size_t k=0; k<vslices; k++){
            for(size_t i=0; i<vrows; i++){
                for(size_t j=0; j<vcols; j++){
                    pscalars(qty)->operator ()(i,j,k) = minRange*(pscalars(qty)->operator ()(i,j,k) - minval)/range;
                }
            }
            qtsignals->EmitFunctionProgress((float)k/(float)vslices);
        }
        desc = "";
        qtsignals->EmitFunctionDesc(desc);
        qtsignals->EmitFunctionProgress(0.0e0);
    }



    retval = UniformVolume<T>::WriteVOLFile(filename,vcols,vrows,vslices,
                                           &pscalars(qty)->operator ()((size_t)0,(size_t)0,(size_t)0));

    if(scalingPerformed){
        /* Undo scaling to data, if scaling was done before
         *
         * Note that IEEE 754 (Standard for Floating Point Arithmetic) requires addition, subtraction,
         * multiplication, and division to be within 0.5 ULP (unit in last place / unit of least precision).
         * This means the operations will produce numbers within 0.5 ULP of the mathematically exact
         * result.  Some operations, such as subtraction, can cause loss-of-significant-digits, though,
         * so the result of the following un-scaling must be treated carefully as some precision may be lost
         * through the scaling and un-scaling processes. */

        desc = "Undoing Scaling for VOL Output";
        qtsignals->EmitFunctionDesc(desc);
        for(size_t k=0; k<vslices; k++){
            for(size_t i=0; i<vrows; i++){
                for(size_t j=0; j<vcols; j++){
                    pscalars(qty)->operator ()(i,j,k) = pscalars(qty)->operator ()(i,j,k)*range/minRange + minval;
                }
            }
            qtsignals->EmitFunctionProgress((float)k/(float)vslices);
        }

        desc = "";
        qtsignals->EmitFunctionDesc(desc);
    }

    if(retval == 0){
        /*
        printf("Data written successfully to %s\n",filename.c_str());
        */
    } else {
        std::cerr << "Data output to VOL file FAILED" << std::endl;
    }

    desc = "";
    qtsignals->EmitFunctionDesc(desc);
    qtsignals->EmitFunctionProgress(0.0e0);
}


template <class T>
void UniformVolume<T>::ReadFile(const std::string filename, const bool isBigEndian)
{
    std::string fileext;
    fileext = StringManip::DetermFileExt(filename);
    if(fileext == "vol" || fileext == "VOL"){
        UniformVolume<T>::ReadVOLFile(filename);
    }
    if(fileext == "vtk" || fileext == "VTK" || fileext == "vti" || fileext == "VTI" || fileext == "pvti" || fileext == "PVTI"){
        UniformVolume<T>::ReadVTKFile(filename,isBigEndian);
    }
    if(fileext == "xmf" || fileext == "XMF" || fileext == "xdmf" || fileext == "XDMF"){
        UniformVolume<T>::ReadXDMFFile(filename);
    }

}


template <class T>
T UniformVolume<T>::Mean(const int scalarqty, const int vectorqty, const int vectorcomp) const
{
    T retval = (T)0.0e0;
    if(scalarqty >= 0){
        retval = pscalars(scalarqty)->Mean();
    }

    if(vectorqty >= 0){
        Array3D<T> tmparray(vrows,vcols,vslices);
        if(vectorcomp >= 0){
            for(int k=0; k<vslices; k++){
                for(int i=0; i<vrows; i++){
                    for(int j=0; j<vcols; j++){
                        tmparray(i,j,k) = pvectors(vectorqty)->operator ()(i,j,k,vectorcomp);
                    }
                }
            }
        } else {
            int n = pvectors(vectorqty)->GetDim(3);
            for(int k=0; k<vslices; k++){
                for(int i=0; i<vrows; i++){
                    for(int j=0; j<vcols; j++){
                        for(int l=0; l<n; l++){
                            tmparray(i,j,k) += pvectors(vectorqty)->operator ()(i,j,k,l)*
                                    pvectors(vectorqty)->operator ()(i,j,k,l);
                        }
                        tmparray(i,j,k) = std::sqrt(tmparray(i,j,k));
                    }
                }
            }
        }
        retval = tmparray.Mean();
    }

    return retval;
}


template <class T>
T UniformVolume<T>::Median(const int scalarqty, const int vectorqty, const int vectorcomp) const
{
    T retval = (T)0.0e0;
    if(scalarqty >= 0){
        retval = pscalars(scalarqty)->MedianVal();
    }

    if(vectorqty >= 0){
        Array3D<T> tmparray(vrows,vcols,vslices);
        if(vectorcomp >= 0){
            for(int k=0; k<vslices; k++){
                for(int i=0; i<vrows; i++){
                    for(int j=0; j<vcols; j++){
                        tmparray(i,j,k) = pvectors(i)->operator ()(i,j,k,vectorcomp);
                    }
                }
            }
        } else {
            int n = pvectors(vectorqty)->GetDim(3);
            for(int k=0; k<vslices; k++){
                for(int i=0; i<vrows; i++){
                    for(int j=0; j<vcols; j++){
                        for(int l=0; l<n; l++){
                            tmparray(i,j,k) += pvectors(vectorqty)->operator ()(i,j,k,l)*
                                    pvectors(vectorqty)->operator ()(i,j,k,l);
                        }
                        tmparray(i,j,k) = std::sqrt(tmparray(i,j,k));
                    }
                }
            }
        }
        retval = tmparray.MedianVal();
    }

    return retval;
}


template <class T>
T UniformVolume<T>::Variance(const int scalarqty, const int vectorqty, const int vectorcomp) const
{
    T retval = (T)0.0e0;
    if(scalarqty >= 0){
        retval = pscalars(scalarqty)->Variance();
    }

    if(vectorqty >= 0){
        Array3D<T> tmparray(vrows,vcols,vslices);
        if(vectorcomp >= 0){
            for(int k=0; k<vslices; k++){
                for(int i=0; i<vrows; i++){
                    for(int j=0; j<vcols; j++){
                        tmparray(i,j,k) = pvectors(vectorqty)->operator ()(i,j,k,vectorcomp);
                    }
                }
            }
        } else {
            int n = pvectors(vectorqty)->GetDim(3);
            for(int k=0; k<vslices; k++){
                for(int i=0; i<vrows; i++){
                    for(int j=0; j<vcols; j++){
                        for(int l=0; l<n; l++){
                            tmparray(i,j,k) += pvectors(vectorqty)->operator ()(i,j,k,l)*
                                    pvectors(vectorqty)->operator ()(i,j,k,l);
                        }
                        tmparray(i,j,k) = std::sqrt(tmparray(i,j,k));
                    }
                }
            }
        }
        retval = tmparray.Variance();
    }

    return retval;
}


template <class T>
T UniformVolume<T>::MinVal(const int scalarqty, const int vectorqty, const int vectorcomp,
         int &loc_x, int &loc_y, int &loc_z) const
{
    T retval = (T)0.0e0;
    if(scalarqty >= 0){
        retval = pscalars(scalarqty)->MinVal(loc_y,loc_x,loc_z);
    }

    if(vectorqty >= 0){
        Array3D<T> tmparray(vrows,vcols,vslices);
        if(vectorcomp >= 0){
            for(int k=0; k<vslices; k++){
                for(int i=0; i<vrows; i++){
                    for(int j=0; j<vcols; j++){
                        tmparray(i,j,k) = pvectors(vectorqty)->operator ()(i,j,k,vectorcomp);
                    }
                }
            }
        } else {
            int n = pvectors(vectorqty)->GetDim(3);
            for(int k=0; k<vslices; k++){
                for(int i=0; i<vrows; i++){
                    for(int j=0; j<vcols; j++){
                        for(int l=0; l<n; l++){
                            tmparray(i,j,k) += pvectors(vectorqty)->operator ()(i,j,k,l)*
                                    pvectors(vectorqty)->operator ()(i,j,k,l);
                        }
                        tmparray(i,j,k) = std::sqrt(tmparray(i,j,k));
                    }
                }
            }
        }
        retval = tmparray.MinVal(loc_y,loc_x,loc_z);
    }

    return retval;
}


template <class T>
T UniformVolume<T>::MaxVal(const int scalarqty, const int vectorqty, const int vectorcomp,
         int &loc_x, int &loc_y, int &loc_z) const
{
    T retval = (T)0.0e0;
    if(scalarqty >= 0){
        retval = pscalars(scalarqty)->MaxVal(loc_y,loc_x,loc_z);
    }

    if(vectorqty >= 0){
        Array3D<T> tmparray(vrows,vcols,vslices);
        if(vectorcomp >= 0){
            for(int k=0; k<vslices; k++){
                for(int i=0; i<vrows; i++){
                    for(int j=0; j<vcols; j++){
                        tmparray(i,j,k) = pvectors(vectorqty)->operator ()(i,j,k,vectorcomp);
                    }
                }
            }
        } else {
            int n = pvectors(vectorqty)->GetDim(3);
            for(int k=0; k<vslices; k++){
                for(int i=0; i<vrows; i++){
                    for(int j=0; j<vcols; j++){
                        for(int l=0; l<n; l++){
                            tmparray(i,j,k) += pvectors(vectorqty)->operator ()(i,j,k,l)*
                                    pvectors(vectorqty)->operator ()(i,j,k,l);
                        }
                        tmparray(i,j,k) = std::sqrt(tmparray(i,j,k));
                    }
                }
            }
        }
        retval = tmparray.MaxVal(loc_y,loc_x,loc_z);
    }

    return retval;
}


template <class T>
size_t UniformVolume<T>::MemoryRequired() const
{
    size_t retval = 0;


    return retval;
}



template <class T>
void UniformVolume<T>::PerformCalculations()
{
    UniformVolume<T>::setScalarDestinationArray(NULL, (size_t)0);
    UniformVolume<T>::ReadFile(conversion_options.input_file,conversion_options.isBigEndian);

    if(conversion_options.output_VTKImageData){
        UniformVolume<T>::setVTKImageDataOutput(true);
        UniformVolume<T>::setVTKRectDataOutput(false);
        UniformVolume<T>::VTKWriteBinaryBigEndian();
    }

    if(conversion_options.output_VTKRectGrid){
        UniformVolume<T>::setVTKImageDataOutput(false);
        UniformVolume<T>::setVTKRectDataOutput(true);
        UniformVolume<T>::VTKWriteBinaryBigEndian();
    }

    if(conversion_options.output_CNDEVOL){
        UniformVolume<T>::VOLWrite(0,conversion_options.output_CNDEVOL_allowScaling,conversion_options.output_CNDEVOL_scalingRange);
    }

    if(conversion_options.output_xdmf){
        UniformVolume<T>::WriteXdmf(0);     /* Hard-code no compression. */
    }

    qtsignals->EmitDone();
    UniformVolume<T>::EmitDone();
}


template <class T>
void UniformVolume<T>::VTKReadLegacyBinary(std::fstream &file, const bool isBigEndian)
{
    char linedata[256];
    std::string sval1("");
    std::string sval2("");
    std::string sval3("");
    long long ival1 = 0;
    long long ival2 = 0;
    long long ival3 = 0;
    T fval1 = (T)0.0e0;
    T fval2 = (T)0.0e0;
    T fval3 = (T)0.0e0;
    bool isImageData = false;
    bool isRectGrid = false;
    bool thisBigEndian = UniformVolume<T>::IsBigEndian();
    std::string desc("");
    std::stringstream sstmp;

    /* Read dataset type. */
    file.getline(linedata,256);
    sstmp.str(linedata);
    sstmp >> sval1 >> sval2;
    if(sval1 != "DATASET"){
        std::cerr << "UniformVolume::VTKReadLegacyBinary() - Dataset keyword not recognized" << std::endl;
        std::cerr << "                                         Keyword read: " << sval1 << std::endl;
        return;
    }

    if(sval2 == "STRUCTURED_POINTS"){
        /* Data is Image Data format. */
        isImageData = true;
        isRectGrid = false;
    }

    if(sval2 == "RECTILINEAR_GRID"){
        /* Data is Rectilinear Grid format. */
        isImageData = false;
        isRectGrid = true;
    }
    sval1 = "";
    sval2 = "";


    /* Read dimensions. */
    file.getline(linedata,256);
    sstmp.clear(); sstmp.str(linedata);
    sstmp >> sval1 >> ival1 >> ival2 >> ival3;
    if(sval1 != "DIMENSIONS"){
        std::cerr << "UniformVolume::VTKReadLegacyBinary() - Dimensions keyword not recognized" << std::endl;
        std::cerr << "                                         Keyword read: " << sval1 << std::endl;
        return;
    }
    if(ival1 > 0 && ival2 > 0 && ival3 > 0){
        UniformVolume<T>::ResetResolution(ival2, ival1, ival3, (T)0);
    } else {
        std::cerr << "UniformVolume::VTKReadLegacyBinary() - Dimensions not read properly" << std::endl;
        std::cerr << "                                         Dimensions read: " << ival1 << " x " << ival2
                                                                          << " x " << ival3 << std::endl;
        return;
    }


    /* Read remainder of header -- Image data. */
    if(isImageData){
        int linecount = 0;
        T o1 = 0.0e0;
        T o2 = 0.0e0;
        T o3 = 0.0e0;
        T s1 = 0.0e0;
        T s2 = 0.0e0;
        T s3 = 0.0e0;
        while(linecount < 2){

            /* Read ORIGIN / SPACING and floating-point values. */
            file.getline(linedata,256);
            sstmp.clear(); sstmp.str(linedata);
            sstmp >> sval1 >> fval1 >> fval2 >> fval3;


            /* Set origin info. */
            if(sval1 == "ORIGIN"){
                o1 = fval1;
                o2 = fval2;
                o3 = fval3;
            }


            /* Set spacing info. */
            if(sval1 == "SPACING"){
                s1 = fval1;
                s2 = fval2;
                s3 = fval3;
            }

            linecount++;


            /* Only set spatial bounds after both origin and spacing have been read-in. */
            if(linecount == 2){
                UniformVolume<T>::setSpatialExtent(1,-1,o1);
                UniformVolume<T>::setSpatialExtent(0,-1,o2);
                UniformVolume<T>::setSpatialExtent(2,-1,o3);
                UniformVolume<T>::setSpatialExtent(1,1,xmin+s1*vcols);
                UniformVolume<T>::setSpatialExtent(0,1,ymin+s2*vrows);
                UniformVolume<T>::setSpatialExtent(2,1,zmin+s3*vslices);
            }
        }
    }



    /* Read remainder of header -- Rectilinear Grid data. */
    if(isRectGrid){

        /* Read next line and set appropriate coordinates.  Only the first and last value are
         * read-in.  Internal points are determined automatically by knowing the min/max values
         * and number of elements along each direction. */
        for(int i=0; i<3; i++){
            file.getline(linedata,256);
            sstmp.clear(); sstmp.str(linedata);
            sstmp >> sval1 >> ival1 >> sval2;

            float ffloat = (float)0.0e0;
            double dfloat = (double)0.0e0;
            bool isfloat = false;
            bool isdouble = false;
            size_t datasize = 0;

            if(sval2 == "float"){
                datasize = sizeof(float);
                isfloat = true;
                isdouble = false;
            }
            if(sval2 == "double"){
                datasize = sizeof(double);
                isfloat = false;
                isdouble = true;
            }


            if(sval1 == "X_COORDINATES"){
                if(isfloat){
                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                    UniformVolume<T>::ByteSwap(&ffloat,datasize);
                    UniformVolume<T>::setSpatialExtent(0,-1,(T)ffloat);

                    for(size_t i=0; i<vcols-2; i++){
                        file.read(reinterpret_cast<char*>(&ffloat),datasize);
                    }

                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                    UniformVolume<T>::ByteSwap(&ffloat,datasize);
                    UniformVolume<T>::setSpatialExtent(0,1,(T)ffloat);
                }
                if(isdouble){
                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                    UniformVolume<T>::ByteSwap(&dfloat,datasize);
                    UniformVolume<T>::setSpatialExtent(0,-1,(T)dfloat);

                    for(size_t i=0; i<vcols-2; i++){
                        file.read(reinterpret_cast<char*>(&dfloat),datasize);
                    }

                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                    UniformVolume<T>::ByteSwap(&dfloat,datasize);
                    UniformVolume<T>::setSpatialExtent(0,1,(T)dfloat);
                }
            }

            if(sval1 == "Y_COORDINATES"){
                if(isfloat){
                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                    UniformVolume<T>::ByteSwap(&ffloat,datasize);
                    UniformVolume<T>::setSpatialExtent(1,-1,(T)ffloat);

                    for(size_t i=0; i<vrows-2; i++){
                        file.read(reinterpret_cast<char*>(&ffloat),datasize);
                    }

                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                    UniformVolume<T>::ByteSwap(&ffloat,datasize);
                    UniformVolume<T>::setSpatialExtent(1,1,(T)ffloat);
                }
                if(isdouble){
                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                    UniformVolume<T>::ByteSwap(&dfloat,datasize);
                    UniformVolume<T>::setSpatialExtent(1,-1,(T)dfloat);

                    for(size_t i=0; i<vrows-2; i++){
                        file.read(reinterpret_cast<char*>(&dfloat),datasize);
                    }

                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                    UniformVolume<T>::ByteSwap(&dfloat,datasize);
                    UniformVolume<T>::setSpatialExtent(1,1,(T)dfloat);
                }
            }

            if(sval1 == "Z_COORDINATES"){
                if(isfloat){
                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                    UniformVolume<T>::ByteSwap(&ffloat,datasize);
                    UniformVolume<T>::setSpatialExtent(2,-1,(T)ffloat);

                    for(size_t i=0; i<vslices-2; i++){
                        file.read(reinterpret_cast<char*>(&ffloat),datasize);
                    }

                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                    UniformVolume<T>::ByteSwap(&ffloat,datasize);
                    UniformVolume<T>::setSpatialExtent(2,1,(T)ffloat);
                }
                if(isdouble){
                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                    UniformVolume<T>::ByteSwap(&dfloat,datasize);
                    UniformVolume<T>::setSpatialExtent(2,-1,(T)dfloat);

                    for(size_t i=0; i<vslices-2; i++){
                        file.read(reinterpret_cast<char*>(&dfloat),datasize);
                    }

                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                    UniformVolume<T>::ByteSwap(&dfloat,datasize);
                    UniformVolume<T>::setSpatialExtent(2,1,(T)dfloat);
                }
            }



        } /* End of loop through reading 3 sets of coordinate values. */

    } /* End of if-block for rectilinear grid format. */





    /* Header-reading complete.  Now read data from file. */


    while(file.eof() == false){


        /* Read point data. */
        file.getline(linedata,256);
        sstmp.clear(); sstmp.str(linedata);
        sstmp >> sval1 >> ival1;
        if(sval1 != "POINT_DATA" && sval1 != ""){
            std::cerr << "UniformVolume::VTKReadLegacyBinary() - Point data keyword not recognized" << std::endl;
            std::cerr << "                                         Keyword read: " << sval1 << std::endl;
            std::cerr << "UniformVolume::VTKReadLegacyBinary() - Only POINT_DATA supported...exiting" << std::endl;
            return;
        }
        sval1 = "";


        /* Read keyword. */
        file.getline(linedata,256);
        sstmp.clear(); sstmp.str(linedata);
        sstmp >> sval1;

        if(sval1 == "SCALARS"){
            sval1 = "";
            sstmp >> sval2 >> sval3 >> ival1;

            if(scalar_data == NULL){
                /* Add scalar quantity to this object. */
                UniformVolume<T>::AddScalarQuantity(sval2);
            } else {
                scalar_names(0) = new std::string;
                scalar_names(0)->reserve(qtysize);
                scalar_names(0)->assign(sval2);
                nscalars = 1;
            }

            /* Read lookup table line. */
            file.getline(linedata,256);
            sstmp.clear(); sstmp.str(linedata);
            sstmp >> sval1 >> sval2;


            /* Read data values. */
            size_t datasize = 0;
            bool isfloat = false;
            bool isdouble = false;
            if(sval3 == "float"){
                datasize = sizeof(float);
                isfloat = true;
                isdouble = false;
            }
            if(sval3 == "double"){
                datasize = sizeof(double);
                isfloat = false;
                isdouble = true;
            }

            if((isBigEndian && !thisBigEndian) || (!isBigEndian && thisBigEndian)){
                /* Read point-by-point and bit-flip each point value. */

                if(false){
                    std::cout << "Big Endian System: " << StringManip::BoolToString(UniformVolume<T>::IsBigEndian()) << std::endl;
                    std::cout << "datasize: " << datasize << "     Big Endian Conversion" << std::endl;
                }

                desc = "Reading binary scalars '" + UniformVolume<T>::ScalarQuantityName(nscalars-1) + "'";
                qtsignals->EmitFunctionDesc(desc);

                if(scalar_data == NULL){
                    /* Read data into Array3D object which is a member of this object.  */
                    float ffloat = (float)0.0e0;
                    double dfloat = (double)0.0e0;
                    for(size_t k=0; k<vslices; k++){
                        for(size_t i=0; i<vrows; i++){
                            for(size_t j=0; j<vcols; j++){
                                if(isfloat){
                                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                                    UniformVolume<T>::ByteSwap(&ffloat,datasize);
                                    pscalars(nscalars-1)->operator ()(i,j,k) = (T)ffloat;
                                }
                                if(isdouble){
                                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                                    UniformVolume<T>::ByteSwap(&dfloat,datasize);
                                    pscalars(nscalars-1)->operator ()(i,j,k) = (T)dfloat;
                                }
                            }
                        }
                        qtsignals->EmitFunctionProgress((float)k/(float)vslices);
                    }
                } else {
                    /* Read data into existing array. */
                    scalar_data_points_read = 0;
                    float ffloat = (float)0.0e0;
                    double dfloat = (double)0.0e0;
                    for(size_t k=0; k<vslices; k++){
                        for(size_t i=0; i<vrows; i++){
                            for(size_t j=0; j<vcols; j++){
                                if(isfloat){
                                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                                    UniformVolume<T>::ByteSwap(&ffloat,datasize);
                                    scalar_data[scalar_data_points_read] = (T)ffloat;
                                }
                                if(isdouble){
                                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                                    UniformVolume<T>::ByteSwap(&dfloat,datasize);
                                    scalar_data[scalar_data_points_read] = (T)dfloat;
                                }

                                /* Update number of points read.  If this value is greater than the size of the
                                 * array, inform the user of the error and return. */
                                scalar_data_points_read++;
                                if(scalar_data_points_read > scalar_data_size){
                                    std::cerr << "UniformVolume::VTKReadLegacyBinary() ERROR: Insufficient array size provided." << std::endl;
                                    return;
                                }
                            }
                        }
                        qtsignals->EmitFunctionProgress((float)k/(float)vslices);
                    }
                }


            } else {
                desc = "Reading binary scalars '" + UniformVolume<T>::ScalarQuantityName(nscalars-1) + "'";
                qtsignals->EmitFunctionDesc(desc);

                /* No byte-swapping necessary.  Type-checking still performed, however, in case
                 * datatype in file does not match datatype of this object. */
                if(scalar_data == NULL){
                    float ffloat = (float)0.0e0;
                    double dfloat = (double)0.0e0;
                    for(size_t k=0; k<vslices; k++){
                        for(size_t i=0; i<vrows; i++){
                            for(size_t j=0; j<vcols; j++){
                                if(isfloat){
                                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                                    pscalars(nscalars-1)->operator ()(i,j,k) = (T)ffloat;
                                }
                                if(isdouble){
                                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                                    pscalars(nscalars-1)->operator ()(i,j,k) = (T)dfloat;
                                }
                            }
                        }
                        qtsignals->EmitFunctionProgress((float)k/(float)vslices);
                    }
                } else {
                    /* Read data into existing array. */
                    scalar_data_points_read = 0;
                    float ffloat = (float)0.0e0;
                    double dfloat = (double)0.0e0;
                    for(size_t k=0; k<vslices; k++){
                        for(size_t i=0; i<vrows; i++){
                            for(size_t j=0; j<vcols; j++){
                                if(isfloat){
                                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                                    scalar_data[scalar_data_points_read] = (T)ffloat;
                                }
                                if(isdouble){
                                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                                    scalar_data[scalar_data_points_read] = (T)dfloat;
                                }

                                /* Update number of points read.  If this value is greater than the size of the
                                 * array, inform the user of the error and return. */
                                scalar_data_points_read++;
                                if(scalar_data_points_read > scalar_data_size){
                                    std::cerr << "UniformVolume::VTKReadLegacyBinary() ERROR: Insufficient array size provided." << std::endl;
                                    return;
                                }
                            }
                        }
                        qtsignals->EmitFunctionProgress((float)k/(float)vslices);
                    }
                }
            }

            /* Set pointer to external array to NULL to prevent overwriting the freshly-read data. */
            if(scalar_data != NULL){
                scalar_data = NULL;
                UniformVolume<T>::RemoveAllData();
            }

            desc = "";
            qtsignals->EmitFunctionDesc(desc);
            qtsignals->EmitFunctionProgress(0.0e0);

        }


        if(sval1 == "VECTORS"){
            sval1 = "";
            sstmp >> sval2 >> sval3 >> ival1;

            /* Add vector quantity to this object, assuming 3 components. */
            UniformVolume<T>::AddVectorQuantity(sval2,3);

            /* Read data values. */
            size_t datasize = 0;
            bool isfloat = false;
            bool isdouble = false;
            if(sval3 == "float"){
                datasize = sizeof(float);
                isfloat = true;
                isdouble = false;
            }
            if(sval3 == "double"){
                datasize = sizeof(double);
                isfloat = false;
                isdouble = true;
            }

            if((isBigEndian && !thisBigEndian) || (!isBigEndian && thisBigEndian)){
                /* Read point-by-point and bit-flip each point value. */

                if(false){
                    std::cout << "Big Endian System: " << StringManip::BoolToString(UniformVolume<T>::IsBigEndian()) << std::endl;
                    std::cout << "datasize: " << datasize << "     Big Endian Conversion" << std::endl;
                }

                desc = "Reading binary vectors '" + UniformVolume<T>::VectorQuantityName(nvectors-1) + "'";
                qtsignals->EmitFunctionDesc(desc);

                float ffloat = (float)0.0e0;
                double dfloat = (double)0.0e0;
                for(size_t k=0; k<vslices; k++){
                    for(size_t i=0; i<vrows; i++){
                        for(size_t j=0; j<vcols; j++){
                            for(size_t l=0; l<3; l++){
                                if(isfloat){
                                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                                    UniformVolume<T>::ByteSwap(&ffloat,datasize);
                                    pvectors(nvectors-1)->operator ()(i,j,k,l) = (T)ffloat;
                                }
                                if(isdouble){
                                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                                    UniformVolume<T>::ByteSwap(&dfloat,datasize);
                                    pvectors(nvectors-1)->operator ()(i,j,k,l) = (T)dfloat;
                            }
                            }
                        }
                    }
                    qtsignals->EmitFunctionProgress((float)k/(float)vslices);
                }


            } else {
                desc = "Reading binary vectors '" + UniformVolume<T>::VectorQuantityName(nvectors-1) + "'";
                qtsignals->EmitFunctionDesc(desc);

                /* No byte-swapping necessary.  Type-checking still performed, however, in case
                 * datatype in file does not match datatype of this object. */
                float ffloat = (float)0.0e0;
                double dfloat = (double)0.0e0;
                for(size_t k=0; k<vslices; k++){
                    for(size_t i=0; i<vrows; i++){
                        for(size_t j=0; j<vcols; j++){
                            for(size_t l=0; l<3; l++){
                                if(isfloat){
                                    file.read(reinterpret_cast<char*>(&ffloat),datasize);
                                    pvectors(nvectors-1)->operator ()(i,j,k,l) = (T)ffloat;
                                }
                                if(isdouble){
                                    file.read(reinterpret_cast<char*>(&dfloat),datasize);
                                    pvectors(nvectors-1)->operator ()(i,j,k,l) = (T)dfloat;
                            }
                            }
                        }
                    }
                    qtsignals->EmitFunctionProgress((float)k/(float)vslices);
                }

            }

            desc = "";
            qtsignals->EmitFunctionDesc(desc);
            qtsignals->EmitFunctionProgress(0.0e0);

        }

        if(sval1 == "FIELD"){
            sval1 = "";
            /* Read number of fields.  Discard field name. */
            int nfields = 0;
            sstmp >> sval2 >> nfields;


            for(int f=0; f<nfields; f++){

                /* Read field name. */
                file.getline(linedata,256);
                sstmp.clear(); sstmp.str(linedata);
                sstmp >> sval2 >> ival1 >> ival2 >> sval3;


                if(scalar_data == NULL){
                    /* Add scalar quantity to this object. */
                    UniformVolume<T>::AddScalarQuantity(sval2);
                } else {
                    scalar_names(0) = new std::string;
                    scalar_names(0)->reserve(qtysize);
                    scalar_names(0)->assign(sval2);
                    nscalars = 1;
                }

                /* Read data values. */
                size_t datasize = 0;
                bool isfloat = false;
                bool isdouble = false;
                if(sval3 == "float"){
                    datasize = sizeof(float);
                    isfloat = true;
                    isdouble = false;
                }
                if(sval3 == "double"){
                    datasize = sizeof(double);
                    isfloat = false;
                    isdouble = true;
                }

                if((isBigEndian && !thisBigEndian) || (!isBigEndian && thisBigEndian)){
                    /* Read point-by-point and bit-flip each point value. */

                    if(false){
                        std::cout << "Big Endian System: " << StringManip::BoolToString(UniformVolume<T>::IsBigEndian()) << std::endl;
                        std::cout << "datasize: " << datasize << "     Big Endian Conversion" << std::endl;
                    }

                    desc = "Reading binary field scalars '" + UniformVolume<T>::ScalarQuantityName(nscalars-1) + "'";
                    qtsignals->EmitFunctionDesc(desc);

                    if(scalar_data == NULL){
                        float ffloat = (float)0.0e0;
                        double dfloat = (double)0.0e0;
                        for(size_t k=0; k<vslices; k++){
                            for(size_t i=0; i<vrows; i++){
                                for(size_t j=0; j<vcols; j++){
                                    if(isfloat){
                                        file.read(reinterpret_cast<char*>(&ffloat),datasize);
                                        UniformVolume<T>::ByteSwap(&ffloat,datasize);
                                        pscalars(nscalars-1)->operator ()(i,j,k) = (T)ffloat;
                                    }
                                    if(isdouble){
                                        file.read(reinterpret_cast<char*>(&dfloat),datasize);
                                        UniformVolume<T>::ByteSwap(&dfloat,datasize);
                                        pscalars(nscalars-1)->operator ()(i,j,k) = (T)dfloat;
                                    }
                                }
                            }
                            qtsignals->EmitFunctionProgress((float)k/(float)vslices);
                        }
                    } else {
                        /* Read data into existing array. */
                        scalar_data_points_read = 0;
                        float ffloat = (float)0.0e0;
                        double dfloat = (double)0.0e0;
                        for(size_t k=0; k<vslices; k++){
                            for(size_t i=0; i<vrows; i++){
                                for(size_t j=0; j<vcols; j++){
                                    if(isfloat){
                                        file.read(reinterpret_cast<char*>(&ffloat),datasize);
                                        UniformVolume<T>::ByteSwap(&ffloat,datasize);
                                        scalar_data[scalar_data_points_read] = (T)ffloat;
                                    }
                                    if(isdouble){
                                        file.read(reinterpret_cast<char*>(&dfloat),datasize);
                                        UniformVolume<T>::ByteSwap(&dfloat,datasize);
                                        scalar_data[scalar_data_points_read] = (T)dfloat;
                                    }

                                    /* Update number of points read.  If this value is greater than the size of the
                                     * array, inform the user of the error and return. */
                                    scalar_data_points_read++;
                                    if(scalar_data_points_read > scalar_data_size){
                                        std::cerr << "UniformVolume::VTKReadLegacyBinary() ERROR: Insufficient array size provided." << std::endl;
                                        return;
                                    }
                                }
                            }
                            qtsignals->EmitFunctionProgress((float)k/(float)vslices);
                        }
                    }


                } else {
                    desc = "Reading binary scalars '" + UniformVolume<T>::ScalarQuantityName(nscalars-1) + "'";
                    qtsignals->EmitFunctionDesc(desc);

                    /* No byte-swapping necessary.  Type-checking still performed, however, in case
                     * datatype in file does not match datatype of this object. */
                    if(scalar_data == NULL){
                        float ffloat = (float)0.0e0;
                        double dfloat = (double)0.0e0;
                        for(size_t k=0; k<vslices; k++){
                            for(size_t i=0; i<vrows; i++){
                                for(size_t j=0; j<vcols; j++){
                                    if(isfloat){
                                        file.read(reinterpret_cast<char*>(&ffloat),datasize);
                                        pscalars(nscalars-1)->operator ()(i,j,k) = (T)ffloat;
                                    }
                                    if(isdouble){
                                        file.read(reinterpret_cast<char*>(&dfloat),datasize);
                                        pscalars(nscalars-1)->operator ()(i,j,k) = (T)dfloat;
                                    }
                                }
                            }
                            qtsignals->EmitFunctionProgress((float)k/(float)vslices);
                        }
                    } else {
                        /* Read data into existing array. */
                        scalar_data_points_read = 0;
                        float ffloat = (float)0.0e0;
                        double dfloat = (double)0.0e0;
                        for(size_t k=0; k<vslices; k++){
                            for(size_t i=0; i<vrows; i++){
                                for(size_t j=0; j<vcols; j++){
                                    if(isfloat){
                                        file.read(reinterpret_cast<char*>(&ffloat),datasize);
                                        scalar_data[scalar_data_points_read] = (T)ffloat;
                                    }
                                    if(isdouble){
                                        file.read(reinterpret_cast<char*>(&dfloat),datasize);
                                        scalar_data[scalar_data_points_read] = (T)dfloat;
                                    }

                                    /* Update number of points read.  If this value is greater than the size of the
                                     * array, inform the user of the error and return. */
                                    scalar_data_points_read++;
                                    if(scalar_data_points_read > scalar_data_size){
                                        std::cerr << "UniformVolume::VTKReadLegacyBinary() ERROR: Insufficient array size provided." << std::endl;
                                        return;
                                    }
                                }
                            }
                            qtsignals->EmitFunctionProgress((float)k/(float)vslices);
                        }
                    }
                }
            }

            /* Set pointer to external array to NULL to prevent overwriting the freshly-read data. */
            if(scalar_data != NULL){
                scalar_data = NULL;
                UniformVolume<T>::RemoveAllData();
            }

            desc = "";
            qtsignals->EmitFunctionDesc(desc);
            qtsignals->EmitFunctionProgress(0.0e0);
        }


    } /* End of while()-loop terminated by EOF. */

    desc = "";
    qtsignals->EmitFunctionDesc(desc);
    qtsignals->EmitFunctionProgress(0.0e0);

} /* UniformVolume::VTKReadLegacyBinary() */


template <class T>
void UniformVolume<T>::VTKReadLegacyASCII(std::fstream &file)
{
    std::string discard;
    std::string discard1;
    std::string discard2;
    std::string datasettype;
    std::stringstream sstmp;
    int itmp1, itmp2, itmp3;
    T d1, d2, d3;
    T o1, o2, o3;
    T Tval;

    file >> discard1 >> discard2;        /* Dataset definition */

    /* Extract dataset definition */
    sstmp << discard2; datasettype = sstmp.str(); sstmp.str("");

    /* Read volume parameters based on dataset type */
    if(datasettype == "STRUCTURED_POINTS"){
        /* Read volume resolution */
        file >> discard1 >> itmp1 >> itmp2 >> itmp3;

        /* Reset resolution of this object */
        UniformVolume<T>::ResetResolution(itmp2,itmp1,itmp3,(T)0.0e0);

        file >> discard1 >> o1 >> o2 >> o3;	/* Volume origin */
        file >> discard1 >> d1 >> d2 >> d3;	/* Voxel spacing */

        /* Determine starting and ending points from origin and spacing info */
        xmin = o1; xmax = xmin + d1*((T)vcols - 1.0e0);
        ymin = o2; ymax = ymin + d2*((T)vrows - 1.0e0);
        zmin = o3; zmax = zmin + d3*((T)vslices - 1.0e0);
        xspacing = d1;
        yspacing = d2;
        zspacing = d3;

    } /* Image Data header */

    if(datasettype == "RECTILINEAR_GRID"){
        /* Read volume resolution */
        file >> discard1 >> itmp1 >> itmp2 >> itmp3;

        /* Reset resolution of this object */
        UniformVolume<T>::ResetResolution(itmp2,itmp1,itmp3,(T)0.0e0);

        /* Read X Coordinates */
        file >> discard1 >> discard1 >> discard1; /* Label line */
        file >> o1;
        UniformVolume<T>::setSpatialExtent(1,-1,o1);
        for(size_t i=0; i<vcols-2; i++){
            file >> discard;     /* Ignore all points except first and last */
        }
        file >> o1;
        UniformVolume<T>::setSpatialExtent(1,1,o1);

        /* Read Y Coordinates */
        file >> discard1 >> discard1 >> discard1; /* Label line */
        file >> o1;
        UniformVolume<T>::setSpatialExtent(0,-1,o1);
        for(size_t i=0; i<vrows-2; i++){
            file >> discard;     /* Ignore all points except first and last */
        }
        file >> o1;
        UniformVolume<T>::setSpatialExtent(0,1,o1);

        /* Read Z Coordinates */
        file >> discard1 >> discard1 >> discard1; /* Label line */
        file >> o1;
        UniformVolume<T>::setSpatialExtent(2,-1,o1);
        for(size_t i=0; i<vslices-2; i++){
            file >> discard;     /* Ignore all points except first and last */
        }
        file >> o1;
        UniformVolume<T>::setSpatialExtent(2,1,o1);
    } /* Rectilinear Grid header */

    /* Read POINT_DATA line */
    file >> discard1 >> discard2;

    /* Read until either "VECTORS" or "SCALARS" is read.  Repeat until end-of-file */
    while(file.eof() == false){
        discard1 = "na";
        while(discard1 != "SCALARS" && discard1 != "VECTORS" && file.eof() == false){
            file >> discard1;
        }

        if(discard1 == "SCALARS" && file.eof() == false){
            file >> discard2;                /* Read name of quantity */
            file >> discard1 >> discard1;    /* Read remainder of line */
            file >> discard1 >> discard1;    /* Read lookup table line */

            /* Add quantity and read-in values */
            T frac = 0.0e0;
            std::string progressstr;

            if(scalar_data == NULL){
                UniformVolume<T>::AddScalarQuantity(discard2);
            } else {
                scalar_names(0) = new std::string;
                scalar_names(0)->reserve(qtysize);
                scalar_names(0)->assign(discard2);
                nscalars = 1;
            }

            progressstr = "Reading ASCII scalars '" + UniformVolume<T>::ScalarQuantityName(nscalars-1) + "'";
            qtsignals->EmitFunctionProgress(frac);
            qtsignals->EmitFunctionProgress(frac,progressstr);

            if(scalar_data == NULL){
                for(size_t s=0; s<vslices; s++){
                    for(size_t r=0; r<vrows; r++){
                        for(size_t c=0; c<vcols; c++){
                            file >> Tval;
                            pscalars(nscalars-1)->operator ()(r,c,s) = (T)Tval;
                        }
                    }
                    frac = ((T)s + 1.0e0)/(T)vslices;
                    qtsignals->EmitFunctionProgress(frac);
                    qtsignals->EmitFunctionProgress(frac,progressstr);
                }
            } else {
                scalar_data_points_read = 0;
                for(size_t s=0; s<vslices; s++){
                    for(size_t r=0; r<vrows; r++){
                        for(size_t c=0; c<vcols; c++){
                            file >> Tval;
                            scalar_data[scalar_data_points_read] = Tval;

                            /* Update number of points read.  If this value is greater than the size of the
                             * array, inform the user of the error and return. */
                            scalar_data_points_read++;
                            if(scalar_data_points_read > scalar_data_size){
                                std::cerr << "UniformVolume::VTKReadLegacyASCII() ERROR: Insufficient array size provided." << std::endl;
                                return;
                            }
                        }
                    }

                    frac = ((T)s + 1.0e0)/(T)vslices;
                    qtsignals->EmitFunctionProgress(frac);
                    qtsignals->EmitFunctionProgress(frac,progressstr);
                }
            }

            /* Set pointer to external array to NULL to prevent overwriting the freshly-read data. */
            scalar_data = NULL;
            UniformVolume<T>::RemoveAllData();

            progressstr = "";
            qtsignals->EmitFunctionDesc(progressstr);
            qtsignals->EmitFunctionProgress(0.0e0);
        }

        if(discard1 == "VECTORS" && file.eof() == false){
            file >> discard2;                /* Read name of quantity */
            file >> discard1;                /* Read remainder of line */

            /*
              Determine number of vector components.  It is assumed that each vector
              is terminated by a newline.
             */
            int n = 0;
            std::stringstream readss;
            std::stringstream countss;
            std::string word;
            getline(file,word);  /* Catches the '\n' remaining on the VECTORS line */
            getline(file,word);  /* Reads the first line of vector data */
            readss.str(word);
            countss.str(word);

            while(getline(countss,word,' ')){
                n++;
            }

            /* Add quantity and read-in values */
            T frac = 0.0e0;
            std::string progressstr;
            UniformVolume<T>::AddVectorQuantity(discard2,n);

            progressstr = "Reading ASCII vectors '" + UniformVolume<T>::VectorQuantityName(nvectors-1) + "'";
            qtsignals->EmitFunctionProgress(frac);
            qtsignals->EmitFunctionProgress(frac,progressstr);

            /* The first line is stored in 'readss' */
            for(int i=0; i<n; i++){
                readss >> pvectors(nvectors-1)->operator ()(0,0,0,i);
            }

            for(size_t s=0; s<vslices; s++){
                for(size_t r=0; r<vrows; r++){
                    for(size_t c=0; c<vcols; c++){
                        if(c == 0 && r == 0 && s == 0){
                            /* First element already read-in, so no action to take.  For all
                               other combinations of (c,r,s), read the file as expected.
                             */
                        } else {
                            for(int l=0; l<n; l++){
                                file >> Tval;
                                pvectors(nvectors-1)->operator ()(r,c,s,l) = (T)Tval;
                            }
                        }
                    }
                }
                frac = ((T)s + 1.0e0)/(T)vslices;
                qtsignals->EmitFunctionProgress(frac);
                qtsignals->EmitFunctionProgress(frac,progressstr);
            }

            progressstr = "";
            qtsignals->EmitFunctionDesc(progressstr);
            qtsignals->EmitFunctionProgress(0.0e0);

        } /* End of conditional on VECTORS. */


        if(discard1 == "FIELD" && file.eof() == false){
            int nfields;
            file >> discard2;                /* Read name of quantity */
            file >> discard1 >> nfields;     /* Read remainder of line */

            for(int f=0; f<nfields; f++){

                /* Read field name. */
                size_t throwaway1, throwaway2;
                file >> discard2 >> throwaway1 >> throwaway2 >> discard1;

                /* Add quantity and read-in values */
                T frac = 0.0e0;
                std::string progressstr;

                if(scalar_data == NULL){
                    UniformVolume<T>::AddScalarQuantity(discard2);
                } else {
                    scalar_names(0) = new std::string;
                    scalar_names(0)->reserve(qtysize);
                    scalar_names(0)->assign(discard2);
                    nscalars = 1;
                }

                progressstr = "Reading ASCII field scalars '" + UniformVolume<T>::ScalarQuantityName(nscalars-1) + "'";
                qtsignals->EmitFunctionProgress(frac);
                qtsignals->EmitFunctionProgress(frac,progressstr);

                if(scalar_data == NULL){
                    for(size_t s=0; s<vslices; s++){
                        for(size_t r=0; r<vrows; r++){
                            for(size_t c=0; c<vcols; c++){
                                file >> Tval;
                                pscalars(nscalars-1)->operator ()(r,c,s) = (T)Tval;
                            }
                        }
                        frac = ((T)s + 1.0e0)/(T)vslices;
                        qtsignals->EmitFunctionProgress(frac);
                        qtsignals->EmitFunctionProgress(frac,progressstr);
                    }
                } else {
                    scalar_data_points_read = 0;
                    for(size_t s=0; s<vslices; s++){
                        for(size_t r=0; r<vrows; r++){
                            for(size_t c=0; c<vcols; c++){
                                file >> Tval;
                                scalar_data[scalar_data_points_read] = Tval;

                                /* Update number of points read.  If this value is greater than the size of the
                                 * array, inform the user of the error and return. */
                                scalar_data_points_read++;
                                if(scalar_data_points_read > scalar_data_size){
                                    std::cerr << "UniformVolume::VTKReadLegacyASCII() ERROR: Insufficient array size provided." << std::endl;
                                    return;
                                }
                            }
                        }

                        frac = ((T)s + 1.0e0)/(T)vslices;
                        qtsignals->EmitFunctionProgress(frac);
                        qtsignals->EmitFunctionProgress(frac,progressstr);
                    }
                }


                /* Set pointer to external array to NULL to prevent overwriting the freshly-read data. */
                scalar_data = NULL;
                UniformVolume<T>::RemoveAllData();

                progressstr = "";
                qtsignals->EmitFunctionDesc(progressstr);
                qtsignals->EmitFunctionProgress(0.0e0);
            }
        }
    }
} /* UniformVolume::VTKReadLegacyASCII() */


template <class T>
void UniformVolume<T>::PointCoordinates(const size_t row, const size_t col, const size_t slice, T &y, T &x, T &z) const
{
    x = xmin + xspacing*(T)col;
    y = ymin + yspacing*(T)row;
    z = zmin + zspacing*(T)slice;
} /* UniformVolume::PointCoordinates() */


template <class T>
void UniformVolume<T>::setScalarDestinationArray(T* array, const size_t array_size)
{
    scalar_data = array;
    scalar_data_size = array_size;
}


template <class T>
size_t UniformVolume<T>::NumberScalarPointsRead()
{
    return scalar_data_points_read;
}


template <class T>
void UniformVolume<T>::WriteXdmf(const int compression)
{
//    size_t narrays = nscalars + nvectors;

//    /* Define a single grid for all data, with multiple data arrays defined to on the grid. */
//    PArray1D<T*> data_ptr(narrays);
//    PArray1D<Array1D<size_t>*> data_dims(1);
//    PArray1D<Array1D<float>*> data_origin(1);
//    PArray1D<Array1D<float>*> data_spacing(1);
//    PArray1D<std::string*> data_name(narrays);


//    data_dims(0) = new Array1D<size_t>;
//    data_dims(0)->ResetSize(3, 0);
//    data_dims(0)->operator ()(0) = vcols;
//    data_dims(0)->operator ()(1) = vrows;
//    data_dims(0)->operator ()(2) = vslices;

//    data_origin(0) = new Array1D<float>;
//    data_origin(0)->ResetSize(3, 0.0f);
//    data_origin(0)->operator ()(0) = (float)xmin;
//    data_origin(0)->operator ()(1) = (float)ymin;
//    data_origin(0)->operator ()(2) = (float)zmin;

//    data_spacing(0) = new Array1D<float>;
//    data_spacing(0)->ResetSize(3, 0.0f);
//    data_spacing(0)->operator ()(0) = (float)xspacing;
//    data_spacing(0)->operator ()(1) = (float)yspacing;
//    data_spacing(0)->operator ()(2) = (float)zspacing;

//    for(size_t i=0; i<narrays; i++){

//        if(i < nscalars){
//            /* Add scalar quantity to output structure. */
//            data_ptr(i) = &pscalars(i)->operator [](0);

//            data_name(i) = new std::string;
//            data_name(i)->assign(scalar_names(i)->substr());

//        } else {
//            /* Add vector quantity to output structure. */
//            data_ptr(i) = &pvectors(i-nscalars)->operator [](0);

//            data_name(i) = new std::string;
//            data_name(i)->assign(vector_names(i-nscalars)->substr());
//        }
//    }

//    XdmfIO::data_info<T> data_info;

//    data_info.data = &data_ptr;
//    data_info.data_name = &data_name;
//    data_info.dims = &data_dims;
//    data_info.origin = &data_origin;
//    data_info.spacing = &data_spacing;

//    std::string filename;
//    filename = outputdir + "/" + filenamestem;

//    qtsignals->EmitFunctionDesc2("Writing XDMF File");

//    /* Call appropriate write function.  Typecasting is used to allow compilation.  During program execution,
//     * the check of type ID should properly handle the calling of the correct function. */
//    if(typeid(T) == typeid(short)){
//        XdmfIO::writeSingleUniformGrid(filename, NULL, NULL, (XdmfIO::data_info<short>*)&data_info, NULL, NULL,
//                                 NULL, NULL, NULL, NULL, compression);
//    }
//    if(typeid(T) == typeid(unsigned short)){
//        XdmfIO::writeSingleUniformGrid(filename, NULL, NULL, NULL, (XdmfIO::data_info<unsigned short>*)&data_info, NULL,
//                                 NULL, NULL, NULL, NULL, compression);
//    }
//    if(typeid(T) == typeid(int)){
//        XdmfIO::writeSingleUniformGrid(filename, NULL, NULL, NULL, NULL, (XdmfIO::data_info<int>*)&data_info,
//                                 NULL, NULL, NULL, NULL, compression);
//    }
//    if(typeid(T) == typeid(unsigned int)){
//        XdmfIO::writeSingleUniformGrid(filename, NULL, NULL, NULL, NULL, NULL,
//                                 (XdmfIO::data_info<unsigned int>*)&data_info, NULL, NULL, NULL, compression);
//    }
//    if(typeid(T) == typeid(long)){
//        XdmfIO::writeSingleUniformGrid(filename, NULL, NULL, NULL, NULL, NULL,
//                                 NULL, (XdmfIO::data_info<long>*)&data_info, NULL, NULL, compression);
//    }
//    if(typeid(T) == typeid(float)){
//        XdmfIO::writeSingleUniformGrid(filename, NULL, NULL, NULL, NULL, NULL,
//                                 NULL, NULL, (XdmfIO::data_info<float>*)&data_info, NULL, compression);
//    }
//    if(typeid(T) == typeid(double)){
//        XdmfIO::writeSingleUniformGrid(filename, NULL, NULL, NULL, NULL, NULL,
//                                 NULL, NULL, NULL, (XdmfIO::data_info<double>*)&data_info, compression);
//    }

//    qtsignals->EmitFunctionDesc2("");


    std::string filename;
    filename = outputdir + "/" + filenamestem;

    qtsignals->EmitFunctionDesc2("Writing XDMF File");

    vtkSmartPointer<vtkImageImport> imageImport = vtkSmartPointer<vtkImageImport>::New();
    imageImport->GlobalWarningDisplayOn();


    imageImport->SetDataSpacing((double)xspacing, (double)yspacing, (double)zspacing);
    imageImport->SetDataOrigin((double)xmin, (double)ymin, (double)zmin);
    imageImport->SetDataExtent(0, (int)vcols-1, 0, (int)vrows-1, 0, (int)vslices-1);
    imageImport->SetWholeExtent(0, (int)vcols-1, 0, (int)vrows-1, 0, (int)vslices-1);
    if(typeid(T) == typeid(float)){
        imageImport->SetDataScalarTypeToFloat();
    }
    if(typeid(T) == typeid(double)){
        imageImport->SetDataScalarTypeToDouble();
    }
    imageImport->SetNumberOfScalarComponents(1);
    imageImport->SetImportVoidPointer(&pscalars(0)->operator [](0), 1);
//    imageImport->CopyImportVoidPointer(&pscalars(0)->operator [](0), sizeof(T)*vcols*vrows*vslices);
    imageImport->Update();


    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData = imageImport->GetOutput();
    imageData->GetPointData()->GetScalars()->SetName(scalar_names(0)->substr().c_str());

    filename = filename + ".xmf";

    vtkSmartPointer<vtkXdmfWriter> writer = vtkSmartPointer<vtkXdmfWriter>::New();
    writer->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(imageData);
#else
    writer->SetInputData(imageData);
#endif
    writer->Write();



    qtsignals->EmitFunctionDesc2("");


} /* UniformVolume<T>::WriteXdmf() */


template <class T>
void UniformVolume<T>::ReadXDMFFile(std::string filename)
{
//#ifdef _DEBUG
//    bool debug = true;          /* Enable/disable debugging output to std::cout. */
//#else
//    bool debug = false;
//#endif

//    UniformVolume<T>::RemoveAllData();


//    /* Read the data into this object.
//     *
//     * Notes:
//     *      - Only 3D or 4D arrays are read in.  Other dimensionalities don't have meaning with a volume object.
//     *      - If the datatype of the data matches the datatype of this object, data is read directly in.  If the
//     *        datatype does not match, a copy operation is required.  A temporary array is created to store the
//     *        read-in data, and then it is copied and type-casted into this object.
//     */


//    /* Determine type of this data object, using the Xdmf integer identifiers. */
//    XdmfInt32 type_this;
//    if(typeid(T) == typeid(char)){
//        type_this = XDMF_INT8_TYPE;
//    }
//    if(typeid(T) == typeid(unsigned char)){
//        type_this = XDMF_UINT8_TYPE;
//    }
//    if(typeid(T) == typeid(short)){
//        type_this = XDMF_INT16_TYPE;
//    }
//    if(typeid(T) == typeid(unsigned short)){
//        type_this = XDMF_UINT16_TYPE;
//    }
//    if(typeid(T) == typeid(int)){
//        type_this = XDMF_INT32_TYPE;
//    }
//    if(typeid(T) == typeid(unsigned int)){
//        type_this = XDMF_UINT32_TYPE;
//    }
//    if(typeid(T) == typeid(float)){
//        type_this = XDMF_FLOAT32_TYPE;
//    }
//    if(typeid(T) == typeid(double)){
//        type_this = XDMF_FLOAT64_TYPE;
//    }


//    /* Create temporary Array1D objects for each datatype which could be used for type-conversion. */
//    Array1D<char> tmp_char(1, 'a');
//    Array1D<unsigned char> tmp_uchar(1, 'a');
//    Array1D<short> tmp_short(1, (short)0);
//    Array1D<unsigned short> tmp_ushort(1, (unsigned short)0);
//    Array1D<int> tmp_int(1, (int)0);
//    Array1D<unsigned int> tmp_uint(1, (unsigned int)0);
//    Array1D<float> tmp_float(1, 0.0f);
//    Array1D<double> tmp_double(1, (double)0.0e0);


//    XdmfDOM *d = new XdmfDOM;
//    XdmfGrid *grid = new XdmfGrid;
//    XdmfGrid *gridsave = new XdmfGrid;
//    XdmfTopology *topo = new XdmfTopology;
//    XdmfGeometry *geo = new XdmfGeometry;



//    d->Parse(filename.c_str());         /* Parse the XML file and allow navigation within it and
//                                         * data extraction. */

//    XdmfXmlNode node;
//    node = d->FindElementByPath("/Xdmf/Domain/Grid");   /* Find first grid element.  This will either be
//                                                         * a grid with attributed data, or a 'tree' or
//                                                         * 'collection' of other grids. */

//    grid->SetDOM(d);
//    grid->SetElement(node);
//    grid->UpdateInformation();                          /* Get information about grid. */
//    grid->Update();
//    int nchildren = grid->GetNumberOfChildren();

//    /* Write top-most grid information. */
//    if(debug){
//        std::cout << "File: " << filename << std::endl;
//        std::cout << "----------------------------------------------------" << std::endl;
//        std::cout << std::endl;
//        std::cout << "Grid" << std::endl;
//        std::cout << "  Name:   " << grid->GetName() << std::endl;
//        std::cout << "  Type:   " << grid->GetGridTypeAsString() << std::endl;
//        std::cout << "  # Child:  " << nchildren << std::endl;
//        std::cout << "  # Attrib: " << grid->GetNumberOfAttributes() << std::endl;
//        std::cout << std::endl;
//    }
//    gridsave = grid;    /* Save pointer to top-most grid. */


//    /* Require that grid is of type UNIFORM. */
//    if(grid->GetGridType() != XDMF_GRID_UNIFORM){
//        std::cerr << "UniformVolume<T>::ReadXDMFFile()  ERROR: Only UNIFORM grid is supported!";
//        std::cerr << "                                         Aborting!" << std::endl;
//        return;
//    }


//    int nloop = nchildren;          /* We want to loop through all the children grids, but if this grid   */
//    if(nloop == 0){                 /* has no children we still want to read its attributes, so force at  */
//        nloop = 1;                  /* least 1 loop iteration. */
//    }
//    for(int c=0; c<nloop; c++){

//        std::string indent("");     /* Used to provide extra indentation, if needed. */
//        int nchild = 0;             /* Number of children for local grid. */

//        if(nchildren == 0){
//            /* Nothing to be done. */
//        } else {

//            /* Set 'grid' to be the child specified by the loop iteration. */
//            grid = grid->GetChild(c);
//            nchild = grid->GetNumberOfChildren();   /* Number of children of the child grid. */
//            indent = "  ";                          /* Provide extra indentation in output. */


//            /* Write info about child grid. */
//            if(debug){
//                std::cout << indent << "-------------------------------" << std::endl;
//                std::cout << indent << "Grid " << c << std::endl;
//                std::cout << indent << "  Name:   " << grid->GetName() << std::endl;
//                std::cout << indent << "  Type:   " << grid->GetGridTypeAsString() << std::endl;
//                std::cout << indent << "  # Child:  " << nchild << std::endl;
//                std::cout << indent << "  # Attrib: " << grid->GetNumberOfAttributes() << std::endl;
//            }
//        }


//        /* Values to store grid size. */
//        XdmfInt64 xres = 0;
//        XdmfInt64 yres = 0;
//        XdmfInt64 zres = 0;
//        XdmfInt64 ncomp = 0;




//        /* Loop through attributes within grid.  Read data into C-style array (i.e., get data out of file
//         * and out of XdmfArray structure and into whatever data structure we want it in), and write
//         * information about the attribute to std::cout. */
//        for(int a=0; a<grid->GetNumberOfAttributes(); a++){

//            XdmfAttribute *at = grid->GetAttribute(a);  /* Get attribute from grid. */

//            at->UpdateInformation();                    /* Data 'metadata' read.  Heavy data not read yet. */

//            XdmfDataDesc *desc = at->GetShapeDesc();    /* Get shape description of data. */
//            int rank = desc->GetRank();                 /* Get data rank (i.e., number of dimensions. */
//            XdmfInt64 *dims = new XdmfInt64[rank];      /* Declare array for array size in each dimension. */

//            rank = desc->GetShape(dims);                /* Get array dimensions. */
//            XdmfInt64 nel = desc->GetNumberOfElements();/* Get number of elements in array. */

//            if(rank < 3){
//                std::cerr << "UniformVolume<T>::ReadXDMFFile()  ERROR: array dimensions less than 3 not supported!" << std::endl;
//            }
//            if(rank > 4){
//                std::cerr << "UniformVolume<T>::ReadXDMFFile()  ERROR: array dimensions greater than 4 not supported!" << std::endl;
//            }


//            if(rank == 3){
//                xres = dims[2];
//                yres = dims[1];
//                zres = dims[0];
//            }
//            if(rank == 4){
//                xres = dims[3];
//                yres = dims[2];
//                zres = dims[1];
//                ncomp = dims[0];
//            }

//            std::string name;
//            name = at->GetName();

//            XdmfInt32 type_file = desc->GetNumberType();

//            size_t num_elements = (size_t)desc->GetNumberOfElements();


//            /* Create scalar/vector array for data.  No memory is meaningfully allocated here since the vrows/vcols/vslices
//             * are not set to the true resolution yet.  True resolution will be set later, and arrays will be properly sized
//             * when they are read-in. */
//            if(rank == 3){
//                UniformVolume<T>::AddScalarQuantity(name);
//            }
//            if(rank == 4){
//                UniformVolume<T>::AddVectorQuantity(name, ncomp);
//            }




//                                                    /* Set to single-element to ensure that a valid pointer
//                                                     * is created.  All we want here is availability of the
//                                                     * pointer.  Allocation of the array here will cause duplicate
//                                                     * memory allocations and thus double memory required to
//                                                     * read-in the data. */

//            XdmfArray *data = new XdmfArray;        /* Create XdmfArray object to receive the data. */

//            /* Setting the XdmfArray data pointer to an existing array
//             * allows the Xdmf library to control memory (e.g., allocate
//             * as-needed), but it flags the memory as not belonging to
//             * the XdmfArray object, so deletion of the object does not
//             * cause the data to be lost. */


//            switch(type_file){

//            case(XDMF_INT8_TYPE):
//                data->SetDataPointer(&tmp_char[0]);
//                break;

//            case(XDMF_UINT8_TYPE):
//                data->SetDataPointer(&tmp_uchar[0]);
//                break;

//            case(XDMF_INT16_TYPE):
//                data->SetDataPointer(&tmp_short[0]);
//                break;

//            case(XDMF_UINT16_TYPE):
//                data->SetDataPointer(&tmp_ushort[0]);
//                break;

//            case(XDMF_INT32_TYPE):
//                data->SetDataPointer(&tmp_int[0]);
//                break;

//            case(XDMF_UINT32_TYPE):
//                data->SetDataPointer(&tmp_uint[0]);
//                break;

//            case(XDMF_FLOAT32_TYPE):
//                data->SetDataPointer(&tmp_float[0]);
//                break;

//            case(XDMF_FLOAT64_TYPE):
//                data->SetDataPointer(&tmp_double[0]);
//                break;

//            default:
//                std::cerr << "UniformVolume3D<T>::ReadXDMFFile()  ERROR: Datatype not recognized!" << std::endl;
//            }


//            at->Update();                           /* Memory allocated for data & data read into memory. */

//            data = at->GetValues(0);                /* Assign the read-in data to the XdmfArray object. */



//            /* If datatypes match, simply move the pointer to the newly-created scalar/vector quantity managed
//             * by this object.
//             *
//             * If datatypes differ, manually loop through the data and copy data into quantity managed by
//             * this object.
//             */
//            if(type_file == type_this && rank == 3){
//                if(rank == 3){
//                    pscalars(nscalars-1)->SetArrayPointer((T*)data->GetDataPointer(),
//                                                                      (size_t)yres, (size_t)xres, (size_t)zres, true);
//                }
//                if(rank == 4){
//                    pvectors(nvectors-1)->SetArrayPointer((T*)data->GetDataPointer(),
//                                                                      (size_t)yres, (size_t)xres, (size_t)zres, (size_t)ncomp,
//                                                                      true);
//                }
//            } else {
//                if(rank == 3){
//                    pscalars(nscalars-1)->ResetSize(yres, xres, zres);

//                    for(size_t k=0; k<num_elements; k++){
//                        pscalars(nscalars-1)->operator [](k) = (T)data->GetValueAsFloat64(k);
//                    }
//                }
//                if(rank == 4){
//                    pvectors(nvectors-1)->ResetSize(yres, xres, zres, ncomp);

//                    for(size_t k=0; k<num_elements; k++){
//                        pvectors(nvectors-1)->operator [](k) = (T)data->GetValueAsFloat64(k);
//                    }
//                }
//            }






//            /* Reset the pointer for data within my Array1D object to be that of the XdmfArray data.  The
//             * Xdmf library allocates memory using 'malloc()', so it must be free'd with 'free()'.  My
//             * array class uses 'new'/'delete', so the third parameter to this function informs the class
//             * that it will need to use 'free()' when releasing the memory occupied by this data.  If
//             * the data is stored in a "normal" c-array (i.e., not in my Array1D container), 'free()'
//             * must still be used to release the memory. */

//            data->Reset();                          /* Resets XdmfArray object to allow clean deletion. */

//            delete data;                            /* Delete XdmfArray.  This is needed because of the 'new'
//                                                     * operation above, and is placed here so that when I
//                                                     * check the data using my Array1D object I can be sure
//                                                     * that the data has survived the destruction of its
//                                                     * original array.  In normal use, data can be deletedtm
//                                                     * at any time after 'Reset()'.  Deletion at this specific
//                                                     * point is to verify proper data management behavior. */

//            tmp_char.SetArrayPointer(NULL, 1, false);
//            tmp_uchar.SetArrayPointer(NULL, 1, false);
//            tmp_short.SetArrayPointer(NULL, 1, false);
//            tmp_ushort.SetArrayPointer(NULL, 1, false);
//            tmp_int.SetArrayPointer(NULL, 1, false);
//            tmp_uint.SetArrayPointer(NULL, 1, false);
//            tmp_float.SetArrayPointer(NULL, 1, false);
//            tmp_double.SetArrayPointer(NULL, 1, false);


//            /* Write attribute information. */
//            if(debug){
//                std::cout << std::endl;
//                std::cout << indent << "    Attrib #" << a << std::endl;
//                std::cout << indent << "      Type: " << at->GetAttributeTypeAsString() << std::endl;

//                std::cout << indent << "      Data" << std::endl;
//                std::cout << indent << "        # Points: " << nel << std::endl;

//                /* Only write extra data information if there are a non-0 number of elements.  */
//                if(nel > 0){
//                    std::cout << indent << "        Shape:    " << desc->GetShapeAsString() << std::endl;
//                    std::cout << indent << "        Datatype: " << desc->GetNumberTypeAsString() << std::endl;
//                }

//            }

//        } /* Loop through attributes of grid. */



//        std::cout << std::endl;

//        /* Get topology information of grid. */
//        topo = grid->GetTopology();


//        /* Get geometry information of grid. */
//        geo = grid->GetGeometry();
//        geo->UpdateInformation();
//        geo->Update();

//        if(debug){
//            std::cout << indent << "Topology" << std::endl;
//            std::cout << indent << "  Type: " << topo->GetTopologyTypeAsString() << std::endl;
//            std::cout << std::endl;
//            std::cout << indent << "Geometry" << std::endl;
//            std::cout << indent << "  Type:    " << geo->GetGeometryTypeAsString() << std::endl;
//            std::cout << indent << "  Origin:  (" << geo->GetOriginX() << " , " << geo->GetOriginY() << " , "
//                      << geo->GetOriginZ() << " )" << std::endl;
//            std::cout << indent << "  Spacing: (" << geo->GetDx() << " , " << geo->GetDy() << " , "
//                      << geo->GetDz() << " )" << std::endl;
//            std::cout << std::endl;
//        }




//        vrows = yres;
//        vcols = xres;
//        vslices = zres;

//        xmin = geo->GetOriginZ();       /* X <--> Z switched due to KIJ ordering in XDMF format. */
//        ymin = geo->GetOriginY();
//        zmin = geo->GetOriginX();

//        xspacing = geo->GetDz();        /* X <--> Z switched due to KIJ ordering in XDMF format. */
//        yspacing = geo->GetDy();
//        zspacing = geo->GetDx();

//        xmax = xmin + (float)vcols*xspacing;
//        ymax = ymin + (float)vrows*yspacing;
//        zmax = zmin + (float)vslices*zspacing;


//        /* Set 'grid' to the head XdmfGrid object so that the GetChild() call at the top of the loop
//         * performs the desired function.  We want to loop over the children of the head node, so
//         * we must reset the pointer to the head grid. */
//        grid = gridsave;


//    } /* Loop through child grids. */


//    delete d;   /* Delete DOM.  As commented in other functions, objects associated with this DOM
//                 * will automatically be deleted. */

//    /* No delete needed for 'gridsave'.  At this point both 'grid' and 'gridsave' point to the same
//     * XdmfGrid object, which is associated with XdmfDOM 'd'.  The actual grid object will be deleted,
//     * which addresses both 'grid' and 'gridsave'. */

//    std::cout << std::endl;



//    if(debug){
//        std::cout << "# scalars: " << nscalars << std::endl;
//        std::cout << "# vectors: " << nvectors << std::endl;
//    }



    UniformVolume<T>::RemoveAllData();


    /* imageExport is used to move read-in data to a C-style array controlled by this object.  If the data cannot
     * be moved to an array controlled by this object, it will be copied.  If copied, we want the VTK library to
     * free the memory it allocates.  This means that the dataset which is returned by the file reading function is
     * not actually assigned to the export object until we know that a move (not a copy) is to be performed. */
    vtkSmartPointer<vtkImageExport> imageExport = vtkSmartPointer<vtkImageExport>::New();


    vtkSmartPointer<vtkXdmfReader> reader = vtkSmartPointer<vtkXdmfReader>::New();
    reader->SetFileName(filename.c_str());
    reader->UpdateInformation();
    reader->UpdateWholeExtent();

    vtkSmartPointer<vtkImageData> dataset = vtkSmartPointer<vtkImageData>::New();
    dataset = vtkImageData::SafeDownCast(reader->GetOutput());


    int dims[3] = {0, 0, 0};
    double *origin;
    double *spacing;

    dataset->GetDimensions(dims);
    origin = dataset->GetOrigin();
    spacing = dataset->GetSpacing();
    char *name = dataset->GetPointData()->GetScalars()->GetName();
    std::string str_name(name);


    int type_this;
    if(typeid(T) == typeid(char)){
        type_this = 2;
    }

    if(typeid(T) == typeid(signed char)){
        type_this = 15;
    }

    if(typeid(T) == typeid(unsigned char)){
        type_this = 3;
    }

    if(typeid(T) == typeid(short)){
        type_this = 4;
    }

    if(typeid(T) == typeid(unsigned short)){
        type_this = 5;
    }

    if(typeid(T) == typeid(int)){
        type_this = 6;
    }

    if(typeid(T) == typeid(unsigned int)){
        type_this = 7;
    }

    if(typeid(T) == typeid(long)){
        type_this = 8;
    }

    if(typeid(T) == typeid(unsigned long)){
        type_this = 9;
    }

    if(typeid(T) == typeid(float)){
        type_this = 10;
    }

    if(typeid(T) == typeid(double)){
        type_this = 11;
    }




    int type_file = dataset->GetScalarType();

    /* If either data is to be read into a user-supplied array (pointer to which is 'scalar_data'), or the datatypes
     * do not match, the read-in data will have to be copied to its destination array. */
    if(scalar_data || type_this != type_file){

        size_t npts = dims[0]*dims[1]*dims[2];

        if(scalar_data){

            /* If data is to be read into pre-existing array, copy values into that array. */

            for(size_t ii=0; ii<scalar_data_size; ii++){

                scalar_data[ii] = (T)dataset->GetPointData()->GetScalars()->GetComponent(ii, 0);
                scalar_data_points_read++;

                if(scalar_data_points_read == npts){
                    std::cerr << "UniformVolume<T>::ReadLegacyVTKFile()  ERROR: Insufficient destination array size." << std::endl;
                    std::cerr << "                                              Aborting." << std::endl;
                    return;
                }
            }

        } else {

            /* If no pre-existing array is to be used, then we must place data into this object, and typecast the
             * values to match the template type. */

            std::string str_name("nothing");

            UniformVolume<T>::ResetResolution((size_t)dims[1], (size_t)dims[0], (size_t)dims[2], (T)0.0e0);
            UniformVolume<T>::AddScalarQuantity(str_name);

            for(size_t i=0; i<npts; i++){
                pscalars(0)->operator [](i) = (T)dataset->GetPointData()->GetScalars()->GetComponent(i, 0);
            }

        }

    } else {

        /* Data destination is this object, and datatypes match, so the data pointer can simply be assigned
         * and no duplication of data is required. */

        reader->Register(reader->GetOutput());          /* Prevents loss of data when 'reader' goes out of scope.  This
                                                         * is desired since we want this object to control the memory
                                                         * rather than the VTK library functions. */
        imageExport->SetInputData(dataset);             /* Assign data to export object for movement into separate
                                                         * C-style array. */

//        UniformVolume<T>::ResetResolution(1, 1, 1, (T)0.0e0);
        UniformVolume<T>::AddScalarQuantity(str_name);

        vrows = dims[1];
        vcols = dims[0];
        vslices = dims[2];

        /* Setting the pointer like this appears to "move" the data into my array structure.  Valgrind
         * does not show any missing deallocation as a result of doing this, so I think this is OK. */
//        T* dataptr = (T*)imageExport->GetPointerToData();
//        imageExport->Export(&pscalars(0)->operator [](0));
//        T* dataptr;
//        imageExport->Export(dataptr);
//        pscalars(0)->SetArraySize(vrows, vcols, vslices, true);

        imageExport->Update();
        pscalars(0)->SetArrayPointer((T*)imageExport->GetPointerToData(), vrows, vcols, vslices, true);
    }


    xmin = origin[0];
    ymin = origin[1];
    zmin = origin[2];

    xspacing = spacing[0];
    yspacing = spacing[1];
    zspacing = spacing[2];

    xmax = xmin + xspacing*(float)vcols;
    ymax = ymin + yspacing*(float)vrows;
    zmax = zmin + zspacing*(float)vslices;


} /* UniformVolume<T>::ReadXDMFFile() */


template <class T>
void UniformVolume<T>::VTKWriteImageData()
{
    std::string filename;
    filename = outputdir + "/" + filenamestem;

    qtsignals->EmitFunctionDesc2("Writing VTK File");

    vtkSmartPointer<vtkImageImport> imageImport = vtkSmartPointer<vtkImageImport>::New();
    imageImport->SetDataSpacing((double)xspacing, (double)yspacing, (double)zspacing);
    imageImport->SetDataOrigin((double)xmin, (double)ymin, (double)zmin);
    imageImport->SetDataExtent(0, (int)vcols-1, 0, (int)vrows-1, 0, (int)vslices-1);
    imageImport->SetWholeExtent(0, (int)vcols-1, 0, (int)vrows-1, 0, (int)vslices-1);
    if(typeid(T) == typeid(float)){
        imageImport->SetDataScalarTypeToFloat();
    }
    if(typeid(T) == typeid(double)){
        imageImport->SetDataScalarTypeToDouble();
    }
    imageImport->SetNumberOfScalarComponents(1);
    imageImport->SetImportVoidPointer(&pscalars(0)->operator [](0), 1);
    imageImport->Update();


    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData = imageImport->GetOutput();

    int npieces = 1;
    if(vrows*vcols*vslices >= 2147483647){
        npieces += (int)((vrows*vcols*vslices)/2147483647);
        filename = filename + ".pvti";

        vtkSmartPointer<vtkXMLPImageDataWriter> writer = vtkSmartPointer<vtkXMLPImageDataWriter>::New();
        writer->SetFileName(filename.c_str());
        writer->SetHeaderTypeToUInt64();
        writer->SetNumberOfPieces(npieces);
        writer->SetStartPiece(0);
        writer->SetEndPiece(npieces-1);
#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(imageData);
#else
        writer->SetInputData(imageData);
#endif
        writer->Write();


    } else {
        filename = filename + ".vti";

        vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        writer->SetFileName(filename.c_str());
        writer->SetHeaderTypeToUInt64();
#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(imageData);
#else
        writer->SetInputData(imageData);
#endif
        writer->Write();
    }

    qtsignals->EmitFunctionDesc2("");
}
