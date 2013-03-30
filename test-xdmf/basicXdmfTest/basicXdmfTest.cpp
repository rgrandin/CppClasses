#include <basicXdmfTest.h>

void basicXdmfTest::writeXDMFFile()
{
    /* This function adapted from:
     *   http://www.xdmf.org/index.php/XDMF_API (see Writing XDMF at bottom of page)
     *
     * This function is very similar to writeXDMFFile_MultiArray(), and detailed comments will be placed
     * in writeXDMFFile_MultiArray() since it is the more-complex function.
     */


    /* Get extra info. */
    std::string arrayname;
    std::string outputname;
    int compression;
    size_t s1, s2, s3;

    std::cout << std::endl;
    std::cout << "Running writeXDMFFile()" << std::endl;
    std::cout << "  Enter 3D array size (rows, cols, slices): ";
    std::cin >> s1 >> s2 >> s3;
    std::cout << "  Enter name for array: ";
    std::cin >> arrayname;
    std::cout << "  Enter name for output file (no extension): ";
    std::cin >> outputname;
    std::cout << "  Enter compression of heavy data (0 - 9, 0 uncompressed): ";
    std::cin >> compression;
    std::cout << std::endl;





    std::cout << "  Creating 3D float array of size ( " << s1 << " , " << s2 << " , "
              << s3 << " )" << std::endl;

    /* Create dummy 3D array of values to be written. */
    Array3D<float> data(s1, s2, s3, 1.0f);

    for(size_t k=0; k<s3; k++){
        for(size_t i=0; i<s1; i++){
            for(size_t j=0; j<s2; j++){
                data(i,j,k) = (float)k*100.0 + (float)i*10.0 + (float)j;

            }
        }
    }

    std::cout << "    Done!" << std::endl;
    std::cout << "    value(i,j,k) = 100*k + 10*i + j" << std::endl;





    std::cout << std::endl;
    std::cout << "  Writing data as XDMF file with heavy data stored as HDF5" << std::endl;
    double t1 = omp_get_wtime();


    /* Represent data as XdmfArray. */
    XdmfArray *xarray = new XdmfArray;
    XdmfInt64 shape[3] = { (XdmfInt64)s3, (XdmfInt64)s1, (XdmfInt64)s2 };

    xarray->SetAllowAllocate(false);
    xarray->SetShape(3, shape);
    xarray->SetDataPointer(&data[0]);
    xarray->SetNumberType(XDMF_FLOAT32_TYPE);


    XdmfDOM *d = new XdmfDOM;
    XdmfRoot *root = new XdmfRoot;
    XdmfDomain *domain = new XdmfDomain;
    XdmfGrid *grid = new XdmfGrid;
    XdmfTopology *topo = new XdmfTopology;
    XdmfGeometry *geo = new XdmfGeometry;
    XdmfAttribute *attrib = new XdmfAttribute;
    XdmfInformation *info = new XdmfInformation;



    root->SetDOM(d);
    root->SetVersion(2.2);
    root->Build();

    info->SetName("InfoName");
    info->SetValue("1.23");
    root->Insert(info);

    root->Insert(domain);

    grid->SetName("Structured Grid");

    topo = grid->GetTopology();
    topo->SetTopologyType(XDMF_3DCORECTMESH);
    topo->GetShapeDesc()->SetShape(3, shape);

    geo = grid->GetGeometry();
    geo->SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ);
    geo->SetOrigin(0.0, 0.0, 0.0);
    geo->SetDxDyDz(1.0, 1.0, 1.0);

    domain->Insert(grid);

    attrib->SetName(arrayname.c_str());
    attrib->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE);
    attrib->SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);
    attrib->SetValues(xarray);

    std::string heavydataname(outputname);
    heavydataname = heavydataname + ".h5:/" + arrayname;

    xarray->SetHeavyDataSetName(heavydataname.c_str());
    xarray->SetCompression(compression);

    grid->Insert(attrib);

    root->Build();

    outputname = outputname + ".xmf";
    d->Write(outputname.c_str());

    double t2 = omp_get_wtime();


    std::cout << "    Done!" << std::endl;
    std::cout << "    XDMF File: " << outputname << std::endl;
    std::cout << "    Heavy Data: " << heavydataname << std::endl;
    std::cout << "    Write time: " << t2 - t1 << " [sec]" << std::endl;
    std::cout << std::endl;



    /* Cleanup pointers. */
    delete xarray;  /* Deletion of XmdfDOM object does not delete array, so 'xarray'
                     * must be manually deleted here. */
    delete d;       /* Deleting XmdfDOM object deletes its children, so other
                     * Xmdf____ objects do not need to be manually deleted here.
                     * This behavior was verified using Valgrind to locate memory
                     * leaks, and this implementation produced no errors. */

}



void basicXdmfTest::writeXDMFFile_MultiArray()
{
    /* This function adapted from:
     *   http://www.xdmf.org/index.php/XDMF_API (see Writing XDMF at bottom of page)
     */


    /* Get extra info. */
    std::string *arrayname;     /* Name of array to be used in heavy data file. */
    std::string outputname;     /* Name of XDMF XML file. */
    int compression;
    int narrays;
    size_t *s1, *s2, *s3;       /* Array sizes in the 1st, 2nd, and 3rd dimensions. */

    std::cout << std::endl;
    std::cout << "Running writeXDMFFile_MultiArray()" << std::endl;
    std::cout << "  Enter number of arrays: ";
    std::cin >> narrays;

    /* Declare arrays for array dimensions and name. */
    s1 = new size_t[narrays];
    s2 = new size_t[narrays];
    s3 = new size_t[narrays];
    arrayname = new std::string[narrays];
    Array3D<float> *data = new Array3D<float>[narrays];

    for(int i=0; i<narrays; i++){
        std::cout << "    Array " << i << ": Enter size (rows, cols, slices): ";
        std::cin >> s1[i] >> s2[i] >> s3[i];

        std::cout << "    Array " << i << ": Enter name for array: ";
        std::cin >> arrayname[i];

        std::cout << std::endl;
    }

    std::cout << "  Enter name for output file (no extension): ";
    std::cin >> outputname;
    std::cout << "  Enter compression of heavy data (0 - 9, 0 uncompressed): ";
    std::cin >> compression;
    std::cout << std::endl;


    /* Create XdmfArray for each array to be written. */
    XdmfArray **xarray = new XdmfArray*[narrays];   /* We want an array of pointers to XdmfArray objects. */

    std::cout << std::endl;
    std::cout << "  Writing data as XDMF file with heavy data stored as HDF5" << std::endl;
    double t1 = omp_get_wtime();


    XdmfDOM *d = new XdmfDOM;
    XdmfRoot *root = new XdmfRoot;
    XdmfDomain *domain = new XdmfDomain;
    XdmfInformation *info = new XdmfInformation;

    XdmfGrid **grid = new XdmfGrid*[narrays];           /* We want an array of pointers to Xdmf____ objects.  */
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
    tree_grid->SetName("MultiArrayTest");
    tree_grid->SetGridType(XDMF_GRID_TREE);

    domain->Insert(tree_grid);


    /* For each data array, create its grid and insert the grid and data into the XML structure. */
    for(int i=0; i<narrays; i++){

        /* Since this function is using cooked-up data, create an array of the specified size and fill
         * it using a known function (allows checking read-in later).
         *
         * Naturally, this part would be omitted if we are attempting to write existing data to the disk. */

        std::cout << "  Creating array " << i << std::endl;

        std::cout << "    Creating 3D float array of size ( " << s1[i] << " , " << s2[i] << " , "
                  << s3[i] << " )" << std::endl;

        /* Create dummy 3D array of values to be written. */
        data[i].ResetSize(s1[i], s2[i], s3[i], 1.0f);

        for(size_t k=0; k<s3[i]; k++){
            for(size_t ii=0; ii<s1[i]; ii++){
                for(size_t j=0; j<s2[i]; j++){
                    data[i](ii,j,k) = (float)k*100.0 + (float)ii*10.0 + (float)j;

                }
            }
        }

        std::cout << "      Done!" << std::endl;
        std::cout << "      value(i,j,k) = 100*k + 10*i + j" << std::endl;



        /* Create an XdmfArray object to represent the data.  No values are actually moved out of the
         * pre-existing array (in this case the Array3D object populated above).  Deletion of the XdmfArray
         * will not affect the actual data, either.  The XdmfArray object only facilitates the interface
         * between the actual data and the data-output. */

        XdmfInt64 shape[3] = { (XdmfInt64)s3[i], (XdmfInt64)s1[i], (XdmfInt64)s2[i] }; /* Array shape is "KIJ", which
                                                                                        * means that the slowest-varying
                                                                                        * index comes first, and the
                                                                                        * fastest-varying (i.e. innermost)
                                                                                        * index comes last. */

        xarray[i] = new XdmfArray;
        xarray[i]->SetAllowAllocate(false);                 /* Prevent the XdmfArray from allocating any memory. */
        xarray[i]->SetShape(3, shape);                      /* Set rank and dimension sizes of array. */
        xarray[i]->SetDataPointer(&data[i][0]);             /* Set pointer to previously-existing data. */
        xarray[i]->SetNumberType(XDMF_FLOAT32_TYPE);        /* Identify the data as 32-bit float ('float' in C/C++). */


        grid[i] = new XdmfGrid;                         /* Create a grid for the data. */
        grid[i]->SetName("Structured Grid");            /* Name the grid.  This is user-choice. */
        grid[i]->SetGridType(XDMF_GRID_UNIFORM);        /* Set grid type.  This example is for image data, which is
                                                         * on a uniform grid.  Grid choice can be other types based
                                                         * on needs. */

        topo[i] = new XdmfTopology;                     /* Create a topology object to help define grid. */
        topo[i] = grid[i]->GetTopology();               /* Set to be grid topology. */
        topo[i]->SetTopologyType(XDMF_3DCORECTMESH);    /* Set grid type, again this is for image data and other
                                                         * options are possible for other applications. */
        topo[i]->GetShapeDesc()->SetShape(3, shape);    /* Set shape to match that of the data array. */

        geo[i] = new XdmfGeometry;                              /* Create a geometry object to complete grid
                                                                  * structure definition. */
        geo[i] = grid[i]->GetGeometry();                        /* Set to be the grid geometry. */
        geo[i]->SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ);   /* Set for image data.  Again, more-complex options exist. */
        geo[i]->SetOrigin((float)i, (float)i, (float)i);        /* Set origin. */
        geo[i]->SetDxDyDz((float)i+1, (float)i+1, (float)i+1);  /* Set point spacing. */

        /* NOTE: more-complex grid/topology/geometry definitions are possible and can require additional information
         *       to be set (e.g., defining an unstructured grid requires specifying grid point coordinates and
         *       connectivity of nodes). */

        /* Insert this data grid into the super-grid. */
        tree_grid->Insert(grid[i]);

        XdmfAttribute *attrib = new XdmfAttribute;              /* Create attribute to describe data. */
        attrib->SetName(arrayname[i].c_str());                  /* Set descriptive name for data. */
        attrib->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE); /* Define data to exist at nodes (rather than within
                                                                 * cells or on cell faces). */
        attrib->SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);   /* Define data to be scalar. */
        attrib->SetValues(xarray[i]);                           /* Define data values to be those stored in the
                                                                 * XdmfArray object created above (which in-turn takes
                                                                 * us to the c-array represented by the XdmfArray object. */

        /* Build name of file/path to be used with HDF5 heavy-data file.  This is set so that the
         * HDF5 filename matches that of the XML file.  The array name is used to identify the data within the
         * HDF5 file.  All arrays will be placed in the same HDF5 file. */
        std::string heavydataname(outputname);
        heavydataname = heavydataname + ".h5:/" + arrayname[i];

        xarray[i]->SetHeavyDataSetName(heavydataname.c_str());  /* Set name used for heavy-data.  This will cause the
                                                                 * actual data values to be written to the specified
                                                                 * HDF5 file/array.  If the array is small-enough,
                                                                 * values will be directly-inserted into the XML.  The
                                                                 * default value for this threshold is 100, but can be
                                                                 * modified by the user.  If data values are located
                                                                 * within the XML file, this function call has no
                                                                 * effect. */

        xarray[i]->SetCompression(compression);                 /* Set compression level to be used by HDF5 file.  This
                                                                 * is also ignored if data is embedded in the XML file.
                                                                 * Compression method will be noted in the HDF5 file and
                                                                 * decompression will happen automatically when the
                                                                 * HDF5 file is opened in a suitable viewer or these
                                                                 * Xdmf____ classes. */

        /* Inform user of file name/array for heavy data. */
        std::cout << "      Heavy Data: " << heavydataname << std::endl;

        /* Insert attribute into grid.  Note that this grid is a child of the "super-grid". */
        grid[i]->Insert(attrib);

    }

    /* Now that all data has been inserted, build the XML document. */
    root->Build();

    /* Add file extension to user-specified file stem and write data to the disk.  This will cause the XML file
     * to be written as well as the HDF5 file, if necessary.  If multiple HDF5 files are specified above, each will
     * be written at this point. */
    outputname = outputname + ".xmf";
    d->Write(outputname.c_str());

    double t2 = omp_get_wtime();


    /* Let user know that file-writing is complete, state the name of the file written, and calculate the time
     * required to build and write the output. */
    std::cout << "    Done!" << std::endl;
    std::cout << "    XDMF File: " << outputname << std::endl;
    std::cout << "    Write time: " << t2 - t1 << " [sec]" << std::endl;
    std::cout << std::endl;



    /* Cleanup pointers. */
    delete [] xarray;  /* Deletion of XmdfDOM object does not delete array, so 'xarray'
                        * must be manually deleted here. */
    delete d;       /* Deleting XmdfDOM object deletes its children, so other
                     * Xmdf____ objects do not need to be manually deleted here.
                     * This behavior was verified using Valgrind to locate memory
                     * leaks, and this implementation produced no errors.  Valgrind shows that there
                     * may still be lost memory, but this appears to be within the Xdmf____ classes
                     * and thus is not under the programmer's control (without digging deep into source
                     * code). */

}



void basicXdmfTest::convertVolume()
{
    /* This function adapted from:
     *   http://www.xdmf.org/index.php/XDMF_API (see Writing XDMF at bottom of page)
     *
     * This function is very similar to writeXDMFFile() and writeXDMFFile_MultiArray(), and detailed comments
     * can be found in writeXDMFFile_MultiArray().
     *
     * Functionality contained within this function that is not found within writeXDMFFile_MultiArray() is
     * the reading of a pre-existing volume file (either VTK or CNDE VOL) to populate a single 3D data array.
     * This data array is then written to the disk as an XDMF/HDF5 file combo.  Generating and writing the
     * XDMF file is identical to writeXDMFFile(), which is a simplification of writeXDMFFile_MultiArray(),
     * and you are referred to those files for detailed commentary regarding the writing of XDMF files.
     */


    /* Get extra info. */
    std::string volumename;
    std::string arrayname;
    std::string outputname;
    size_t s1, s2, s3;
    int compression;

    std::cout << std::endl;
    std::cout << "Running convertVolume()" << std::endl;
    std::cout << "  Enter name of volume file to be converted: ";
    std::cin >> volumename;
    std::cout << "  Enter name for array: ";
    std::cin >> arrayname;
    std::cout << "  Enter name for output file (no extension): ";
    std::cin >> outputname;
    std::cout << "  Enter compression of heavy data (0 - 9, 0 uncompressed): ";
    std::cin >> compression;
    std::cout << std::endl;






    std::cout << "  Reading volume file" << std::endl;

    /* If file is VTK, assume it to be Big Endian byte ordering (automatically ignored if specified file is ASCII). */
    bool isBigEndian = false;
    if(StringManip::DetermFileExt(volumename) == "vtk" || StringManip::DetermFileExt(volumename) == "VTK"){
        isBigEndian = true;
    }

    /* Create volume object and read specified file. */
    UniformVolume<float> uv;
    uv.ReadFile(volumename, isBigEndian);


    std::cout << "    Done!" << std::endl;



    /* Get pointer to start of data array (want actual data, not Array3D which contains the data. */
    float *data_start = &uv(0,0,0,0);

    /* Get volume size, point spacing, and origin. */
    s1 = uv.GetResolution(1);
    s2 = uv.GetResolution(0);
    s3 = uv.GetResolution(2);

    double sp1 = (double)uv.PointSpacing(0);
    double sp2 = (double)uv.PointSpacing(1);
    double sp3 = (double)uv.PointSpacing(2);

    double o1 = (double)uv.SpatialExtent(1, -1);
    double o2 = (double)uv.SpatialExtent(0, -1);
    double o3 = (double)uv.SpatialExtent(2, -1);




    /* Write data to disk as an XDMF file, with heavy data stored as HDF5.  From here forward the functionality
     * is very similar to writeXDMFFile() and writeXDMFFile_MultiArray() and you are referred to those functions
     * for detailed comments. */


    std::cout << std::endl;
    std::cout << "  Writing data as XDMF file with heavy data stored as HDF5" << std::endl;
    double t1 = omp_get_wtime();


    /* Represent data as XdmfArray. */
    XdmfArray *xarray = new XdmfArray;
    XdmfInt64 shape[3] = { (XdmfInt64)s3, (XdmfInt64)s2, (XdmfInt64)s1 };

    xarray->SetAllowAllocate(false);
    xarray->SetShape(3, shape);
    xarray->SetDataPointer(data_start);
    xarray->SetNumberType(XDMF_FLOAT32_TYPE);


    XdmfDOM *d = new XdmfDOM;
    XdmfRoot *root = new XdmfRoot;
    XdmfDomain *domain = new XdmfDomain;
    XdmfGrid *grid = new XdmfGrid;
    XdmfTopology *topo = new XdmfTopology;
    XdmfGeometry *geo = new XdmfGeometry;
    XdmfAttribute *attrib = new XdmfAttribute;
    XdmfInformation *info = new XdmfInformation;



    root->SetDOM(d);
    root->SetVersion(2.2);
    root->Build();

    info->SetName("InfoName");
    info->SetValue("1.23");
    root->Insert(info);

    root->Insert(domain);

    grid->SetName("Structured Grid");

    topo = grid->GetTopology();
    topo->SetTopologyType(XDMF_3DCORECTMESH);
    topo->GetShapeDesc()->SetShape(3, shape);

    geo = grid->GetGeometry();
    geo->SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ);
    geo->SetOrigin(o1, o2, o3);
    geo->SetDxDyDz(sp1, sp2, sp3);

    domain->Insert(grid);

    attrib->SetName(arrayname.c_str());
    attrib->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE);
    attrib->SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);
    attrib->SetValues(xarray);

    std::string heavydataname(outputname);
    heavydataname = heavydataname + ".h5:/" + arrayname;

    xarray->SetHeavyDataSetName(heavydataname.c_str());
    xarray->SetCompression(compression);

    grid->Insert(attrib);

    root->Build();

    outputname = outputname + ".xmf";
    d->Write(outputname.c_str());

    double t2 = omp_get_wtime();


    std::cout << "    Done!" << std::endl;
    std::cout << "    XDMF File: " << outputname << std::endl;
    std::cout << "    Heavy Data: " << heavydataname << std::endl;
    std::cout << "    Write time: " << t2 - t1 << " [sec]" << std::endl;
    std::cout << std::endl;



    /* Cleanup pointers. */
    delete xarray;  /* Deletion of XmdfDOM object does not delete array, so 'xarray'
                     * must be manually deleted here. */
    delete d;       /* Deleting XmdfDOM object deletes its children, so other
                     * Xmdf____ objects do not need to be manually deleted here.
                     * This behavior was verified using Valgrind to locate memory
                     * leaks, and this implementation produced no errors. */
}



void basicXdmfTest::readXDMFFile()
{
    /* This function adapted from:
     *   http://www.xdmf.org/index.php/XDMF_API (see Reading XDMF at bottom of page)
     *
     * This function is very basic.  It can handle a collection of uniform grids (or a single uniform grid),
     * but no checks are made to prevent erros resulting from trying to read other grids.  I expect functionality
     * for other grids will be similar, and much of the code required to read files should remain the same.
     */


    /* Get extra info. */
    std::string filename("c.xmf");

    std::cout << std::endl;
    std::cout << "Running readXDMFFile()" << std::endl;
    std::cout << "  Enter name of file to be read: ";
    std::cin >> filename;
    std::cout << std::endl;

    XdmfDOM *d = new XdmfDOM;
    XdmfGrid *grid = new XdmfGrid;
    XdmfGrid *gridsave = new XdmfGrid;
    XdmfTopology *topo = new XdmfTopology;
    XdmfGeometry *geo = new XdmfGeometry;


    Array1D<float> datavals(1, 0.0f);



    d->Parse(filename.c_str());         /* Parse the XML file and allow navigation within it and
                                         * data extraction. */

    XdmfXmlNode node;
    node = d->FindElementByPath("/Xdmf/Domain/Grid");   /* Find first grid element.  This will either be
                                                         * a grid with attributed data, or a 'tree' or
                                                         * 'collection' of other grids. */

    grid->SetDOM(d);
    grid->SetElement(node);
    grid->UpdateInformation();                          /* Get information about grid. */
    grid->Update();
    int nchildren = grid->GetNumberOfChildren();

    /* Write top-most grid information. */
    std::cout << "Grid" << std::endl;
    std::cout << "  Name:   " << grid->GetName() << std::endl;
    std::cout << "  Type:   " << grid->GetGridTypeAsString() << std::endl;
    std::cout << "  # Child:  " << nchildren << std::endl;
    std::cout << "  # Attrib: " << grid->GetNumberOfAttributes() << std::endl;
    std::cout << std::endl;

    gridsave = grid;    /* Save pointer to top-most grid. */


    int nloop = nchildren;          /* We want to loop through all the children grids, but if this grid   */
    if(nloop == 0){                 /* has no children we still want to read its attributes, so force at  */
        nloop = 1;                  /* least 1 loop iteration. */
    }
    for(int c=0; c<nloop; c++){

        std::string indent("");     /* Used to provide extra indentation, if needed. */
        int nchild = 0;             /* Number of children for local grid. */

        if(nchildren == 0){
            /* Nothing to be done. */
        } else {

            /* Set 'grid' to be the child specified by the loop iteration. */
            grid = grid->GetChild(c);
            nchild = grid->GetNumberOfChildren();   /* Number of children of the child grid. */
            indent = "  ";                          /* Provide extra indentation in output. */


            /* Write info about child grid. */
            std::cout << indent << "-------------------------------" << std::endl;
            std::cout << indent << "Grid " << c << std::endl;
            std::cout << indent << "  Name:   " << grid->GetName() << std::endl;
            std::cout << indent << "  Type:   " << grid->GetGridTypeAsString() << std::endl;
            std::cout << indent << "  # Child:  " << nchild << std::endl;
            std::cout << indent << "  # Attrib: " << grid->GetNumberOfAttributes() << std::endl;

        }




        /* Loop through attributes within grid.  Read data into C-style array (i.e., get data out of file
         * and out of XdmfArray structure and into whatever data structure we want it in), and write
         * information about the attribute to std::cout. */
        for(int a=0; a<grid->GetNumberOfAttributes(); a++){

            XdmfAttribute *at = grid->GetAttribute(a);  /* Get attribute from grid. */

            at->UpdateInformation();                    /* Data 'metadata' read.  Heavy data not read yet. */

            XdmfDataDesc *desc = at->GetShapeDesc();    /* Get shape description of data. */
            int rank = desc->GetRank();                 /* Get data rank (i.e., number of dimensions. */
            XdmfInt64 dims[rank];                       /* Declare array for array size in each dimension. */

            rank = desc->GetShape(dims);                /* Get array dimensions. */
            XdmfInt64 nel = desc->GetNumberOfElements();/* Get number of elements in array. */

            datavals.ResetSize(1);                  /* Set to single-element to ensure that a valid pointer
                                                     * is created.  All we want here is availability of the
                                                     * pointer.  Allocation of the array here will cause duplicate
                                                     * memory allocations and thus double memory required to
                                                     * read-in the data. */

            XdmfArray *data = new XdmfArray;        /* Create XdmfArray object to receive the data. */

            data->SetDataPointer(&datavals[0]);     /* Setting the XdmfArray data pointer to an existing array
                                                     * allows the Xdmf library to control memory (e.g., allocate
                                                     * as-needed), but it flags the memory as not belonging to
                                                     * the XdmfArray object, so deletion of the object does not
                                                     * cause the data to be lost. */

            at->Update();                           /* Memory allocated for data & data read into memory. */

            data = at->GetValues(0);                /* Assign the read-in data to the XdmfArray object. */


            /* Reset the pointer for data within my Array1D object to be that of the XdmfArray data.  The
             * Xdmf library allocates memory using 'malloc()', so it must be free'd with 'free()'.  My
             * array class uses 'new'/'delete', so the third parameter to this function informs the class
             * that it will need to use 'free()' when releasing the memory occupied by this data.  If
             * the data is stored in a "normal" c-array (i.e., not in my Array1D container), 'free()'
             * must still be used to release the memory. */
            datavals.SetArrayPointer((float*)data->GetDataPointer(), (size_t)nel, true);

            data->Reset();                          /* Resets XdmfArray object to allow clean deletion. */

            delete data;                            /* Delete XdmfArray.  This is needed because of the 'new'
                                                     * operation above, and is placed here so that when I
                                                     * check the data using my Array1D object I can be sure
                                                     * that the data has survived the destruction of its
                                                     * original array.  In normal use, data can be deleted
                                                     * at any time after 'Reset()'.  Deletion at this specific
                                                     * point is to verify proper data management behavior. */


            /* Write attribute information. */
            std::cout << std::endl;
            std::cout << indent << "    Attrib #" << a << std::endl;
            std::cout << indent << "      Type: " << at->GetAttributeTypeAsString() << std::endl;

            std::cout << indent << "      Data" << std::endl;
            std::cout << indent << "        # Points: " << nel << std::endl;

            /* Only write extra data information if there are a non-0 number of elements.  */
            if(nel > 0){
                std::cout << indent << "        Shape:    " << desc->GetShapeAsString() << std::endl;
                std::cout << indent << "        Datatype: " << desc->GetNumberTypeAsString() << std::endl;

                std::cout << indent << "          test val:  " << datavals[nel-1] << std::endl;
                std::cout << indent << "          should be: " << (dims[0]-1)*100 + (dims[1]-1)*10 + dims[2]-1 << std::endl;
            }


        } /* Loop through attributes of grid. */



        std::cout << std::endl;

        /* Get topology information of grid. */
        topo = grid->GetTopology();
        std::cout << indent << "Topology" << std::endl;
        std::cout << indent << "  Type: " << topo->GetTopologyTypeAsString() << std::endl;
        std::cout << std::endl;


        /* Get geometry information of grid. */
        geo = grid->GetGeometry();
        geo->UpdateInformation();
        geo->Update();
        std::cout << indent << "Geometry" << std::endl;
        std::cout << indent << "  Type:    " << geo->GetGeometryTypeAsString() << std::endl;
        std::cout << indent << "  Origin:  (" << geo->GetOriginX() << " , " << geo->GetOriginY() << " , "
                  << geo->GetOriginZ() << " )" << std::endl;
        std::cout << indent << "  Spacing: (" << geo->GetDx() << " , " << geo->GetDy() << " , "
                  << geo->GetDz() << " )" << std::endl;
        std::cout << std::endl;


        /* Set 'grid' to the head XdmfGrid object so that the GetChild() call at the top of the loop
         * performs the desired function.  We want to loop over the children of the head node, so
         * we must reset the pointer to the head grid. */
        grid = gridsave;


    } /* Loop through child grids. */


    delete d;   /* Delete DOM.  As commented in other functions, objects associated with this DOM
                 * will automatically be deleted. */

    /* No delete needed for 'gridsave'.  At this point both 'grid' and 'gridsave' point to the same
     * XdmfGrid object, which is associated with XdmfDOM 'd'.  The actual grid object will be deleted,
     * which addresses both 'grid' and 'gridsave'. */

    std::cout << std::endl;

}


void basicXdmfTest::writeXDMFFile_MultiTest()
{
    std::string filename("multi-test-file");

    /* Create 1D char array. */
    Array1D<char> char_array(101, 'b');
    char_array.setName("char_test");

    /* Create 1D short array. */
    Array1D<short> short_1d(203, (short)2);
    short_1d.setName("short_1d_test");

    /* Create 2D short array. */
    Array2D<short> short_array(203, 31, (short)3112);
    short_array.setName("short_2d_test");

    /* Create 4D double array. */
    Array4D<double> double_array(11, 23, 2, 7, (double)2.374e0);
    double_array.setName("double_4d_test");


    XdmfIO::data_info<char> char_data;
    XdmfIO::data_info<short> short_data;
    XdmfIO::data_info<double> double_data;



    /* Create data-struct for char data. */
    PArray1D<char*> char_ptr(1);
    PArray1D<Array1D<size_t>*> char_dims(1);
    PArray1D<Array1D<float>*> char_origin(1);
    PArray1D<Array1D<float>*> char_spacing(1);
    PArray1D<std::string*> char_name(1);


    char_dims(0) = new Array1D<size_t>;
    char_origin(0) = new Array1D<float>;
    char_spacing(0) = new Array1D<float>;
    char_name(0) = new std::string;

    char_ptr(0) = &char_array[0];
    char_dims(0)->ResetSize(1);    char_dims(0)->operator ()(0) = 101;
    char_origin(0)->ResetSize(1);  char_origin(0)->operator ()(0) = 0.0e0;
    char_spacing(0)->ResetSize(1); char_spacing(0)->operator ()(0) = 1.0e0;
    //char_name(0)->assign(char_array.Name());
    char_name(0)->assign("char_array");

    char_data.data = &char_ptr;
    char_data.dims = &char_dims;
    char_data.origin = &char_origin;
    char_data.spacing = &char_spacing;
    char_data.data_name = &char_name;


    /* Create data struct for short data. */
    PArray1D<short*> short_ptr(2);
    PArray1D<Array1D<size_t>*> short_dims(2);
    PArray1D<Array1D<float>*> short_origin(2);
    PArray1D<Array1D<float>*> short_spacing(2);
    PArray1D<std::string*> short_name(2);

    short_dims(0) = new Array1D<size_t>;
    short_origin(0) = new Array1D<float>;
    short_spacing(0) = new Array1D<float>;
    short_name(0) = new std::string;

    short_dims(1) = new Array1D<size_t>;
    short_origin(1) = new Array1D<float>;
    short_spacing(1) = new Array1D<float>;
    short_name(1) = new std::string;

    short_ptr(0) = &short_array[0];
    short_dims(0)->ResetSize(2);    short_dims(0)->operator ()(0) = 203;      short_dims(0)->operator ()(1) = 31;
    short_origin(0)->ResetSize(2);  short_origin(0)->operator ()(0) = 0.1e0;  short_origin(0)->operator ()(1) = 0.15e0;
    short_spacing(0)->ResetSize(2); short_spacing(0)->operator ()(0) = 0.5e0; short_spacing(0)->operator ()(1) = 0.5e0;
    short_name(0)->assign(short_array.Name());

    short_ptr(1) = &short_1d[0];
    short_dims(1)->ResetSize(1);    short_dims(1)->operator ()(0) = 203;
    short_origin(1)->ResetSize(1);  short_origin(1)->operator ()(0) = 0.25e0;
    short_spacing(1)->ResetSize(1); short_spacing(1)->operator ()(0) = 0.75e0;
    short_name(1)->assign(short_1d.Name());

    short_data.data = &short_ptr;
    short_data.dims = &short_dims;
    short_data.origin = &short_origin;
    short_data.spacing = &short_spacing;
    short_data.data_name = &short_name;


    /* Create data struct for double data. */
    PArray1D<double*> double_ptr(1);
    PArray1D<Array1D<size_t>*> double_dims(1);
    PArray1D<Array1D<float>*> double_origin(1);
    PArray1D<Array1D<float>*> double_spacing(1);
    PArray1D<std::string*> double_name(1);

    double_dims(0) = new Array1D<size_t>;
    double_origin(0) = new Array1D<float>;
    double_spacing(0) = new Array1D<float>;
    double_name(0) = new std::string;

    double_ptr(0) = &double_array[0];
    double_dims(0)->ResetSize(4);
    for(int i=0; i<4; i++){
        double_dims(0)->operator ()(i) = double_array.GetDim(i+1);
    }
    double_origin(0)->ResetSize(4);  double_origin(0)->ResetVal(2.0e0);
    double_spacing(0)->ResetSize(4); double_spacing(0)->ResetVal(1.7e0);
    double_name(0)->assign(double_array.Name());

    double_data.data = &double_ptr;
    double_data.dims = &double_dims;
    double_data.origin = &double_origin;
    double_data.spacing = &double_spacing;
    double_data.data_name = &double_name;


    XdmfIO::writeUniformGrid(filename, &char_data, NULL, &short_data, NULL, NULL, NULL, NULL, NULL, &double_data, 0);


    /* Set pointers to data to be NULL to prevent attemptd double-deletion. */
    char_ptr(0) = NULL;
    short_ptr(0) = NULL;
    short_ptr(1) = NULL;
    double_ptr(0) = NULL;

}
