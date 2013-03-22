#include <basicXdmfTest.h>

void basicXdmfTest::writeXDMFFile()
{
    /* This function adapted from:
     *   http://www.xdmf.org/index.php/XDMF_API (see Writing XDMF at bottom of page)
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
    std::string *arrayname;
    std::string outputname;
    int compression;
    int narrays;
    size_t *s1, *s2, *s3;

    std::cout << std::endl;
    std::cout << "Running writeXDMFFile_MultiArray()" << std::endl;
    std::cout << "  Enter number of arrays: ";
    std::cin >> narrays;

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
    XdmfArray **xarray = new XdmfArray*[narrays];

    std::cout << std::endl;
    std::cout << "  Writing data as XDMF file with heavy data stored as HDF5" << std::endl;
    double t1 = omp_get_wtime();


    XdmfDOM *d = new XdmfDOM;
    XdmfRoot *root = new XdmfRoot;
    XdmfDomain *domain = new XdmfDomain;
    XdmfInformation *info = new XdmfInformation;

    XdmfGrid **grid = new XdmfGrid*[narrays];
    XdmfTopology **topo = new XdmfTopology*[narrays];
    XdmfGeometry **geo = new XdmfGeometry*[narrays];


    root->SetDOM(d);
    root->SetVersion(2.2);
    root->Build();

    info->SetName("InfoName");
    info->SetValue("1.23");
    root->Insert(info);

    root->Insert(domain);

    XdmfGrid *tree_grid = new XdmfGrid;
    tree_grid->SetName("MultiArrayTest");
    tree_grid->SetGridType(XDMF_GRID_TREE);

    domain->Insert(tree_grid);

    for(int i=0; i<narrays; i++){

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


        XdmfInt64 shape[3] = { (XdmfInt64)s3[i], (XdmfInt64)s1[i], (XdmfInt64)s2[i] };
        xarray[i] = new XdmfArray;
        xarray[i]->SetAllowAllocate(false);
        xarray[i]->SetShape(3, shape);
        xarray[i]->SetDataPointer(&data[i][0]);
        xarray[i]->SetNumberType(XDMF_FLOAT32_TYPE);


        grid[i] = new XdmfGrid;
        grid[i]->SetName("Structured Grid");
        grid[i]->SetGridType(XDMF_GRID_UNIFORM);

        topo[i] = new XdmfTopology;
        topo[i] = grid[i]->GetTopology();
        topo[i]->SetTopologyType(XDMF_3DCORECTMESH);
        topo[i]->GetShapeDesc()->SetShape(3, shape);

        geo[i] = new XdmfGeometry;
        geo[i] = grid[i]->GetGeometry();
        geo[i]->SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ);
        geo[i]->SetOrigin((float)i, (float)i, (float)i);
        geo[i]->SetDxDyDz((float)i+1, (float)i+1, (float)i+1);

        tree_grid->Insert(grid[i]);

        XdmfAttribute *attrib = new XdmfAttribute;
        attrib->SetName(arrayname[i].c_str());
        attrib->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE);
        attrib->SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);
        attrib->SetValues(xarray[i]);

        std::string heavydataname(outputname);
        heavydataname = heavydataname + ".h5:/" + arrayname[i];

        xarray[i]->SetHeavyDataSetName(heavydataname.c_str());
        xarray[i]->SetCompression(compression);

        std::cout << "      Heavy Data: " << heavydataname << std::endl;

        grid[i]->Insert(attrib);

    }

    root->Build();

    outputname = outputname + ".xmf";
    d->Write(outputname.c_str());

    double t2 = omp_get_wtime();


    std::cout << "    Done!" << std::endl;
    std::cout << "    XDMF File: " << outputname << std::endl;
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



void basicXdmfTest::convertVolume()
{
    /* This function adapted from:
     *   http://www.xdmf.org/index.php/XDMF_API (see Writing XDMF at bottom of page)
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

    bool isBigEndian = false;
    if(StringManip::DetermFileExt(volumename) == "vtk" || StringManip::DetermFileExt(volumename) == "VTK"){
        isBigEndian = true;
    }

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
     */


    /* Get extra info. */
    //std::string filename("test01.xmf");   /* 800 MiB, good for verifying non-duplication of data. */
    std::string filename("c.xmf");          /* ~32 kiB, good for quick testing. */

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



    d->Parse(filename.c_str());

    XdmfXmlNode node;
    node = d->FindElementByPath("/Xdmf/Domain/Grid");

    grid->SetDOM(d);
    grid->SetElement(node);
    grid->UpdateInformation();
    grid->Update();
    int nchildren = grid->GetNumberOfChildren();
    std::cout << "Grid" << std::endl;
    std::cout << "  Name:   " << grid->GetName() << std::endl;
    std::cout << "  Type:   " << grid->GetGridTypeAsString() << std::endl;
    std::cout << "  # Child:  " << nchildren << std::endl;
    std::cout << "  # Attrib: " << grid->GetNumberOfAttributes() << std::endl;
    std::cout << std::endl;

    gridsave = grid;


    int nloop = nchildren;
    if(nloop == 0){
        nloop = 1;
    }
    for(int c=0; c<nloop; c++){

        std::string indent("");
        int nchild = 0;

        if(nchildren == 0){
            nchildren = 1;
        } else {
            grid = grid->GetChild(c);
            nchild = grid->GetNumberOfChildren();
            indent = "  ";

            std::cout << indent << "-------------------------------" << std::endl;
            std::cout << indent << "Grid " << c << std::endl;
            std::cout << indent << "  Name:   " << grid->GetName() << std::endl;
            std::cout << indent << "  Type:   " << grid->GetGridTypeAsString() << std::endl;
            std::cout << indent << "  # Child:  " << nchild << std::endl;
            std::cout << indent << "  # Attrib: " << grid->GetNumberOfAttributes() << std::endl;

        }





        for(int a=0; a<grid->GetNumberOfAttributes(); a++){
            std::cout << std::endl;
            std::cout << indent << "    Attrib #" << a << std::endl;
            XdmfAttribute *at = grid->GetAttribute(a);

            std::cout << indent << "      Type: " << at->GetAttributeTypeAsString() << std::endl;

            at->UpdateInformation();    /* Data 'metadata' read.  Heavy data not read yet. */

            XdmfDataDesc *desc = at->GetShapeDesc();
            int rank = desc->GetRank();
            XdmfInt64 dims[rank];

            rank = desc->GetShape(dims);
            XdmfInt64 nel = desc->GetNumberOfElements();

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
                                                     * original array. */




            std::cout << indent << "      Data" << std::endl;
            std::cout << indent << "        # Points: " << nel << std::endl;
            if(nel > 0){
                std::cout << indent << "        Shape:    " << desc->GetShapeAsString() << std::endl;
                std::cout << indent << "        Datatype: " << desc->GetNumberTypeAsString() << std::endl;

                std::cout << indent << "          test val:  " << datavals[nel-1] << std::endl;
                std::cout << indent << "          should be: " << (dims[0]-1)*100 + (dims[1]-1)*10 + dims[2]-1 << std::endl;
            }


        } /* Loop through attributes of grid. */



        std::cout << std::endl;

        topo = grid->GetTopology();
        std::cout << indent << "Topology" << std::endl;
        std::cout << indent << "  Type: " << topo->GetTopologyTypeAsString() << std::endl;
        std::cout << std::endl;

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


        grid = gridsave;


    } /* Loop through grid child grids. */






    delete d;

    std::cout << std::endl;

}
