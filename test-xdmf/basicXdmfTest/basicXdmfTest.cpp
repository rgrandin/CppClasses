#include <basicXdmfTest.h>

void basicXdmfTest::writeXDMFFile(size_t size1, size_t size2, size_t size3)
{
    /* This function adapted from:
     *   http://www.xdmf.org/index.php/XDMF_API (see Writing XDMF at bottom of page)
     */


    /* Create dummy 3D array of values to be written. */
    Array3D<float> data(size1, size2, size3, 1.0f);

    for(size_t k=0; k<size3; k++){
        for(size_t i=0; i<size1; i++){
            for(size_t j=0; j<size2; j++){
                data(i,j,k) = (float)k*100.0 + (float)i*10.0 + (float)j;

            }
        }
    }



    /* Represent data as XdmfArray. */
    XdmfArray *xarray = new XdmfArray;
    XdmfInt64 shape[3] = { (XdmfInt64)size3, (XdmfInt64)size1, (XdmfInt64)size2 };

    xarray->SetAllowAllocate(false);
    xarray->SetShape(3, shape);
    xarray->SetDataPointer(&data[0]);
    xarray->SetNumberType(XDMF_FLOAT32_TYPE);


    XdmfHDF *xhdf = new XdmfHDF;
    xhdf->CopyShape(xarray);
    xhdf->CopyType(xarray);
    xhdf->Open("h5test.h5:/Sample01", "w");
    xhdf->Write(xarray);
    xhdf->Close();

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

    attrib->SetName("TestData");
    attrib->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE);
    attrib->SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);
    attrib->SetValues(xarray);

    xarray->SetHeavyDataSetName("h5test.h5:/Sample01");

    grid->Insert(attrib);

    root->Build();

    d->Write("test.xmf");


//    std::cout << "Debugging:" << std::endl;
//    std::cout << "  Topology type: " << topo->GetTopologyType() << std::endl;
//    std::cout << "  Geometry type: " << geo->GetGeometryType() << std::endl;
//    std::cout << std::endl;


    /* Need to include this cleanup to prevent memory leaks. */
    delete xarray;
    delete xhdf;
    delete d;
    delete root;
    delete domain;
    delete grid;
//    delete topo;      /* Commented-out to avoid runtime seg-fault. */
//    delete geo;       /* Commented-out to avoid runtime seg-fault. */
    delete attrib;
    delete info;
}
