#ifndef VTKXMLIOTEST_H
#define VTKXMLIOTEST_H

#include <Array3D.h>

//#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <vtkStructuredPoints.h>

#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPImageDataWriter.h>
#include <vtkImageImport.h>

#include <vtkXMLImageDataReader.h>
#include <vtkXMLPImageDataReader.h>
#include <vtkImageExport.h>




template<class TReader> vtkDataSet *ReadAnXMLFile(const char*fileName)
{
    vtkSmartPointer<TReader> reader = vtkSmartPointer<TReader>::New();
    reader->SetFileName(fileName);
    reader->Update();
    reader->GetOutput()->Register(reader);
    return vtkDataSet::SafeDownCast(reader->GetOutput());
}



int writeVTKFile(size_t size1, size_t size2, size_t size3)
{
    if(false){
        /* Copied from http://www.vtk.org/Wiki/VTK/Examples/Cxx/IO/XMLStructuredGridWriter */
        vtkSmartPointer<vtkStructuredGrid> structuredGrid =
                vtkSmartPointer<vtkStructuredGrid>::New();

        vtkSmartPointer<vtkPoints> points =
                vtkSmartPointer<vtkPoints>::New();

        points->InsertNextPoint(0, 0, 0);
        points->InsertNextPoint(1, 0, 0);
        points->InsertNextPoint(0, 1, 0);
        points->InsertNextPoint(1, 1, 0);
        points->InsertNextPoint(0, 2, 0);
        points->InsertNextPoint(1, 2, 1);

        // Specify the dimensions of the grid
        structuredGrid->SetDimensions(2,3,1);
        structuredGrid->SetPoints(points);

        // Write file
        vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
                vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName("output.vts");
#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(structuredGrid);
#else
        writer->SetInputData(structuredGrid);
#endif
        writer->Write();
    }




    /* Test actually writing and reading image data using the VTK file IO classes. */

    float datavalue = 3.0f;
    Array3D<float> *data = new Array3D<float>(size2, size1, size3, datavalue);

    int width = (int)size1;
    int height = (int)size2;
    int depth = (int)size3;

    vtkSmartPointer<vtkImageImport> imageImport = vtkSmartPointer<vtkImageImport>::New();
    imageImport->SetDataSpacing(1, 1, 1);
    imageImport->SetDataOrigin(0, 0, 0);
    imageImport->SetDataExtent(0, width-1, 0, height-1, 0, depth-1);
    imageImport->SetWholeExtent(0, width-1, 0, height-1, 0, depth-1);
    imageImport->SetDataScalarTypeToFloat();
    imageImport->SetNumberOfScalarComponents(1);
    imageImport->SetImportVoidPointer(&data->operator [](0), 1);
            imageImport->Update();


    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData = imageImport->GetOutput();


    //int npieces = 16;
    //std::string filename("output.pvti");
    //vtkSmartPointer<vtkXMLPImageDataWriter> writer2 = vtkSmartPointer<vtkXMLPImageDataWriter>::New();
    std::string filename("output_serial.vti");
    vtkSmartPointer<vtkXMLImageDataWriter> writer2 = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer2->SetFileName(filename.c_str());
    //writer2->SetNumberOfPieces(npieces);
    //writer2->SetStartPiece(0);
    //writer2->SetEndPiece(npieces-1);
#if VTK_MAJOR_VERSION <= 5
    writer2->SetInput(imageData);
#else
    writer2->SetInputData(imageData);
#endif
    writer2->Write();

    std::cout << "File written" << std::endl;


    data->ResetSize(1,1,1);


    vtkSmartPointer<vtkImageExport> imageExport = vtkSmartPointer<vtkImageExport>::New();
    vtkDataSet *dataset;

//    dataset = ReadAnXMLFile<vtkXMLPImageDataReader>(filename.c_str());
    dataset = ReadAnXMLFile<vtkXMLImageDataReader>(filename.c_str());

    int dims[3] = {0, 0, 0};
    imageExport->SetInput(dataset);
    imageExport->GetDataDimensions(dims);

    data->ResetSize((size_t)dims[1], (size_t)dims[2], (size_t)dims[0], datavalue-1.0f);

    imageExport->SetExportVoidPointer(&data->operator [](0));

    imageExport->Export();


    std::cout << std::endl;
    std::cout << "Test values: " << std::endl;
    std::cout << "  " << data->operator ()(0,0,0) << " " << data->operator ()(0,1,0)
              << " " << data->operator ()(0,2,0) << std::endl;
    std::cout << std::endl;
    std::cout << "Correct value: " << datavalue << std::endl;



    delete data;


    return EXIT_SUCCESS;
}


#endif // VTKXMLIOTEST_H
