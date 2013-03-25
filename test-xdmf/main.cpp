#include <iostream>

#include <XdmfDemoConfig.h>

#include <basicXdmfTest.h>

#include <XdmfIO.h>

int main()
{
    int selection = 2;

    std::cout << "XdmfDemo -- version " << XdmfDemo_VERSION_MAJOR << "." << XdmfDemo_VERSION_MINOR << std::endl;
    std::cout << std::endl;
    std::cout << "Enter demo selection: " << std::endl;
    std::cout << "   0 - Write 3D array, user specified size, value = 100*k + 10*i + j" << std::endl;
    std::cout << "   1 - Convert volume file to XDMF" << std::endl;
    std::cout << "   2 - Read XDMF file and write information to std::cout" << std::endl;
    std::cout << "   3 - Write multiple 3D arrays, user specified size, value = 100*k + 10*i + j" << std::endl;
    std::cout << "   4 - Write multiple arrays of multiple dimensionality, size, and datatype" << std::endl;
    std::cout << std::endl;

    std::cout << "Selection: ";
    std::cin >> selection;
    std::cout << std::endl;
    std::cout << std::endl;

    switch(selection){
    case(0):
        basicXdmfTest::writeXDMFFile();
        break;

    case(1):
        basicXdmfTest::convertVolume();
        break;

    case(2):
        basicXdmfTest::readXDMFFile();
        break;

    case(3):
        basicXdmfTest::writeXDMFFile_MultiArray();
        break;

    case(4):
        basicXdmfTest::writeXDMFFile_MultiTest();
        break;

    case(999):
        /* Include this function call force compilation of this function for testing purposes. */
        XdmfIO::writeUniformGrid("nofile", NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
        break;

    default:
        std::cerr << "ERROR: Invalid demo option specified!" << std::endl;
        break;
    }



    return 0;
}

