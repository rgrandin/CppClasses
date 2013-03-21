#include <iostream>

#include <XdmfDemoConfig.h>

#include <basicXdmfTest.h>

int main()
{
    int selection;

    std::cout << "XdmfDemo -- version " << XdmfDemo_VERSION_MAJOR << "." << XdmfDemo_VERSION_MINOR << std::endl;
    std::cout << std::endl;
    std::cout << "Enter demo selection: " << std::endl;
    std::cout << "   0 - Write 3D array, user specified size, value = 100*k + 10*i + j" << std::endl;
    std::cout << std::endl;

    std::cout << "Selection: ";
    std::cin >> selection;
    std::cout << std::endl;
    std::cout << std::endl;

    switch(selection){
    case(0):
        size_t s1, s2, s3;
        std::cout << "Enter 3D array size (rows, cols, slices): ";
        std::cin >> s1 >> s2 >> s3;
        basicXdmfTest::writeXDMFFile(s1,s2,s3);
        break;

    default:
        std::cerr << "ERROR: Invalid demo option specified!" << std::endl;
        break;
    }



    return 0;
}

