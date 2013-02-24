#include <ArrayBase.h>
#include <DataFilters.h>

int main()
{
    std::string result;

    ArrayBase<float> abf;
    abf.Test(result);

    if(result == "SUCCESS"){
        std::cout << "ArrayBase result: SUCCESS" << std::endl;
    } else {
        std::cout << "ArrayBase result: FAILED" << std::endl;
        std::cout << result << std::endl;
    }


    result = "";
    DataFilters<float> dff;
    dff.Test(result);

    if(result == "SUCCESS"){
        std::cout << "DataFilters result: SUCCESS" << std::endl;
    } else {
        std::cout << "DataFilters result: FAILED" << std::endl;
        std::cout << result << std::endl;
    }


}
