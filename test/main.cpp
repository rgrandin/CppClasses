#include <ArrayBase.h>
#include <DataFilters.h>

int main()
{
    std::string result;

    ArrayBase<float> abf;
    abf.Test(result);
    std::cout << "ArrayBase result: " << result << std::endl;


    result = "";
    DataFilters<float> dff;
    dff.Test(result);
    std::cout << "DataFilters result: " << result << std::endl;


}
