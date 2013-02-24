#include <ArrayBase.h>

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


}
