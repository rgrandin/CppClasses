/**
 * @file ProgressInfo.cpp
 * @author Robert Grandin
 * @brief Implementation of ProgressInfo class.
 */


#include <ProgressInfo.h>


ProgressInfo::ProgressInfo()
{
    completion_progress = (float)0.0e0;

    function_description = "NA";

    function_description2 = "NA";

}


ProgressInfo::ProgressInfo(const ProgressInfo &a) :
    completion_progress(a.completion_progress),
    function_description(a.function_description),
    function_description2(a.function_description2)
{
}


#ifdef CXX11
ProgressInfo::ProgressInfo(ProgressInfo &&a) : ProgressInfo()
{
    ProgressInfoSwap(*this, a);
}
#endif



ProgressInfo::~ProgressInfo()
{
}


ProgressInfo& ProgressInfo::operator=(ProgressInfo &a)
{
    ProgressInfoSwap(*this, a);
    return *this;
}
