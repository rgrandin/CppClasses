/**
 * @file QtIntermediary.cxx
 * @author Robert Grandin
 * @brief Implementation of QtIntermediary class.
 */


#include "QtIntermediary.h"

#ifdef USEQT
QtIntermediary::QtIntermediary(QObject *parent) :
    QtIntermediaryBase(parent)
{
}


#else

QtIntermediary::QtIntermediary()
{
}
#endif


QtIntermediary::~QtIntermediary()
{
}



