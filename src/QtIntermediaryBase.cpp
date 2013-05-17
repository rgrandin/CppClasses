/**
 * @file QtIntermediaryBase.cpp
 * @author Robert Grandin
 * @brief Implementation of QtIntermediaryBase class.
 */

#include "QtIntermediaryBase.h"

#ifdef USEQT
QtIntermediaryBase::QtIntermediaryBase(QObject *parent) :
    QObject(parent)
{
}


#else

QtIntermediaryBase::QtIntermediaryBase()
{
}
#endif


QtIntermediaryBase::~QtIntermediaryBase()
{
}


void QtIntermediaryBase::PerformCalculations()
{
    /* Nothing to do here.  This function is to be overloaded. */
}


void QtIntermediaryBase::ExitThread()
{
    /* Nothing to do here.  This function is to be overloaded. */
}


void QtIntermediaryBase::EnableThread()
{
    /* Nothing to do here.  This function is to be overloaded. */
}


void QtIntermediaryBase::EmitFunctionProgress(const int frac)
{
#ifdef USEQT
    emit FunctionProgress(frac);
#else
    /* Dummy code to suppress 'unused parameter' compiler warnings. */
    int frac2 = frac;
    frac2++;
#endif
}


void QtIntermediaryBase::EmitFunctionProgress(const float frac)
{
#ifdef USEQT
    int qfrac = (int)(frac*100.0e0);
    emit FunctionProgress(qfrac);
#else
    /* Dummy code to suppress 'unused parameter' compiler warnings. */
    float frac2 = frac;
    frac2 += 1.0;
#endif
}


void QtIntermediaryBase::EmitFunctionProgress(const double frac)
{
#ifdef USEQT
    int qfrac = (int)(frac*100.0e0);
    emit FunctionProgress(qfrac);
#else
    /* Dummy code to suppress 'unused parameter' compiler warnings. */
    double frac2 = frac;
    frac2 += 1.0;
#endif
}


void QtIntermediaryBase::EmitFunctionProgress(const float frac, const std::string desc)
{
#ifdef USEQT
    int qfrac = (int)(frac*100.0e0);
    QString qdesc = QString::fromStdString(desc);
    emit FunctionProgress(qfrac,qdesc);
#else
    /* Dummy code to suppress 'unused parameter' compiler warnings. */
    float frac2 = frac;
    frac2 += 1.0;

    std::string desc2;
    desc2 = desc + "2";
#endif
}


void QtIntermediaryBase::EmitFunctionEstimatedTimeRemaining(const double time_seconds)
{
#ifdef USEQT
    emit EstimatedTimeRemaining(time_seconds);
#else
    /* Dummy code to suppress 'unused parameter' compiler warnings. */
    double dummyval = time_seconds;
    dummyval += 1.0e0;
#endif
}


void QtIntermediaryBase::EmitFunctionEstimatedTimeRemaining2(const double time_seconds, const std::string &desc)
{
#ifdef USEQT
    emit EstimatedTimeRemaining2(time_seconds, QString::fromStdString(desc));
#else
    /* Dummy code to suppress 'unused parameter' compiler warnings. */
    double dummyval = time_seconds;
    dummyval += 1.0e0;
    std::string dummydesc = "NA";
    dummydesc = dummydesc + desc;
#endif
}


#ifdef USEQT
void QtIntermediaryBase::EmitFunctionDesc(QString desc)
{
    emit FunctionLabel(desc);
}
#endif


void QtIntermediaryBase::EmitFunctionDesc(std::string desc)
{
#ifdef USEQT
    QString qdesc = QString::fromStdString(desc);
    emit FunctionLabel(qdesc);
#else
    /* Dummy code to suppress 'unused parameter' compiler warnings. */
    std::string desc2;
    desc2 = desc + "2";
#endif
}


void QtIntermediaryBase::EmitFunctionDesc2(std::string desc)
{
#ifdef USEQT
    QString qdesc = QString::fromStdString(desc);
    emit FunctionLabel2(qdesc);
#else
    /* Dummy code to suppress 'unused parameter' compiler warnings. */
    std::string desc2;
    desc2 = desc + "2";
#endif
}


void QtIntermediaryBase::EmitDone()
{
#ifdef USEQT
    emit done();
#endif
}
