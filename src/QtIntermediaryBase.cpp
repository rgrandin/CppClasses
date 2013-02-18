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
#endif
}


void QtIntermediaryBase::EmitFunctionProgress(const float frac)
{
#ifdef USEQT
    int qfrac = (int)(frac*100.0e0);
    emit FunctionProgress(qfrac);
#endif
}


void QtIntermediaryBase::EmitFunctionProgress(const double frac)
{
#ifdef USEQT
    int qfrac = (int)(frac*100.0e0);
    emit FunctionProgress(qfrac);
#endif
}


void QtIntermediaryBase::EmitFunctionProgress(const float frac, const std::string desc)
{
#ifdef USEQT
    int qfrac = (int)(frac*100.0e0);
    QString qdesc = QString::fromStdString(desc);
    emit FunctionProgress(qfrac,qdesc);
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
#endif
}


void QtIntermediaryBase::EmitFunctionDesc2(std::string desc)
{
#ifdef USEQT
    QString qdesc = QString::fromStdString(desc);
    emit FunctionLabel2(qdesc);
#endif
}


void QtIntermediaryBase::EmitDone()
{
#ifdef USEQT
    emit done();
#endif
}
