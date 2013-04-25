/**
 * @file QtIntermediaryBase.h
 * @author Robert Grandin
 * @date 16 Nov 2012
 * @brief Definition of QtIntermediaryBase class.
 *
 * @section Class Description & Notes
 *
 * This class serves as an intermediary between QObject and other classes.  It is
 * intended that other classes which want to make use of Qt signals will be derived
 * from this class.  If Qt is enabled, by defining USEQT at compile-time, Qt
 * functionality will be provided.  If Qt is disabled, by defining NOQT at compile-
 * time, the signal-emitting wrapper functions in this class may still be called.
 * However, no signals will actually be emitted and this class will not induce any
 * dependence upon Qt.  The use of the wrapper functions allows the classes which
 * inherit this class to be unmodified when using with either Qt or non-Qt applications.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 *
 * @section Revisions
 *
 * @date 16 November 2012
 *	- Creation date.
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2012, Robert Grandin
 * All rights reserved.
 *
 * Redistribution and use of this file is permitted provided that the following
 * conditions are met:
 * 	-# 	Redistributions must produce the above copyright notice, this list of
 * 		conditions, and the following disclaimer in the documentation and/or
 * 		other materials provided with the distribution.
 * 	-#	Neither the name of the organization nor the names of its contributors
 * 		may be used to endorse or promote products derived from this software
 * 		without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
 * EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING BUT NOT
 * LIMITING TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 */



#ifndef QtIntermediaryBase_
#define QtIntermediaryBase_

#include <stdio.h>
#include <string>

#ifdef USEQT

#include <QObject>
#include <QString>


/**
  @brief Class to provide Qt interoperability with other classes.  This class contains NO signal
        functionality.
 */
class QtIntermediaryBase : public QObject
{
    /** @brief Requried macro for QObject-derived classes providing Qt functionality. */
    Q_OBJECT

public:
    /**
      @brief Default constructor with Qt enabled.
      @return None.
    */
    explicit QtIntermediaryBase(QObject *parent = 0);

#else
/**
  @brief Class to provide Qt interoperability with other classes.
 */
class QtIntermediaryBase
{
public:
    /**
      @brief Default constructor with Qt disabled.
    */
    QtIntermediaryBase();
#endif


    /**
      @brief Destructor.
      @pre QtIntermediaryBase object exists.
      @post Object destoryed.
      @return None.
     */
    virtual ~QtIntermediaryBase();

#ifdef USEQT
public slots:
#endif
    /**
     * @brief PerformCalculations is to be overloaded in a class derived from
     *      QtIntermediaryBase.
     *
     * When using QThreads it is suggested that the thread
     *  signal 'started()' is connected to this slot so that the desired
     *  calculations are started when the thread is launched.
     */
    virtual void PerformCalculations();


    /**
     * @brief ExitThread provides a slot for derived classes which allows
     *      the setting of an execution-control flag which will cease calculations.
     */
    virtual void ExitThread();


    /**
     * @brief EnableThread provides a slot for derived classes which allows
     *      the setting of an execution-control flag which will allow calculations.
     */
    virtual void EnableThread();


    /**
      @brief Function to cause the FuctionProgress signal to be emitted.

      The emitted signal contains an integer which indicates the progress on a 0-100 scale.
      @pre QtIntermediary object exists.
      @param frac Fraction of progress, expressed as a decimal between 0.0 and 1.0.
      @post No changes to object.
      @return None.
     */
    void EmitFunctionProgress(const int frac);


    /**
      @brief Function to cause the FuctionProgress signal to be emitted.

      The emitted signal contains an integer which indicates the progress on a 0-100 scale.
      @pre QtIntermediary object exists.
      @param frac Fraction of progress, expressed as a decimal between 0.0 and 1.0.
      @post No changes to object.
      @return None.
     */
    void EmitFunctionProgress(const float frac);


    /**
      @brief Function to cause the FuctionProgress signal to be emitted.

      The emitted signal contains an integer which indicates the progress on a 0-100 scale.
      @pre QtIntermediary object exists.
      @param frac Fraction of progress, expressed as a decimal between 0.0 and 1.0.
      @post No changes to object.
      @return None.
     */
    void EmitFunctionProgress(const double frac);


    /**
      @brief Function to cause the FuctionProgress signal to be emitted.

      The emitted signal contains an integer which indicates the progress on a 0-100 scale and a
        string description.
      @pre QtIntermediary object exists.
      @param frac Fraction of progress, expressed as a decimal between 0.0 and 1.0.
      @param desc Description of the function being executed.
      @post No changes to object.
      @return None.
     */
    void EmitFunctionProgress(const float frac, const std::string desc);


    /**
     * @brief Function to cause the ElapsedTime signal to be emitted.
     *
     *  The emitted signal contains a double, which indicates the number of seconds elapsed since
     *  the function began.
     * @param time_seconds Number of seconds elapsed.
     */
    void EmitFunctionElapsedTime(const double time_seconds);


    /**
     * @brief Function to cause the ElapsedTime signal to be emitted.
     *
     *  The emitted signal contains a double, which indicates the number of seconds elapsed since
     *  the function began.
     * @param time_seconds Number of seconds elapsed.
     * @param desc Descriptive text of function.
     */
    void EmitFunctionElapsedTime2(const double time_seconds, const std::string &desc);


#ifdef USEQT
    /**
      @brief Function to cause FunctionLabel signal to be emitted.
      @pre QtIntermediary object exists.
      @param desc Description of function being executed.
      @post No changes to object.
      @return None.
    */
    void EmitFunctionDesc(QString desc);
#endif


    /**
      @brief Function to cause FunctionLabel to be emitted.
      @pre QtIntermediary object exists.
      @param desc Description of function being executed.
      @post No changes to object.
      @return None.
    */
    void EmitFunctionDesc(std::string desc);


    /**
      @brief Function to cause FunctionLabel to be emitted.

      This provides the same functionality as EmitFunctionDesc(), but provides a
        means to have multiple function descriptions emitted and connected to
        different GUI widgets.
      @pre QtIntermediary object exists.
      @param desc Description of function being executed.
      @post No changes to object.
      @return None.
    */
    void EmitFunctionDesc2(std::string desc);


    /**
     * @brief EmitDone emits signal 'done', which is intended to be used in an
     *      analogous way to 'finished'.
     */
    void EmitDone();



#ifdef USEQT


signals:
    /**
      @brief Signal containing function-progress information, intended
            for a progressbar.
      @param frac Fraction of progress, expressed as a decimal between
            0.0 and 1.0.
     */
    void FunctionProgress(int frac);


    /**
      @brief Signal containing function-progress information, intended
            for a progress dialog box.
      @param frac Fraction of progress, expressed as a decimal between
            0.0 and 1.0.
      @param desc Description of the function being executed.
     */
    void FunctionProgress(int frac, QString desc);


    /**
      @brief Signal containing a QString label describing a function.
      @param proglabel Label describing function.
    */
    void FunctionLabel(QString proglabel);


    /**
      @brief Signal containing a QString label describing a function.
      @param proglabel Label describing function.
    */
    void FunctionLabel2(QString proglabel);


    /**
      @brief Signal containing a std::string label describing a function.
      @param proglabel Label describing function.
    */
    void FunctionLabel(std::string);


    /**
     * @brief Signal containing number of seconds elapsed since function began.
     * @param time_seconds Elapsed time, in seconds.
     */
    void ElapsedTime(double time_seconds);


    /**
     * @brief Signal containing number of seconds elapsed since function began.
     * @param time_seconds Elapsed time, in seconds.
     */
    void ElapsedTime2(double time_seconds, QString desc);


    /**
     * @brief done signals that the desired action(s) have been completed.
     */
    void done();

#endif


};


#endif /* QtIntermediaryBase_ */
