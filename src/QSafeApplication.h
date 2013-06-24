/**
 * @file QSafeApplication.h
 * @author 	Robert Grandin
 * @date 17 Nov 2012
 * @brief Definition @e and implementation of QSafeApplication class.
 *
 * @section Class Description & Notes
 *
 * This class reimplements QApplication::notify to enable exception, catching.
 * If an exception is caught, a descriptive message is written to std::cerr.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 17 November 2012
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


#ifndef QSAFEAPPLICATION_H
#define QSAFEAPPLICATION_H

#include <QtWidgets/QApplication>
#include <iostream>
#include <stdio.h>


/**
 * @brief The SafeApplication class reimplements QApplication.
 *
 * This is done to redefine QApplication::notify to enable exception-catching.
 *  If an exception is caught, a descriptive message is written to std::cerr.
 */
class SafeApplication : public QApplication {
public:
    /**
     * @brief SafeApplication constructor.
     * @param argc Number of parameters passed-in on the command-line.
     * @param argv Command-line arguments.
     */
    SafeApplication(int &argc, char *argv[]) : QApplication(argc, argv) { }


    /**
     * @brief ~SafeApplication destructor.
     */
    virtual ~SafeApplication() { }


    /**
     * @brief notify reimplements QApplication::notify in order to allow handling of
     *      thrown exceptions from slots.
     * @param receiver QObject which was to receive the event.
     * @param event Event to be received.
     * @return True if exception is not caught.  False otherwise.  If false, exception is printed
     *      to std::cerr.
     */
    virtual bool notify(QObject *receiver, QEvent *event) {
        try {
            return QApplication::notify(receiver, event);
        } catch(std::exception &e) {
            std::cerr << "std::exception was caught: " << e.what() << std::endl;
        }
        return false;
    }
};

#endif // QSAFEAPPLICATION_H
