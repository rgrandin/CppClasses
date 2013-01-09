/**
 * @file QtIntermediary.h
 * @author Robert Grandin
 * @date 5 Feb 2012
 * @brief Definition of QtIntermediary class.
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
 * @date 5 February 2012
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



#ifndef QtIntermediary_
#define QtIntermediary_

#include <QtIntermediaryBase.h>

#ifdef USEQT

/**
  @brief Class to provide Qt interoperability with other classes.

   This class provides no
    additional features with respect to QtIntermediaryBase.  As such, it is an attempt to
    further-hide/abstract Qt functionality.
 */
class QtIntermediary : public QtIntermediaryBase
{

public:
    /**
      @brief Default constructor with Qt disabled.
      @pre Sufficient memory exists for object.
      @post QtIntermediary object created.  Default values are
          - Spatial extents: -1.0 to 1.0 in each direction
          - 1 scalar quantity, 0 vector quantities
      @return None.
    */
    explicit QtIntermediary(QObject *parent = 0);

#else
/**
  @brief Class to provide Qt interoperability with other classes.
 */
class QtIntermediary
{
public:
    /**
      @brief Default constructor with Qt disabled.
      @pre Sufficient memory exists for object.
      @post QtIntermediary object created.  Default values are
          - Spatial extents: -1.0 to 1.0 in each direction
          - 1 scalar quantity, 0 vector quantities
      @return None.
    */
    QtIntermediary();
#endif


    /**
      @brief Destructor.
      @pre QtIntermediary object exists.
      @post Object destoryed.
      @return None.
     */
    virtual ~QtIntermediary();


};


#endif /* QtIntermediary_ */
