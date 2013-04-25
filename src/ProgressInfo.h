/**
 * @file ProgressInfo.h
 * @author Robert Grandin
 * @date 19 Nov 2012
 * @brief Definition of ProgressInfo class.
 *
 * @section Class Description & Notes
 *
 * This class contains information about function execution progress.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 *
 * @section Revisions
 *
 * @date 19 November 2012
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



#ifndef ProgressInfo_
#define ProgressInfo_

#include <iostream>


/**
  @brief Class to contain function progress information.
  @warning This class contains no checks for thread-safety.  It is left to the programmer to ensure
    that the accessing and setting of member variables is done correctly.
  @warning C++11 functionality, such as the move-constructor and move-assignment, requires the symbol "CXX11"
    to be defined.
  @todo Implement thread-safety checks.
 */
class ProgressInfo
{
public:
    /**
      @brief Default constructor.

      Member variables initialized as follows:
        - completion_progress: 0.0
        - function_description: "NA"
        - function_description2: "NA"
        - continue_calculations: true
      */
    ProgressInfo();


    /**
     * @brief Copy constructor.
     * @param a Reference to existing ProgressInfo object to be copied.
     */
    ProgressInfo(const ProgressInfo &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Rvalue to existing ProgressInfo object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    ProgressInfo(ProgressInfo &&a);
#endif


    /**
      @brief Destructor.
     */
    virtual ~ProgressInfo();


    /**
     * @brief Copy-assignment operator.
     * @param a Reference to ProgressInfo object being assigned.
     * @return Reference to instance of ProgressInfo.
     */
    ProgressInfo& operator=(ProgressInfo &a);


    /** @brief Function completion progress, expressed as a decimal between 0.0 and 1.0. */
    float completion_progress;

    /** @brief Function description. */
    std::string function_description;

    /** @brief Alternate/additional function description. */
    std::string function_description2;

    /** @brief Elapsed time since function began. */
    double elapsed_time;



private:

    /**
     * @brief ProgressInfoSwap swaps member information between two ProgressInfo objects.
     * @param first First ProgressInfo object.
     * @param second Second ProgressInfo object.
     */
    friend void ProgressInfoSwap(ProgressInfo &first, ProgressInfo &second)
    {
        std::swap(first.completion_progress, second.completion_progress);
        std::swap(first.function_description, second.function_description);
        std::swap(first.function_description2, second.function_description2);
    }



};


#endif /* ProgressInfo_ */
