/**
 * @file StringManip.h
 * @author 	Robert Grandin
 * @date 24 November 2010
 * @brief Definition @e and implementation of StringManip namespace.
 *
 * @section Description & Notes
 *
 * This header file contains functions useful in manipulating string data.
 *
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 24 November 2010
 *	- Creation date.
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2010-2013, Robert Grandin
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

#ifndef StringManip_
#define StringManip_

#include <algorithm>
#include <cctype>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <string.h>





/**
 * @brief Namespace for functions used to manipulate strings.
 */
namespace StringManip {

/** 
 * @brief Convert boolean variable state to string for output.
 * @pre Boolean variable is set to true or false.
 * @param b Boolean variable to be converted.
 * @post No changes to data.
 * @return String of "true" or "false".
 */
std::string BoolToString(bool b);


/**
 * @brief Convert boolean variable state to integer for output.
 * @pre Boolean variable is set to true or false.
 * @param b Boolean variable to be converted.
 * @post No changes to data.
 * @return Integer of 0 or 1 corresponding to "False" or "True".
 */
int BoolToInt(bool b);


/**
 * @brief Convert integer to boolean.
 * @pre Integer is 1 or 0.
 * @param i Integer to be converted.  0:False, 1:True
 * @post No changes to data.
 * @return Boolean of "True" or "False".
 */
bool IntToBool(int i);


/**
 * @brief Determine the unsigned integer represented by a string (i.e., convert
 * 			the string "561" to the integer 561).
 * @pre Input string has been set.
 * @param str String to be converted to an integer.
 * @post No changes to data.
 * @return Integer corresponding to the input string.
 * @warning If the returned value is equal to the maximum value of an integer, an
 *  error may have occurred.
 */
int StrToUInt(std::string str);


/**
 * @brief Generate a string expressing a time in a nice HH:MM:SS.sss format.
 * @pre Input value, in seconds, has been determined.
 * @param time Time, in seconds, to be processed.
 * @post Formatted output generated.  Input time is not modified.
 * @return String expressing time in HH:MM:SS.sss format.
 */
std::string FormatTime(double time);



/**
 * @brief Determine the extension of a file.
 * @pre Input variable has been assigned a value.
 * @param file Name of file for which extension is to be determined.
 * @post No variables modified.
 * @return String containing the extension, without the ".".
 */
std::string DetermFileExt(std::string file);


/**
 * @brief Estimate the memory required to store a number of points of a given
 * 		datatype.
 * @pre None.
 * @param npts Number of points to be stored in memory.  This is double-precision
 * 		to reduce possible memory overflow errors.
 * @param size Size of each data point to be stored.  This is the value retruned
 * 		by the sizeof() command.
 * @post No variables modified.
 * @return String containing the estimated memory usage expressed in common sizes.
 * 		Note that the SI-esque prefixes used are binary (1024-based) and not the
 * 		usual (1000-based).
 */
std::string DetermMemReq(double npts, size_t size);



/**
 * @brief Estimate the memory required to store a number of points of a given
 * 		datatype.
 * @pre None.
 * @param npts Number of points to be stored in memory.  This is double-precision
 * 		to reduce possible memory overflow errors.
 * @param size Size of each data point to be stored.  This is the value retruned
 * 		by the sizeof() command.
 * @post No variables modified.
 * @return String containing the estimated memory usage expressed in common sizes.
 * 		Note that the SI-esque prefixes used are binary (1024-based) and not the
 * 		usual (1000-based).
 */
std::string DetermMemReq(size_t npts, size_t size);


/**
 * @brief Express the number of Bytes in a human-readable form using SI prefixes.
 * @pre None.
 * @param nbytes Number of Bytes which is to be converted to a human-friendly format.
 * @post No variables modified.
 * @return String containing the Bytes expressed in common sizes.
 * 		Note that the SI-esque prefixes used are binary (1024-based) and not
 * 		decimal (1000-based).
 */
std::string DetermBytesLabel(double nbytes);


/**
  @brief Convert integer to string.
  @pre None.
  @param val Integer value to be converted.
  @post No variables modified.
  @return String containing integer.
*/
template <class T>
std::string NumToStr(const T val)
{
    std::stringstream tmpss;
    tmpss.clear(); tmpss.str("");
    tmpss << val;
    return tmpss.str();
}


/**
  @brief Convert integer to string with specified minimum number of digits.
  @pre None.
  @param val Integer to be converted.
  @param mindigits Minimum number of digits to be used.  Leading zeros used as-necessary.
  @return String containing integer.
*/
std::string NumToStr(const int val, const int mindigits);



/**
 * @brief FormattedNumber returns the string representation of a formatted number.
 * @param val Value to be formated.
 * @param width Minimum width of field.  If val contains more characters than width, the full
 *      val will be shown.  If val contains fewer characters, the difference will be filled
 *      with fill.
 * @param precision Number of decimal points to be displayed.
 * @param fill Fill character to be used, if necessary.
 * @param showpos Show '+' in front of positive numbers.
 * @return String of val, formatted as specified.
 */
template <class T>
std::string FormattedNumber(const T val, const int width, const int precision, const char fill,
                            const bool showpos)
{
    std::stringstream tmpss;
    tmpss << std::setw(width);
    tmpss << std::setprecision(precision);
    tmpss << std::setfill(fill);
    if(showpos == true){
        tmpss << std::setiosflags(std::ios::fixed | std::ios::showpos);
    } else {
        tmpss << std::setiosflags(std::ios::fixed);
    }
    tmpss << val;

    if(width == 0){
        tmpss.str("");
    }

    return tmpss.str();
}


/**
 * @brief DetermFileStem determines the stem of a filename, size of numeric index used,
 *  the numeric index of the provided filename, and extension of the file.
 * @param filename Name of file for which the stem is to be determined.  It is expected that
 *  this is just a file name, without a path.
 * @param stem Stem of 'filename'.
 * @param ndigits Number of digits used for the index number.
 * @param ext Extension of filename.
 * @param index_value Numerical value of index.
 */
void DetermFileStem(const std::string &filename, std::string &stem, size_t &ndigits,
                    std::string &ext, size_t &index_value);


/**
 * @brief Convert string to all uppercase letters.
 * @param s String to be converted.
 */
void uppercaseString(std::string &s);


/**
 * @brief Trim leading whitespace from string.
 *
 *  Copied from http://stackoverflow.com/a/217605 on 8 March 2013.
 * @param s String to be trimmed.
 * @return Trimmed version of string.
 */
std::string& ltrim(std::string &s);


/**
 * @brief Trim trailing whitespace from string.
 *
 *  Copied from http://stackoverflow.com/a/217605 on 8 March 2013.
 * @param s String to be trimmed.
 * @return Trimmed version of string.
 */
std::string& rtrim(std::string &s);


/**
 * @brief Trim whitespace from front <b>and</b> back of string.
 *
 *  Copied from http://stackoverflow.com/a/217605 on 8 March 2013.
 * @param s String to be trimmed.
 * @return Trimmed version of string.
 */
std::string& trim(std::string &s);



/**
 * @brief DetermNumElements counts the number of elements in a std::stringstream.
 * @param stream std::stringstream for which the number of elements is to be determined.
 * @return Number of elements in stream.
 * @warning It is assumed that whitespace is the separating character between elements.
 */
int DetermNumElements(std::stringstream &stream);


/**
 * @brief DetermNumElements counts the number of elements in a std::stringstream.
 * @param stream std::stringstream for which the number of elements is to be determined.
 * @param Delimiting character to be used.
 * @return Number of elements in stream.
 */
int DetermNumElements(std::stringstream &stream, const char delim);


/**
 * @brief Compare two strings for equality.  Case insensitive.  This is a wrapper
 *  to allow cleaner cross-platform code.
 * @param a String 1.
 * @param b String 2.
 * @return Comparison result.  '0' indicates strings are equal.  Any other value
 *  indicates strings were not equal.
 */
int str_compare(const char *a, const char *b);


/**
 * @brief Sanitize a string by replacing "\" with "/" and enclosing the entire string
 *  in quotes to allow spaces.
 * @param input String to be sanitized.
 * @return Sanitized string.
 */
std::string SanitizeString(std::string &input);



} /* StringManip Namespace */

#endif /* StringManip_ */
