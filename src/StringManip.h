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
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <sstream>
#include <math.h>
#include <limits>




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
inline
std::string BoolToString(bool b)
{
	std::string retval("ERR");
	if(b == true){retval = "true";}
	if(b == false){retval = "false";}
	return retval;
}


/**
 * @brief Convert boolean variable state to integer for output.
 * @pre Boolean variable is set to true or false.
 * @param b Boolean variable to be converted.
 * @post No changes to data.
 * @return Integer of 0 or 1 corresponding to "False" or "True".
 */
inline
int BoolToInt(bool b)
{
	int retval = -1;
	if(b == true){retval = 1;}
	if(b == false){retval = 0;}
	return retval;
}


/**
 * @brief Convert integer to boolean.
 * @pre Integer is 1 or 0.
 * @param i Integer to be converted.  0:False, 1:True
 * @post No changes to data.
 * @return Boolean of "True" or "False".
 */
inline
bool IntToBool(int i)
{
	bool retval = false;
	if(i == 1){retval = true;}
	if(i == 0){retval = false;}
	return retval;
}


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
inline
int StrToUInt(std::string str)
{
    int n = (int)str.length();                      // GET THE LENGTH OF THE INPUT STRING
    int retval = std::numeric_limits<int>::max();   // VALUE TO BE RETURNED TO CALLING FUNCTION

    int powerten = 0;
    for(int i=n-1; i>=0; i--){
		// DETERMINE POWER OF TEN ASSOCIATED WITH THIS POSITION
        powerten *= 10;

        if(i == n-1){
            powerten = 1;
        }

		// LOOP THROUGH THE NUMBER, FROM RIGHT TO LEFT
		std::string sub;
		sub = str.substr(i,1);	// GET CHARACTER

		// COMPARE TO DIGITS TO DETERMINE MATCH
		bool match = false;
		for(int j=0; j<10; j++){	// LOOP THROUGH POSSIBLE DIGITS
			if(match == false){	// CONDITIONAL ON NO PREVIOUS MATCH
				std::stringstream ssj;
				std::string sj;
				int comparevalue;

				// PUT j INTO STRINGSTREAM FOR COVERSION TO STRING
				ssj << j;
				sj = ssj.str();

				// COMPARE TO SUBSTRING FROM INPUT STRING
				comparevalue = sub.compare(sj);

				if(comparevalue == 0){	// MATCH FOUND
					// CHANGE FLAG STATUS
					match = true;

                    if(retval == std::numeric_limits<int>::max()){
                        retval = 0;
                    }

					// ADD APPROPRIATE VALUE TO retval.  pow FUNCTION REQUIRES
					// INPUT ARGUMENTS TO BE float
                    retval += j*powerten;
				}
			}
		}
	}

	return(retval);
}


/**
 * @brief Generate a string expressing a time in a nice HH:MM:SS.sss format.
 * @pre Input value, in seconds, has been determined.
 * @param time Time, in seconds, to be processed.
 * @post Formatted output generated.  Input time is not modified.
 * @return String expressing time in HH:MM:SS.sss format.
 */
inline
std::string FormatTime(double time)
{
	int hours = 0;
	int mins = 0;
	double secs = 0.0;
	double remainder;
	std::string nicetime("ERROR");
	std::stringstream ss;
	char ctime[50];

	hours = (int)(floor(time/3600.0));
	remainder = time - (double)hours*3600.0;
	mins = (int)(floor(remainder/60.0));
	secs = remainder - (double)mins*60.0;

	sprintf(ctime,"%02d:%02d:%06.3lf",hours,mins,secs);
	ss << ctime;
	nicetime = ss.str();

	return nicetime;

}



/**
 * @brief Determine the extension of a file.
 * @pre Input variable has been assigned a value.
 * @param file Name of file for which extension is to be determined.
 * @post No variables modified.
 * @return String containing the extension, without the ".".
 */
inline
std::string DetermFileExt(std::string file)
{
    size_t index1 = file.find_last_of(".");
	std::string filenameext;
	filenameext = file.substr(index1+1);

	return filenameext;
}


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
inline
std::string DetermMemReq(double npts, size_t size)
{
	double nbytes = 0;

	nbytes = npts*(double)size;

	std::string output;
	std::stringstream ss;
	output = "ERR";

	/*
	 * DETERMINE APPROPRIATE POWER OF 10 TO MAKE nbytes EASILY READABLE.  THE
	 * VALUE WILL BE DISPLAYED IN THE MAXIMUM ALLOWABLE POWER OF 10 SUBJECT TO
	 * THE FOLLOWING CONSTRAINTS:
	 * 		-THE LEADING DIGIT IS NON-ZERO (I.E., 1.054 IS VALID, 0.1054 IS NOT)
	 * 		-THE INCREMENT BETWEEN POWERS OF 10 IS 3 TO MATCH THE STANDARD
	 * 		 QUASI-SCIENTIFIC PREFIXES (1024-BASED, NOT 1000-BASED).
	 *
	 */
	double lowerlim;
	double upperlim;
	for(int i=0; i<6; i++){
		lowerlim = 1.0e0;
		upperlim = 1024.0e0;
		for(int j=0; j<i; j++){
			lowerlim *= 1024.0e0;
			upperlim *= 1024.0e0;
		}

		// CHECK IF f_spacing IS WITHIN BOUNDS
		if(nbytes >= lowerlim && nbytes < upperlim){
			// RESET OUTPUT STRING
            output = "---";

			// DETERMINE "CORRECTED" VALUE
			double f_corr = nbytes/lowerlim;

			// PLACE "CORRECTED" VALUE INTO STRINGSTREAM WITH APPROPRIATE LABEL
			ss << f_corr;

			if(i == 0){
				ss << " B";
			}
			if(i == 1){
                ss << " kiB";
			}
			if(i == 2){
                ss << " MiB";
			}
			if(i == 3){
                ss << " GiB";
			}
			if(i == 4){
                ss << " TiB";
			}
			if(i == 5){
                ss << " PiB";
			}
			if(i == 6){
                ss << " EiB";
			}
		}
	}

	// CONVERT STRINGSTREAM TO STRING FOR OUTPUT
	output = ss.str();

	return output;
}



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
inline
std::string DetermMemReq(size_t npts, size_t size)
{
    return StringManip::DetermMemReq((double)npts,size);
}


/**
 * @brief Express the number of Bytes in a human-readable form using SI prefixes.
 * @pre None.
 * @param nbytes Number of Bytes which is to be converted to a human-friendly format.
 * @post No variables modified.
 * @return String containing the Bytes expressed in common sizes.
 * 		Note that the SI-esque prefixes used are binary (1024-based) and not
 * 		decimal (1000-based).
 */
inline
std::string DetermBytesLabel(double nbytes)
{
	std::string output;
	std::stringstream ss;
	output = "ERR";

	/*
	 * DETERMINE APPROPRIATE POWER OF 10 TO MAKE nbytes EASILY READABLE.  THE
	 * VALUE WILL BE DISPLAYED IN THE MAXIMUM ALLOWABLE POWER OF 10 SUBJECT TO
	 * THE FOLLOWING CONSTRAINTS:
	 * 		-THE LEADING DIGIT IS NON-ZERO (I.E., 1.054 IS VALID, 0.1054 IS NOT)
	 * 		-THE INCREMENT BETWEEN POWERS OF 10 IS 3 TO MATCH THE STANDARD
	 * 		 QUASI-SCIENTIFIC PREFIXES (1024-BASED, NOT 1000-BASED).
	 *
	 */
	double lowerlim;
	double upperlim;
	for(int i=0; i<6; i++){
		lowerlim = 1.0e0;
		upperlim = 1024.0e0;
		for(int j=0; j<i; j++){
			lowerlim *= 1024.0e0;
			upperlim *= 1024.0e0;
		}

		// CHECK IF f_spacing IS WITHIN BOUNDS
		if(nbytes >= lowerlim && nbytes < upperlim){
			// RESET OUTPUT STRING
			output = "";

			// DETERMINE "CORRECTED" VALUE
			double f_corr = nbytes/lowerlim;

			// PLACE "CORRECTED" VALUE INTO STRINGSTREAM WITH APPROPRIATE LABEL
			ss << f_corr;

			if(i == 0){
				ss << " B";
			}
			if(i == 1){
				ss << " kB";
			}
			if(i == 2){
				ss << " MB";
			}
			if(i == 3){
				ss << " GB";
			}
			if(i == 4){
				ss << " TB";
			}
			if(i == 5){
				ss << " PB";
			}
			if(i == 6){
				ss << " EB";
			}
		}
	}

	// CONVERT STRINGSTREAM TO STRING FOR OUTPUT
	output = ss.str();

	return output;
}


/**
  @brief Convert integer to string.
  @pre None.
  @param val Integer value to be converted.
  @post No variables modified.
  @return String containing integer.
*/
inline
std::string NumToStr(const int val)
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
inline
std::string NumToStr(const int val, const int mindigits)
{
    std::stringstream tmpss;
    tmpss.clear(); tmpss.str("");

    int digitsreqd = 1 + (int)log10((double)val);

    for(int i=digitsreqd; i<mindigits; i++){
        tmpss << "0";
    }

    tmpss << val;
    return tmpss.str();
}


/**
  @brief Convert float to string.
  @pre None.
  @param val Float value to be converted.
  @post No variables modified.
  @return String containing integer.
*/
inline
std::string NumToStr(const float val)
{
    std::stringstream tmpss;
    tmpss.clear(); tmpss.str("");
    tmpss << val;
    return tmpss.str();
}


/**
  @brief Convert double to string.
  @pre None.
  @param val Double value to be converted.
  @post No variables modified.
  @return String containing integer.
*/
inline
std::string NumToStr(const double val)
{
    std::stringstream tmpss;
    tmpss.clear(); tmpss.str("");
    tmpss << val;
    return tmpss.str();
}


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
inline
std::string FormattedNumber(const int val, const int width, const int precision, const char fill,
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

    return tmpss.str();
}


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
inline
std::string FormattedNumber(const float val, const int width, const int precision, const char fill,
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

    return tmpss.str();
}


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
inline
void DetermFileStem(const std::string &filename, std::string &stem, size_t &ndigits,
                    std::string &ext, size_t &index_value)
{
    size_t index_ext = filename.find_last_of(".");
    ext = filename.substr(index_ext+1);

    std::string tmpstem;
    tmpstem = filename.substr(0, index_ext);

    if(ext == "VTK" || ext == "vtk"){
        /* Filename may contain "_BigEndian" if it was created by UniformVolume class. */
        std::string tmp2;
        size_t index_us = filename.find_last_of("_");

        if(index_ext != filename.npos){
            tmp2 = filename.substr(index_us,(index_ext-index_us));
            if(tmp2 == "_BigEndian"){
                tmpstem = filename.substr(0, index_us);
            }
        }
    }


    size_t tmpsize = tmpstem.length();

    bool character_found = false;
    size_t idx = 0;
    size_t maxint = (size_t)std::numeric_limits<int>::max();
    while(!character_found){
        std::string sub;
        sub = tmpstem.substr(tmpsize-idx-1, 1);

        size_t digit = (size_t)StringManip::StrToUInt(sub);

        if(digit == maxint){
            character_found = true;
        } else {
            idx++;
        }
    }

    ndigits = idx;
    stem = filename.substr(0,(tmpsize-ndigits));

    index_value = (size_t)StringManip::StrToUInt(tmpstem.substr(tmpsize-ndigits));
}


/**
 * @brief Convert string to all uppercase letters.
 * @param s String to be converted.
 */
inline
void uppercaseString(std::string &s)
{
    for(int i=0; s[i]!='\0'; i++)
    {
        s[i] = toupper(s[i]);
    }
}


/**
 * @brief Trim leading whitespace from string.
 *
 *  Copied from http://stackoverflow.com/a/217605 on 8 March 2013.
 * @param s String to be trimmed.
 * @return Trimmed version of string.
 */
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}


/**
 * @brief Trim trailing whitespace from string.
 *
 *  Copied from http://stackoverflow.com/a/217605 on 8 March 2013.
 * @param s String to be trimmed.
 * @return Trimmed version of string.
 */
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}


/**
 * @brief Trim whitespace from front <b>and</b> back of string.
 *
 *  Copied from http://stackoverflow.com/a/217605 on 8 March 2013.
 * @param s String to be trimmed.
 * @return Trimmed version of string.
 */
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}



/**
 * @brief DetermNumElements counts the number of elements in a std::stringstream.
 * @param stream std::stringstream for which the number of elements is to be determined.
 * @return Number of elements in stream.
 */
static int DetermNumElements(std::stringstream &stream)
{
    /* Save current state of stream into string. */
    std::string originalstream;
    originalstream = stream.str();

    /* Count elements in stream. */
    std::string word;
    int count = 0;

    while(getline(stream,word,' ')){
        count++;
    }

    /* Set stream to original value using saved string. */
    stream.str(originalstream);

    /* Return number of elements in stream. */
    return count;
}



} /* StringManip Namespace */

#endif /* StringManip_ */
