/**
 * @file UnitConvert.h
 * @author Robert Grandin
 * @date 4 June 2011
 * @brief Definition of UnitConvert class.
 *
 * @section Class Description & Notes
 *
 * This class contains routines for converting values between unit systems.
 * Case-sensitive strings are used to identify incoming and outgoing units.  The
 * following unit identifiers are valid:
 * 	- Metric Units
 * 		- "m", meters
 * 		- "g", grams
 * 		- "N", Newtons
 * 		- "dyn", dynes
 * 		- "C", degrees Celcius
 * 		- "A", Ampere
 * 		- "K", Kelvin
 * 		- "J", Joules
 * 		- "erg", ergs
 * 		- "cal", small (gram) calories
 * 		- "Cal", large (dietary) calories
 * 		- "eV", electron-volts
 * 		- "W", Watts
 * 		- "N_m", Newton-meters (moments/torque)
 * 			- Note the _ between the "N" and the "m"
 * 		- "ha", hectare
 * 		- "mol", mole
 * 		- "L", liter
 * 		- "Pa", Pascals
 * 		- "AU", astronomical-unit
 * 		- "ly", light-year
 * 		- "pc", parsec
 * 		- "R", Roentgen
 * 		- "rad", radiation dose
 * 		- "rem", rad-equivalent-man (rad <--> rem assumes Q = 1)
 * 		- "Gy", Grays
 * 		- "Sv", Sieverts
 * 		- "Ci", Curies
 * 		- "Bq", becquerels
 *
 * 	- English Units
 * 		- "in", inches
 * 		- "mil", 1/1000th of an inch
 * 		- "ft", feet
 * 		- "mile", miles
 * 		- "slug", slugs
 * 		- "lbf", pounds-force
 * 		- "hp", horsepower (mechanical, ~746 W)
 * 		- "ft_lbf", Foot-pounds-force (moments/torque)
 * 			- Note the space between the "ft" and the "lbf"
 * 		- "BTU", British Thermal Units
 * 		- "F", degrees Fahrenheight
 * 		- "psi", pounds per square inch
 * 		- "ksi", kilopounds per square inch
 * 		- "atm", atmospheres
 * 		- "acre", acres
 * 		- "gee", g's (acceleration)
 * 		- "bu", bushel
 * 		- "bbl", oil barrel
 * 		- "gal", US gallons
 * 		- "imp_gal", imperial gallons
 * 			- Note the space between the "imp" and the "gal"
 *
 * 	- Not associated with Metric/English
 * 		- "deg", degrees
 * 		- "rad", radians
 * 		- "sec", seconds
 * 		- "min", minutes
 * 		- "hr", hours
 * 		- "day", days
 *
 * 	- Biblical (from http://www.biblestudy.org/beginner/bible-weights-and-measures.html)
 * 		- "span", span
 *		- "cubit", cubits (definitions vary.  18 inches assumed here)
 *		- "sdj", Sabbath day's journey
 *	 	- "ephah", ephahs
 *		- "omer", omer
 *		- "homer", homers
 *		- "hin", hin
 *		- "bath", baths
 *		- "shek", sheckels
 *		- "gerah", gerahs
 *		- "mina", minas
 *		- "bpound", pounds (Biblical)
 *		- "tal", talents
 *
 * Additionally, the standard SI prefixes can be combined with any unit (metric
 * and English):
 * 	- "y", yocto: 10^-24
 * 	- "z", zepto: 10^-21
 * 	- "a", atto: 10^-18
 * 	- "f", femto: 10^-15
 * 	- "p", pico: 10^-12
 * 	- "n", nano: 10^-9
 * 	- "u", micro: 10^-6 (standard notation is the Greek letter @f$ \mu @f$, approximated
 * 	  here as "u")
 * 	- "m", milli: 10^-3
 * 	- "c", centi: 10^-2
 * 	- "d", deci: 10^-1
 * 	- "da", deca: 10^1
 * 	- "h", hecto: 10^2
 * 	- "k", kilo: 10^3
 * 	- "M", mega: 10^6
 * 	- "G", giga: 10^9
 * 	- "T", tera: 10^12
 * 	- "P", peta: 10^15
 * 	- "E", exa: 10^18
 * 	- "Z", zetta: 10^21
 * 	- "Y", yotta: 10^24
 *
 * Combining these prefixes with unit labels should be done as "prefix:unit", with
 * a colon separating the two sets of characters.  For example:
 * 	- "n:m:1" is nanometers
 * 	- "m:in:1" is milli-inches (also known as "mils")
 * 	- "k:g:1" is kilograms
 *
 * For cases where a prefix is not required, the period must still preceed the
 * unit label, althogh the prefix may be ommitted.  Additionally, a prefix of "0"
 * (or any other characters not matching SI prefixes) may be used.  For example:
 * 	- ":m:1" is meters
 * 	- "0:m:1" is also meters
 * 	- ":min:1" is minutes
 * 	- ":ft:1" is feet
 *
 * A power can be designated with a unit by appending the "prefix:unit" with
 * ":power".  The power designation is required, and a value of 1 is used to
 * indicate that just the base-unit is desired.  For example:
 * 	- ":ft:0" is feet
 * 	- ":ft:2" is feet-squared
 * 	- "c:m:2" is centimeters-squared
 * 	- ":ft:2" is feet-squared
 * 	- ":ft:-3" is 1/(feet-cubed)
 *
 * Multiple units can be combined, with each separated by a pipe: "|".  Multiplication
 * is assumed and division must be noted with negative power.  For example:
 * 	- ":N:1|c:m:-2" is Newtons per square centimeter
 * 	- ":N:1|:m:-2" is Newtons per square meter (also known as Pascals)
 *
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 4 June 2011
 *	- Creation date.
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2011, Robert Grandin
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

#ifndef PI_
#define PI_

/** @brief Define mathematical constant @f$ \pi @f$ */
#define PI (T)3.14159265358979e0

#endif /* PI_ */


#ifndef UnitConvert_

/** @brief Define the unit-conversion classes. */
#define UnitConvert_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <Array1D.h>


/**
 * @brief Class defining a particular unit.  This class is intended to simply
 * hold the data and only be accessed from the UnitConvert object.  All members
 * are public and can be directly-accessed without Set/Get functions.
 */
template <class T> class Unit {

public:
	/** @brief Unit label. */
	std::string label;

	/** @brief Power on unit label. */
	T power;

	/** @brief Conversion factor. */
	T conversionfactor;

	/** @brief Prefix value. */
	T prefix;

	Unit()
	{
		label = "unit";
		power = 1.0e0;
		conversionfactor = 1.0e0;
		prefix = 1.0e0;
	}
};



/**
 * @brief Class for converting units.
 *
 * This class contains routines for converting values between unit systems.
 * Case-sensitive strings are used to identify incoming and outgoing units.  The
 * following unit identifiers are valid:
 * 	- Metric Units
 * 		- "m", meters
 * 		- "g", grams
 * 		- "N", Newtons
 * 		- "dyn", dynes
 * 		- "C", degrees Celcius
 * 		- "A", Ampere
 * 		- "K", Kelvin
 * 		- "J", Joules
 * 		- "erg", ergs
 * 		- "cal", small (gram) calories
 * 		- "Cal", large (dietary) calories
 * 		- "eV", electron-volts
 * 		- "W", Watts
 * 		- "N_m", Newton-meters (moments/torque)
 * 			- Note the _ between the "N" and the "m"
 * 		- "ha", hectare
 * 		- "mol", mole
 * 		- "L", liter
 * 		- "Pa", Pascals
 * 		- "AU", astronomical-unit
 * 		- "ly", light-year
 * 		- "pc", parsec
 * 		- "R", Roentgen
 * 		- "rad", radiation dose
 * 		- "rem", rad-equivalent-man (rad <--> rem assumes Q = 1)
 * 		- "Gy", Grays
 * 		- "Sv", Sieverts
 * 		- "Ci", Curies
 * 		- "Bq", becquerels
 *
 * 	- English Units
 * 		- "in", inches
 * 		- "mil", 1/1000th of an inch
 * 		- "ft", feet
 * 		- "mile", miles
 * 		- "slug", slugs
 * 		- "lbf", pounds-force
 * 		- "hp", horsepower (mechanical, ~746 W)
 * 		- "ft_lbf", Foot-pounds-force (moments/torque)
 * 			- Note the space between the "ft" and the "lbf"
 * 		- "BTU", British Thermal Units
 * 		- "F", degrees Fahrenheight
 * 		- "psi", pounds per square inch
 * 		- "ksi", kilopounds per square inch
 * 		- "atm", atmospheres
 * 		- "acre", acres
 * 		- "gee", g's (acceleration)
 * 		- "bu", bushel
 * 		- "bbl", oil barrel
 * 		- "gal", US gallons
 * 		- "imp_gal", imperial gallons
 * 			- Note the space between the "imp" and the "gal"
 *
 * 	- Not associated with Metric/English
 * 		- "deg", degrees
 * 		- "rad", radians
 * 		- "sec", seconds
 * 		- "min", minutes
 * 		- "hr", hours
 * 		- "day", days
 *
 * 	- Biblical (from http://www.biblestudy.org/beginner/bible-weights-and-measures.html)
 * 		- "span", span
 *		- "cubit", cubits (definitions vary.  18 inches assumed here)
 *		- "sdj", Sabbath day's journey
 *	 	- "ephah", ephahs
 *		- "omer", omer
 *		- "homer", homers
 *		- "hin", hin
 *		- "bath", baths
 *		- "shek", sheckels
 *		- "gerah", gerahs
 *		- "mina", minas
 *		- "bpound", pounds (Biblical)
 *		- "tal", talents
 *
 * Additionally, the standard SI prefixes can be combined with any unit (metric
 * and English):
 * 	- "y", yocto: 10^-24
 * 	- "z", zepto: 10^-21
 * 	- "a", atto: 10^-18
 * 	- "f", femto: 10^-15
 * 	- "p", pico: 10^-12
 * 	- "n", nano: 10^-9
 * 	- "u", micro: 10^-6 (standard notation is the Greek letter @f$ \mu @f$, approximated
 * 	  here as "u")
 * 	- "m", milli: 10^-3
 * 	- "c", centi: 10^-2
 * 	- "d", deci: 10^-1
 * 	- "da", deca: 10^1
 * 	- "h", hecto: 10^2
 * 	- "k", kilo: 10^3
 * 	- "M", mega: 10^6
 * 	- "G", giga: 10^9
 * 	- "T", tera: 10^12
 * 	- "P", peta: 10^15
 * 	- "E", exa: 10^18
 * 	- "Z", zetta: 10^21
 * 	- "Y", yotta: 10^24
 *
 * Combining these prefixes with unit labels should be done as "prefix:unit", with
 * a colon separating the two sets of characters.  For example:
 * 	- "n:m:1" is nanometers
 * 	- "m:in:1" is milli-inches (also known as "mils")
 * 	- "k:g:1" is kilograms
 *
 * For cases where a prefix is not required, the period must still preceed the
 * unit label, althogh the prefix may be ommitted.  Additionally, a prefix of "0"
 * (or any other characters not matching SI prefixes) may be used.  For example:
 * 	- ":m:1" is meters
 * 	- "0:m:1" is also meters
 * 	- ":min:1" is minutes
 * 	- ":ft:1" is feet
 *
 * A power can be designated with a unit by appending the "prefix:unit" with
 * ":power".  The power designation is required, and a value of 1 is used to
 * indicate that just the base-unit is desired.  For example:
 * 	- ":ft:0" is feet
 * 	- ":ft:2" is feet-squared
 * 	- "c:m:2" is centimeters-squared
 * 	- ":ft:2" is feet-squared
 * 	- ":ft:-3" is 1/(feet-cubed)
 *
 * Multiple units can be combined, with each separated by a pipe: "|".  Multiplication
 * is assumed and division must be noted with negative power.  For example:
 * 	- ":N:1|c:m:-2" is Newtons per square centimeter
 * 	- ":N:1|:m:-2" is Newtons per square meter (also known as Pascals)
 */
template <class T> class UnitConvert {

public:
	/**
	 * @brief Constructor.
	 * @pre Sufficient memory exists.
	 * @post UnitConvert object created.
	 * @return None.
	 */
	UnitConvert();


	/**
	 * @brief Deconstructor.
	 * @pre UnitConvert object exists.
	 * @post UnitConvert object destroyed.
	 * @return None.
	 */
    virtual ~UnitConvert();


	/**
	 * @brief Convert value from one unit to another.
	 * @pre UnitConvert object exists.
	 * @param val Value to be converted.
	 * @param unit1 Units of input value.
	 * @param unit2 Units of output value.
	 * @warning A value of 0.0 is returned if the specified unit conversion
	 * 		is not physically possible (e.g., "m" to "kg").
	 * @post No object attributes changed.
	 * @return Value transformed into output units.
	 */
	T ConvertUnits(T val, std::string unit1, std::string unit2);


	/**
	 * @brief Print valid units and prefixes for reference.
	 * @pre UnitConvert object exists.
	 * @post No object attributes changed.
	 * @return String containing the units and keywords.
	 */
	std::string PrintUnits() const;




protected:



private:
	/**
	 * @brief Initialization function to set conversion factors.
	 * @pre UnitConvert object exists.
	 * @post Conversion factors set.
	 * @return None.
	 */
	void UnitConvertInitialize();


	/**
	 * @brief Count the number of occurances of a particular character in a string.
	 * @pre UnitConvert object exists.
	 * @param fullstr Full string.
	 * @param specstr String for which the number of occurances is desired.
	 * @post No object attributes changed.
	 * @return Number of occurances of specified character.
	 */
	int GetNumOccurences(std::string fullstr, std::string specstr) const;


	/**
	 * @brief Determine the prefix value associated with the string label.
	 * @pre UnitConvert object exists.
	 * @param prefixlabel String identifying the SI prefix.
	 * @post No object attributes changed.
	 * @return Numerical prefix value (e.g., "1.0e-3 for 'milli' prefix).
	 */
	T GetPrefixValue(std::string prefixlabel) const;


	/**
	 * @brief Determine the fundamental unit represented by a unit string.
	 * @pre UnitConvert object exists.
	 * @param str Unit to be tested.
	 * @param unitresult Reference to Array1D object containing the powers of
	 * 		the fundamental units corresponding to the unit identified in str.
	 * @warning It is expected that unitresult is initialized to all 0's prior
	 * 		to this function's being called.
	 * @post No object attributes changed.
	 * @return None.
	 */
	void GetFundamentalType(std::string str, Array1D<T> &unitresult) const;


	/**
	 * @brief Determine the conversion factor required to convert a given unit
	 * 		to/from the fundamental units used in this class.
	 * @pre UnitConvert object exists.
	 * @param dir Direction of conversion.
	 * 		- 1: Supplied unit is to be converted to fundamental unit.
	 * 		- 2: Fundamental unit is to be converted to supplied unit.
	 * @param unitObj Unit object to be converted.
	 * @post No UnitConvert object parameters changed.  Unit object has its
	 * 		conversionfactor member changed.
	 * @return None.
	 */
	void SetConversionFactor(int dir, Unit<T> &unitObj);


	/*
	 * UNIT-TRACKING
	 */
	/** @brief Powers of fundamental units for input,
	 * 	[mass,length,time,temperature,quantity,current,luminosity]
	 */
	Array1D<T> unitsin;

	/** @brief Powers of fundamental units for output,
	 * 	[mass,length,time,temperature,quantity,current,luminosity]
	 */
	Array1D<T> unitsout;

    /** @brief Internally-stored value of converted quantity.  This is stored as a
     *      single number with the following fundamental units (powers of these units
     *      identified using the unitsin and unitsout arrays).
     *      - Mass: grams
     *      - Length: meters
     *      - Time: seconds
     *      - Temperature: Kelvin
     *      - Quantity: gram-moles
     *      - Current: amperes
     *      - Luminosity: candellas
     *      - Angle: radians
     */
	T internalvalue;




	/*
	 * CONVERSION FACTORS
	 */

	/** @brief Inches per meter. */
	T inchespermeter;

	/** @brief Meters per inch. */
	T metersperinch;

	/** @brief Inches per foot. */
	T inchesperfoot;

	/** @brief Feet per inch. */
	T feetperinch;

	/** @brief Feet per mile. */
	T feetpermile;

	/** @brief Miles per feet. */
	T milesperfoot;

	/** @brief Meters per light-year. */
	T metersperlightyear;

	/** @brief Light-years per meter. */
	T lightyearspermeter;

	/** @brief Meters per astronomical-unit. */
	T metersperau;

	/** @brief Astronomical-units per meter. */
	T auspermeter;

	/** @brief Meters per parsec. */
	T metersperparsec;

	/** @brief Parsecs per meter. */
	T parsecspermeter;

	/** @brief Grams per slug. */
	T gramsperslug;

	/** @brief Slugs per gram. */
	T slugspergram;

	/** @brief Pounds-force per Newton. */
	T lbfspernewton;

	/** @brief Newtons per pound-force. */
	T newtonsperlbf;

	/** @brief Seconds per minute. */
	T secspermin;

	/** @brief Minutes per second. */
	T minspersec;

	/** @brief Minutes per hour. */
	T minsperhour;

	/** @brief Hours per minute. */
	T hourspermin;

	/** @brief Hours per day. */
	T hoursperday;

	/** @brief Days per hour. */
	T daysperhour;

	/** @brief Ratio of Fahrenheit degree to Kelvin (relative size of degrees). */
	T ratioFtoK;

	/** @brief Ratio of Kelvin to Fahrenheit degree (relative size of degrees). */
	T ratioKtoF;

	/** @brief Avagadro's number. */
	T itemspermol;

	/** @brief Inverse of Avagadro's number. */
	T molsperitem;

	/** @brief Gee's in [m/s^2]. */
	T geeinms2;

	/** @brief Degrees per radian. */
	T degperrad;

	/** @brief Radians per degree. */
	T radperdeg;


	/*
	 * PREFIX VALUES
	 */
	/** @brief yocto: 10^-24 */
	T prefix_n24;

	/** @brief zepto: 10^-21 */
	T prefix_n21;

	/** @brief atto: 10^-18 */
	T prefix_n18;

	/** @brief femto: 10^-15 */
	T prefix_n15;

	/** @brief pico: 10^-12 */
	T prefix_n12;

	/** @brief nano: 10^-9 */
	T prefix_n9;

	/** @brief micro: 10^-6 */
	T prefix_n6;

	/** @brief milli: 10^-3 */
	T prefix_n3;

	/** @brief centi: 10^-2 */
	T prefix_n2;

	/** @brief deci: 10^-1 */
	T prefix_n1;

	/** @brief deca: 10^1 */
	T prefix_p1;

	/** @brief hecto: 10^2 */
	T prefix_p2;

	/** @brief kilo: 10^3 */
	T prefix_p3;

	/** @brief Mega: 10^6 */
	T prefix_p6;

	/** @brief Giga: 10^9 */
	T prefix_p9;

	/** @brief Tera: 10^12 */
	T prefix_p12;

	/** @brief Peta: 10^15 */
	T prefix_p15;

	/** @brief Exa: 10^18 */
	T prefix_p18;

	/** @brief Zetta: 10^21 */
	T prefix_p21;

	/** @brief Yotta: 10^24 */
	T prefix_p24;

};
 
#include "UnitConvert.cpp"


#endif /* UnitConvert_ */
