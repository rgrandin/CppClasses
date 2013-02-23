/**
 * @file UnitConvert.cpp
 * @author Robert Grandin
 * @brief Implementation of UnitConvert class.
 */


#include <UnitConvert.h>


/* ==================================================================
 * ================
 * ================    PRIVATE FUNCTIONS
 * ================
 */

template <class T>
void UnitConvert<T>::UnitConvertInitialize()
{
	/*
	 * UNIT TRACKING
	 */
	unitsin.ResetSize(8,0);
	unitsout.ResetSize(8,0);

	/*
	 * INTERNAL VALUE
	 */
    internalvalue = (T)0.0;

	/*
	 * CONVERSION FACTORS
	 */
	// LENGTH CONVERSIONS
    inchespermeter = (T)39.37e0;
    metersperinch = (T)0.0254e0;
    inchesperfoot = (T)12.0e0;
    feetperinch = (T)1.0e0/(T)12.0e0;
    feetpermile = (T)5280.0e0;
    milesperfoot = (T)1.0e0/(T)5280.0e0;
    metersperlightyear = (T)9.4605284e15;
    lightyearspermeter = (T)1.0/(T)9.4605284e15;
    metersperau = (T)149.598000e9;
    auspermeter = (T)1.0/(T)149.598000e9;
    metersperparsec = (T)3.08568025e16;
    parsecspermeter = (T)1.0/(T)3.08568025e16;

	// MASS CONVERSIONS
    slugspergram = (T)6.852e-5;
    gramsperslug = (T)14594.279e0;

	// FORCE CONVERSIONS
    lbfspernewton = (T)0.224820e0;
    newtonsperlbf = (T)4.448e0;

	// TIME CONVERSIONS
    secspermin = (T)60.0e0;
    minspersec = (T)1.0e0/(T)60.0e0;
    minsperhour = (T)60.0e0;
    hourspermin = (T)1.0e0/(T)60.0e0;
    hoursperday = (T)24.0e0;
    daysperhour = (T)1.0e0/(T)24.0e0;

	// TEMPERATURE CONVERSIONS
    ratioFtoK = (T)5.0e0/(T)9.0e0;
    ratioKtoF = (T)1.8e0;

	// AVAGADRO'S NUMBER
    itemspermol = (T)6.02214179e23;
    molsperitem = (T)1.0e0/itemspermol;

	// ACCELERATION
    geeinms2 = (T)9.80665e0;

	// ANGLES
    degperrad = (T)180.0e0/PI;
    radperdeg = PI/(T)180.0e0;


	/*
	 * PREFIXES
	 */
    prefix_n24 = (T)1.0e-24;
    prefix_n21 = (T)1.0e-21;
    prefix_n18 = (T)1.0e-18;
    prefix_n15 = (T)1.0e-15;
    prefix_n12 = (T)1.0e-12;
    prefix_n9 = (T)1.0e-9;
    prefix_n6 = (T)1.0e-6;
    prefix_n3 = (T)1.0e-3;
    prefix_n2 = (T)1.0e-2;
    prefix_n1 = (T)1.0e-1;
    prefix_p1 = (T)1.0e1;
    prefix_p2 = (T)1.0e2;
    prefix_p3 = (T)1.0e3;
    prefix_p6 = (T)1.0e6;
    prefix_p9 = (T)1.0e9;
    prefix_p12 = (T)1.0e12;
    prefix_p15 = (T)1.0e15;
    prefix_p18 = (T)1.0e18;
    prefix_p21 = (T)1.0e21;
    prefix_p24 = (T)1.0e24;
}


template <class T>
int UnitConvert<T>::GetNumOccurences(std::string fullstr, std::string specstr) const
{
	int count = 0;
	int stringsize = (int)fullstr.size();
	for(int i=0; i<stringsize; i++){
		if(fullstr.substr(i,1) == specstr){
			count++;
		}
	}
	//printf("In %s, %s appears %d times.\n",fullstr.c_str(),specstr.c_str(),count);
	return count;
}


template <class T>
T UnitConvert<T>::GetPrefixValue(std::string prefix) const
{
	T retval = 0.0;
    if(prefix == "y"){ retval = (T)1.0e-24; }
    if(prefix == "z"){ retval = (T)1.0e-21; }
    if(prefix == "a"){ retval = (T)1.0e-18; }
    if(prefix == "f"){ retval = (T)1.0e-15; }
    if(prefix == "p"){ retval = (T)1.0e-12; }
    if(prefix == "n"){ retval = (T)1.0e-9; }
    if(prefix == "u"){ retval = (T)1.0e-6; }
    if(prefix == "m"){ retval = (T)1.0e-3; }
    if(prefix == "c"){ retval = (T)1.0e-2; }
    if(prefix == "d"){ retval = (T)1.0e-1; }
    if(prefix == "da"){ retval =(T) 1.0e1; }
    if(prefix == "h"){ retval = (T)1.0e2; }
    if(prefix == "k"){ retval = (T)1.0e3; }
    if(prefix == "M"){ retval = (T)1.0e6; }
    if(prefix == "G"){ retval = (T)1.0e9; }
    if(prefix == "T"){ retval = (T)1.0e12; }
    if(prefix == "P"){ retval = (T)1.0e15; }
    if(prefix == "E"){ retval = (T)1.0e18; }
    if(prefix == "Z"){ retval = (T)1.0e21; }
    if(prefix == "Y"){ retval = (T)1.0e24; }

	return retval;
}


template <class T>
void UnitConvert<T>::GetFundamentalType(std::string str, Array1D<T> &unitresult) const
{
	/*
	 * CHECK EACH POSSIBLE VALID VALUE FOR str AND SET unitresult ACCORDINGLY
	 */

	bool valset = false;

	/*
	 * SINGLE-FUNDAMENTAL UNITS
	 */
	// LENGTH UNITS
	if((str == "m" || str == "mil" || str == "in" || str == "ft" || str == "mile" ||
	   str == "AU" || str == "ly" || str == "pc" || str == "cubit" || str == "sdj" ||
	   str == "span" || str == "fur" || str == "rod" || str == "chain" || str == "p" ||
	   str == "P" || str == "li" || str == "lea" || str == "ftm" || str == "cb" ||
	   str == "nmi" || str == "hand") && valset == false){
        unitresult(1) = (T)1.0e0;
		valset = true;
	}

	// MASS UNITS
	if((str == "g" || str == "slug" || str == "gr" || str == "dr" || str == "lbm")
		&& valset == false){
        unitresult(0) = (T)1.0e0;
		valset = true;
	}

	// TIME UNITS
	if((str == "sec" || str == "min" || str == "hr" || str == "day") && valset == false){
        unitresult(2) = (T)1.0e0;
		valset = true;
	}

	// TEMPERATURE UNITS
	if((str == "K" || str == "F" || str == "C") && valset == false){
        unitresult(3) = (T)1.0e0;
		valset = true;
	}

	// QUANTITY UNITS
	if((str == "mol") && valset == false){
        unitresult(4) = (T)1.0e0;
		valset = true;
	}

	// CURRENT UNITS
	if((str == "A") && valset == false){
        unitresult(5) = (T)1.0e0;
		valset = true;
	}

	// LUMINOSITY UNITS


	// ANGLE UNITS
	if((str == "deg" || str == "radian" || str == "grad") && valset == false){
        unitresult(7) = (T)1.0e0;
		valset = true;
	}


	/*
	 * COMPOSITE-FUNDAMENTAL UNITS
	 */

	// ACCELERATION
	if((str == "gee") && valset == false){
        unitresult(1) = (T)1.0e0;
        unitresult(2) = (T)-2.0e0;
		valset = true;
	}

	// FORCE
	if((str == "N" || str == "lbf" || str == "mite" || str == "farth" ||
	   str == "gerah" || str == "mina" || str == "penny" || str == "bpound" ||
	   str == "shek" || str == "tal" || str == "oz" || str == "dyn") && valset == false){
        unitresult(0) = (T)1.0e0;
        unitresult(1) = (T)1.0e0;
        unitresult(2) = (T)-2.0e0;
		valset = true;
	}

	// PRESSURE
	if((str == "Pa" || str == "psi" || str == "ksi" || str == "atm" || str == "at" ||
		str == "torr" || str == "bar") && valset == false){
        unitresult(0) = (T)1.0e0;
        unitresult(1) = (T)-1.0e0;
        unitresult(2) = (T)-2.0e0;
		valset = true;
	}

	// WORK/ENERGY AND MOMENT/TORQUE
	if((str == "J" || str == "N_m" || str == "ft_lbf" || str == "erg" || str == "cal" ||
		str == "Cal" || str == "BTU" || str == "eV") && valset == false){
        unitresult(0) = (T)1.0e0;
        unitresult(1) = (T)2.0e0;
        unitresult(2) = (T)-2.0e0;
		valset = true;
	}

	// POWER
	if((str == "W" || str == "hp") && valset == false){
        unitresult(0) = (T)1.0e0;
        unitresult(1) = (T)2.0e0;
        unitresult(2) = (T)-3.0e0;
		valset = true;
	}

	// AREA
	if((str == "ha" || str == "acre") && valset == false){
        unitresult(1) = (T)2.0e0;
		valset = true;
	}

	// VOLUME
	if((str == "L" || str == "bu" || str == "bbl" || str == "gal" || str == "imp_gal" ||
	   str == "bath" || str == "omer" || str == "homer" || str == "hin" ||
	   str == "ephah" || str == "lqt" || str == "lpt" || str == "cup" || str == "minim" ||
	   str == "fl_dr" || str == "tsp" || str == "Tbsp" || "fl_oz" || str == "jig" ||
	   str == "gi" || str == "lbbl" || str == "hghd" || str == "dpt" || str == "dqt" ||
	   str == "pk" || str == "dbbl") && valset == false){
        unitresult(1) = (T)3.0e0;
		valset = true;
	}

	// RADIATION (charge/mass)
	if((str == "R") && valset == false){
        unitresult(5) = (T)1.0e0;
        unitresult(2) = (T)1.0e0;
        unitresult(0) = (T)-1.0e0;
		valset = true;
	}

	// RADIATION DOSE (energy/mass)
	if((str == "rad" || str == "Gy" || str == "Sv" || str == "rem") && valset == false){
        unitresult(0) = (T)1.0e0;
        unitresult(1) = (T)1.0e0;
        unitresult(2) = (T)-2.0e0;
		valset = true;
	}

	// RADIOACTIVE DECAY (decays/time)
	if((str == "Ci" || str == "Bq") && valset == false){
        unitresult(2) = (T)-1.0e0;
		valset = true;
	}

}


template <class T>
void UnitConvert<T>::SetConversionFactor(int dir, Unit<T> &unitObj)
{

	/*
	 * DETERMINE CONVERSION FACTOR ASSUMING INPUT UNIT IS BEING CONVERTED TO
	 * FUNDAMENTAL UNITS.  THEN, IF THE CONVERSION IS SUPPOSED TO BE THE OTHER
	 * WAY USE THE INVERSE OF THE CONVERSION FACTOR.
	 */
    T convfactor = (T)1.0e0;
	std::string str(unitObj.label);
	bool valset = false;

	/*
	 * LENGTH UNITS
	 */
	if(str == "m" && valset == false){
        convfactor = (T)1.0e0;
		valset = true;
	}

	if(str == "ft" && valset == false){
		convfactor = metersperinch*inchesperfoot;
		valset = true;
	}

	if(str == "in" && valset == false){
		convfactor = metersperinch;
		valset = true;
	}

	if(str == "mil" && valset == false){
        convfactor = metersperinch*(T)1.0e-3;
		valset = true;
	}

	if(str == "mile" && valset == false){
		convfactor = metersperinch*inchesperfoot*feetpermile;
		valset = true;
	}

	if(str == "AU" && valset == false){
		convfactor = metersperau;
		valset = true;
	}

	if(str == "ly" && valset == false){
		convfactor = metersperlightyear;
		valset = true;
	}

	if(str == "pc" && valset == false){
		convfactor = metersperparsec;
		valset = true;
	}
	if(str == "fur" && valset == false){
        convfactor = metersperinch*inchesperfoot*(T)660.0e0;
		valset = true;
	}
	if(str == "rod" && valset == false){
        convfactor = metersperinch*inchesperfoot*(T)16.5e0;
		valset = true;
	}
	if(str == "chain" && valset == false){
        convfactor = metersperinch*inchesperfoot*(T)66.0e0;
		valset = true;
	}
	if(str == "p" && valset == false){
        convfactor = (T)3.528e-4;
		valset = true;
	}
	if(str == "P" && valset == false){
        convfactor = (T)4.233e-3;
		valset = true;
	}
	if(str == "li" && valset == false){
        convfactor = metersperinch*inchesperfoot*(T)33.0e0/(T)50.0e0;
		valset = true;
	}
	if(str == "lea" && valset == false){
        convfactor = metersperinch*inchesperfoot*feetpermile*(T)3.0e0;
		valset = true;
	}
	if(str == "ftm" && valset == false){
        convfactor = metersperinch*inchesperfoot*(T)6.0e0;
		valset = true;
	}
	if(str == "cb" && valset == false){
        convfactor = metersperinch*inchesperfoot*(T)6.0e0*(T)120.0e0;
		valset = true;
	}
	if(str == "nmi" && valset == false){
        convfactor = metersperinch*inchesperfoot*feetpermile*(T)1.151e0;
		valset = true;
	}
	if(str == "hand" && valset == false){
        convfactor = (T)10.16e-2;
		valset = true;
	}

	/*
	 * MASS UNITS
	 */
	if(str == "g" && valset == false){
        convfactor = (T)1.0e0;
		valset = true;
	}

	if(str == "slug" && valset == false){
		convfactor = gramsperslug;
		valset = true;
	}
	if(str == "lbm" && valset == false){
        convfactor = (T)453.59237e0;
		valset = true;
	}
	if(str == "gr" && valset == false){
        convfactor = (T)64.79891e-3;
		valset = true;
	}
	if(str == "dr" && valset == false){
        convfactor = (T)(27.0e0 + (T)11.0e0/(T)32.0e0)*(T)64.79891e-3;
		valset = true;
	}
	if(str == "dwt" && valset == false){
        convfactor = (T)24.0e0*(T)64.79891e-3;
		valset = true;
	}

	/*
	 * TIME UNITS
	 */
	if(str == "sec" && valset == false){
        convfactor = (T)1.0e0;
		valset = true;
	}
	if(str == "min" && valset == false){
		convfactor = secspermin;
		valset = true;
	}
	if(str == "hr" && valset == false){
		convfactor = secspermin*minsperhour;
		valset = true;
	}
	if(str == "day" && valset == false){
		convfactor = secspermin*minsperhour*hoursperday;
		valset = true;
	}


	/*
	 * QUANTITY UNITS
	 */
	if(str == "mol" && valset == false){
        convfactor = (T)1.0e0;
		valset = true;
	}


	/*
	 * CURRENT UNITS
	 */
	if(str == "A" && valset == false){
        convfactor = (T)1.0e0;
		valset = true;
	}


	/*
	 * LUMINOSITY UNITS
	 */



	/*
	 * ANGULAR UNITS
	 */
	if(str == "deg" && valset == false){
		convfactor = radperdeg;
		valset = true;
	}
	if(str == "radian" && valset == false){
        convfactor = (T)1.0e0;
		valset = true;
	}
	if(str == "grad" && valset == false){
        convfactor = PI/(T)200.0e0;
		valset = true;
	}

	/*
	 * TEMPERATURE UNITS
	 */
	if(str == "F" && valset == false){
		convfactor = ratioFtoK;
		valset = true;
	}
	if(str == "K" && valset == false){
        convfactor = (T)1.0e0;
		valset = true;
	}
	if(str == "C" && valset == false){
        convfactor = (T)1.0e0;
		valset = true;
	}

	/*
	 * ACCELERATION UNITS
	 */
	if(str == "gee" && valset == false){
		convfactor = geeinms2;
		valset = true;
	}


	/*
	 * FORCE UNITS
	 */
	if(str == "N" && valset == false){
        convfactor = (T)1.0e3;
		valset = true;
	}
	if(str == "lbf" && valset == false){
        convfactor = (T)1.0e3*newtonsperlbf;
		valset = true;
	}
	if(str == "oz" && valset == false){
        convfactor = (T)1.0e3*newtonsperlbf*(T)16.0e0;
		valset = true;
	}
	if(str == "dyn" && valset == false){
        convfactor = (T)1.0e3*(T)1.0e-5;
		valset = true;
	}
	if(str == "cwt" && valset == false){
        convfactor = (T)1.0e3*newtonsperlbf*(T)1.0e2;
		valset = true;
	}



	/*
	 * PRESSURE UNITS
	 */
	if(str == "Pa" && valset == false){
        convfactor = (T)1.0e3;
		valset = true;
	}
	if(str == "psi" && valset == false){
        convfactor = (T)1.0e3*inchespermeter*inchespermeter*newtonsperlbf;
		valset = true;
	}
	if(str == "ksi" && valset == false){
        convfactor = (T)1.0e6*inchespermeter*inchespermeter*newtonsperlbf;
		valset = true;
	}
	if(str == "atm" && valset == false){
        convfactor = (T)101325.0e0*(T)1.0e3;	// Defined to be 101,325 Pa
		valset = true;
	}
	if(str == "at" && valset == false){
        convfactor = (T)1.0e3*(T)98066.5e0;
		valset = true;
	}
	if(str == "torr" && valset == false){
        convfactor = (T)101325.0e0*(T)1.0e3/(T)760.0e0;
		valset = true;
	}
	if(str == "bar" && valset == false){
        convfactor = (T)1.0e3*(T)100.0e3;
		valset = true;
	}


	/*
	 * WORK/ENERGY/MOMENT/TORQUE UNITS
	 */
	if(str == "J" && valset == false){
        convfactor = (T)1.0e3;
		valset = true;
	}
	if(str == "N_m" && valset == false){
        convfactor = (T)1.0e3;
		valset = true;
	}
	if(str == "ft_lbf" && valset == false){
        convfactor = (T)1.0e3*inchesperfoot*metersperinch*newtonsperlbf;
		valset = true;
	}
	if(str == "erg" && valset == false){
        convfactor = (T)1.0e3*(T)1.0e-7;
		valset = true;
	}
	if(str == "cal" && valset == false){
        convfactor = (T)1.0e3*(T)4.184e0;
		valset = true;
	}
	if(str == "Cal" && valset == false){
        convfactor = (T)1.0e3*(T)4.184e3;
		valset = true;
	}
	if(str == "BTU" && valset == false){
        convfactor = (T)1.0e3*(T)1055.05585262e0;
		valset = true;
	}
	if(str == "eV" && valset == false){
        convfactor = (T)1.0e3*(T)1.602176565e-19;
		valset = true;
	}


	/*
	 * POWER UNITS
	 */
	if(str == "W" && valset == false){
        convfactor = (T)1.0e3;
		valset = true;
	}
	if(str == "hp" && valset == false){
        convfactor = (T)1.0e3*(T)745.699872e0;
		valset = true;
	}


	/*
	 * AREA UNITS
	 */
	if(str == "ha" && valset == false){
        convfactor = (T)1.0e4;
		valset = true;
	}
	if(str == "acre" && valset == false){
        convfactor = (T)4046.85642e0;
		valset = true;
	}


	/*
	 * VOLUME UNITS
	 */
	if(str == "L" && valset == false){
        convfactor = (T)1.0e-3;
		valset = true;
	}
	if(str == "bu" && valset == false){
        convfactor = (T)0.035239072e0;
		valset = true;
	}
	if(str == "bbl" && valset == false){
        convfactor = (T)0.158987295e0;
		valset = true;
	}
	if(str == "gal" && valset == false){
        convfactor = (T)3.78541178e-3;
		valset = true;
	}
	if(str == "imp_gal" && valset == false){
        convfactor = (T)4.54609188e-3;
		valset = true;
	}
	if(str == "lqt" && valset == false){
        convfactor = (T)3.78541178e-3*(T)0.25e0;
		valset = true;
	}
	if(str == "lpt" && valset == false){
        convfactor = (T)3.78541178e-3*(T)0.125e0;
		valset = true;
	}
	if(str == "cup" && valset == false){
        convfactor = (T)3.78541178e-3*(T)0.0625e0;
		valset = true;
	}
	if(str == "minim" && valset == false){
        convfactor = (T)61.61152e-9;
		valset = true;
	}
	if(str == "fl_dr" && valset == false){
        convfactor = (T)60.0e0*(T)61.61152e-9;
		valset = true;
	}
	if(str == "tsp" && valset == false){
        convfactor = (T)80.0e0*(T)61.61152e-9;
		valset = true;
	}
	if(str == "Tbsp" && valset == false){
        convfactor = (T)3.0e0*(T)80.0e0*(T)61.61152e-9;
		valset = true;
	}
	if(str == "fl_oz" && valset == false){
        convfactor = (T)2.0e0*(T)3.0e0*(T)80.0e0*(T)61.61152e-9;
		valset = true;
	}
	if(str == "jig" && valset == false){
        convfactor = (T)3.0e0*(T)3.0e0*(T)80.0e0*(T)61.61152e-9;
		valset = true;
	}
	if(str == "gi" && valset == false){
        convfactor = (T)4.0e0*(T)2.0e0*(T)3.0e0*(T)80.0e0*(T)61.61152e-9;
		valset = true;
	}
	if(str == "lbbl" && valset == false){
        convfactor = (T)31.5e0*(T)3.78541178e-3;
		valset = true;
	}
	if(str == "hghd" && valset == false){
        convfactor = (T)63.0e0*(T)3.78541178e-3;
		valset = true;
	}
	if(str == "dpt" && valset == false){
        convfactor = (T)5.50605e-4;
		valset = true;
	}
	if(str == "dqt" && valset == false){
        convfactor = (T)2.0e0*(T)5.50605e-4;
		valset = true;
	}
	if(str == "pk" && valset == false){
        convfactor = (T)2.0e0*(T)3.78541178e-3;
		valset = true;
	}
	if(str == "dbbl" && valset == false){
        convfactor = (T)0.115627e0;
		valset = true;
	}


	/*
	 * BIBLICAL UNITS
	 */
	if(str == "bath" && valset == false){
        convfactor = (T)1.0e-3*(T)22.0e0;			// 22 Liters
		valset = true;
	}
	if(str == "cubit" && valset == false){
        convfactor = (T)18.0e0*metersperinch;	// 18 inches
		valset = true;
	}
	if(str == "ephah" && valset == false){
        convfactor = (T)1.0e-3*(T)22.0e0;			// 22 Liters
		valset = true;
	}
	if(str == "gerah" && valset == false){
        convfactor = (T)1.0e3*newtonsperlbf*(T)16.0e0*(T)2.0e-2;	// 0.02 oz per gerah
		valset = true;
	}
	if(str == "hin" && valset == false){
        convfactor = (T)1.0e-3*(T)22.0e0/(T)6.0e0;	// 1/6 of 22 Liters
		valset = true;
	}
	if(str == "homer" && valset == false){
        convfactor = (T)1.0e-3*(T)22.0e0*(T)10.0e0;	// 10x 22 Liters
		valset = true;
	}
	if(str == "mina" && valset == false){
        convfactor = (T)1.0e3*newtonsperlbf*(T)1.26e0;	// 1.26 lbf per mina
		valset = true;
	}
	if(str == "omer" && valset == false){
        convfactor = (T)1.0e-3*(T)22.0e0*(T)0.1e0;	// 1/10 of 22 Liters
		valset = true;
	}
	if(str == "bpound" && valset == false){
        convfactor = (T)1.0e3*newtonsperlbf*(T)16.0e0*(T)40.0e2;	// 40.0 oz per bpound
		valset = true;
	}
	if(str == "sdj" && valset == false){
        convfactor = (T)2.0e3*(T)18.0e0*metersperinch;	// 2000x 18-inch cubits
		valset = true;
	}
	if(str == "shek" && valset == false){
        convfactor = (T)1.0e3*newtonsperlbf*(T)16.0e0*(T)0.4e0;	// 0.4 oz per talent
		valset = true;
	}
	if(str == "span" && valset == false){
        convfactor = (T)9.0e0*metersperinch;	// 18-inch cubit
		valset = true;
	}
	if(str == "tal" && valset == false){
        convfactor = (T)1.0e3*newtonsperlbf*(T)75.0e0;// 75 lbf per talent
		valset = true;
	}

	/*
	 * RADIATION UNITS
	 */
	if(str == "R" && valset == false){
        convfactor = (T)2.58e-7;
		valset = true;
	}
	if(str == "rad" && valset == false){
        convfactor = (T)1.0e3*(T)1.0e-5;
		valset = true;
	}
	if(str == "rem" && valset == false){
        convfactor = (T)1.0e3*(T)1.0e-5;
		valset = true;
	}
	if(str == "Gy" && valset == false){
        convfactor = (T)1.0e3*(T)1.0e-3;
		valset = true;
	}
	if(str == "Sv" && valset == false){
        convfactor = (T)1.0e3*(T)1.0e-3;
		valset = true;
	}
	if(str == "Ci" && valset == false){
        convfactor = (T)3.7e10;
		valset = true;
	}
	if(str == "Bq" && valset == false){
        convfactor = (T)1.0e0;
		valset = true;
	}


	if(dir == 1){
		/*
		 * CONVERSION FROM INPUT UNIT TO FUNDAMENTAL UNITS
		 */
        if(unitObj.prefix > (T)1.0e-30){
			unitObj.conversionfactor = convfactor*unitObj.prefix;
		} else {
			unitObj.conversionfactor = convfactor;
		}
	}
	if(dir == 2){
		/*
		 * CONVERSION FROM FUNDAMENTAL UNITS TO INPUT UNIT
		 */
        if(unitObj.prefix > (T)1.0e-30){
            unitObj.conversionfactor = (T)1.0e0/convfactor/unitObj.prefix;
		} else {
            unitObj.conversionfactor = (T)1.0e0/convfactor;
		}
	}
}




/* ==================================================================
 * ================
 * ================    PUBLIC FUNCTIONS
 * ================
 */

template <class T>
UnitConvert<T>::UnitConvert()
{
	UnitConvert<T>::UnitConvertInitialize();
}


template <class T>
UnitConvert<T>::UnitConvert(UnitConvert<T> &a) : UnitConvert()
{
    UnitConvertSwap(*this, a);
}


#ifdef CXX11
template <class T>
UnitConvert<T>::UnitConvert(UnitConvert<T> &&a) : UnitConvert()
{
    UnitConvertSwap(*this, a);
}
#endif


template <class T>
UnitConvert<T>::~UnitConvert()
{
	// NOTHING TO DO
}


template <class T>
UnitConvert<T>& UnitConvert<T>::operator=(UnitConvert<T> a)
{
    UnitConvertSwap(*this, a);
    return *this;
}


#ifdef CXX11
template <class T>
UnitConvert<T>& UnitConvert<T>::operator=(UnitConvert<T> &&a)
{
    UnitConvertSwap(*this, a);
    return *this;
}
#endif


template <class T>
T UnitConvert<T>::ConvertUnits(T val, std::string unit1, std::string unit2)
{
	/*
	 * DEFINE A VALID RETURN VALUE.  IF UNIT-CONVERSION CAN BE PERFORMED THIS
	 * WILL BE UPDATED ACCORDINGLY.  IF UNIT-CONVERSION IS NOT PHYSICALLY
	 * POSSIBLE (E.G., "m" TO "kg") THIS VALUE OF 0.0 WILL BE RETURNED.  FOR A
	 * NON-ZERO INPUT VALUE, A NON-ZERO OUTPUT IS EXPECTED FOR VALID CONVERSIONS.
	 * ZERO-INPUT IS TRIVIAL SINCE ZERO IS ZERO REGARDLESS OF UNITS (E.G.,
	 * 0 INCHES EQUALS 0 METERS).
	 */
	T retval = 0.0;

	// GET THE NUMBER OF UNITS IN THE FIRST UNIT STRING (ONE MORE THAN THE
	// NUMBER OF PIPES PRESENT)
	int nunits1 = UnitConvert<T>::GetNumOccurences(unit1,"|");
	nunits1++;

	// GET THE NUMBER OF UNITS IN THE SECOND UNIT STRING (ONE MORE THAN THE
	// NUMBER OF PIPES PRESENT)
	int nunits2 = UnitConvert<T>::GetNumOccurences(unit2,"|");
	nunits2++;

	// CREATE Unit OBJECTS FOR EACH INPUT AND OUTPUT UNIT
	Unit<T> *inputUnits = new Unit<T>[nunits1];
	Unit<T> *outputUnits = new Unit<T>[nunits2];

	// DETERMINE THE POWER, PREFIX, AND LABEL FOR EACH UNIT
	std::string fulllabel,word;
	std::stringstream stream,stream2,stream3;
	T powerval;
	int count = 0;
	int item = 0;

	/*
	 * SPLIT INPUT LABEL INTO EACH UNIT AND SET THE APPROPRIATE OBJECT PARAMETER.
	 * NOTE THAT THE clear() METHOD FOR THE STRINGSTREAMS IS REQUIRED TO CLEAR
	 * THE SETTING OF eofbits GENERATED BY getline().  THE FOLLOWING LOOPS DO NOT
	 * BOUNDS-CHECK THE inputUnits AND outputUnits OBJECTS.  IT IS ASSUMED THAT
	 * BASING THE while() LOOPS ON THE SUB-STRINGS WILL BE ACCURATE SINCE THE
	 * NUMBER OF Unit OBJECTS IS DETERMINED BASED ON THE NUMBER OF SUB-STRINGS.
	 */
	stream.clear();
	stream.str(unit1);
	while(getline(stream,fulllabel,'|')){
		stream2.clear();
		stream2.str(fulllabel);
		count = 0;
		word = "";
		while(getline(stream2,word,':')){
			if(count == 0){
				inputUnits[item].prefix = UnitConvert<T>::GetPrefixValue(word);
			}
			if(count == 1){
				inputUnits[item].label = word;
			}
			if(count == 2){
				stream3.clear();
				stream3.str(word); stream3 >> powerval; stream3.str("");
				inputUnits[item].power = powerval;
			}
			count++;
		}
		item++;
	}

	item = 0;
	count = 0;
	stream.clear();
	stream.str(unit2);
	while(getline(stream,fulllabel,'|')){
		stream2.clear();
		stream2.str(fulllabel);
		count = 0;
		word = "";
		while(getline(stream2,word,':')){
			if(count == 0){
				outputUnits[item].prefix = UnitConvert<T>::GetPrefixValue(word);
			}
			if(count == 1){
				outputUnits[item].label = word;
			}
			if(count == 2){
				stream3.clear();
				stream3.str(""); stream3 << word; stream3 >> powerval; stream3.str("");
				outputUnits[item].power = powerval;
			}
			count++;
		}
		item++;
	}

	/*
	 * BEFORE ATTEMPTING ANY CONVERSIONS, FIRST DETERMINE THE POWERS ON THE
	 * FUNDAMENTAL UNITS FOR BOTH THE USER-SUPPLIED INPUT AND OUTPUT UNITS.
	 * THEN, CHECK THAT THE USER-SUPPLIED UNITS ARE CONSISTENT.  IF THEY ARE,
	 * PROCEED TO CALCULATE THE REQUIRED CONVERSION.  IF THEY ARE NOT, SKIP THE
	 * CONVERSION AND RETURN A VALUE OF 0.0.
	 */
	bool unitscompatible = false;
	std::string unitstr;
	T unitpower;
	Array1D<T> unitstmp;
	unitstmp.ResetSize(8,0.0);

	// DETERMINE FUNDAMENTAL POWERS OF INPUT UNITS
	for(int i=0; i<nunits1; i++){
		unitstr = inputUnits[i].label;
		unitpower = inputUnits[i].power;
		unitstmp.ResetVal(0.0);

		UnitConvert<T>::GetFundamentalType(unitstr,unitstmp);
		UnitConvert<T>::SetConversionFactor(1,inputUnits[i]);
		for(int j=0; j<8; j++){
			if(unitpower != 0){
				unitsin(j) += unitpower*unitstmp(j);
			} else {
				unitsin(j) += unitstmp(j);
			}
		}
		/*
		printf("Input Units: ");
		for(int i=0; i<8; i++){
			cout << unitsin(i) << " ";
		}
		cout << endl;
		*/
	}

	// DETERMINE FUNDAMENTAL POWERS OF OUTPUT UNITS
	for(int i=0; i<nunits2; i++){
		unitstr = outputUnits[i].label;
		unitpower = outputUnits[i].power;
		unitstmp.ResetVal(0.0);

		UnitConvert<T>::GetFundamentalType(unitstr,unitstmp);
		/*
		cout << "  Temp Units: ";
		for(int k=0; k<8; k++){
			cout << unitstmp(k) << " ";
		}
		cout << endl;
		*/
		UnitConvert<T>::SetConversionFactor(2,outputUnits[i]);
		for(int j=0; j<8; j++){
			if(unitpower != 0){
				unitsout(j) += unitpower*unitstmp(j);
			} else {
				unitsout(j) += unitstmp(j);
			}
		}
		/*
		cout << "Output Units: ";
		for(int k=0; k<8; k++){
			cout << unitsout(k) << " ";
		}
		cout << endl;
		*/
	}

	// CHECK TO SEE IF INPUT AND OUTPUT ARE COMPATIBLE
	bool tmpbool = true;
	bool forcemass = false;
	bool massforce = false;
	for(int i=0; i<8; i++){
		if(unitsin(i) != unitsout(i)){
			// SET TO FALSE IF *ANY* COMPONENT DIFFERS
			tmpbool = false;
		}
	}

	// CHECK FOR FORCE <--> MASS CONVERSION
    T eps = (T)0.1e0;
	/*
	cout << "Input Units: ";
	for(int i=0; i<8; i++){
		cout << unitsin(i) << " ";
	}
	cout << endl;
	cout << "Output Units: ";
	for(int i=0; i<8; i++){
		cout << unitsout(i) << " ";
	}
	cout << endl;
	cout << endl;
	cout << "Input Unit Test Results (mass): ";
	cout << abs(unitsin(0) - 1.0e0) << " " << abs(unitsin(1) - 0.0e0) << " " <<
			abs(unitsin(2) - 0.0e0) << " " << abs(unitsin(3) - 0.0e0) << " " <<
			abs(unitsin(4) - 0.0e0) << " " << abs(unitsin(5) - 0.0e0) << " " <<
			abs(unitsin(6) - 0.0e0) << " " << abs(unitsin(7) - 0.0e0) << endl;
	cout << "Input Tolerance Comparison (mass): ";
	if(abs(unitsin(0) - 1.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsin(1) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsin(2) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsin(3) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsin(4) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsin(5) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsin(6) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsin(7) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	cout << endl;
	cout << "Output Unit Test Results (force): ";
	cout << abs(unitsout(0) - 1.0e0) << " " << abs(unitsout(1) - 1.0e0) << " " <<
			abs(unitsout(2) + 2.0e0) << " " << abs(unitsout(3) - 0.0e0) << " " <<
			abs(unitsout(4) - 0.0e0) << " " << abs(unitsout(5) - 0.0e0) << " " <<
			abs(unitsout(6) - 0.0e0) << " " << abs(unitsout(7) - 0.0e0) << endl;
	cout << "Output Tolerance Comparison (mass): ";
	if(abs(unitsout(0) - 1.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsout(1) - 1.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsout(2) + 2.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsout(3) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsout(4) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsout(5) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsout(6) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	if(abs(unitsout(7) - 0.0e0) < eps){ cout << "True ";} else { cout << "False "; }
	cout << endl;
	cout << "Tolerance: " << eps << endl;
	*/
	if(abs(unitsin(0) - 1.0e0) < eps && abs(unitsin(1) - 0.0e0) < eps &&
		abs(unitsin(2) - 0.0e0) < eps && abs(unitsin(3) - 0.0e0) < eps &&
		abs(unitsin(4) - 0.0e0) < eps && abs(unitsin(5) - 0.0e0) < eps &&
		abs(unitsin(6) - 0.0e0) < eps && abs(unitsin(7) - 0.0e0) < eps &&
		abs(unitsout(0) - 1.0e0) < eps && abs(unitsout(1) - 1.0e0) < eps &&
		abs(unitsout(2) + 2.0e0) < eps && abs(unitsout(3) - 0.0e0) < eps &&
		abs(unitsout(4) - 0.0e0) < eps && abs(unitsout(5) - 0.0e0) < eps &&
		abs(unitsout(6) - 0.0e0) < eps && abs(unitsout(7) - 0.0e0) < eps){
		massforce = true;
		//cout << "Converting mass to force" << endl;
	}
	if(abs(unitsout(0) - 1.0e0) < eps && abs(unitsout(1) - 0.0e0) < eps &&
		abs(unitsout(2) - 0.0e0) < eps && abs(unitsout(3) - 0.0e0) < eps &&
		abs(unitsout(4) - 0.0e0) < eps && abs(unitsout(5) - 0.0e0) < eps &&
		abs(unitsout(6) - 0.0e0) < eps && abs(unitsout(7) - 0.0e0) < eps &&
		abs(unitsin(0) - 1.0e0) < eps && abs(unitsin(1) - 1.0e0) < eps &&
		abs(unitsin(2) + 2.0e0) < eps && abs(unitsin(3) - 0.0e0) < eps &&
		abs(unitsin(4) - 0.0e0) < eps && abs(unitsin(5) - 0.0e0) < eps &&
		abs(unitsin(6) - 0.0e0) < eps && abs(unitsin(7) - 0.0e0) < eps){
		forcemass = true;
		//cout << "Converting force to mass" << endl;
	}


	if(tmpbool == true){
		// RESET unitscompatible ONLY IF tmpbool REMAINS TRUE AFTER CHECKING
		// EACH COMPONENT OF THE ARRAYS
		unitscompatible = true;
	}


	/*
	 * IF UNITS ARE COMPATIBLE, DETERMINE THE CONVERSION FACTORS.  FOR EACH UNIT
	 * IN THE INPUT AND OUTPUT UNITS, DETERMINE A CONVERSION FACTOR.  THIS
	 * CONVERSION FACTOR IS DEFINED SUCH THAT THE PRODUCT OF CONVERSION FACTORS
	 * FOR ALL UNITS IN A GROUP MULTIPLIED BY THE INPUT VALUE WILL RESULT IN THE
	 * VALUE'S BEING EXPRESSED IN THE FUNDAMENTAL UNITS OF THIS CLASS.
	 *
	 * SIMILARLY, MULTIPLYING THE INTERNAL VALUE (WHICH IS IN TERMS OF THE
	 * FUNDAMENTAL UNITS) BY THE PRODUCT OF THE CONVERSION FACTORS OF THE OUTPUT
	 * UNITS WILL PRODUCE A NUMERICAL RESULT EXPRESSED IN THE USER-SPECIFIED
	 * OUTPUT UNITS.
	 */
	T convfactor = 1.0;
	if(unitscompatible == true || forcemass == true || massforce == true){
		/*
		 * NOW, SET THE internalvalue VARIABLE BASED ON THE USER-SUPPLIED INPUT
		 * VALUE AND INPUT UNITS.
		 */
		for(int i=0; i<nunits1; i++){
			if(abs(inputUnits[i].power) > 1.0e-5){
				convfactor *= pow(inputUnits[i].conversionfactor,inputUnits[i].power);
			} else {
				convfactor *= inputUnits[i].conversionfactor;
			}
		}

		internalvalue = val*convfactor;

		// HANDLE SPECIAL CASE OF TEMPERATURE CONVERSIONS
		if(nunits1 == 1){
			if(inputUnits[0].label == "K"){ internalvalue = val; }
            if(inputUnits[0].label == "F"){ internalvalue = (val-(T)32.0e0)/(T)1.8e0 + (T)273.15e0; }
            if(inputUnits[0].label == "C"){ internalvalue = val + (T)273.15e0; }
		}

		/*
		 * HANDLE SPECIAL CASE OF MASS <--> FORCE CONVERSIONS.  IF PERFORMING SUCH
		 * A CONVERSION IS TO BE DONE, EITHER MULTIPLY OR DIVIDE BY STANDARD
		 * GRAVITY (THUS VALUES ONLY VALID FOR ASSUMED VALUE OF GRAVITY)
		 */
		if(massforce == true){
			// INPUT: MASS, OUTPUT: FORCE
			internalvalue *= geeinms2;
		}
		if(forcemass == true){
			// INPUT: FORCE, OUTPUT: MASS
			internalvalue /= geeinms2;
		}


		/*
		 * NEXT, DETERMINE THE MODIFICATIONS REQUIRED TO EXPRESS internalvalue
		 * IN TERMS OF THE USER-SUPPLIED OUTPUT UNITS.
		 */
		convfactor = 1.0;
		for(int i=0; i<nunits2; i++){
			if(abs(outputUnits[i].power) > 1.0e-5){
				convfactor *= pow(outputUnits[i].conversionfactor,outputUnits[i].power);
			} else {
				convfactor *= outputUnits[i].conversionfactor;
			}
		}

		retval = internalvalue*convfactor;

		// HANDLE SPECIAL CASE OF TEMPERATURE CONVERSIONS
		if(nunits1 == 1){
			if(outputUnits[0].label == "K"){ retval = internalvalue; }
            if(outputUnits[0].label == "F"){ retval = (internalvalue-(T)273.15e0)*(T)1.8e0 + (T)32.0e0; }
            if(outputUnits[0].label == "C"){ retval = internalvalue - (T)273.15e0; }
		}
	} else {
		/*
		 * NOTHING TO DO SINCE retval IS DEFINED TO BE 0.0 AT THE TOP OF THIS
		 * ROUTINE.
		 */
	}


	/*
	 * CLEAN-UP Unit OBJECTS PRIOR TO EXITING ROUTINE
	 */
	/*
	cout << "Input Units: ";
	for(int i=0; i<8; i++){
		cout << unitsin(i) << " ";
	}
	cout << endl;
	cout << "Output Units: ";
	for(int i=0; i<8; i++){
		cout << unitsout(i) << " ";
	}
	cout << endl;
	*/

	unitsin.ResetVal(0.0);
	unitsout.ResetVal(0.0);

	delete [] inputUnits;
	delete [] outputUnits;

	return retval;
}


template <class T>
std::string UnitConvert<T>::PrintUnits() const
{
	std::string retval;
	retval = "Error";

	/*
	 * GENERATE STRING CONTAINING VALID UNITS AND PREFIXES
	 */

	retval =  "UNIT ENTRY FORMAT: \n";
	retval += "    [prefix]:unitlabel:power \n";
	retval += "  \n";
	retval += "  Where: \n";
	retval += "    [prefix] is an optional standard-SI prefix \n";
	retval += "    unitlabel is the base unit \n";
	retval += "    power is a decimal power associated with the unit \n";
	retval += "  \n";
	retval += "  Note that both colons (:) are required, although the prefix \n";
	retval += "  can be omitted.  If no prefix is required, a '0' can be supplied \n";
	retval += "  as a placeholder. \n";
	retval += "  \n";
	retval += "  Also, multiple units can be combined using the pipe (|) character. \n";
	retval += "  Such combinations imply multiplication of units, thus division \n";
	retval += "  requires negative powers. \n";
	retval += "  \n";
	retval += "  All prefixes and labels are case-sensitive \n";
	retval += "  \n";
	retval += "  Notice that some units use an underscore (_).  This is required and \n";
	retval += "  cannot be replaced with a space. \n";
	retval += "  \n";
	retval += "  Examples: \n";
	retval += "    :ft:1 - feet \n";
	retval += "    0.ft:1 - feet \n";
	retval += "    :ft:2 - feet-squared \n";
	retval += "    :lbf:1|:ft:-2 - pounds-force per square foot \n";
	retval += "  \n";
	retval += "  \n";
	retval += "VALID UNIT LABELS: (denoted by quantity within single quotes) \n";
	retval += "  Metric Units \n";
	retval += "    - 'A':       Amperes \n";
	retval += "    - 'AU':      astronomical-units \n";
	retval += "    - 'Bq':      becquerels \n";
	retval += "    - 'C':       degrees Celcius \n";
	retval += "    - 'cal':     small (gram) calorie \n";
	retval += "    - 'Cal':     large (dietary) calorie \n";
	retval += "    - 'Ci':      Curies \n";
	retval += "    - 'dyn':     dynes \n";
	retval += "    - 'erg':     ergs \n";
	retval += "    - 'eV':      electron-volts \n";
	retval += "    - 'g':       grams \n";
	retval += "    - 'Gy':      Grays \n";
	retval += "    - 'ha':      hectares \n";
	retval += "    - 'J':       Joules \n";
	retval += "    - 'K':       Kelvin \n";
	retval += "    - 'L':       liters \n";
	retval += "    - 'ly':      light-years \n";
	retval += "    - 'm':       meters \n";
	retval += "    - 'mol':     mols \n";
	retval += "    - 'N':       Newtons \n";
	retval += "    - 'N_m':     Newton-meters \n";
	retval += "    - 'Pa':      Pascals \n";
	retval += "    - 'pc':      parsecs \n";
	retval += "    - 'rad':     radiation dose \n";
	retval += "    - 'radian':  radians \n";
	retval += "    - 'rem':     rad-equivalent-man (rad <--> rem assumes Q = 1) \n";
	retval += "    - 'R':       Roentgens \n";
	retval += "    - 'Sv':      Sieverts \n";
	retval += "    - 'W':       Watts \n";
	retval += "  \n";
	retval += "  English Units \n";
	retval += "    - 'acre':    acres \n";
	retval += "    - 'atm':     atmospheres \n";
	retval += "    - 'bbl':     oil barrels \n";
	retval += "    - 'BTU':     British Thermal Units \n";
	retval += "    - 'bu':      bushels \n";
	retval += "    - 'cb':      cables \n";
	retval += "    - 'chain':   chains \n";
	retval += "    - 'cup':     US cups (8 fluid ounces ea.) \n";
	retval += "    - 'cwt':     US hundredweight \n";
	retval += "    - 'dbbl':    dry barrels \n";
	retval += "    - 'dpt':     US dry-pint \n";
	retval += "    - 'dqt':     US dry-quart \n";
	retval += "    - 'dr':      drams \n";
	retval += "    - 'dwt':     pennyweight \n";
	retval += "    - 'F':       degrees Fahrenheit \n";
	retval += "    - 'fl_dr':   US fluid-drams \n";
	retval += "    - 'fl_oz':   US fluid-ounces \n";
	retval += "    - 'ft':      feet \n";
	retval += "    - 'ft_lbf':  foot-pounds-force \n";
	retval += "    - 'ftm':     fathoms \n";
	retval += "    - 'fur':     furlongs \n";
	retval += "    - 'gal':     US gallons \n";
	retval += "    - 'gee':     g's (1 gee = acceleration of gravity at Earth's surface) \n";
	retval += "    - 'gi':      gills \n";
	retval += "    - 'gr':      grains \n";
	retval += "    - 'hand':    hands \n";
	retval += "    - 'hghd':    hogsheads \n";
	retval += "    - 'hp':      horse-power \n";
	retval += "    - 'imp_gal': imperial gallons \n";
	retval += "    - 'in':      inches \n";
	retval += "    - 'jig':     jiggers \n";
	retval += "    - 'ksi':     1000 pounds-force-per-square-inch \n";
	retval += "    - 'lbbl':    liquid barrels \n";
	retval += "    - 'lbf':     pounds-force \n";
	retval += "    - 'lbm':     pounds-mass \n";
	retval += "    - 'lea':     leagues \n";
	retval += "    - 'li':      links \n";
	retval += "    - 'lpt':     US liquid pints \n";
	retval += "    - 'lqt':     US liquid quarts \n";
	retval += "    - 'mil;      1/1000th of an inch \n";
	retval += "    - 'mile':    miles \n";
	retval += "    - 'minim':   minims \n";
	retval += "    - 'nmi':     nautical miles \n";
	retval += "    - 'oz':      ounces \n";
	retval += "    - 'p':       points \n";
	retval += "    - 'P':       picas \n";
	retval += "    - 'pk':      pecks \n";
	retval += "    - 'psi':     pounds-force-per-square-inch \n";
	retval += "    - 'rod':     rods \n";
	retval += "    - 'slug':    slug \n";
	retval += "    - 'Tbsp':    tablespoons \n";
	retval += "    - 'tsp':     teaspoon \n";
	retval += "  \n";
	retval += "  Biblical Units \n";
	retval += "    - 'bath':    baths \n";
	retval += "    - 'cubit':   cubits (definitions vary.  18 inches assumed here) \n";
	retval += "    - 'ephah':   ephahs \n";
	retval += "    - 'gerah':   gerahs \n";
	retval += "    - 'hin':     hin \n";
	retval += "    - 'homer':   homers \n";
	retval += "    - 'mina':    minas \n";
	retval += "    - 'omer':    omer \n";
	retval += "    - 'bpound':  pounds (Biblical) \n";
	retval += "    - 'sdj':     Sabbath day's journey \n";
	retval += "    - 'shek':    sheckels \n";
	retval += "    - 'span':    span \n";
	retval += "    - 'tal':     talents \n";
	retval += "  \n";
	retval += "  Time \n";
	retval += "    - 'sec':     seconds \n";
	retval += "    - 'min':     minutes \n";
	retval += "    - 'hr':      hours \n";
	retval += "    - 'day':     days \n";
	retval += "  \n";
	retval += "  \n";
	retval += "VALID SI PREFIXES: (can be applied to any unit) \n";
	retval += "    - 'y':  yocto (10^-24) \n";
	retval += "    - 'z':  zepto (10^-21) \n";
	retval += "    - 'a':  atto  (10^-18) \n";
	retval += "    - 'f':  femto (10^-15) \n";
	retval += "    - 'p':  pico  (10^-12) \n";
	retval += "    - 'n':  nano  (10^-09) \n";
	retval += "    - 'u':  micro (10^-06) (this is the english 'u', not a greek 'mu')\n";
	retval += "    - 'm':  milli (10^-03) \n";
	retval += "    - 'c': centi (10^-02) \n";
	retval += "    - 'd':  deci  (10^-01) \n";
	retval += "    - '0':  None  (10^+00) (this is notation here, not a standard prefix)\n";
	retval += "    - 'da': deca  (10^+01) \n";
	retval += "    - 'h':  hecto (10^+02) \n";
	retval += "    - 'k':  kilo  (10^+03) \n";
	retval += "    - 'M':  mega  (10^+06) \n";
	retval += "    - 'G':  giga  (10^+09) \n";
	retval += "    - 'T':  tera  (10^+12) \n";
	retval += "    - 'P':  peta  (10^+15) \n";
	retval += "    - 'E':  exa   (10^+18) \n";
	retval += "    - 'Z':  zetta (10^+21) \n";
	retval += "    - 'Y':  yotta (10^+24) \n";


	return retval;
}

