/**
 * @file StringManip.cxx
 * @author Robert Grandin
 * @brief Implementation of non-templated functions in StringManip namespace.
 */


#include <StringManip.h>


std::string StringManip::BoolToString(bool b)
{
    std::string retval("ERR");
    if(b == true){retval = "true";}
    if(b == false){retval = "false";}
    return retval;
}


int StringManip::BoolToInt(bool b)
{
    int retval = -1;
    if(b == true){retval = 1;}
    if(b == false){retval = 0;}
    return retval;
}


bool StringManip::IntToBool(int i)
{
    bool retval = false;
    if(i == 1){retval = true;}
    if(i == 0){retval = false;}
    return retval;
}


int StringManip::StrToUInt(std::string str)
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


std::string StringManip::FormatTime(double time)
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


std::string StringManip::DetermFileExt(std::string file)
{
    size_t index1 = file.find_last_of(".");
    std::string filenameext;
    filenameext = file.substr(index1+1);

    return filenameext;
}


std::string StringManip::DetermMemReq(double npts, size_t size)
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


std::string StringManip::DetermMemReq(size_t npts, size_t size)
{
    return StringManip::DetermMemReq((double)npts,size);
}


std::string StringManip::DetermBytesLabel(double nbytes)
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



std::string StringManip::NumToStr(const int val, const int mindigits)
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


void StringManip::DetermFileStem(const std::string &filename, std::string &stem, size_t &ndigits,
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

        if(index_us != filename.npos){
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


void StringManip::uppercaseString(std::string &s)
{
    for(int i=0; s[i]!='\0'; i++){
        s[i] = toupper(s[i]);
    }
}


std::string& StringManip::ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}


std::string& StringManip::rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}


std::string& StringManip::trim(std::string &s) {
    return ltrim(rtrim(s));
}


int StringManip::DetermNumElements(std::stringstream &stream)
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


int StringManip::DetermNumElements(std::stringstream &stream, const char delim)
{
    /* Save current state of stream into string. */
    std::string originalstream;
    originalstream = stream.str();

    /* Count elements in stream. */
    std::string word;
    int count = 0;

    while(getline(stream,word,delim)){
        count++;
    }

    /* Set stream to original value using saved string. */
    stream.str(originalstream);

    /* Return number of elements in stream. */
    return count;
}


int StringManip::str_compare(const char *a, const char *b)
{
    int retval = 0;

#ifdef COMPILELINUX
    retval = strcasecmp(a, b);
#endif
#ifdef WIN_MSVC
    retval = _strcmpi(a, b);
#endif

    return retval;
}


std::string StringManip::SanitizeString(std::string &input)
{
    std::string retval(input);

    for(int i=0; retval[i]!='\0'; i++){
        if(retval[i] == '\\'){
            retval[i] = '/';
        }
    }

    //retval = "\"" + retval + "\"";

    return retval;
}
