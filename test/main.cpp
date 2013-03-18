#include <ArrayBase.h>
#include <Array2D.h>
#include <Array3D.h>
#include <DataFilters.h>
#include <UniformVolume.h>

#include <omp.h>


/* Define symbols to control implementation of code. */
#define COPY_CTOR 1     /* Cases:
                         *  1   CORRECT METHOD
                         *          Base class:     - Member variables set in initialization list.  This
                         *                            makes use of the copy-constructor which is guaranteed
                         *                            to exist for all variables, and does not require the
                         *                            programmer to manually copy values.
                         *          Derived class:  - Initialization list calls base copy-constructor and
                         *                            then sets derived member variables.  Prevents need to
                         *                            manually copy data.
                         *
                         *          All objects in C++ have copy constructors.  If no copy constructor is
                         *          explicitly defined by the programmer, the compiler creates one and performs
                         *          and element-copy of member variables.
                         *
                         *
                         *  2   Correct, but less-efficient, method
                         *          Base class:     - Member values set manually within the function.
                         *          Derived class:  - Member values set manually within the function.
                         *
                         *          This method introduces opportunity for error since the user must manually
                         *          copy data to initialize member variables.  Using the variables' copy
                         *          constructors, as in method 1, removes this additional potential source
                         *          of error and makes use of already-existing code.
                         *
                         *  3   PARTIALLY-INCORRECT METHOD
                         *          Base class:     - Initialization list used.  No explicit call to
                         *                            the default constructor (this is OK and matches case 1).
                         *          Derived class:  - Initialization list used.  No explicit call to
                         *                            initialize the base class (this is a problem).
                         *
                         *  4   INCORRECT METHOD
                         *          Base class:     - Initialization list not used.  No explicit call to
                         *                            the default constructor.  This is problem because the
                         *                            member variables are never initialized until they are
                         *                            manually managed within the function.
                         *          Derived class:  - Initialization list not used.  No explicit call to
                         *                            initialize the base class (this is a problem).
                         *
                         *  Cases 3 and 4 may compile, depending on the compiler, but due to the lack of explicit
                         *  initialization (either of the member variables and/or the base class), unexpected errors may occur.
                         */

#define MOVE_CTOR 1     /* Cases:
                         *  1   CORRECT METHOD      - Derived move constructor calls the base move constructor from the
                         *                            initialization list.  From the code, note the use of "std::move(obj)"
                         *                            to force the calling of the base move constructor rather than the base
                         *                            copy constructor.
                         *
                         *  2   INCORRECT METHOD    - Derived move constructor calls base copy constructor rather than the
                         *                            base move constructor.  The example in this code will still compile and
                         *                            run, but the output will show that the move constructor is not properly
                         *                            implemented.
                         */





/**
 * @brief Basic class demonstrating the functionality of a stand-alone class as well as
 *  a base class.
 */
class TestBase {

public:
    /** @brief Member variable. */
    std::string name;

    /**
     * @brief Default constructor.
     */
    TestBase() : name("defaut_name") {
        std::cout << "TestBase default constructor" << std::endl;
        std::cout << "  - Name: " << name << std::endl;
    }

    /**
     * @brief Constructor which sets member variable to user-specified value.
     * @param username Value to be assigned to member variable.
     */
    TestBase(std::string username) : name(username) {
        std::cout << "TestBase name-specified constructor" << std::endl;
        std::cout << "  - Name: " << name << std::endl;
    }

#if COPY_CTOR == 1 || COPY_CTOR == 3
    /**
     * @brief Copy constructor using initialization list to set member variables.
     * @param tb Instance of TestBase class to be copied.
     */
    TestBase(const TestBase &tb) : name(tb.name) {
        std::cout << "TestBase copy constructor" << std::endl;
        std::cout << "  - Name (this):  " << this->name << std::endl;
        std::cout << "  - Name (input): " << tb.name << std::endl;
    }
#endif
#if COPY_CTOR == 2
    /**
     * @brief Copy constructor using default constructor to initialize object, and then
     *  manually copy values to set member variables.
     * @param tb Instance of TestBase class to be copied.
     */
    TestBase(const TestBase &tb) : TestBase() {
        std::cout << "TestBase copy constructor" << std::endl;
        std::cout << "  - Name (this):  " << this->name << std::endl;
        std::cout << "  - Name (input): " << tb.name << std::endl;

        std::cout << "    Copying name" << std::endl;
        name = tb.name;

        std::cout << "  - Name (this):  " << this->name << std::endl;
    }
#endif
#if COPY_CTOR == 4
    /**
     * @brief Copy constructor which performs no pre-function initialization of member variables.
     * @param tb Instance of TestBase class to be copied.
     */
    TestBase(const TestBase &tb) {
        std::cout << "TestBase copy constructor" << std::endl;
        std::cout << "  - Name (this):  " << this->name << std::endl;
        std::cout << "  - Name (input): " << tb.name << std::endl;

        std::cout << "    Copying name" << std::endl;
        name = tb.name;

        std::cout << "  - Name (this):  " << this->name << std::endl;
    }
#endif

    /**
     * @brief Move constructor.
     * @param tb Instance of TestBase class to be moved.
     */
    TestBase(TestBase &&tb) {
        /* Use default constructor to create a default instance of this object. Then,
         * use std::swap to replace the default instance values with those of the
         * passed-in instance. */

        std::cout << "TestBase move constructor" << std::endl;
        std::cout << "  - Name (this):  " << this->name << std::endl;
        std::cout << "  - Name (input): " << tb.name << std::endl;

        std::cout << "    Moving name" << std::endl;
        std::swap(name, tb.name);

        std::cout << "  - Name (this):  " << this->name << std::endl;
    }

    /**
     * @brief Destructor.
     */
    virtual ~TestBase(){
        std::cout << "TestBase destructor" << std::endl;
        std::cout << "  - Name: " << name << std::endl;
    }


    /**
     * @brief Copy-assignment operator.
     *
     *  This is written using the copy-and-swap idiom, making it also function as the
     *  move-assignment operator.
     * @param tb Instance of TestBase class to be assigned.
     * @return Instance of TestBase class.
     */
    TestBase& operator=(TestBase tb){
        /* Copy-assignment operator, using the copy-and-swap idiom.  Note that
         * argument is passed by-value, not by-reference.  This makes use of the
         * copy constructor.  Using the copy-and-swap idiom means that this also
         * works as the move constructor. */
        std::cout << "TestBase copy-assignment (pre)" << std::endl;
        std::cout << "  - Name:         " << name << std::endl;

        std::cout << "TestBase copy-assignment (rvalue)" << std::endl;
        std::cout << "  - Name:         " << tb.name << std::endl;

        /* Swap member variables.  This function can be replaced with a friend
         * function which performs the swap for all member variables. */
        std::swap(this->name, tb.name);

        std::cout << "TestBase copy-assignment (post)" << std::endl;
        std::cout << "  - Name:         " << name << std::endl;

        /* Return reference to this object. */
        return *this;
    }

};



/**
 * @brief Basic class demonstrating derivation.
 */
class TestDerived : public TestBase {

public:
    /** @brief Member variable used to track changes. */
    std::string derived_name;


    /**
     * @brief Default constructor.
     */
    TestDerived() : derived_name("defaut_derived_name") {
        std::cout << "TestDerived default constructor" << std::endl;
        std::cout << "  - Name:         " << name << std::endl;
        std::cout << "  - Derived Name: " << derived_name << std::endl;
    }


    /**
     * @brief Constructor which sets member variable to user-specified value.
     * @param username Value to be assigned to member variable.
     */
    TestDerived(std::string username) : derived_name(username) {
        std::cout << "TestDerived name-specified constructor" << std::endl;
        std::cout << "  - Name:         " << name << std::endl;
        std::cout << "  - Derived Name: " << derived_name << std::endl;
    }

#if COPY_CTOR == 1
    /**
     * @brief Copy constructor using copy constructor of base class and initialization list to
     *  initialize base class and set member variables.
     * @param td Instance of TestDerived class to be copied.
     */
    TestDerived(const TestDerived &td) : TestBase(td), derived_name(td.derived_name) {
        std::cout << "TestDerived copy constructor" << std::endl;
        std::cout << "  - Name (this):          " << name << std::endl;
        std::cout << "  - Derived Name (this):  " << derived_name << std::endl;
        std::cout << "  - Name (input):         " << td.name << std::endl;
        std::cout << "  - Derived Name (input): " << td.derived_name << std::endl;
    }
#endif
#if COPY_CTOR == 2
    /**
     * @brief Copy constructor using copy constructor of base class to initialize object, and then
     *  manually copy values to set member variables.
     * @param tb Instance of TestDerived class to be copied.
     */
    TestDerived(const TestDerived &td) : TestBase(td) {
        std::cout << "TestDerived copy constructor" << std::endl;
        std::cout << "  - Name (this):          " << name << std::endl;
        std::cout << "  - Derived Name (this):  " << derived_name << std::endl;
        std::cout << "  - Name (input):         " << name << std::endl;
        std::cout << "  - Derived Name (input): " << derived_name << std::endl;

        std::cout << "    Copying derived_name" << std::endl;
        derived_name = td.derived_name;

        std::cout << "  - Name (this):          " << name << std::endl;
        std::cout << "  - Derived Name (this):  " << derived_name << std::endl;
    }
#endif
#if COPY_CTOR == 3
    /**
     * @brief Copy constructor with no explicit initialization of base class and member variables
     *  set via initialization list.
     * @param td Instance of TestDerived class to be copied.
     */
    TestDerived(const TestDerived &td) : derived_name(td.derived_name) {
        std::cout << "TestDerived copy constructor" << std::endl;
        std::cout << "  - Name (this):          " << name << std::endl;
        std::cout << "  - Derived Name (this):  " << derived_name << std::endl;
        std::cout << "  - Name (input):         " << name << std::endl;
        std::cout << "  - Derived Name (input): " << derived_name << std::endl;
    }
#endif
#if COPY_CTOR == 4
    /**
     * @brief Copy constructor with no explicit initialization.
     * @param td Instance of TestDerived class to be copied.
     */
    TestDerived(const TestDerived &td) {
        std::cout << "TestDerived copy constructor" << std::endl;
        std::cout << "  - Name (this):          " << name << std::endl;
        std::cout << "  - Derived Name (this):  " << derived_name << std::endl;
        std::cout << "  - Name (input):         " << name << std::endl;
        std::cout << "  - Derived Name (input): " << derived_name << std::endl;

        std::cout << "    Copying derived_name" << std::endl;
        derived_name = td.derived_name;

        std::cout << "  - Name (this):          " << name << std::endl;
        std::cout << "  - Derived Name (this):  " << derived_name << std::endl;
    }
#endif

#if MOVE_CTOR == 1
    /**
     * @brief Move constructor which initializes the base class via the base class move constructor.
     * @param td Instance of TestDerived class to be moved.
     */
    TestDerived(TestDerived &&td) : TestBase(std::move(td)) {

        /* Note the use of "std::move(td)" in the initialization list.  This is what forces the
         * base class move constructor to be called.
         *
         * Also, the member variable "derived_name" for this instance will not be initialized.  However,
         * the swap operation will set it.  The swap will result in there being no value defined in 'td'
         * post-swap, but this is not a problem since 'td' is a temporary. */

        std::cout << "TestDerived move constructor" << std::endl;
        std::cout << "  - Name (this):          " << name << std::endl;
        std::cout << "  - Derived Name (this):  " << derived_name << std::endl;
        std::cout << "  - Name (input):         " << td.name << std::endl;
        std::cout << "  - Derived Name (input): " << td.derived_name << std::endl;

        std::cout << "    Moving derived_name" << std::endl;
        std::swap(derived_name, td.derived_name);

        std::cout << "  - Name (this):          " << name << std::endl;
        std::cout << "  - Derived Name (this):  " << derived_name << std::endl;
    }
#endif
#if MOVE_CTOR == 2
    /**
     * @brief Move constructor which initializes the base class via the base class copy constructor.
     * @param td Instance of TestDerived class to be moved.
     */
    TestDerived(TestDerived &&td) : TestBase(td) {

        /* The member variable "derived_name" for this instance will not be initialized.  However,
         * the swap operation will set it.  The swap will result in there being no value defined in 'td'
         * post-swap, but this is not a problem since 'td' is a temporary. */

        std::cout << "TestDerived move constructor" << std::endl;
        std::cout << "  - Name (this):          " << name << std::endl;
        std::cout << "  - Derived Name (this):  " << derived_name << std::endl;
        std::cout << "  - Name (input):         " << td.name << std::endl;
        std::cout << "  - Derived Name (input): " << td.derived_name << std::endl;

        std::cout << "    Moving derived_name" << std::endl;
        std::swap(derived_name, td.derived_name);

        std::cout << "  - Name (this):          " << name << std::endl;
        std::cout << "  - Derived Name (this):  " << derived_name << std::endl;
    }
#endif

    /**
     * @brief Destructor.
     */
    virtual ~TestDerived(){
        std::cout << "TestDerived destructor" << std::endl;
        std::cout << "  - Name:         " << name << std::endl;
        std::cout << "  - Derived Name: " << derived_name << std::endl;
    }


    /**
     * @brief Copy-assignment operator.
     *
     *  This is written using the copy-and-swap idiom, making it also function as the
     *  move-assignment operator.
     * @param td Instance of TestDerived class to be assigned.
     * @return Instance of TestDerived class.
     */
    TestDerived& operator=(TestDerived td){

        std::cout << "TestDerived copy-assignment (pre)" << std::endl;
        std::cout << "  - Name:         " << name << std::endl;
        std::cout << "  - Derived Name: " << derived_name << std::endl;

        std::cout << "TestDerived copy-assignment (rvalue)" << std::endl;
        std::cout << "  - Name:         " << td.name << std::endl;
        std::cout << "  - Derived Name: " << td.derived_name << std::endl;

        /* Call assignment operator for base class.  This handles the swap of base class member variables.
         * For more info: http://stackoverflow.com/a/8867477. */
        TestBase::operator=(static_cast<TestBase>(td));

        /* Swap derived class member variables.  This function can be replaced with a friend
         * function which performs the swap for all member variables. */
        std::swap(this->derived_name, td.derived_name);

        std::cout << "TestDerived copy-assignment (post)" << std::endl;
        std::cout << "  - Name:         " << name << std::endl;
        std::cout << "  - Derived Name: " << derived_name << std::endl;

        /* Return reference to this object. */
        return *this;
    }

};














int main()
{
    int testcase = 2;   /* 0: Illustrative sample classes defined above. */
                        /* 1: Run Test() functions for general-use C++ classes. */
                        /* 2: Sandbox area for general testing/debugging. */

    if(testcase == 0){

        std::cout << "======================================================" << std::endl;
        std::cout << "====" << std::endl;
        std::cout << "====    Testing Base Class" << std::endl;
        std::cout << "====" << std::endl;
        std::cout << "====" << std::endl;

        std::cout << std::endl << std::endl;
        std::cout << "Create default instance of TestBase" << std::endl;
        TestBase tb1;

        std::cout << std::endl << std::endl;
        std::cout << "Create alternative instance of TestBase" << std::endl;
        TestBase tb2("test2");

        std::cout << std::endl << std::endl;
        std::cout << "Create copy instance of TestBase" << std::endl;
        TestBase tb3(tb2);
        tb3.name = "test3";

        std::cout << std::endl << std::endl;
        std::cout << "Create move instance of TestBase" << std::endl;
        TestBase tb4 = std::move(tb3);   /* 'std::move()' used to force move-constructor. */


        /* Create two objects to test assignment operators. */
        std::cout << std::endl << std::endl;
        TestBase tb5("assign1");
        std::cout << std::endl;
        TestBase tb6("assign2");
        std::cout << std::endl;
        TestBase tb7("assign3");

        std::cout << std::endl << std::endl;
        std::cout << "Copy-assignment operator" << std::endl;
        tb5 = tb6;

        std::cout << std::endl << std::endl;
        std::cout << "Move-assignment operator" << std::endl;
        tb7 = std::move(tb6);       /* 'std::move()' used to force move-assignment. */



        std::cout << std::endl << std::endl;
        std::cout << std::endl << std::endl;
        std::cout << "======================================================" << std::endl;
        std::cout << "====" << std::endl;
        std::cout << "====    Testing Derived Class" << std::endl;
        std::cout << "====" << std::endl;
        std::cout << "====" << std::endl;

        std::cout << std::endl << std::endl;
        std::cout << "Create default instance of TestDerived" << std::endl;
        TestDerived td1;    td1.name = "td1";


        std::cout << std::endl << std::endl;
        std::cout << "Create alternative instance of TestDerived" << std::endl;
        TestDerived td2("test4");   td2.name = "td2";

        std::cout << std::endl << std::endl;
        std::cout << "Create copy instance of TestDerived" << std::endl;
        TestDerived td3(td2);       td3.name = "td3";
        td3.derived_name = "test5";

        std::cout << std::endl << std::endl;
        std::cout << "Create move instance of TestDerived" << std::endl;
        TestDerived td4 = std::move(td3);   /* 'std::move()' used to force move-constructor. */


        /* Create two objects to test assignment operators. */
        std::cout << std::endl << std::endl;
        TestDerived td5("assign4");     td5.name = "assign1";
        std::cout << std::endl;
        TestDerived td6("assign5");     td6.name = "assign2";
        std::cout << std::endl;
        TestDerived td7("assign6");     td7.name = "assign3";

        std::cout << std::endl << std::endl;
        std::cout << "Copy-assignment operator" << std::endl;
        td5 = td6;

        std::cout << std::endl << std::endl;
        std::cout << "Move-assignment operator" << std::endl;
        td7 = std::move(td6);       /* 'std::move()' used to force move-assignment. */



        std::cout << std::endl << std::endl;
        std::cout << "Cleanup" << std::endl;

    } /* if(testcase == 0) */









    if(testcase == 1){
        std::string result;

        ArrayBase<float> abf;
        abf.Test(result);
        std::cout << "ArrayBase result: " << result << std::endl;


        result = "";
        DataFilters<float> dff;
        dff.Test(result);
        std::cout << "DataFilters result: " << result << std::endl;


    } /* if(testcase == 1) */







    if(testcase == 2){

//        UniformVolume<float> *uv = new UniformVolume<float>;
//        uv->AddScalarQuantity("test");

//        delete uv;

        /* Test VTK XML. */
        if(true){
//            size_t size1 = 2;
//            size_t size2 = 2147483647;      /* Max-value of 'int' */

            size_t size1 = 2147483647/200;
            size_t size2 = 199;

//            while(size1*size2 <= 2147483647){
//                size2++;
//            }

            size1 = 16384;
            size2 = 16384;

            size_t npieces = 9;
            size_t start_extent = 0;
            size_t end_extent = size2 - 1;

            Array2D<float> *data = new Array2D<float>(size1, size2*npieces, (float)3);

            std::ofstream file;
            file.open("TestFile_BINARY.vti", std::ios::binary);
//            file.open("TestFile_ASCII.vti");
            char buffer[512];
            int size;

            if(file.is_open()){
                size = sprintf(buffer,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
                file.write(buffer, size);

                size = sprintf(buffer,"  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n", (int)size1-1, (int)size2*(int)npieces-1);
                file.write(buffer, size);

                for(size_t n=0; n<npieces; n++){

                    size = sprintf(buffer,"    <Piece Extent=\"0 %d %d %d 0 0\">\n", (int)size1-1, (int)start_extent, (int)end_extent);
                    file.write(buffer, size);

                    start_extent += size2;
                    end_extent += size2;

                    size = sprintf(buffer,"      <PointData Scalars=\"intensity\">\n");
                    file.write(buffer, size);

    //                size = sprintf(buffer,"        <DataArray type=\"Int16\" Name=\"array_name\" format=\"appended\" RangeMin=\"3\" RangeMax=\"3\" offset=\"0\" />\n");

                    std::stringstream tmpss;
                    tmpss << "        <DataArray type=\"Float32\" Name=\"array_name\" format=\"appended\" RangeMin=\"3\" RangeMax=\"3\" offset=\"" << size1*size2*n << "\" />\n";
                    size = sprintf(buffer,"%s", tmpss.str().c_str());
                    file.write(buffer, size);

                    size = sprintf(buffer,"      </PointData>\n");
                    file.write(buffer, size);

                    size = sprintf(buffer,"    </Piece>\n");
                    file.write(buffer, size);

                }

                size = sprintf(buffer,"  </ImageData>\n");
                file.write(buffer, size);

                size = sprintf(buffer,"  <AppendedData encoding=\"raw\">\n");
                file.write(buffer, size);

                size = sprintf(buffer,"   _");
                file.write(buffer, size);

                file.write((char*)&data->operator [](0), sizeof(float)*size1*size2*npieces);

                size = sprintf(buffer,"\n  </AppendedData>\n");
                file.write(buffer, size);

                size = sprintf(buffer,"</VTKFile>\n");
                file.write(buffer, size);

//                file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
//                file << "  <ImageData WholeExtent=\"0 " << size1-1 << " 0 " << size2-1 << " 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
//                file << "    <Piece Extent=\"0 " << size1-1 << " 0 " << size2-1 << " 0 0\">" << std::endl;
//                file << "      <PointData Scalars=\"intensity\">" << std::endl;
//                file << "        <DataArray type=\"Int16\" Name=\"array_name\" format=\"ascii\">" << std::endl;
//                for(size_t i=0; i<size1; i++){
//                    for(size_t j=0; j<size2; j++){
//                        file << data->operator ()(i,j) << " ";
//                    }
//                    file << std::endl;
//                }
//                file << "        </DataArray>" << std::endl;
//                file << "      </PointData>" << std::endl;
//                file << "    </Piece>" << std::endl;
//                file << "  </ImageData>" << std::endl;
//                file << "</VTKFile>" << std::endl;

//                file.close();

            } else {
                std::cerr << "ERROR: Could not open file." << std::endl;
            }


        }





        /* Test array-transposition. */
        if(false){

#ifdef FFTW_TRANSPOSE
            std::cout << "Transposition tests using FFTW (2D and 3D for 0 <--> 1)" << std::endl;
#else
            std::cout << "Transposition tests using copy arrays" << std::endl;
#endif

            size_t size1 = 2;
            size_t size2 = 3;
            size_t size3 = 4;

            if(true){
                size1 = 450;
                size2 = 1000;
                size3 = 360;
            } else {
                size1 = 2;
                size2 = 2;
                size3 = 2;
            }
            Array2D<int> a2d(size1, size2, 0.0f);
            Array3D<float> a3d(size1, size2, size3, 0.0f);

            float count = 0.0f;
            for(size_t k=0; k<size3; k++){
                for(size_t i=0; i<size1; i++){
                    for(size_t j=0; j<size2; j++){
                        if(k == 0){
                            a2d(i,j) = (int)count;
                        }

                        a3d(i,j,k) = count;

                        count += 1.0f;
                    }
                }
            }


            Array2D<int> a2d_2(a2d);
            Array3D<float> a3d_2(a3d);
            Array3D<float> a3d_original(a3d);

            double t1 = omp_get_wtime();
            a2d.Transpose();
            double t2 = omp_get_wtime();

            float eps = (float)1.0e-5;
            bool pass = true;
            for(size_t i=0; i<size1; i++){
                for(size_t j=0; j<size2; j++){

                    float val1 = (float)a2d(j,i);
                    float val2 = (float)a2d_2(i,j);
                    float diff = val2 - val1;

                    if(fabs(diff) > eps){
                        pass = false;
                    }

                }
            }

            std::cout << std::endl;
            std::cout << "2D transpose check passed: " << StringManip::BoolToString(pass) << std::endl;
            std::cout << "   Time: " << t2-t1 << " [sec]" << std::endl;



            t1 = omp_get_wtime();
            a3d.Transpose(1,0);
            t2 = omp_get_wtime();
            eps = (float)1.0e-5;
            pass = true;
            for(size_t k=0; k<size3; k++){
                for(size_t i=0; i<size1; i++){
                    for(size_t j=0; j<size2; j++){

                        float val1 = a3d(j,i,k);
                        float val2 = a3d_2(i,j,k);
                        float diff = val2 - val1;

                        if(fabs(diff) > eps){
                            pass = false;
                        }

                    }
                }
            }
            a3d = a3d_original;

            std::cout << std::endl;
            std::cout << "3D transpose (0 <--> 1) check passed: " << StringManip::BoolToString(pass) << std::endl;
            std::cout << "   Time: " << t2-t1 << " [sec]" << std::endl;


            t1 = omp_get_wtime();
            a3d.Transpose(0,2);
            t2 = omp_get_wtime();
            pass = true;
            for(size_t k=0; k<size3; k++){
                for(size_t i=0; i<size1; i++){
                    for(size_t j=0; j<size2; j++){

                        float val1 = a3d(k,j,i);
                        float val2 = a3d_2(i,j,k);
                        float diff = val2 - val1;

                        if(fabs(diff) > eps){
                            pass = false;
                        }

                    }
                }
            }
            a3d = a3d_original;

            std::cout << std::endl;
            std::cout << "3D transpose (0 <--> 2) check passed: " << StringManip::BoolToString(pass) << std::endl;
            std::cout << "   Time: " << t2-t1 << " [sec]" << std::endl;

            t1 = omp_get_wtime();
            a3d.Transpose(1,2);
            t2 = omp_get_wtime();
            pass = true;
            for(size_t k=0; k<size3; k++){
                for(size_t i=0; i<size1; i++){
                    for(size_t j=0; j<size2; j++){

                        float val1 = a3d(i,k,j);
                        float val2 = a3d_2(i,j,k);
                        float diff = val2 - val1;

                        if(fabs(diff) > eps){
                            pass = false;
                        }

                    }
                }
            }
            a3d = a3d_original;

            std::cout << std::endl;
            std::cout << "3D transpose (1 <--> 2) check passed: " << StringManip::BoolToString(pass) << std::endl;
            std::cout << "   Time: " << t2-t1 << " [sec]" << std::endl;

        } /* Control testing of array transposition. */

    } /* if(testcase == 2) */

} /* main() */
