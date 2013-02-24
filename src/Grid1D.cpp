/**
 * @file Grid1D.cpp
 * @author Robert Grandin
 * @brief Implementation of Grid1D class.
 */



#include <Grid1D.h>


/* ==================================================================
 * ================
 * ================    PRIVATE FUNCTIONS
 * ================
 */







/* ==================================================================
 * ================
 * ================    PUBLIC FUNCTIONS
 * ================
 */

template <class T>
Grid1D<T>::Grid1D()
{
    pi = 4.0e0*atan(1.0e0);
    nscalars = 0;
    nvectors = 0;
    gridname = "Gridname_Not_Set";
    qtysize =256;
    inc_snapshots = -1;
    count_snapshots = 0;
    name_snapshots = "snapshot";
}


template <class T>
Grid1D<T>::Grid1D(const Grid1D<T> &a) : Grid1D()
{
    Grid1DSwap(*this, a);
}


#ifdef CXX11
template <class T>
Grid1D<T>::Grid1D(Grid1D<T> &&a) : Grid1D<T>()
{
    Grid1DSwap(*this, a);
}
#endif



template <class T>
Grid1D<T>::~Grid1D()
{
    /*
      Delete arrays of 2D and 3D quantities not required here due to
      the destructor of PArray1D doing this.
     */

}


template <class T>
void Grid1D<T>::WriteGrid(std::string file)
{
	/*
	 * VTK STRUCTURED GRID FILE CREATED TO SHOW THE GRID.  ALL GRID POINTS ARE
	 * ASSOCIATED WITH A SCALAR VALUE OF ZERO.  THE INTENT IS FOR THE GRID TO
	 * BE VIEWED IN PARAVIEW.
	 */

    std::fstream ssfile;
    ssfile.open(file.c_str(),std::ios::out);
    ssfile << "# vtk DataFile Version 2.0" << std::endl;
    ssfile << gridname << std::endl;
    ssfile << "ASCII" << std::endl;
    ssfile << "DATASET STRUCTURED_GRID" << std::endl;
    ssfile << "DIMENSIONS " << nx << " 1 1" << std::endl;
    ssfile << "POINTS " << nx << " float" << std::endl;
    for(int j=0; j<nx; j++){
        ssfile << xcoords(i,j) << " 0.0 0.0" << std::endl;
    }
    ssfile << "POINT_DATA " << nx << std::endl;
    ssfile << "SCALARS grid float 1" << std::endl;
    ssfile << "LOOKUP_TABLE default" << std::endl;
    for(int j=0; j<nx; j++){
        ssfile << "0" << std::endl;
    }
    ssfile.close();
}


template <class T>
void Grid1D<T>::WriteSolution(std::string file)
{
    /*
      Write VTK Structured Grid file to show the solution.  Each scalar and
      vector quantity will be included in the file.
     */

    std::fstream ssfile;
    ssfile.open(file.c_str(),std::ios::out);
    ssfile << "# vtk DataFile Version 2.0" << std::endl;
    ssfile << gridname << std::endl;
    ssfile << "ASCII" << std::endl;
    ssfile << "DATASET STRUCTURED_GRID" << std::endl;
    ssfile << "DIMENSIONS " << nx << " 1 1" << std::endl;
    ssfile << "POINTS " << nx << " float" << std::endl;
    for(int j=0; j<nx; j++){
        ssfile << xcoords(i,j) << " " << ycoords(i,j) << " 0.0" << std::endl;
    }
    ssfile << "POINT_DATA " << nx << std::endl;
    for(int s=0; s<nscalars; s++){
        std::string name(scalar_names(s)->substr());
        ssfile << "SCALARS " << name << " float 1" << std::endl;
        ssfile << "LOOKUP_TABLE default" << std::endl;
        for(int j=0; j<nx; j++){
            ssfile << pscalars(s)->operator ()(i) << std::endl;
        }
    }
    for(int v=0; v<nvectors; v++){
        std::string name(vector_names(v)->substr());
        ssfile << "VECTORS " << name << " float" << std::endl;
        for(int j=0; j<nx; j++){
            ssfile << pvectors(v)->operator ()(i,0) << " " << pvectors(v)->operator ()(i,1) << " " <<
                   pvectors(v)->operator ()(i,2) << std::endl;
        }
    }
    ssfile.close();

}


template <class T>
int Grid1D<T>::GetSize()
{
    return nx;
}


template <class T>
void Grid1D<T>::AddScalarQuantity(const std::string name)
{
    // Require that a grid has been generated
    if(isogrid == true || ishgrid == true){
        // Increase count of number of scalar quantities
        nscalars++;

        /*
          Increase the size of the PArray1D object containing pointers to
          Array2D objects containing data.  The data must be copied into temporary
          arrays since the current Array2D objects are destroyed when the
          PArray1D object is resized.  The data is then copied back into the
          resized PArray1D object.

          During the copy process, the source object is deleted immediately after
          its data has been copied.  This is to reduce memory requirements so that
          only one extra object's worth of memory is required during the
          copy process.  Note that setting the resulting pointer to NULL is
          needed to avoid freeing the memory more than once upon deletion of
          the temporary PArray1D object.

          After the copy is complete, add a new Array2D object to the final
          entry in the resized PArray1D object, set its size equal to the grid
          size, and all points initialized to 0.0.
         */
        PArray1D<Array1D<T>*> tmparray(nscalars);
        for(int i=0; i<nscalars-1; i++){
            tmparray(i) = new Array1D<T>;
            tmparray(i)->ResetSize(ny,nx);
            for(int ii=0; ii<ny; ii++){
                for(int jj=0; jj<nx; jj++){
                    tmparray(i)->operator ()(ii,jj) = pscalars(i)->operator ()(ii,jj);
                }
            }
            delete pscalars(i);
            pscalars(i) = NULL;
        }
        pscalars.ResetSize(nscalars);
        for(int i=0; i<nscalars-1; i++){
            pscalars(i) = new Array1D<T>;
            pscalars(i)->ResetSize(ny,nx);
            for(int ii=0; ii<ny; ii++){
                for(int jj=0; jj<nx; jj++){
                    pscalars(i)->operator ()(ii,jj) = tmparray(i)->operator ()(ii,jj);
                }
            }
            delete tmparray(i);
            tmparray(i) = NULL;
        }
        pscalars(nscalars-1) = new Array1D<T>;
        pscalars(nscalars-1)->ResetSize(ny,nx);
        pscalars(nscalars-1)->ResetVal((T)0.0e0);


        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nscalars);
        for(int i=0; i<nscalars-1; i++){
            tmparray2(i) = new std::string;
            tmparray2(i)->reserve(qtysize);
            tmparray2(i)->assign(scalar_names(i)->substr());
            delete scalar_names(i);
            scalar_names(i) = NULL;
        }
        scalar_names.ResetSize(nscalars);
        for(int i=0; i<nscalars-1; i++){
            scalar_names(i) = new std::string;
            scalar_names(i)->reserve(qtysize);
            scalar_names(i)->assign(tmparray2(i)->substr());
            delete tmparray2(i);
            tmparray2(i) = NULL;
        }
        scalar_names(nscalars-1) = new std::string;
        scalar_names(nscalars-1)->reserve(qtysize);
        scalar_names(nscalars-1)->assign(name);
    }
}


template <class T>
void Grid1D<T>::RemoveScalarQuantity(const int qty)
{
    // Require that a grid has been generated
    if((isogrid == true || ishgrid == true) && qty < nscalars){
        // Decrease count of number of scalar quantities
        nscalars--;

        /*
          Decrease the size of the PArray1D object containing pointers to
          Array2D objects containing data.  The data must be copied into temporary
          arrays since the current Array2D objects are destroyed when the
          PArray1D object is resized.  The data is then copied back into the
          resized PArray1D object.

          During the copy process, the source object is deleted immediately after
          its data has been copied.  This is to reduce memory requirements so that
          only one extra object's worth of memory is required during the
          copy process.  Note that setting the resulting pointer to NULL is
          needed to avoid freeing the memory more than once upon deletion of
          the temporary PArray1D object.

          After the copy is complete, add a new Array2D object to the final
          entry in the resized PArray1D object, set its size equal to the grid
          size, and all points initialized to 0.0.
         */
        PArray1D<Array1D<T>*> tmparray(nscalars+1);
        for(int i=0; i<nscalars+1; i++){
            tmparray(i) = new Array1D<T>;
            tmparray(i)->ResetSize(ny,nx);
            for(int ii=0; ii<ny; ii++){
                for(int jj=0; jj<nx; jj++){
                    tmparray(i)->operator ()(ii,jj) = pscalars(i)->operator ()(ii,jj);
                }
            }
            delete pscalars(i);
            pscalars(i) = NULL;
        }
        pscalars.ResetSize(nscalars);
        for(int i=0; i<nscalars+1; i++){
            int itmp = i;
            if(i != qty){
                if(i > qty){
                    itmp = i - 1;   // Reduce index by 1 since loop is past removed quantity
                }
                pscalars(itmp) = new Array1D<T>;
                pscalars(itmp)->ResetSize(ny,nx);
                for(int ii=0; ii<ny; ii++){
                    for(int jj=0; jj<nx; jj++){
                        pscalars(itmp)->operator ()(ii,jj) = tmparray(i)->operator ()(ii,jj);
                    }
                }
            }
            delete tmparray(i);
            tmparray(i) = NULL;
        }


        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nscalars+1);
        for(int i=0; i<nscalars+1; i++){
            tmparray2(i) = new std::string;
            tmparray2(i)->reserve(qtysize);
            tmparray2(i)->assign(scalar_names(i)->substr());
            delete scalar_names(i);
            scalar_names(i) = NULL;
        }
        scalar_names.ResetSize(nscalars);
        for(int i=0; i<nscalars+1; i++){
            int itmp = i;
            if(i != qty){
                if(i > qty){
                    itmp = i - 1;
                }
                scalar_names(itmp) = new std::string;
                scalar_names(itmp)->reserve(qtysize);
                scalar_names(itmp)->assign(tmparray2(i)->substr());
                delete tmparray2(i);
                tmparray2(i) = NULL;
            }
        }
    }
}


template <class T>
void Grid1D<T>::AddVectorQuantity(const std::string name)
{
    // Require that a grid has been generated
    if(isogrid == true || ishgrid == true){
        // Increase count of number of vector quantities
        nvectors++;

        /*
          Increase the size of the PArray1D object containing pointers to
          Array3D objects containing data.  The data must be copied into temporary
          arrays since the current Array3D objects are destroyed when the
          PArray1D object is resized.  The data is then copied back into the
          resized PArray1D object.

          During the copy process, the source object is deleted immediately after
          its data has been copied.  This is to reduce memory requirements so that
          only one extra object's worth of memory is required during the
          copy process.  Note that setting the resulting pointer to NULL is
          needed to avoid freeing the memory more than once upon deletion of
          the temporary PArray1D object.

          After the copy is complete, add a new Array3D object to the final
          entry in the resized PArray1D object, set its size equal to the grid
          size, and all points initialized to 0.0.
         */
        PArray1D<Array2D<T>*> tmparray(nvectors);
        for(int i=0; i<nvectors-1; i++){
            tmparray(i) = new Array2D<T>;
            tmparray(i)->ResetSize(ny,nx,3);
            for(int kk=00; kk<3; kk++){
                for(int ii=0; ii<ny; ii++){
                    for(int jj=0; jj<nx; jj++){
                        tmparray(i)->operator ()(ii,jj,kk) = pvectors(i)->operator ()(ii,jj,kk);
                    }
                }
            }
            delete pvectors(i);
            pvectors(i) = NULL;
        }
        pvectors.ResetSize(nvectors);
        for(int i=0; i<nvectors-1; i++){
            pvectors(i) = new Array2D<T>;
            pvectors(i)->ResetSize(ny,nx,3);
            for(int kk=0; kk<3; kk++){
                for(int ii=0; ii<ny; ii++){
                    for(int jj=0; jj<nx; jj++){
                        pvectors(i)->operator ()(ii,jj,kk) = tmparray(i)->operator ()(ii,jj,kk);
                    }
                }
            }
            delete tmparray(i);
            tmparray(i) = NULL;
        }
        pvectors(nvectors-1) = new Array2D<T>;
        pvectors(nvectors-1)->ResetSize(ny,nx,3);
        pvectors(nvectors-1)->ResetVal((T)0.0e0);


        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nvectors);
        for(int i=0; i<nvectors-1; i++){
            tmparray2(i) = new std::string;
            tmparray2(i)->reserve(qtysize);
            tmparray2(i)->assign(vector_names(i)->substr());
            delete vector_names(i);
            vector_names(i) = NULL;
        }
        vector_names.ResetSize(nvectors);
        for(int i=0; i<nvectors-1; i++){
            vector_names(i) = new std::string;
            vector_names(i)->reserve(qtysize);
            vector_names(i)->assign(tmparray2(i)->substr());
            delete tmparray2(i);
            tmparray2(i) = NULL;
        }
        vector_names(nvectors-1) = new std::string;
        vector_names(nvectors-1)->reserve(qtysize);
        vector_names(nvectors-1)->assign(name);
    }
}


template <class T>
void Grid1D<T>::RemoveVectorQuantity(const int qty)
{
    // Require that a grid has been generated
    if((isogrid == true || ishgrid == true) && qty < nvectors){
        // Decrease count of number of scalar quantities
        nvectors--;

        /*
          Decrease the size of the PArray1D object containing pointers to
          Array3D objects containing data.  The data must be copied into temporary
          arrays since the current Array3D objects are destroyed when the
          PArray1D object is resized.  The data is then copied back into the
          resized PArray1D object.

          During the copy process, the source object is deleted immediately after
          its data has been copied.  This is to reduce memory requirements so that
          only one extra object's worth of memory is required during the
          copy process.  Note that setting the resulting pointer to NULL is
          needed to avoid freeing the memory more than once upon deletion of
          the temporary PArray1D object.

          After the copy is complete, add a new Array3D object to the final
          entry in the resized PArray1D object, set its size equal to the grid
          size, and all points initialized to 0.0.
         */
        PArray1D<Array2D<T>*> tmparray(nvectors+1);
        for(int i=0; i<nvectors+1; i++){
            tmparray(i) = new Array2D<T>;
            tmparray(i)->ResetSize(ny,nx,3);
            for(int kk=0; kk<3; kk++){
                for(int ii=0; ii<ny; ii++){
                    for(int jj=0; jj<nx; jj++){
                        tmparray(i)->operator ()(ii,jj,kk) = pvectors(i)->operator ()(ii,jj,kk);
                    }
                }
            }
            delete pvectors(i);
            pvectors(i) = NULL;
        }
        pvectors.ResetSize(nvectors);
        for(int i=0; i<nvectors+1; i++){
            int itmp = i;
            if(i != qty){
                if(i > qty){
                    itmp = i - 1;   // Reduce index by 1 since loop is past removed quantity
                }
                pvectors(itmp) = new Array2D<T>;
                pvectors(itmp)->ResetSize(ny,nx,3);
                for(int kk=0; kk<3; kk++){
                    for(int ii=0; ii<ny; ii++){
                        for(int jj=0; jj<nx; jj++){
                            pvectors(itmp)->operator ()(ii,jj,kk) = tmparray(i)->operator ()(ii,jj,kk);
                        }
                    }
                }
            }
            delete tmparray(i);
            tmparray(i) = NULL;
        }


        /*
          Like above, expand the arrays related to quantity labels.
         */
        PArray1D<std::string*> tmparray2(nvectors+1);
        for(int i=0; i<nvectors+1; i++){
            tmparray2(i) = new std::string;
            tmparray2(i)->reserve(qtysize);
            tmparray2(i)->assign(vector_names(i)->substr());
            delete vector_names(i);
            vector_names(i) = NULL;
        }
        vector_names.ResetSize(nvectors);
        for(int i=0; i<nvectors+1; i++){
            int itmp = i;
            if(i != qty){
                if(i > qty){
                    itmp = i - 1;
                }
                vector_names(itmp) = new std::string;
                vector_names(itmp)->reserve(qtysize);
                vector_names(itmp)->assign(tmparray2(i)->substr());
                delete tmparray2(i);
                tmparray2(i) = NULL;
            }
        }
    }
}



// () OPERATOR
template < class T > inline
T& Grid1D<T>::operator()(int ind1, int qty)
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 >= 0 && ind1 < ny);
        assert(ind2 >= 0 && ind2 < nx);
        return pscalars(qty)->operator()(ind1);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pscalars(qty)->operator()(ind1);
    #endif
}


template < class T > inline
const T& Grid1D<T>::operator()(int ind1, int qty) const
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 >= 0 && ind1 < ny);
        assert(ind2 >= 0 && ind2 < nx);
        return pscalars(qty)->operator ()(ind1);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pscalars(qty)->operator()(ind1);
    #endif
}


template < class T > inline
T& Grid1D<T>::operator()(int ind1, int ind2, int qty)
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 >= 0 && ind1 < ny);
        assert(ind2 >= 0 && ind2 < ny);
        assert(ind3 >= 0 && ind3 < 3);
        return pvectors(qty)->operator ()(ind1,ind2);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pvectors(qty)->operator ()(ind1,ind2);
    #endif
}


template < class T > inline
const T& Grid1D<T>::operator()(int ind1, int ind2, int qty) const
{
    #ifndef RELEASE
        /*
         * "RELEASE" NOT DEFINED, SO DEFAULT TO BOUNDS-CHECKING ON ARRAY ACCESS
         */
        assert(ind1 >= 0 && ind1 < ny);
        assert(ind2 >= 0 && ind2 < nx);
        assert(ind3 >= 0 && ind3 < 3);
        return pvectors(qty)->operator ()(ind1,ind2,);
    #else
        /*
         * "RELEASE" DEFINED, SO DISABLE BOUNDS-CHECKING
         */
        return pvectors(qty)->operator ()(ind1,ind2,);
    #endif
}


template <class T>
Grid1D<T>& Grid1D<T>::operator=(Grid1D<T> a)
{
    Grid1DSwap(*this, a);
    return *this;
}


template <class T>
void Grid1D<T>::setGridName(const std::string name)
{
    gridname = name;
}


template <class T>
std::string Grid1D<T>::GridName() const
{
    return gridname;
}


template <class T>
void Grid1D<T>::setScalarQuantityName(const int qty, const std::string name)
{
    scalar_names(qty)->assign(name);
}


template <class T>
std::string Grid1D<T>::ScalarQuantityName(const int qty)
{
    return scalar_names(qty)->substr();
}


template <class T>
void Grid1D<T>::setVectorQuantityName(const int qty, const std::string name)
{
    vector_names(qty)->assign(name);
}


template <class T>
std::string Grid1D<T>::VectorQuantityName(const int qty)
{
    return vector_names(qty)->substr();
}
