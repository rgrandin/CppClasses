/**
 * @file MPI_Custom.cpp
 * @author Robert Grandin
 * @brief Implementation of MPI_Custom namespace.
 */


#include "MPI_Custom.h"

void MPI_Custom::LocalWorkUnits(const int rank, const int np, const int nwork, int &qstart, int &qend)
{
    int basecount = (int)(nwork/np);
    int remainder = nwork - basecount*np;
    int localwork = 0;

    if(rank < remainder){
        localwork = basecount + 1;
        qstart = rank*localwork;
        qend = qstart + localwork - 1;
    } else {
        localwork = basecount;
        qstart = remainder*(basecount + 1) + (rank - remainder)*localwork;
        qend = qstart + localwork - 1;
    }
}


void MPI_Custom::Idx_Transfer_2D(const int singleval, const int size1, int &ival, int &jval)
{
    ival = (int)singleval/size1;
    jval = singleval - size1*ival;
}


void MPI_Custom::Idx_Transfer_To_1D(int &singleval, const int size1, const int ival, const int jval)
{
    singleval = ival*size1 + jval;
}
