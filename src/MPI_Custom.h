/**
 * @file MPI_Custom.h
 * @author Robert Grandin
 * @date 14 April 2012
 * @brief Definition of MPI_Custom namespace.
 *
 *
 * @section Revisions
 *
 * @date 14 April 2012
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

#ifndef MPI_CUSTOM_
#define MPI_CUSTOM_


/**
 * @brief Collection of functions to assist MPI programming.
 *
 * This namespace contains utility functions related to MPI programming.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
*/
namespace MPI_Custom {


/**
  @brief Calculate starting and ending indices for a single MPI process.

  This assumes that the work units to be divided are indexed by a single parameter.
    Additionally, this assumes that the work will be divided evenly among processes
    and that work assigned to a single process is contiguous in the index identifying
    the work units.  If the work units are not able to be evenly divided among the
    number of MPI processes (i.e., nwork % np != 0), a single additional unit of
    work is assigned to the lower-rank processes until all remaining work units
    have been assigned.
  @param rank Rank of local MPI process.
  @param np Number of MPI processes.
  @param nwork Number of work units to be divided among MPI processes.
  @param qstart Index of first work-unit to be performed by local MPI process.
  @param qend Index of final work-unit to be performed by local MPI process.
  @return None.
*/
void LocalWorkUnits(const int rank, const int np, const int nwork, int &qstart, int &qend);


/**
  @brief Convert from single index to double-index, such as when a single index
    is used to define elements on a two-dimensional grid.

  For example, 2D
    indices can be mapped to a single index by @f$ q = j + n_x*i @f$.
  @param singleval Value of single-index ('q' in example equation).
  @param size1 Size of first dimension of two-dimensional indices ('nx' in example equation).
  @param ival Corresponding value for outer dimension ('i' in example equation).
  @param jval Corresponding value for inner dimension ('j' in example equation).
  @return None.
*/
void Idx_Transfer_2D(const int singleval, const int size1, int &ival, int &jval);


/**
  @brief Convert from double index to single index, such when a single index
    is used to define elements on a two-dimensional grid.

  For example, 2D
    indices can be mapped to a single index by @f$ q = j + n_x*i @f$.
  @param singleval Value of single-index ('q' in example equation).
  @param size1 Size of first dimension of two-dimensional indices ('nx' in example equation).
  @param ival Corresponding value for outer dimension ('i' in example equation).
  @param jval Corresponding value for inner dimension ('j' in example equation).
  @return None.
*/
void Idx_Transfer_To_1D(int &singleval, const int size1, const int ival, const int jval);


}


#endif
