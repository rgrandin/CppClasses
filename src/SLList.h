/**
 * @file SLList.h
 * @author Robert Grandin
 * @date 28 September 2011
 * @brief Definition of SLList class.
 *
 * @section Class Description & Notes
 *
 * This class is intended to assist with the management of singly-linked lists.
 * To use this class, the structure for each node must be previously declared
 * using a 'typedef' statement.  This newly-defined type is then used as the
 * template parameter for this class.  This class is heavily inspired by sample
 * code found at http://www.inversereality.org/tutorials/c++/linkedlists.html.
 *
 * Each node in this class stores data in an array, allowing multiple values to
 * be associated with each node.  Each value must have the same datatype since
 * the all reside in the same array.
 *
 * Any instance of this class created with the copy constructor can access all
 * data in the list and traverse the nodes.  These copied instances cannot add,
 * delete, or modify nodes.  Only the original instance can modify nodes.
 *
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 28 September 2011
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

#ifndef SLList_
#define SLList_

#include <stdlib.h>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <string>
#include <assert.h>
#include <limits>


/*
 * CLASS DEFINITION HERE
 *
 */

/**
 * @brief Class for the management of singly-linked lists.
 *
 * This is templated by
 * 	the datatype stored at each node and the number of data elements at each
 * 	node.  Multiple datatypes at each node is not supported.
 * @param T Datatype used to store the data at each node.
 * @param I Number of elements of datatype T stored at each node.
 */
template <class T, int I>
class SLList {

public:
	/** @brief Node definition */
	struct node
	{
		/** @brief Array of values stored at node */
		T value[I];

		/** @brief Pointer to node following this node */
		node *pNext;
	};


	/**
     * @brief Default constructor.
     *
     * The resulting list contains a single node with
     * 	data initialized to 0.
	 * @pre Memory exists to create this list.
	 * @post SLList object created.
	 * @return None.
	 */
	SLList();


	/**
     * @brief Copy constructor.
     *
     * The newly-created list contains the same nodes
     * 	as the previously-existing list.  The current position in the newly-
     * 	created list is the head node.
	 * @pre Memory exists to create this list and the input list exists.
	 * @param list List copied to  make a new list.
	 * @post SLList object created.
	 * @return None.
	 * @warning Copying a SLList object only copies the pointers to the head and
	 * 		current nodes.  The node data itself is not copied.  Use this copy
	 * 		functionality carefully as the addition or deletion of nodes will be
	 * 		seen by all SLList objects that share the same head node.
	 */
	SLList(SLList<T,I> &list);


	/**
	 * @brief Destructor.  All nodes and member data deleted.
	 * @pre SLList object exists.
	 * @post SLList object destroyed, along with all nodes.
	 * @return None.
	 * @warning This destructor only deletes pointers/data if the list being
	 * 		destroyed is not a copy of another list (i.e., created with the copy
	 * 		constructor).  Since copied lists share the same nodes, only the
	 * 		original list for a set of nodes will perform the necessary delete
	 * 		operations.
	 */
    virtual ~SLList();


	/**
	 * @brief Add a single node to the end of the list.
	 * @pre SLList object exists.
	 * @post Node appended to end of list.
	 * @return None.
	 * @warning Only the original instance of an SLList object can add a node.  Any
	 * 		instances created with the copy constructor will have no effect when this
	 * 		function is called.
	 */
	void AddNode();


	/**
     * @brief Advance forward in the list one node.
     *
     * If current node is the end
     * 	of the list, the pointer will remain on the last node.
	 * @pre SLList object exists.
	 * @post Current node advanced to next node.
	 * @return None.
	 */
	void Advance();


	/**
     * @brief Move backward in the list one node.
     *
     * If the current node is the head
     * 	of the list, the pointer will remain on the first node.
	 * @pre SLList object exists.
	 * @post Current node moved to previous node.
	 * @return None.
	 */
	void Rewind();


	/**
     * @brief Delete current node from the list.
     *
     * After deletion, the current position
     *  within the list is set to the first node.
	 * @pre SLList object exists.
	 * @post Current node removed from list
	 * @return None.
	 * @warning Data stored in the current node will be lost.
	 * @warning Only the original instance of an SLList object can delete a node.  Any
	 * 		instances created with the copy constructor will have no effect when this
	 * 		function is called.
	 */
	void DeleteNode();


	/**
     * @brief Delete all nodes from the list.
     *
     * After deletion a single node will
     * 	be created with its data initialized to 0.
	 * @pre SLList object exists.
	 * @post All nodes deleted.
	 * @return None.
	 * @warning All data stored in the nodes will be lost.
	 * @warning Only the original instance of an SLList object can delete all nodes.  Any
	 * 		instances created with the copy constructor will have no effect when this
	 * 		function is called.
	 */
	void DeleteAllNodes();


	/**
	 * @brief Get data from current node.
	 * @pre SLList object exists.
	 * @param idx 0-based index of the data to be retrieved.
	 * @post No changes to object.
	 * @return Value stored in current node at index 'idx'.
	 */
    T GetDataValue(int idx) const;


	/**
	 * @brief Set data in current node.
	 * @pre SLList object exists.
	 * @param val Value to be placed in node.
	 * @param idx 0-based index of the data to be set.
	 * @post Data value set.
	 * @return None.
	 * @warning Only the original instance of an SLList object can modify a node.  Any
	 * 		instances created with the copy constructor will have no effect when this
	 * 		function is called.
	 */
	void SetDataValue(const T val, int idx);


	/**
	 * @brief Get the number of nodes in the list.
	 * @pre SLList object exists.
	 * @post No changes to object.
	 * @return Number of nodes in the list.
	 */
    int GetNumNodes() const;


	/**
	 * @brief Print the node values to the screen.
	 * @pre SLList object exists.
	 * @param wraplim Number of nodes to print on a single line before wrapping
	 * 		output to another line.  Default value is 10.
	 * @param blankline Specify if a blank line should be printed between each
	 * 		node's output.  Default is false.
	 * @post Current position within the list set to first node.
	 * @return None.
	 */
	void PrintNodes(int wraplim, bool blankline);


	/**
	 * @brief Set position in list to specified node.
	 * @pre SLList object exists.
	 * @param nodenum 1-based index node number corresponding to the position desired.
	 * @post Current position at specified node.
	 * @return None.
	 * @warning If the node number specified is greater than the number of nodes
	 * 		in the list, the position when this function exits will be the final
	 * 		node in the list.
	 */
	void GoToNode(const int nodenum);


	/**
	 * @brief Set position in list to final node.
	 * @pre SLList object exists.
	 * @post Current position at final node.
	 * @return None.
	 */
	void GoToFinalNode();


	/**
	 * @brief Get pointer to the first data node in the list.
	 * @pre SLList object exists.
	 * @post No changes to object.
	 * @return Pointer to first data node.
	 */
	node* GetPointerToHead();


	/**
	 * @brief Get pointer to the current data node in the list.
	 * @pre SLList object exists.
	 * @post No changes to object.
	 * @return Pointer to current data node.
	 */
	node* GetPointerToCurrent();


	/**
	 * @brief Copy values from one SLList object to another.  This is a true copy
	 * 		and not a reference.
	 * @pre SLList object exists.
	 * @param sourcelist List from which values are to be copied.
	 * @post SLList pointer addresses to head and current nodes copied.
	 * @return None.
	 */
	SLList<T,I>& operator=(const SLList<T,I>& sourcelist);



protected:
	/** @brief Pointer to head element in list */
	node* pHead;

	/** @brief Pointer to current element in list */
	node* pCurrent;

	/** @brief Number of nodes in the list */
	int nodecount;

	/** @brief Track if instance of list was created as a copy of another list.
	 * 		If so, the destructor doesn't perform any action and the original
	 * 		list's destructor will handle the pointer-deletion.	 */
	bool iscopy;
};

#include "SLList.cpp"

#endif /* SLList_ */
