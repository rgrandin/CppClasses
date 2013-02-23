/**
 * @file Tree2.h
 * @author Robert Grandin
 * @date 15 March 2012
 * @brief Definition of Tree2 class.
 *
 * @section Class Description & Notes
 *
 * This class representds data structered as a tree with each node having 0, 1, or 2
 * children.  Further, the link between each node may be quantified (e.g., a distance
 * between nodes may be specified).
 *
 * The tree can be visualized in the following image (image provided by Dr. Karin Dorman,
 * Iowa State University).  Note that in this particular image the dashed line has a length
 * of 0 so that S1, S2, and S11 are all children of S0 while still only requiring the use
 * of two children per node.
 * @image html tree2.png "Example tree structure"
 * @image latex tree2.png "Example tree structure" width=7.5cm
 *
 * Any instance of this class created with the copy constructor can access all
 * data in the tree and traverse the nodes.  These copied instances cannot add,
 * delete, or modify nodes.  Only the original instance can modify nodes.  Due to a
 * large number of the member functions requiring the ability to modify nodal properties,
 * the copy constructor is not currently provided.
 *
 * All functions contained within this class are intended for use with the GNU
 * C++ compiler (g++).  Use with other compilers may produce unexpected results
 * and such use is at the users' own risk.
 *
 *
 * @section Revisions
 *
 * @date 15 March 2012
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

#ifndef Tree2_
#define Tree2_

#include <stdlib.h>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <string>
#include <assert.h>
#include <limits>

#include <Array1D.h>
#include <Array2D.h>
#include <StringManip.h>


/**
 * @brief Class for the management of tree structures allowing 0, 1, or 2 children per
 *  per node.
 *
 * This is templated by the datatype used to quantify the link between
 *  nodes.  Multiple datatypes for these links within the same tree is not permitted.
 *  @image html tree2.png
 * @param T Datatype used to quantify the inter-node links.
 * @warning C++11 features, such as move-constructor and move-assignment, require the symbol
 *  "CXX11" to be defined.
 */
template <class T>
class Tree2 {

public:
	/** @brief Node definition */
	struct node
	{
        /** @brief Link quantity to child 1 */
        T qty1;

        /** @brief Link quantity to child 2 */
        T qty2;

        /** @brief Link quantity to parent */
        T qtyparent;

        /** @brief Pointer to child node 1 */
        node *pchild1;

        /** @brief Pointer to child node 2 */
        node *pchild2;

        /** @brief Pointer to parent node */
        node *pparent;

        /** @brief Count number of times node has been visited while traversing the tree */
        int nvisits;

        /** @brief Count number of children of the node */
        int nchildren;

        /** @brief ID number for the node */
        int nodeid;
	};


	/**
     * @brief Default constructor.  The resulting tree contains a single root node with no
     *  children.  The root node is automatically assigned the ID '0'.
     * @pre Memory exists to create this tree.
     * @post Tree2 object created.
	 * @return None.
	 */
    Tree2();


    /**
     * @brief Copy constructor.
     * @param a Reference to existing Tree2 object to be copied.
     */
    Tree2(Tree2<T> &a);


#ifdef CXX11
    /**
     * @brief Move constructor (C++11).
     * @param a Reference to existing Tree2 object to be copied.
     * @warning This function requires C++11 compiler support.
     */
    Tree2(Tree2<T> &&a);
#endif


	/**
	 * @brief Destructor.  All nodes and member data deleted.
     * @pre Tree2 object exists.
     * @post Tree2 object destroyed, along with all nodes.
	 * @return None.
	 * @warning This destructor only deletes pointers/data if the list being
     * 		destroyed is not a copy of another tree (i.e., created with the copy
	 * 		constructor).  Since copied lists share the same nodes, only the
	 * 		original list for a set of nodes will perform the necessary delete
	 * 		operations.
	 */
    virtual ~Tree2();


    /**
     * @brief Copy-assignment operator.
     * @param a Tree2 object being assigned.
     * @return Reference to instance of Tree2.
     */
    Tree2& operator=(Tree2<T> a);


#ifdef CXX11
    /**
     * @brief Move-assignment operator (C++11).
     * @param a Reference to Tree2 object being assigned.
     * @return Reference to instance of Tree2.
     * @warning This function requires C++11 compiler support.
     */
    Tree2& operator=(Tree2<T> &&a);
#endif


	/**
     * @brief Add a child node to the current node.
     *
     * If the current node already has two
     *  children, no action is taken.  If current node has no children, the added node is
     *  considered "child 1".  If one child node is already present, the added node is considered
     *  "child 2".
     * @pre Tree2 object exists.
     * @param id ID number to be assigned to this node.  This must be unique to this node.
     * @param qty Quantity associated with the added node.
     * @post Child node added.  Current node changed to added node.
	 * @return None.
     * @warning Only the original instance of an Tree2 object can add a node.  Any
	 * 		instances created with the copy constructor will have no effect when this
	 * 		function is called.
	 */
    void AddNode(const int id, const T qty);


    /**
     * @brief Add a child node to the current node.
     *
     * If the current node already has two
     *  children, no action is taken.  If current node has no children, the added node is
     *  considered "child 1".  If one child node is already present, the added node is considered
     *  "child 2".
     * @pre Tree2 object exists.
     * @param id ID number to be assigned to this node.  This must be unique to this node.
     * @param qty Quantity associated with the added node.
     * @post Child node added.  Current node changed to added node.
     * @return Integer value identifying if the node was successfully added.  Possible values:
     *      - -1: General failure to add the node.
     *      - 0: Success
     *      - 1: Both children already exist for current node, so node not added.
     *      - 2: User-supplied ID is not unique, so node not added.
     * @warning Only the original instance of an Tree2 object can add a node.  Any
     * 		instances created with the copy constructor will have no effect when this
     * 		function is called.
     */
    int AddNode2(const int id, const T qty);


    /**
      @brief Move to parent of current node.
      @pre Tree2 object exists.
      @post Current node updated to parent of current node.
      @return None.
    */
    void GoToParentNode();


    /**
      @brief Move to root node.
      @pre Tree2 object exists.
      @post Current node updated to root node.
      @return None.
    */
    void GoToRootNode();


    /**
      @brief Move to child of current node.
      @pre Tree2 object exists.
      @param idx 1-based index identifying which child to move to.  Any value other than 1 or 2
        will cause no action to be taken.
      @post Current node updated to the specified child of current node.
      @return None.
      @warning If the specified child does not exist, the current node remains unchanged.
    */
    void GoToChildNode(const int idx);


	/**
     * @brief Delete current node from the tree.
     *
     * After deletion, the current position
     * 		within the list is set to the root node.
     * @pre Tree2 object exists.
	 * @post Current node removed from list
	 * @return None.
	 * @warning Data stored in the current node will be lost.
     * @warning Only the original instance of an Tree2 object can delete a node.  Any
	 * 		instances created with the copy constructor will have no effect when this
	 * 		function is called.
     * @warning Deleting a node will cause all its children to be deleted.
	 */
	void DeleteNode();


	/**
     * @brief Delete all nodes from the list.
     *
     * After deletion a single root node will
     * 		be created with no children.
     * @pre Tree2 object exists.
	 * @post All nodes deleted.
	 * @return None.
	 * @warning All data stored in the nodes will be lost.
     * @warning Only the original instance of an Tree2 object can delete all nodes.  Any
	 * 		instances created with the copy constructor will have no effect when this
	 * 		function is called.
	 */
	void DeleteAllNodes();


	/**
	 * @brief Get data from current node.
     * @pre Tree2 object exists.
     * @param idx Index of quantity to be returned.
     *      - -1: Quantity for link to parent node.
     *      - 1: Quantity for link to child 1.
     *      - 2: Quantity for link to child 2.
	 * @post No changes to object.
     * @return Quantity associated with child 'idx'.  If 'idx' is not -1, 1, or 2, nan is returned.
	 */
    T GetDataValue(int idx) const;


	/**
	 * @brief Set data in current node.
     * @pre Tree2 object exists.
     * @param val Value to be assigned to quantity.
     * @param idx 1-based index of quantity to be set.
     *      - -1: Quantity for link to parent node.
     *      - 1: Quantity for link to child 1.
     *      - 2: Quantity for link to child 2.
	 * @post Data value set.
	 * @return None.
     * @warning Only the original instance of an Tree2 object can modify a node.  Any
	 * 		instances created with the copy constructor will have no effect when this
	 * 		function is called.
	 */
	void SetDataValue(const T val, int idx);


	/**
     * @brief Get the number of nodes in the tree.
     * @pre Tree2 object exists.
	 * @post No changes to object.
     * @return Number of nodes in the tree.
	 */
    int GetNumNodes() const;


	/**
	 * @brief Get pointer to the first data node in the list.
     * @pre Tree2 object exists.
	 * @post No changes to object.
	 * @return Pointer to first data node.
	 */
	node* GetPointerToHead();


	/**
	 * @brief Get pointer to the current data node in the list.
     * @pre Tree2 object exists.
	 * @post No changes to object.
	 * @return Pointer to current data node.
	 */
	node* GetPointerToCurrent();


    /**
      @brief Traverse through the tree, searching for node with the specified ID.

      When intending to simply write the tree structure to std::cout choose an ID which doesn't
        exist in the list.
      @pre Tree2 object exists.
      @param id ID which is being searched for.
      @param output Controls if node information is to be written to std::cout as the
        tree is traversed.  The letter "S" is prepended to each nodal ID and output uses
        the Newick format (http://en.wikipedia.org/wiki/Newick_format).
      @post If specified node is found, the current node will be set to the desired node.
        If specified node is not found, the current node will be set to the root node.
      @return None.
    */
    void TraverseTree(const int id, const bool output);


    /**
      @brief Calculate likelihood of evolved data using the tree structure.
      @pre Tree2 object exists.
      @param params Array1D containing parameters needed by the likelihood calculations.
      @param observed_vals Array1D containing observed data at the tree nodes.  This assumes
        single char data is observed at each node.  If data was not observed, simply use an
        invalid value (e.g., if observed data must be 'a', 'c', 'g', or 't', use any other
        letter for unobserved data).
      @param possible_vals Array1D containing allowable values for observed data.
      @param f Pointer to function which evaluates the evolution transition probabilities.
      @post Current node set to tree root at the end of this function.
      @return Likelihood of the data evolving to its observed state given the current parameter values.
      @warning It is assumed that the first values of params correspond to the unconditional probabilities
        associated with each value in possible_vals (i.e., params(0) is the unconditional probability of
        possible_vals(0)).
    */
    T CalculateLikelihood(const Array1D<T> &params, const Array1D<char> &observed_vals,
                          const Array1D<char> &possible_vals,
                          T(*f)(const char, const char, const T, const Array1D<T>&));


protected:
    /**
      @brief Initialize tree to have a single root node with no children.
      @pre Tree2 object exists.
      @post Tree initialized.
      @return None.
    */
    void Initialize();


    /**
      @brief Check if the user-supplied ID is unique (i.e., hasn't been previously used).
      @pre Tree2 object exists.
      @param id User-supplied ID.
      @post No changes to object.
      @return Boolean value specifying if input value is unique.
    */
    bool IDIsUnique(const int id) const;


    /**
      @brief Reset the visit counter for each node.
      @pre Tree2 object exists.
      @post Visit counter for each node reset to 0.
      @return None.
    */
    void ResetVisitCounter();


    /** @brief Pointer to head element in tree */
	node* pHead;

    /** @brief Pointer to current element in tree */
	node* pCurrent;

    /** @brief Number of nodes in the tree */
	int nodecount;

    /** @brief Track if instance of tree was created as a copy of another tree.
	 * 		If so, the destructor doesn't perform any action and the original
     * 		tree's destructor will handle the pointer-deletion.	 */
	bool iscopy;

    /** @brief Array of node ID numbers which have been set */
    Array1D<int> usedids;

    /** @brief NaN */
    T nanval;


private:
    /**
     * @brief Tree2Swap swaps member information between two Tree2 objects.
     * @param first First Tree2 object.
     * @param second Second Tree2 object.
     */
    friend void Tree2Swap(Tree2<T> &first, Tree2<T> &second)
    {
        std::swap(first.pHead, second.pHead);
        std::swap(first.pCurrent, second.pCurrent);
        std::swap(first.nodecount, second.nodecount);
        std::swap(first.iscopy, second.iscopy);
        std::swap(first.usedids, second.usedids);
        std:;swap(first.nanval, second.nanval);
    }

};

#include "Tree2.cpp"

#endif /* Tree2_ */
