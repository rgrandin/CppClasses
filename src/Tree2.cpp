/**
 * @file Tree2.cpp
 * @author Robert Grandin
 * @brief Implementation of Tree2 class.
 */


#include <Tree2.h>

/*
 * CONSTRUCTORS
 */

template <class T>
Tree2<T>::Tree2()
{
    /* Identify this tree as not being a copy (i.e., it's an original tree) */
	iscopy = false;

    /* Initialize nanval */
    nanval = std::numeric_limits<T>::quiet_NaN();

    /* Initialize root node */
    Initialize();
}


template <class T>
Tree2<T>::Tree2(const Tree2<T> &a) : Tree2()
{
    Tree2Swap(*this, a);
}


#ifdef CXX11
template <class T>
Tree2<T>::Tree2(Tree2<T> &&a) : Tree2<T>()
{
    Tree2Swap(*this, a);
}
#endif



/*
 * DESTRUCTOR
 */
template <class T>
Tree2<T>::~Tree2()
{
	/*
     * ONLY PERFORM THE DESTRUCTION IF THIS INSTANCE OF THE TREE IS *NOT* A COPY.
	 */

    if(iscopy == false){
        /* Go to root node and then delete it, taking all child nodes with it */
        pCurrent = pHead;
        DeleteNode();
    }

}




/*
 * PUBLIC MEMBER FUNCTIONS
 */


template <class T>
Tree2<T>& Tree2<T>::operator=(Tree2<T> a)
{
    Tree2Swap(*this, a);
    return *this;
}



template <class T>
void Tree2<T>::Initialize()
{
    /* Create root node */
    pHead = new node;

    /* Set current node pointer to root node */
    pCurrent = pHead;

    /* Set pointer to parent node equal to NULL */
    pCurrent->pparent = NULL;

    /* Set pointers to children nodes to NULL */
    pCurrent->pchild1 = NULL;
    pCurrent->pchild2 = NULL;

    /* Initialize quatities for children to 0 */
    pCurrent->qty1 = nanval;
    pCurrent->qty2 = nanval;
    pCurrent->qtyparent = nanval;

    nodecount = 1;

    /* Initialize 'usedids' to 20 elements. */
    usedids.ResetSize(20,(int)nanval);

    /* Set and record the automatically-assigned ID number of '0' for this root node */
    pCurrent->nodeid = 0;
    usedids(0) = 0;

    /* Initialize nvisits */
    pCurrent->nvisits = 0;

    /* Initialize nchildren */
    pCurrent->nchildren = 0;
}


template <class T>
void Tree2<T>::AddNode(const int id, const T qty)
{
    int tmpval = AddNode2(id,qty);
    tmpval++;                           /* Added to remove 'unused variable' compiler warnings */
}


template <class T>
int Tree2<T>::AddNode2(const int id, const T qty)
{
    /* If number of nodes is within 5 of the size of the 'usedids' array, add 20 elements to 'usedids'  */
    if((usedids.GetDim() - nodecount) < 5){
        Array1D<int> tmp(usedids.GetDim());
        for(int i=0; i<usedids.GetDim(); i++){
            tmp(i) = usedids(i);
        }

        usedids.ResetSize(usedids.GetDim()+20,(int)nanval);
        for(int i=0; i<nodecount; i++){
            usedids(i) = tmp(i);
        }
    }

    int retval = -1;
	if(iscopy == false){
        /* Check if first child has been created.  If not, add first child. */
        if(pCurrent->nchildren == 0){
            /* Check that user-supplied ID is valid */
            if(IDIsUnique(id) == true){
                pCurrent->pchild1 = new node;
                pCurrent->qty1 = qty;
                pCurrent->pchild1->pparent = pCurrent;
                pCurrent->pchild1->pchild1 = NULL;
                pCurrent->pchild1->pchild2 = NULL;
                pCurrent->pchild1->qty1 = nanval;
                pCurrent->pchild1->qty2 = nanval;
                pCurrent->pchild1->qtyparent = qty;
                pCurrent->pchild1->nvisits = 0;
                pCurrent->pchild1->nchildren = 0;
                pCurrent->pchild1->nodeid = id;
                pCurrent->nchildren += 1;
                pCurrent = pCurrent->pchild1;
                usedids(nodecount) = id;
                nodecount++;
                retval = 0;
            } else {
                retval = 2;
            }
        } else {
            /* If first child has been created, check second child */
            if(pCurrent->nchildren == 1){
                if(IDIsUnique(id) == true){
                    pCurrent->pchild2 = new node;
                    pCurrent->qty2 = qty;
                    pCurrent->pchild2->pparent = pCurrent;
                    pCurrent->pchild2->pchild1 = NULL;
                    pCurrent->pchild2->pchild2 = NULL;
                    pCurrent->pchild2->qty1 = nanval;
                    pCurrent->pchild2->qty2 = nanval;
                    pCurrent->pchild2->qtyparent = qty;
                    pCurrent->pchild2->nvisits = 0;
                    pCurrent->pchild2->nchildren = 0;
                    pCurrent->pchild2->nodeid = id;
                    pCurrent->nchildren += 1;
                    pCurrent = pCurrent->pchild2;
                    usedids(nodecount) = id;
                    nodecount++;
                    retval = 0;
                } else {
                    retval = 2;
                }
            } else {
                /* If this point is reached, both children already exist, so the node
                  cannot be added. */
                retval = 1;
            }
        }
	}

    return retval;
}


template <class T>
void Tree2<T>::GoToParentNode()
{
    if(pCurrent->pparent != NULL){
        pCurrent = pCurrent->pparent;
    }
}


template <class T>
void Tree2<T>::GoToRootNode()
{
    pCurrent = pHead;
}


template <class T>
void Tree2<T>::GoToChildNode(const int idx)
{
    switch(idx){
    case 1:
        if(pCurrent->pchild1 != NULL){
            pCurrent = pCurrent->pchild1;
        }
        break;
    case 2:
        if(pCurrent->pchild2 != NULL){
            pCurrent = pCurrent->pchild2;
        }
        break;
    default:
        break;
    }
}


template <class T>
void Tree2<T>::DeleteNode()
{
	if(iscopy == false){
        /*
          NOTE: The terminology "parent" and "child" used when referring to the relationships between nodes.  As
                a result, this genealogical reference scheme will be used to describe the deletion algorithm.  A
                node is termed "older" if it is closer to the root node, and "younger" if it is farther from the
                root node.

                The general approach for this function is to find the youngest node, and begin deleting the tree
                working from youngest to oldest.  Failure to delete in this fashion will result in memory leaks.
                Once a node is deleted, access to its descendants is lost, requiring deletion to begin at the
                youngest generation and work "up".
        */

        /* Save location of current node. */
        node *psave = pCurrent;

        /* Set flag signifying that I'm working with descendant nodes. */
        bool isdescendant = true;

        /* Create variable for last-visited node. */
        node *plast = pCurrent;

        /* Loop through descendants until the youngest node is reached.  Such a node will not have any children. */
        while(isdescendant == true){
            while(pCurrent->pchild1 != NULL || pCurrent->pchild2 != NULL){

                /* If child1 exists, go to child 1 */
                if(pCurrent->pchild1 != NULL){
                    GoToChildNode(1);
                } else {
                    /* If child2 exists, go to child 2 */
                    if(pCurrent->pchild2 != NULL){
                        GoToChildNode(2);
                    }
                }
            }

            /* When this point is reached, the current node has no children and thus can be deleted. */
            plast = pCurrent;                   /* Pointer to node to be deleted */
            pCurrent = pCurrent->pparent;       /* Move current pointer to parent node */
            if(plast == pCurrent->pchild1){
                pCurrent->pchild1 = NULL;       /* Identify child1 as the node to be deleted and set it to NULL */
                pCurrent->qty1 = nanval;        /* Reset quantity */
                pCurrent->nchildren--;
            }
            if(plast == pCurrent->pchild2){
                pCurrent->pchild2 = NULL;       /* Identify child2 as the node to be deleted and set it to NULL */
                pCurrent->qty2 = nanval;        /* Reset quantity */
                pCurrent->nchildren--;
            }
            delete plast;                       /* Delete the node 'plast' */
            plast = NULL;
            nodecount--;                        /* Decrease nodecount */

            /* Check that if the current node only has one child that the child is child1 */
            if(pCurrent != NULL){
                if(pCurrent->pchild1 == NULL && pCurrent->pchild2 != NULL){
                    pCurrent->pchild1 = pCurrent->pchild2;
                    pCurrent->qty1 = pCurrent->qty2;
                    pCurrent->pchild2 = NULL;
                    pCurrent->qty2 = nanval;
                }
            }

            /* If current node is the saved node, change 'isdescendant'. */
            if(pCurrent == psave && pCurrent->pchild1 == NULL && pCurrent->pchild2 == NULL){
                isdescendant = false;
            }
        }

        /* When this point is reached, all descendents of 'psave' have been deleted.  Finally, 'psave' is deleted.
           Note that this step requires checking that the node to be deleted isn't the root node.  If it is the root
           node, it must be handled differently as the root node has no parent.
        */
        if(pCurrent != pHead){  /* Case where node to be deleted is *not* the root node */
            plast = pCurrent;                   /* Pointer to node to be deleted */
            pCurrent = pCurrent->pparent;       /* Move current pointer to parent node */
            if(plast == pCurrent->pchild1){
                pCurrent->pchild1 = NULL;       /* Identify child1 as the node to be deleted and set it to NULL */
            }
            if(plast == pCurrent->pchild2){
                pCurrent->pchild2 = NULL;       /* Identify child2 as the node to be deleted and set it to NULL */
            }
            delete plast;                       /* Delete the node 'plast' */
            plast = NULL;
            nodecount--;                        /* Decrease nodecount */
        } else {
            /* Case where node to be deleted is the root node */
            delete pCurrent;
            pCurrent = NULL;
            pHead = NULL;
        }


        /* Cleanup temporary pointers to nodes. */
        psave = NULL;
        plast = NULL;
	}
}


template <class T>
void Tree2<T>::DeleteAllNodes()
{
	if(iscopy == false){
        /* Go to root node and then delete it, taking all child nodes with it */
        pCurrent = pHead;
        DeleteNode();

        /* Create new root node with no children */
        Initialize();
	}
}


template <class T>
T Tree2<T>::GetDataValue(int idx) const
{
    T retval = 0;
    switch(idx){
    case -1:
        retval = pCurrent->qtyparent;
        break;
    case 1:
        retval = pCurrent->qty1;
        break;
    case 2:
        retval = pCurrent->qty2;
        break;
    default:
        retval = nanval;
        break;
    }
    return retval;
}


template <class T>
void Tree2<T>::SetDataValue(const T val, int idx)
{
	if(iscopy == false){
        switch(idx){
        case -1:
            pCurrent->qtyparent = val;
            if(pCurrent != pHead && pCurrent == pCurrent->pparent->pchild1){
                pCurrent->pparent->qty1 = val;
            } else {
                pCurrent->pparent->qty2 = val;
            }
            break;
        case 1:
            pCurrent->qty1 = val;
            pCurrent->pchild1->qtyparent = val;
            break;
        case 2:
            pCurrent->qty2 = val;
            pCurrent->pchild2->qtyparent = val;
            break;
        default:
            break;
        }
	}
}


template <class T>
int Tree2<T>::GetNumNodes() const
{
	return nodecount;
}


template <class T>
typename Tree2<T>::node* Tree2<T>::GetPointerToHead()
{
	return pHead;
}


template <class T>
typename Tree2<T>::node* Tree2<T>::GetPointerToCurrent()
{
	return pCurrent;
}


template <class T>
void Tree2<T>::TraverseTree(const int id, const bool output)
{
    std::string newickstr;
    int nodestraversed = 0;
    int nodeid = std::numeric_limits<int>::quiet_NaN();

    /* Start at root node */
    pCurrent = pHead;

    while(nodeid != id && nodestraversed < nodecount){
        pCurrent->nvisits += 1;
        nodeid = pCurrent->nodeid;
        if(nodeid == id){ break; }

        if(pCurrent->nvisits == pCurrent->nchildren-1 &&
           pCurrent->pchild1 != NULL){
            if(pCurrent->pchild1->nvisits != pCurrent->pchild1->nchildren+1){
                newickstr += "(";
                GoToChildNode(1);
            }
        }
        if(pCurrent->nvisits == pCurrent->nchildren &&
           pCurrent->pchild1 != NULL){
            if(pCurrent->pchild1->nvisits != pCurrent->pchild1->nchildren+1){
                newickstr += "(";
                GoToChildNode(1);
            }
        }
        if(pCurrent->nvisits == pCurrent->nchildren &&
           pCurrent->pchild2 != NULL){
            if(pCurrent->pchild2->nvisits != pCurrent->pchild2->nchildren+1){
                newickstr += ",";
                GoToChildNode(2);
            }
        }
        if(pCurrent->nvisits == pCurrent->nchildren+1 &&
          (pCurrent->pchild1 != NULL || pCurrent->pchild2 != NULL)){
            if(pCurrent != pHead){
                newickstr += ")S" + NumToStr(pCurrent->nodeid) + ":" + NumToStr(pCurrent->qtyparent);
                GoToParentNode();
                nodeid = pCurrent->nodeid;
                nodestraversed += 1;
            } else {
                newickstr += ")S" + NumToStr(pCurrent->nodeid);
                nodestraversed += 1;
            }
        }
        if(pCurrent->nvisits == pCurrent->nchildren+1 &&
          (pCurrent->pchild1 == NULL && pCurrent->pchild2 == NULL)){
            if(pCurrent != pHead){
                newickstr += "S" + NumToStr(pCurrent->nodeid) + ":" + NumToStr(pCurrent->qtyparent);
                GoToParentNode();
                nodeid = pCurrent->nodeid;
                nodestraversed += 1;
            } else {
                newickstr += "S" + NumToStr(pCurrent->nodeid);
                nodestraversed += 1;
            }
        }

        nodestraversed = nodestraversed;

    }

    if(output == true){
        std::cout << "Newick String: " << newickstr << std::endl;
    }

    /* Set current node to root if specified ID not found, otherwise save ID */
    node *psave = NULL;
    if(nodeid != id){
        pCurrent = pHead;
    } else {
        psave = pCurrent;
    }

    /* Reset visit counter */
    ResetVisitCounter();

    /* Set current pointer to saved location */
    pCurrent = psave;
    psave = NULL;
}


template <class T>
T Tree2<T>::CalculateLikelihood(const Array1D<T> &params, const Array1D<char> &observed_vals,
                         const Array1D<char> &possible_vals, T(*f)(const char, const char, const T, const Array1D<T>&))
{
    T retval = 0.0e0;
    int nodestraversed = 0;

    int nodeid = std::numeric_limits<int>::quiet_NaN();

    int npossible = possible_vals.GetDim();
    Array2D<T> B(nodecount,npossible,0.0e0);

    /* Start at root node */
    pCurrent = pHead;

    while(nodestraversed < nodecount){
        pCurrent->nvisits += 1;
        nodeid = pCurrent->nodeid;

        if(pCurrent->nvisits == pCurrent->nchildren-1 &&
           pCurrent->pchild1 != NULL){
            if(pCurrent->pchild1->nvisits != pCurrent->pchild1->nchildren+1){
                //newickstr += "(";
                GoToChildNode(1);
            }
        }
        if(pCurrent->nvisits == pCurrent->nchildren &&
           pCurrent->pchild1 != NULL){
            if(pCurrent->pchild1->nvisits != pCurrent->pchild1->nchildren+1){
                //newickstr += "(";
                GoToChildNode(1);
            }
        }
        if(pCurrent->nvisits == pCurrent->nchildren &&
           pCurrent->pchild2 != NULL){
            if(pCurrent->pchild2->nvisits != pCurrent->pchild2->nchildren+1){
                //newickstr += ",";
                GoToChildNode(2);
            }
        }
        if(pCurrent->nvisits == pCurrent->nchildren+1/* &&
          (pCurrent->pchild1 != NULL || pCurrent->pchild2 != NULL)*/){
            /*
                Check type of node.  Note that this DOES NOT account for a node with a single
                child!
            */
            int i = pCurrent->nodeid;
            if(pCurrent->nchildren == 0){
                /* Calculate B(:,:) values for leaf */
                for(int j=0; j<npossible; j++){
                    char sp = possible_vals(j);
                    char si = observed_vals(i-1);
                    T t = pCurrent->qtyparent;
                    T tmp = f(sp,si,t,params);
                    //B(i,j) = f(sp,si,t,params);
                    B(i,j) = tmp;
                }
            }
            if(pCurrent->nchildren == 2 && pCurrent != pHead){
                /* Calculate B(:,:) values for non-root branch */
                for(int j=0; j<npossible; j++){
                    char sp = possible_vals(j);
                    T t = pCurrent->qtyparent;
                    int c1 = pCurrent->pchild1->nodeid;
                    int c2 = pCurrent->pchild2->nodeid;
                    for(int k=0; k<npossible; k++){
                        char si = possible_vals(k);
                        T tmp = f(sp,si,t,params)*B(c1,k)*B(c2,k);
                        //B(i,j) += f(sp,si,t,params)*B(c1,k)*B(c2,k);
                        B(i,j) += tmp;
                    }
                }
            }
            if(pCurrent->nchildren == 2 && pCurrent == pHead){
                /* Calculate total likelihood */
                for(int j=0; j<npossible; j++){
                    int c1 = pCurrent->pchild1->nodeid;
                    int c2 = pCurrent->pchild2->nodeid;
                    T tmp = params(i)*B(c1,j)*B(c2,j);
                    //retval += params(i)*B(c1,j)*B(c2,j);
                    retval += tmp;
                }
            }
            GoToParentNode();
            nodestraversed += 1;
        }
//        if(pCurrent->nvisits == pCurrent->nchildren+1 &&
//          (pCurrent->pchild1 == NULL && pCurrent->pchild2 == NULL)){
//            /*
//                Check type of node.  Note that this DOES NOT account for a node with a single
//                child!
//            */
//            int i = pCurrent->nodeid;
//            if(pCurrent->nchildren == 0){
//                /* Calculate B(:,:) values for leaf */
//                for(int j=0; j<npossible; j++){
//                    char sp = possible_vals(j);
//                    char si = observed_vals(i);
//                    T t = pCurrent->qtyparent;
//                    T tmp = f(sp,si,t,params);
//                    //B(i,j) = f(sp,si,t,params);
//                    B(i,j) = tmp;
//                }
//            }
//            if(pCurrent->nchildren == 2 && pCurrent != pHead){
//                /* Calculate B(:,:) values for non-root branch */
//                for(int j=0; j<npossible; j++){
//                    char sp = possible_vals(j);
//                    T t = pCurrent->qtyparent;
//                    int c1 = pCurrent->pchild1->nodeid;
//                    int c2 = pCurrent->pchild2->nodeid;
//                    for(int k=0; k<npossible; k++){
//                        char si = observed_vals(i);
//                        T tmp = f(sp,si,t,params)*B(c1,k)*B(c2,k);
//                        //B(i,j) += f(sp,si,t,params)*B(c1,k)*B(c2,k);
//                        B(i,j) += tmp;
//                    }
//                }
//            }
//            if(pCurrent->nchildren == 2 && pCurrent == pHead){
//                /* Calculate total likelihood */
//                for(int j=0; j<npossible; j++){
//                    int c1 = pCurrent->pchild1->nodeid;
//                    int c2 = pCurrent->pchild2->nodeid;
//                    T tmp = params(i)*B(c1,j)*B(c2,j);
//                    //retval += params(i)*B(c1,j)*B(c2,j);
//                    retval += tmp;
//                }
//            }
//            GoToParentNode();
//            nodestraversed += 1;
//        }
    }

    /* Reset visit counter */
    ResetVisitCounter();

    return retval;
}


template <class T>
bool Tree2<T>::IDIsUnique(const int id) const
{
    bool retval = true;
    for(int i=0; i<nodecount; i++){
        if(id == usedids(i)){
            retval = false;
        }
    }
    return retval;
}


template <class T>
void Tree2<T>::ResetVisitCounter()
{
    int nodestraversed = 0;

    /* Start at root node */
    pCurrent = pHead;

    while(nodestraversed < nodecount){
        while((pCurrent->pchild1 != NULL && pCurrent->pchild1->nvisits != 0) ||
              (pCurrent->pchild2 != NULL && pCurrent->pchild2->nvisits != 0)){

            if(pCurrent->pchild1 != NULL && pCurrent->pchild1->nvisits != 0){
                GoToChildNode(1);
            } else {
                if(pCurrent->pchild2 != NULL && pCurrent->pchild2->nvisits != 0){
                    GoToChildNode(2);
                }
            }

        }

        pCurrent->nvisits = 0;
        GoToParentNode();
        nodestraversed++;
    }
}
