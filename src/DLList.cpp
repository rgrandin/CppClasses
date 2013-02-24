/**
 * @file DLList.cpp
 * @author Robert Grandin
 * @brief Implementation of DLList class.
 */



#include <DLList.h>     /* Included for syntax highlighting in Qt Creator */

/*
 * CONSTRUCTORS
 */

template <class T, int I>
DLList<T,I>::DLList()
{
	// IDENTIFY THIS INSTANCE AS *NOT* A COPY
	iscopy = false;

	// CREATE FIRST NODE
	pHead = new node;

	// IDENTIFY NODE AS BOTH THE TAIL AND CURRENT LOCATIONS
	pCurrent = pHead;

	// SET pNext TO NULL SINCE THIS IS THE TAIL
	pCurrent->pNext = NULL;
        pCurrent->pPrev = NULL;

	nodecount = 1;
}


template <class T, int I>
DLList<T,I>::DLList(const DLList<T,I> &list) :
    pHead(list.pHead), pCurrent(list.pCurrent), nodecount(list.nodecount)
{
    /* Identify this instance as a copy. */
	iscopy = true;
}


#ifdef CXX11
template <class T, int I>
DLList<T,I>::DLList(DLList<T, I> &&a) : DLList<T,I>()
{
    DLListSwap(*this, a);
}
#endif



/*
 * DESTRUCTOR
 */
template <class T, int I>
DLList<T,I>::~DLList()
{
	/*
	 * ONLY PERFORM THE DESTRUCTION IF THIS INSTANCE OF THE LIST IS *NOT* A COPY.
	 */
	if(iscopy == false){
		// DELETE POINTERS TO EACH NODE, BEGINNING WITH HEAD NODE
		node* pTemp = pHead;
		pCurrent = pHead;

		while(pCurrent != NULL){
			pCurrent = pCurrent->pNext;
			delete pTemp;
			pTemp = pCurrent;
		}
	}
}




/*
 * PUBLIC MEMBER FUNCTIONS
 */

template <class T, int I>
DLList<T,I>& DLList<T,I>::operator=(DLList<T,I> sourcelist)
{
    DLListSwap(*this, sourcelist);

    /* Identify this instance as a copy. */
	iscopy = true;

	return *this;
}


template <class T, int I>
void DLList<T,I>::AddNode()
{
	if(iscopy == false){
        node *pTemp;
		GoToFinalNode();
		pCurrent->pNext = new node;
        pTemp = pCurrent;
		Advance();
		pCurrent->pNext = NULL;
        pCurrent->pPrev = pTemp;

		nodecount++;
	}
}


template <class T, int I>
void DLList<T,I>::AddNodeAfter()
{
    if(iscopy == false){
        node *rightneighbor;
        node *leftneighbor;
        node *middleneighbor;

        rightneighbor = pCurrent->pNext;
        leftneighbor = pCurrent;
        pCurrent->pNext = new node;

        Advance();
        middleneighbor = pCurrent;
        pCurrent->pPrev = leftneighbor;
        pCurrent->pNext = rightneighbor;

        if(leftneighbor != NULL){
            leftneighbor->pNext = middleneighbor;
        }
        if(rightneighbor != NULL){
            rightneighbor->pPrev = middleneighbor;
        }

        nodecount++;
    }
}


template <class T, int I>
void DLList<T,I>::AddNodeBefore()
{
    if(iscopy == false){
        node *rightneighbor;
        node *leftneighbor;
        node *middleneighbor;

        rightneighbor = pCurrent;
        leftneighbor = pCurrent->pPrev;
        pCurrent->pPrev = new node;

        Rewind();
        middleneighbor = pCurrent;
        pCurrent->pPrev = leftneighbor;
        pCurrent->pNext = rightneighbor;

        if(leftneighbor != NULL){
            leftneighbor->pNext = middleneighbor;
        }
        if(rightneighbor != NULL){
            rightneighbor->pPrev = middleneighbor;
        }

        if(leftneighbor == NULL){
            pHead = pCurrent;
        }

        nodecount++;
    }
}


template <class T, int I>
void DLList<T,I>::Advance()
{
	if(pCurrent->pNext != NULL){
		pCurrent = pCurrent->pNext;
	}
}


template <class T, int I>
void DLList<T,I>::Rewind()
{
    if(pCurrent->pPrev != NULL){
        pCurrent = pCurrent->pPrev;
    }
}


template <class T, int I>
void DLList<T,I>::DeleteNode()
{
	if(iscopy == false){
		node* pTemp;
        node* pTemp2;

		// CHECK IF NODE TO BE DELETED IS HEAD NODE
		if(pCurrent == pHead)
		{
			pTemp = pHead;
			pHead = pHead->pNext;
			delete pTemp;
		}

		// CHECK IF NODE TO BE DELETED IS TAIL NODE
        if(pCurrent->pNext == NULL)
		{
			pTemp = pCurrent;
            pTemp2 = pCurrent->pPrev->pPrev;
			Rewind();
			pCurrent->pNext = NULL;
            pCurrent->pPrev = pTemp2;
			delete pTemp;
		}

		// IF PREVIOUS CONDITIONS FAIL, NODE IS BETWEEN ENDPOINTS OF LIST
        if(pCurrent != pHead && pCurrent->pNext != NULL)
		{
			pTemp = pCurrent;
            pTemp2 = pCurrent->pPrev->pPrev;
			Rewind();
			pCurrent->pNext = pTemp->pNext;
            pCurrent->pPrev = pTemp2;
			delete pTemp;
		}

		// SET LOCATION TO FIRST NODE IN THE LIST
		pCurrent = pHead;

		nodecount--;
	}
}


template <class T, int I>
void DLList<T,I>::DeleteAllNodes()
{
	if(iscopy == false){
		// DELETE POINTERS TO EACH NODE, BEGINNING WITH HEAD NODE
		node* pTemp = pHead;
		pCurrent = pHead;

		while(pCurrent != NULL){
			pCurrent = pCurrent->pNext;
			delete pTemp;
			pTemp = pCurrent;
		}

		// CREATE FIRST NODE
		pHead = new node;

		// IDENTIFY NODE AS BOTH THE TAIL AND CURRENT LOCATIONS
		pCurrent = pHead;

		// SET pNext TO NULL SINCE THIS IS THE TAIL
		pCurrent->pNext = NULL;

		// INITIALIZE DATA TO 0
                //for(int i=0; i<I; i++){
                //	pCurrent->value[i] = (T)0.0e0;
                //}

		nodecount = 1;
	}
}


template <class T, int I>
T DLList<T,I>::GetDataValue(int idx) const
{
	return pCurrent->value[idx];
}


template <class T, int I>
void DLList<T,I>::SetDataValue(const T val, int idx)
{
	if(iscopy == false){
		pCurrent->value[idx] = val;
	}
}


template <class T, int I>
int DLList<T,I>::GetNumNodes() const
{
	return nodecount;
}


template <class T, int I>
void DLList<T,I>::PrintNodes(int wraplim=10, bool blankline=false)
{
	pCurrent = pHead;
	int currpos = 1;
	int wrapcount = 0;

	// PRINT ALL NODES PRIOR TO FINAL NODE
        while(pCurrent->pNext != NULL){
		wrapcount = 0;
		std::cout << "Node " << currpos << ": ";
		for(int i=0; i<I; i++){
			if(wrapcount < wraplim-1){
				std::cout << pCurrent->value[i] << " ";
				wrapcount++;
			} else {
				std::cout << pCurrent->value[i] << std::endl;
				std::cout << "          ";
				wrapcount = 0;
			}
		}
		if(blankline){
			std::cout << std::endl;
		}

		std::cout << std::endl;
		currpos++;
		pCurrent = pCurrent->pNext;
	}

	// PRINT FINAL NODE

	wrapcount = 0;
	std::cout << "Node " << currpos << ": ";
	for(int i=0; i<I; i++){
		if(wrapcount < wraplim-1){
			std::cout << pCurrent->value[i] << " ";
			wrapcount++;
		} else {
			std::cout << pCurrent->value[i] << std::endl;
			std::cout << "          ";
			wrapcount = 0;
		}
        }
	std::cout << std::endl;
}


template <class T, int I>
void DLList<T,I>::GoToNode(const int nodenum)
{
	// GO TO FIRST NODE
	pCurrent = pHead;

	// ADVANCE THROUGH LIST UNTIL DESIRED NODE IS REACHED
	int count = 1;
	while(pCurrent->pNext != NULL && count != nodenum){
		Advance();
		count++;
	}
}


template <class T, int I>
void DLList<T,I>::GoToFinalNode()
{
	pCurrent = pHead;
	while(pCurrent->pNext != NULL){
		Advance();
	}
}


template <class T, int I>
typename DLList<T,I>::node* DLList<T,I>::GetPointerToHead()
{
	return pHead;
}


template <class T, int I>
typename DLList<T,I>::node* DLList<T,I>::GetPointerToCurrent()
{
	return pCurrent;
}


template <class T, int I>
bool DLList<T,I>::IsHeadNode()
{
    if(pCurrent == pHead){
        return true;
    } else {
        return false;
    }
}
