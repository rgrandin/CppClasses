/**
 * @file SLList.cpp
 * @author Robert Grandin
 * @brief Implementation of SLList class.
 */


#include <SLList.h>



/*
 * CONSTRUCTORS
 */

template <class T, int I>
SLList<T,I>::SLList()
{
	// IDENTIFY THIS INSTANCE AS *NOT* A COPY
	iscopy = false;

	// CREATE FIRST NODE
	pHead = new node;

	// IDENTIFY NODE AS BOTH THE TAIL AND CURRENT LOCATIONS
	pCurrent = pHead;

	// SET pNext TO NULL SINCE THIS IS THE TAIL
	pCurrent->pNext = NULL;

	// INITIALIZE DATA TO 0
	for(int i=0; i<I; i++){
		pCurrent->value[i] = (T)0.0e0;
	}

	nodecount = 1;
}


template <class T, int I>
SLList<T, I>::SLList(const SLList<T,I> &list) :
    pHead(list.pHead), pCurrent(list.pCurrent), nodecount(list.nodecount)
{
    iscopy = true;
}


#ifdef CXX11
template <class T, int I>
SLList<T, I>::SLList(SLList<T,I> &&a) : SLList<T,I>()
{
    SLListSwap(*this, a);
}
#endif




/*
 * DESTRUCTOR
 */
template <class T, int I>
SLList<T,I>::~SLList()
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
SLList<T>& SLList<T,I>::operator=(SLList<T,I> a)
{
    SLListSwap(*this, a);
    iscopy = true;
    return *this;
}



template <class T, int I>
SLList<T,I>& SLList<T,I>::operator=(const SLList<T,I>& sourcelist)
{
	// IDENTIFY THIS INSTANCE AS A COPY
	iscopy = true;

	// APPLY POINTER VALUES TO LOCAL MEMBER VARIABLES
	pHead = sourcelist.GetPointerToHead();
	pCurrent = pHead;

	// NOTE: NO NODE-INITIALIZATION SINCE THAT WOULD OVERWRITE THE DATA STORED
	//		 IN list.

	// SET THE NODECOUNT
	nodecount = sourcelist.GetNumNodes();

	return *this;
}



template <class T, int I>
void SLList<T,I>::AddNode()
{
	if(iscopy == false){
		GoToFinalNode();
		pCurrent->pNext = new node;
		Advance();
		pCurrent->pNext = NULL;

		nodecount++;
	}
}


template <class T, int I>
void SLList<T,I>::Advance()
{
	if(pCurrent->pNext != NULL){
		pCurrent = pCurrent->pNext;
	}
}


template <class T, int I>
void SLList<T,I>::Rewind()
{
	if(pCurrent != pHead){
		node* pTemp = pHead;

		if(pCurrent == pHead)
		{
			pCurrent = pHead;
		}

		while(pTemp->pNext != pCurrent)
		{
			pTemp = pTemp->pNext;
		}
		pCurrent = pTemp;
	}
}


template <class T, int I>
void SLList<T,I>::DeleteNode()
{
	if(iscopy == false){
		node* pTemp;

		// CHECK IF NODE TO BE DELETED IS HEAD NODE
		if(pCurrent == pHead)
		{
			pTemp = pHead;
			pHead = pHead->pNext;
			delete pTemp;
		}

		// CHECK IF NODE TO BE DELETED IS TAIL NODE
		else if(pCurrent->pNext == NULL)
		{
			pTemp = pCurrent;
			Rewind();
			pCurrent->pNext = NULL;
			delete pTemp;
		}

		// IF PREVIOUS CONDITIONS FAIL, NODE IS BETWEEN ENDPOINTS OF LIST
		else
		{
			pTemp = pCurrent;
			Rewind();
			pCurrent->pNext = pTemp->pNext;
			delete pTemp;
		}

		// SET LOCATION TO FIRST NODE IN THE LIST
		pCurrent = pHead;

		nodecount--;
	}
}


template <class T, int I>
void SLList<T,I>::DeleteAllNodes()
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
		for(int i=0; i<I; i++){
			pCurrent->value[i] = (T)0.0e0;
		}

		nodecount = 1;
	}
}


template <class T, int I>
T SLList<T,I>::GetDataValue(int idx) const
{
	return pCurrent->value[idx];
}


template <class T, int I>
void SLList<T,I>::SetDataValue(const T val, int idx)
{
	if(iscopy == false){
		pCurrent->value[idx] = val;
	}
}


template <class T, int I>
int SLList<T,I>::GetNumNodes() const
{
	return nodecount;
}


template <class T, int I>
void SLList<T,I>::PrintNodes(int wraplim=10, bool blankline=false)
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
void SLList<T,I>::GoToNode(const int nodenum)
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
void SLList<T,I>::GoToFinalNode()
{
	pCurrent = pHead;
	while(pCurrent->pNext != NULL){
		Advance();
	}
}


template <class T, int I>
typename SLList<T,I>::node* SLList<T,I>::GetPointerToHead()
{
	return pHead;
}


template <class T, int I>
typename SLList<T,I>::node* SLList<T,I>::GetPointerToCurrent()
{
	return pCurrent;
}
