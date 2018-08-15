#include "GeneDecoder.h"

namespace	SeqAnsis
{

CGeneDecoder::CGeneDecoder():m_pPath(NULL),m_iPathLength(0)
{
}

CGeneDecoder::~CGeneDecoder()
{
	//destroy the path link list
	LinkList *p(NULL), *q(NULL); 
	p = m_pPath->pNext;
	SAFE_DELETE(m_pPath);
	while( NULL != p )
	{
		q = p;
		p = p->pElem;
		SAFE_DELETE(q);
	}
}

void CGeneDecoder::Decode( const char* pGene, int geneLen )
{
	int iValidGeneLen = CheckValidGeneLen( pGene, geneLen );

	if ( NULL != m_pPath )  
	{
		//destroy the path link list
		LinkList *p(NULL), *q(NULL); 
		p = m_pPath->pNext;
		SAFE_DELETE(m_pPath);
		while( NULL != p )
		{
			q = p;
			p = p->pElem;
			SAFE_DELETE(q);
		}
	}

	m_pPath = new LinkList;//the first elem of Path has nothing.
	LinkList *pTail = m_pPath;
	LinkList *p=m_pPath, *q=NULL;
	for( int j=0; j<iValidGeneLen; ++j )
	{
		q = new LinkList;
		pTail = q;
		q->elem = pGene[j];
		q->pPrior = p;
		p->pNext = q;
		p = q;
	}

	m_iPathLength = iValidGeneLen;
	if ( m_pPath != pTail )//the list is NOT empty
	{
		LinkList *pFuncPointer = pTail;
		while( pFuncPointer != m_pPath )
		{
			while( 'A' != pFuncPointer->elem && pFuncPointer != m_pPath )  pFuncPointer = pFuncPointer->pPrior;
			if ( pFuncPointer != m_pPath )//found one function character
			{
				//calculate the path length;
				--m_iPathLength;

				//do the operation of A function
				LinkList *pOpnd2 = pTail;					
				LinkList *pOpnd1 = pTail->pPrior;
				q = pOpnd1;
				p = pOpnd1->pElem;
				while( p )  
				{
					q = p;
					p = p->pElem;
				}
				q->pElem = pOpnd2;//add Opnd2 to the tail of Opnd1 list

				//move the tail pointer to the tail of list after the calculation
				pTail = pOpnd1->pPrior;
				//if the tail just the position to be replaced from function character to the result string of the calculation
				if( pFuncPointer == pTail )  pTail = pOpnd1;
				//replace the function character with the result string of the calculation.
				p = pFuncPointer->pPrior;	
				if(pFuncPointer->pNext == pOpnd1)
				{
					p->pNext = pOpnd1;
					pOpnd1->pPrior = p;
				}
				else
				{										
					q = pFuncPointer->pNext;
					p->pNext = pOpnd1;
					pOpnd1->pPrior = p;
					q->pPrior = pOpnd1;
					pOpnd1->pNext = q;
				}

				//prepare for the next round of function search
				q = pFuncPointer;
				pFuncPointer = pFuncPointer->pPrior;	

				//release the space of function character have been used for string calculation.
				SAFE_DELETE(q);
			}
		}
	}
}

int CGeneDecoder::CheckValidGeneLen(const char *pGene, int geneLen)
{
	int Nt=0, Nf=0;
	for ( int i=0; i<geneLen; ++i )
	{
		switch( pGene[i] )
		{
		case 'A'://function
			{
				++Nf;
			}
			break;
		case 'R':
		case 'B':
		case 'D':
			{
				++Nt;
			}
			break;
		}
		if ( Nt == (Nf+1) )//The parameter number of 'A' is just 2. 
		{
			break;
		}
	}
	return Nt+Nf;
}

bool    CGeneDecoder::GetPath( std::string& voPath )
{
	LinkList *p = m_pPath->pNext;
	for( int i=0; i<m_iPathLength; ++i )
	{
		voPath.push_back( p->elem );
		p = p->pElem;
	}
	return true;
}

}