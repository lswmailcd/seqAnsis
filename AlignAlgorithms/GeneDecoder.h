//GeneDecoder is used to translate the linear gene code to link list expressed linear path.
#pragma once

#include <string>
#include "Sequence.h"

namespace	SeqAnsis
{

class CGeneDecoder
{
public:
	CGeneDecoder();
	~CGeneDecoder();
	void	Decode( const char* pGene, int geneLen );
	int		GetPathLength() { return m_iPathLength; }
	bool    GetPath( std::string& voPath );
private:
	int		CheckValidGeneLen(const char *pGene, int geneLen);
private:
	struct LinkList
	{
		char			elem;
		struct LinkList	*pNext;
		struct LinkList *pPrior;
		struct LinkList *pElem;
		LinkList::LinkList()
		{
			elem = GENESPACE;
			pNext = pPrior = pElem = NULL;
		}
	};

	int			m_iPathLength;
	LinkList	*m_pPath;
};

}