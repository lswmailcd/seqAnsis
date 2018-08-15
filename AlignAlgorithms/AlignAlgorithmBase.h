#pragma once

#include<vector>
#include"Sequence.h"

namespace SeqAnsis
{

#define POS(r,c,w)  (r)*(w)+(c)		

class CAlignAlgorithmBase
{
public:
	virtual void    Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& sequences )=0;
};

}