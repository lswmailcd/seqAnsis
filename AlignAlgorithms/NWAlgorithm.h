//This is the implementation of needleman_wunsch algorithm
#pragma once

#include"AlignAlgorithmBase.h"
#include"Sequence.h"

namespace SeqAnsis
{

#define	NW_MAT_TYPE   BLOSUM62

class CNWAlgorithm:public CAlignAlgorithmBase
{
public:
	virtual void    Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences );
};

}