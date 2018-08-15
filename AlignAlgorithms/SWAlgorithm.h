//This is the implementation of smith-waterman algorithm
#pragma once

#include"AlignAlgorithmBase.h"
#include"Sequence.h"

namespace SeqAnsis
{

#define	SW_MAT_TYPE   PAM250

class CSWAlgorithm:public CAlignAlgorithmBase
{
public:
	virtual void    Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences );
};

}