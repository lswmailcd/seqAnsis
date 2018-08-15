//This is the implementation of blast algorithm
#pragma once

#include <memory>
#include"AlignAlgorithmBase.h"
#include"Sequence.h"

namespace SeqAnsis
{

const int WORD_LENGTH = 4;
const int WORD_EXTEND_SCORE_DISTANCE = 2;
const int MSP_SCORE_CUTOFF = 17;//not used, it will be used later

class CHit
{
public:
	int idxWord1;
	int idxWord2;
	CHit(){}
	CHit( const CHit& vHit )
	{
		idxWord1 = vHit.idxWord1;
		idxWord2 = vHit.idxWord2;
	}
	CHit& operator=( const CHit& vRight )
	{
		if(this==&vRight)		return *this;

		this->idxWord1 = vRight.idxWord1;
		this->idxWord2 = vRight.idxWord2;
		return *this;
	}
};

class CMSP
{
public:
	int idxStart1;
	int idxEnd1;
	int idxStart2;
	int idxEnd2;
	int len;
	int score;
	CMSP(){}
	CMSP( const CMSP& vHit )
	{
		idxStart1 = vHit.idxStart1;
		idxStart2 = vHit.idxStart2;
		idxEnd1 = vHit.idxEnd1;
		idxEnd2 = vHit.idxEnd2;		
		len = vHit.len;		
		score = vHit.score;
	}
	CMSP& operator=( const CMSP& vRight )
	{
		if(this==&vRight)		return *this;

		this->idxStart1 = vRight.idxStart1;
		this->idxStart2 = vRight.idxStart2;
		this->idxEnd1 = vRight.idxEnd1;
		this->idxEnd2 = vRight.idxEnd2;		
		this->len = vRight.len;		
		this->score = vRight.score;
		return *this;
	}
};

typedef struct 
{
	char *pReSeq1;
	char *pReSeq2;
	int   len;
	float score;
}StruReSeqPair;

class CBlastAlgorithm:public CAlignAlgorithmBase
{
public:
	virtual void    Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences );
private:
	int		CompareWordPair( const std::string& vWord1, const std::string& vWord2 );
	int		CompareResidue( const char& vCh1, const char& vCh2 );
	std::vector<std::shared_ptr<CHit>>		m_Hits;
	int         m_iMSPCount;
	std::vector<std::shared_ptr<CMSP>>		m_MSP;	
};

}