#include<hash_map>
#include"Common.h"
#include "NWAlgorithm.h"
#include"BlastAlgorithm.h"
#include"GlobalSpace.h"


namespace SeqAnsis
{

typedef std::hash_map<std::string, int>	 CStringNumPair;

void	CBlastAlgorithm::Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences )
{
	if( vSequences.size()!=2 )// || vAlignedSequences.size()!=2 )
	{
		CGlobalSpace::m_sEventLog.writeEvent("Wrong sequences size to align by blast algorithm!");
		return;
	}

	int nSeq1 = vSequences[0].getSequenceContext().size();
	int nSeq2 = vSequences[1].getSequenceContext().size();
	const CSeqData& Seq1 = vSequences[0].getSequenceContext();
	const CSeqData& Seq2 = vSequences[1].getSequenceContext();

//Step1: Compiling a list of high-scoring words
	if ( nSeq1<WORD_LENGTH || nSeq2<WORD_LENGTH )
	{
		CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "Sequence 1 or/and 2 is/are shorter than the WORD_LENGTH!" );
	    return;
	}

	CStringNumPair	WordList;
	//generate the hash table for the word list of sequence1
	for ( int i=0; i<(nSeq1-WORD_LENGTH+1); ++i )
	{
		std::string str;
		for(int j=0; j<WORD_LENGTH; ++j)
		{
			str.push_back( Seq1[i+j].m_char );
		}
		WordList.insert(std::make_pair(str, i));
	}

//Step2: Scanning the database for hits
	bool bFoundHit = false;
	for ( int i=0; i<nSeq2-WORD_LENGTH+1; ++i )
	{
		std::string str;
		for(int j=0; j<WORD_LENGTH; ++j)
		{
			str.push_back( Seq2[i+j].m_char );
		}
		if( WordList.find(str)!=WordList.end() )
		{//found a hit
			std::auto_ptr<CHit> pHit(new CHit);
			pHit->idxWord1 = WordList[str];
			pHit->idxWord2 = i;
			m_Hits.push_back(pHit);
			bFoundHit = true;
		}
	}
	if ( !bFoundHit )
	{
		CGlobalSpace::m_sEventLog.writeEvent("No hit found in database!");
		return;
	}

	//Step3: Extending hits
	for ( int i=0; i<(int)(m_Hits.size()); ++i )
	{
		//extending left direction, extend one residue pair of the hit 
		//until the score is <highScore-Distance or arrive the left boundary
		std::auto_ptr<CMSP> pMsp(new CMSP);
		pMsp->idxStart1 = m_Hits[i]->idxWord1;
		pMsp->idxStart2 = m_Hits[i]->idxWord2;
		pMsp->len = WORD_LENGTH;
		int maxScore = 0;
		//extending toward left
		while ( (pMsp->idxStart1-1)>=0 && (pMsp->idxStart2-1)>=0 )
		{
			int nLen = m_Hits[i]->idxWord1-(pMsp->idxStart1-1)+WORD_LENGTH;
			std::string str1, str2;
			for(int j=0; j<nLen; ++j)
			{
				str1.push_back( Seq1[pMsp->idxStart1-1+j].m_char );
				str2.push_back( Seq2[pMsp->idxStart2-1+j].m_char );
			}
			int score = CompareWordPair( str1, str2 );

			if ( score>maxScore )
			{
				maxScore = score;
			}
			else
			{
				if( (maxScore-score)>WORD_EXTEND_SCORE_DISTANCE ) 
				{
					break;
				}
			}
			--(pMsp->idxStart1);
			--(pMsp->idxStart2);
			++(pMsp->len);
		}
		//extending right direction, extend one residue pair of the hit 
		//until the score is <highScore-Distance or arrive the right boundary
		pMsp->idxEnd1 = m_Hits[i]->idxWord1+WORD_LENGTH-1;
		pMsp->idxEnd2 = m_Hits[i]->idxWord2+WORD_LENGTH-1;
		maxScore = 0;
		while ( (pMsp->idxEnd1+1)<(int)(Seq1.size()) && (pMsp->idxEnd2+1)<(int)(Seq2.size()) )
		{
			int nLen = pMsp->idxEnd1+1-(m_Hits[i]->idxWord1)+1;
			std::string str1, str2;
			for(int j=0; j<nLen; ++j)
			{
				str1.push_back( Seq1[m_Hits[i]->idxWord1+j].m_char );
				str2.push_back( Seq2[m_Hits[i]->idxWord2+j].m_char );
			}
			int score = CompareWordPair( str1, str2 );


			if ( score>maxScore )
			{
				maxScore = score;
			}
			else
			{
				if( (maxScore-score)>WORD_EXTEND_SCORE_DISTANCE ) 
				{
					break;
				}
			}
			++(pMsp->idxEnd1);
			++(pMsp->idxEnd2);
			++(pMsp->len);
		}
		m_MSP.push_back(pMsp);
	}

	//use SW algorithm to find the final match of the MSPs.
	CNWAlgorithm	NWAlgorithm;
	std::vector<CSequence> AlignedSeqsSW;
	for ( int i=0; i<(int)(m_MSP.size()); ++i )
	{
		std::vector<CSequence> MSPSeqs;
		CStringNumPair sn;
		CSeqData data1,data2;
		for( int j=0; j<m_MSP[i]->len; ++j )
		{	
			StruSeqElem elem1( Seq1[m_MSP[i]->idxStart1+j].m_char, m_MSP[i]->idxStart1+j );
			data1.push_back(elem1);
				
			StruSeqElem elem2( Seq2[m_MSP[i]->idxStart2+j].m_char, m_MSP[i]->idxStart2+j );
			data2.push_back(elem2);
		}
		CSequence  MSPSeq1(data1, vSequences[0].getName(), "");
		CSequence  MSPSeq2(data2, vSequences[1].getName(), "");
		MSPSeqs.push_back(MSPSeq1);
		MSPSeqs.push_back(MSPSeq2);

		AlignedSeqsSW.clear();
		NWAlgorithm.Align( AlignedSeqsSW, MSPSeqs );

		vAlignedSequences.push_back( AlignedSeqsSW[0] );
		vAlignedSequences.push_back( AlignedSeqsSW[1] );
	}
}

int		CBlastAlgorithm::CompareResidue( const char& vCh1, const char& vCh2 )
{
	int score = 0;
	if (vCh1==vCh2 )
	{
		score = 5;
	}
	else
	{
		score = -4;
	}
	//else if( r1=='A'&&r2=='G' || r1=='G'&&r2=='A' 
	//	||r1=='C'&&r2=='T' || r1=='T'&&r2=='C' )
	//{
	//	score = -5;
	//}
	//else
	//{
	//	score = -7;
	//}
	return score;
}

int		CBlastAlgorithm::CompareWordPair( const std::string& vWord1, const std::string& vWord2 )
{
	if( vWord1.size()!=vWord2.size() )
	{
		throw	CAppException( DEF_EXCEPTION_INVALID_PARAMETER, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "word length is not equal!" );
	}

	int score = 0;
	for ( int i=0; i<(int)(vWord1.size()); ++i )
	{
		score += CompareResidue( vWord1[i], vWord2[i] );
	}

	return score;
}

}