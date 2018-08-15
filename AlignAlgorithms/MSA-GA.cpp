#include <ctime>
#include"Common.h"
#include"MSA-GA.h"
#include"GlobalSpace.h"
#include"Timer.h"
#include "NWAlgorithm.h"

namespace SeqAnsis
{

int  _cdecl compScore( const void* s1, const void* s2 )
{
	return *(CCompSPS *)s2 - *(CCompSPS *)s1; 
}

CMSAGAAlgorithm::CMSAGAAlgorithm() :m_pPopulation(NULL), m_fMaxScore(-INT_MAX), m_iRun(0), MSA_GA_POPULATION_SIZE(100)
{

}

CMSAGAAlgorithm::~CMSAGAAlgorithm()
{
	SAFE_DELETE_ARRAY(m_pPopulation);
}

void		CMSAGAAlgorithm::dbgCheckSequences(const std::vector<CSequence>& vSequences)
{
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		for ( int j=0; j<m_pPopulation[i].nSeqSize; ++j )
		{
			const CSeqData& seq = m_pPopulation[i].pSequence[j].sequence.getSequenceContext();
			int nSeq=0;
			for( int k=0; k<seq.size(); ++k )
			{
				if ( seq[k].m_char!='-' )
				{
					++nSeq;
				}
			}
			if ( nSeq!=vSequences[j].getLen() )
			{
				char ch[128];
				sprintf_s(ch,  "The seq %i of the population %i is wrong! The evolution is %i" , j,  i, m_iRun);
				CGlobalSpace::m_sEventLog.writeEvent(ch);
			}
		}
	}
}

void CMSAGAAlgorithm::SetAlignParams(int nPopulationNum, const std::vector<float>& SeqWeight)
{
	MSA_GA_POPULATION_SIZE = nPopulationNum;

	m_SeqWeight.resize(SeqWeight.size());
#ifdef CLUSTALW_SP_SCORE
	for (int i = 0; i < (int)SeqWeight.size(); ++i)
	{
		m_SeqWeight[i] = SeqWeight[i];
	}
#else
	for (int i = 0; i < (int)SeqWeight.size(); ++i)
	{
		m_SeqWeight[i] = 1.0f;
	}
#endif
}

void CMSAGAAlgorithm::Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences )
{
	m_bdbgFirstRun = true;
	////find the longest sequence in the sequence vector
	int iLongestSeqSize = -INT_MAX;
	for ( int i=0; i<(int)vSequences.size(); ++i )
	{
		m_iLongestSeqSize = iLongestSeqSize<vSequences[i].getLen()?vSequences[i].getLen():iLongestSeqSize;

	}

	//define the column of the sequence matrix
	int nCol = (int)( m_iLongestSeqSize*(MSA_GA_SCALING_FACTOR_K+MSA_GA_OFFSET_FACTOR_X) );

	//population initialization
	InitPopulation2009( vSequences, nCol );

	//calculate the fitness of every organism. if not meet the requirement ,then evolution begin.
	char ch[100];
	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA: Population evolution begin!" );
	m_iRun = 0;
	bool  bStop=false;
	int iCount=0;
	float maxScore=INT_MIN;

	CTimer  time;
	double dStartTime = time.getCurrentTime();

	while(!Fitness2009()  && m_iRun<MSA_GA_GENERATION_NUM && iCount<MSA_GA_NO_ADV_GENERATION_NUM)
	{
		if( m_fMaxScore>maxScore )
		{
			iCount=0;
			maxScore = m_fMaxScore;
		}
		else
		{
			++iCount;
		}
		if ( m_bdbgFirstRun )
		{
			sprintf_s(ch, "generation=%i, maxscore=%i", m_iRun, maxScore);
			CGlobalSpace::m_sEventLog.writeEvent(ch);
			m_bdbgFirstRun = false;
		}
		++m_iRun;
		Evolution();
	}
	double dInterval = time.getCurrentTime() - dStartTime;
	sprintf_s(ch, "the process time of MSA-GA algorithm is %f", dInterval);
	SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent(ch);
	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA: Population evolution end!" );

	sprintf_s(ch, "generation=%i, maxscore=%f, iCount=%i", m_iRun, maxScore, iCount);
	CGlobalSpace::m_sEventLog.writeEvent(ch);

	//return the best organism as the aligned sequences
	int bestScore = 0;
	int idxBestOrganism = 0;
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		//sprintf_s(ch,  "seq No=%i, spScore=%i" , i, m_pPopulation[i].score );
		//CGlobalSpace::m_sEventLog.writeEvent(ch);

		if (  bestScore<m_pPopulation[i].score )
		{
			bestScore = m_pPopulation[i].score;
			idxBestOrganism = i;
		}
	}

	//if( !CheckSumCharNum( m_pPopulation[idxBestOrganism] ) )
	//{
	//	CGlobalSpace::m_sEventLog.writeEvent("CheckSumCharNum  ERROR!");
	//}

	for ( int i=0; i<m_pPopulation[idxBestOrganism].nSeqSize; ++i )
	{
		vAlignedSequences.push_back( m_pPopulation[idxBestOrganism].pSequence[i].sequence );
	}
}

bool		CMSAGAAlgorithm::CheckSumCharNum( const COrganism& vOG )
{
	for ( int i=0; i<vOG.nSeqSize; ++i )
	{
		if (  vOG.pSequence[i].charNum != CountSeqCharNum( vOG.pSequence[i].sequence ) )
		{
			return false;
		}
	}
	return true;
}

void		CMSAGAAlgorithm::Evolution()
{
	Selection();
	Mutation();	
	Recombination();	
}

void		CMSAGAAlgorithm::Recombination2004()
{
	VerticalRecombination();
}

void		CMSAGAAlgorithm::Mutation2004()
{
	GapInsertMutation2004();
	GapExtensionMutation();
	GapReductionMutation();
}

void		CMSAGAAlgorithm::GapInsertMutation2004()
{
	for ( int i=1; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < MSA_GA_BLOCK_MUTATION )
		{
			int iSeqNo = CGlobalSpace::m_sUtility.getRandomNumber(m_pPopulation[i].nSeqSize-1);
			//only one gap insert
			COrganism  org;
			org = m_pPopulation[i];

			int iBlockNum = 1;
			CSeqData& Seq = org.pSequence[iSeqNo].sequence.getSequenceContext();				
			int idxSelectedPos = CGlobalSpace::m_sUtility.getRandomNumber(org.pSequence[iSeqNo].sequence.getLen()-1);
			CSeqData::iterator&	itrSeqData = Seq.begin();
			for ( int j=0; j<idxSelectedPos; ++j )
			{
				++itrSeqData;
			}
			if ( itrSeqData==Seq.end() )
			{
				throw CAppException( DEF_EXCEPTION_INDEX_OUT_OF_RANGE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "false in the pointer!" );
			}
			//insert gaps
			for ( int j=0; j<iBlockNum; ++j )
			{
				itrSeqData = Seq.insert( itrSeqData, StruSeqElem(GENESPACE,-1) );
			}
			ResizeOrgan( org );
			if ( SPScore(org)>SPScore(m_pPopulation[i]) )
			{
				m_pPopulation[i] = org;
			}
		}
	}
}

void		CMSAGAAlgorithm::Mutation()
{
	GapInsertMutation();
	GapExtensionMutation();
	GapReductionMutation();
}

void		CMSAGAAlgorithm::GapInsertMutation()
{
	for ( int i=1; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < MSA_GA_BLOCK_MUTATION )
		{
				int iSeqNo = CGlobalSpace::m_sUtility.getRandomNumber(m_pPopulation[i].nSeqSize);
				//only one gap insert?
				int iBlockNum = CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_MAX_INSERT_GAP)+1;
				CSeqData& Seq = m_pPopulation[i].pSequence[iSeqNo].sequence.getSequenceContext();				
				int idxSelectedPos = CGlobalSpace::m_sUtility.getRandomNumber(m_pPopulation[i].pSequence[iSeqNo].sequence.getLen()-1);
				CSeqData::iterator&	itrSeqData = Seq.begin();
				for ( int j=0; j<idxSelectedPos; ++j )
				{
					++itrSeqData;
				}
				if ( itrSeqData==Seq.end() )
				{
					throw CAppException( DEF_EXCEPTION_INDEX_OUT_OF_RANGE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "false in the pointer!" );
				}
				//insert gaps
				for ( int j=0; j<iBlockNum; ++j )
				{
					itrSeqData = Seq.insert( itrSeqData, StruSeqElem(GENESPACE,-1) );
				}
				ResizeOrgan( m_pPopulation[i] );
			}
	}
}

int	    CMSAGAAlgorithm::SelectRandomGap( CSeqData& vSeqSelected )
{
	int idxPos = 0;
	std::vector<int>	VectIdxGap;
	//find all gap blocks in the sequence
	while ( idxPos<(int)vSeqSelected.size() )
	{
		while( idxPos<(int)vSeqSelected.size() && vSeqSelected[idxPos].m_iCode!=CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )	 ++idxPos;
		if ( idxPos<(int)vSeqSelected.size() )
		{
			VectIdxGap.push_back(idxPos);
			//find continuous gaps in current position.
			while( idxPos<(int)vSeqSelected.size() && vSeqSelected[idxPos].m_iCode==CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) ) ++idxPos;
		}
	}
	//select one gap block randomly and return the first gap position of this gap block.
	switch( (int)VectIdxGap.size() )
	{
	case 0://no gap in the sequence
		{
			idxPos = -1;
			break;
		}
	case 1:// only one gap, we use it directly.
		{
			idxPos = VectIdxGap[0];
			break;
		}
	default:
		{
			idxPos = VectIdxGap[CGlobalSpace::m_sUtility.getRandomNumber(VectIdxGap.size()-1)];
		}
	}

	return idxPos;
}

void		CMSAGAAlgorithm::GapExtensionMutation2004()
{
	for ( int i=1; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < MSA_GA_EXTENSION_MUTATION )
		{
			COrganism  org;
			org = m_pPopulation[i];
			int iSeqNo = CGlobalSpace::m_sUtility.getRandomNumber(org.nSeqSize-1);
			CSeqData& Seq = org.pSequence[iSeqNo].sequence.getSequenceContext();
			int iSelectedPos = SelectRandomGap(Seq);
			if ( iSelectedPos<0 )
			{// no gap, we can not extend it.
				return;
			}
			int j=0;
			std::vector<StruSeqElem>::iterator  itrSelectPos = Seq.begin();
			while( j<iSelectedPos )
			{
				++j;
				++itrSelectPos;
			}
			Seq.insert(itrSelectPos , StruSeqElem(GENESPACE,-1) );

			ResizeOrgan( org );
			if ( SPScore(org)>SPScore(m_pPopulation[i]) )
			{
				m_pPopulation[i] = org;
			}
		}
	}
}

void		CMSAGAAlgorithm::GapExtensionMutation()
{
	for ( int i=1; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < MSA_GA_EXTENSION_MUTATION )
		{
				int iSeqNo = CGlobalSpace::m_sUtility.getRandomNumber(m_pPopulation[i].nSeqSize);
				CSeqData& Seq = m_pPopulation[i].pSequence[iSeqNo].sequence.getSequenceContext();
				int iSelectedPos = SelectRandomGap(Seq);
				if ( iSelectedPos<0 )
				{// no gap, we can not extend it.
					return;
				}
				int j=0;
				std::vector<StruSeqElem>::iterator  itrSelectPos = Seq.begin();
				while( j<iSelectedPos )
				{
					++j;
					++itrSelectPos;
				}
				Seq.insert(itrSelectPos , StruSeqElem(GENESPACE,-1) );

				ResizeOrgan( m_pPopulation[i] );
		}
	}
}

void		CMSAGAAlgorithm::GapReductionMutation2004()
{
	for ( int i=1; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < MSA_GA_REDUCTION_MUTATION )
		{
			COrganism  org;
			org = m_pPopulation[i];
			int iSeqNo = CGlobalSpace::m_sUtility.getRandomNumber(org.nSeqSize-1);
			CSeqData& Seq = org.pSequence[iSeqNo].sequence.getSequenceContext();
			//select one gap block randomly and remove one gap of it
			int idxPos = SelectRandomGap(Seq);
			if ( idxPos<0 )
			{// no gap, we can not extend it.
				return;
			}
			int j=0;
			CSeqData::iterator itrSelectPos = Seq.begin();
			while( j<idxPos )
			{
				++j;
				++itrSelectPos;
			}

			int iGapBlockLen = 0;
			CSeqData::iterator	itrPos=itrSelectPos;
			while( Seq[idxPos].m_char == GENESPACE && itrPos!=Seq.end() )
			{
				++itrPos;
				++iGapBlockLen;
			}
			if ( CGlobalSpace::m_sUtility.getRandomNumber() < 1.0f/iGapBlockLen )
			{
				Seq.erase( itrSelectPos );	
			}
			ResizeOrgan( org );
			if ( SPScore(org)>SPScore(m_pPopulation[i]) )
			{
				m_pPopulation[i] = org;
			}
		}
	}
}

void		CMSAGAAlgorithm::GapReductionMutation()
{
	for ( int i=1; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < MSA_GA_REDUCTION_MUTATION )
		{
				int iSeqNo = CGlobalSpace::m_sUtility.getRandomNumber(m_pPopulation[i].nSeqSize);
				CSeqData& Seq = m_pPopulation[i].pSequence[iSeqNo].sequence.getSequenceContext();
				//select one gap block randomly and remove one gap of it
				int idxPos = SelectRandomGap(Seq);
				if ( idxPos<0 )
				{// no gap, we can not extend it.
					return;
				}
				int j=0;
				CSeqData::iterator itrSelectPos = Seq.begin();
				while( j<idxPos )
				{
					++j;
					++itrSelectPos;
				}

				int iGapBlockLen = 0;
				CSeqData::iterator	itrPos=itrSelectPos;
				while( Seq[idxPos].m_char == GENESPACE && itrPos!=Seq.end() )
				{
					++itrPos;
					++iGapBlockLen;
				}
				if ( CGlobalSpace::m_sUtility.getRandomNumber() < 1.0f/iGapBlockLen )
				{
					Seq.erase( itrSelectPos );	
				}
				ResizeOrgan( m_pPopulation[i] );
		}
	}
}

void		CMSAGAAlgorithm::Recombination2009()
{
	HorizentalRecombination2009();
	VerticalRecombination2009();
}

void		CMSAGAAlgorithm::Recombination()
{
	HorizentalRecombination();
	VerticalRecombination();
}

void		CMSAGAAlgorithm::HorizentalRecombination2009()
{
	COrganism  *pPopulation( new COrganism[MSA_GA_POPULATION_SIZE] );
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		pPopulation[i] =m_pPopulation[i];
	}
	for ( int i=1; i<MSA_GA_POPULATION_SIZE-1; i+=2 )
	{
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < MSA_GA_HORIZENTAL_RECOMB_RATIO )
		{
			//randomly choose the organism to recombine with.
			int iOrgan1 = CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_POPULATION_SIZE-1);
			int iOrgan2 = CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_POPULATION_SIZE-1);

			if( iOrgan1==iOrgan2) 
			{
				m_pPopulation[i] = pPopulation[iOrgan1];
				m_pPopulation[i+1] = pPopulation[iOrgan2];
			}
			else
			{
				//randomly choose the sequence to recombine with.
				int iIndex =  CGlobalSpace::m_sUtility.getRandomNumber(pPopulation[iOrgan1].nSeqSize-1);

				//horizontal recombination of the selected two organism.
				m_pPopulation[i] = pPopulation[iOrgan1];
				m_pPopulation[i].pSequence[iIndex] = pPopulation[iOrgan2].pSequence[iIndex];

				m_pPopulation[i+1] = pPopulation[iOrgan2];
				m_pPopulation[i+1].pSequence[iIndex] = pPopulation[iOrgan1].pSequence[iIndex];

				//find the best two as the children
				CCompSPS sps[4];
				ResizeOrgan( m_pPopulation[i] );
				sps[0].sps = SPScore( m_pPopulation[i] );
				sps[0].pOrgan = &(m_pPopulation[i]);

				ResizeOrgan( m_pPopulation[i+1] );
				sps[1].sps = SPScore( m_pPopulation[i+1] );
				sps[1].pOrgan = &(m_pPopulation[i+1]);

				sps[2].sps = pPopulation[iOrgan1].score;
				sps[2].pOrgan = &(pPopulation[iOrgan1]);

				sps[3].sps = pPopulation[iOrgan2].score;
				sps[3].pOrgan = &(pPopulation[iOrgan2]);

				qsort( sps, 4, sizeof(CCompSPS), compScore );
				m_pPopulation[i] = *sps[0].pOrgan;
				m_pPopulation[i+1] = *sps[1].pOrgan;
			}
		}
	}
	SAFE_DELETE_ARRAY(pPopulation);
}

void		CMSAGAAlgorithm::HorizentalRecombination()
{
	COrganism  *pPopulation( new COrganism[MSA_GA_POPULATION_SIZE] );
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		pPopulation[i] =m_pPopulation[i];
	}
	for ( int i=1; i<MSA_GA_POPULATION_SIZE-1; i+=2 )
	{
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < MSA_GA_HORIZENTAL_RECOMB_RATIO )
		{
			//randomly choose the organism to recombine with.
			int iOrgan1 = CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_POPULATION_SIZE);
			int iOrgan2 = CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_POPULATION_SIZE);

			if( iOrgan1==iOrgan2) 
			{
				m_pPopulation[i] = pPopulation[iOrgan1];
				m_pPopulation[i+1] = pPopulation[iOrgan2];
			}
			else
			{
				//randomly choose the sequence to recombine with.
				int iIndex =  CGlobalSpace::m_sUtility.getRandomNumber(pPopulation[iOrgan1].nSeqSize);

				//horizontal recombination of the selected two organism.
				m_pPopulation[i] = pPopulation[iOrgan1];
				m_pPopulation[i].pSequence[iIndex] = pPopulation[iOrgan2].pSequence[iIndex];
				ResizeOrgan( m_pPopulation[i] );

				m_pPopulation[i+1] = pPopulation[iOrgan2];
				m_pPopulation[i+1].pSequence[iIndex] = pPopulation[iOrgan1].pSequence[iIndex];
				ResizeOrgan( m_pPopulation[i+1] );
			}
		}
	}
	SAFE_DELETE_ARRAY(pPopulation);
}

void		CMSAGAAlgorithm::VerticalRecombination2009()
{
	COrganism  *pPopulation( new COrganism[MSA_GA_POPULATION_SIZE] );
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		pPopulation[i] =m_pPopulation[i];
	}

	for ( int i=1; i<MSA_GA_POPULATION_SIZE-1; i+=2 )
	{
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < MSA_GA_VERTICAL_RECOMB_RATIO )
		{
			//randomly choose the organism to recombine with.
			int iOrgan1 = CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_POPULATION_SIZE-1);
			int iOrgan2 = CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_POPULATION_SIZE-1);

			if( iOrgan1==iOrgan2) 
			{
				m_pPopulation[i] = pPopulation[iOrgan1];
				m_pPopulation[i+1] = pPopulation[iOrgan2];
			}
			else
			{
				//randomly choose the sequence position to recombine with.
				int *pCutPos2 = new int[pPopulation[iOrgan1].nSeqSize];
				int idxCutPos1 = CGlobalSpace::m_sUtility.getRandomNumber(pPopulation[iOrgan1].pSequence[0].sequence.getLen());
				for ( int k=0; k<pPopulation[iOrgan1].nSeqSize; ++k )
				{
					//find the cut point. we count the none GENSPACE character number to locate the same cut point of the sequence in another organism
					int	chNum=0;
					for ( int j=0; j<(idxCutPos1+1); ++j )
					{
						if ( pPopulation[iOrgan1].pSequence[k].sequence.getSequenceContext().at(j).m_iCode !=CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
						{
							++chNum;
						}
					}
					//find the cut point of the sequence in another organism
					int chNum2=0;
					int idxCutPos2=0;
					while( chNum2<chNum )
					{
						if ( pPopulation[iOrgan2].pSequence[k].sequence.getSequenceContext().at(idxCutPos2).m_iCode !=CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
						{
							++chNum2;
						}
						++idxCutPos2;
					}
					pCutPos2[k]= idxCutPos2;
				}
				//vertical recombination at the two sequence
				int idxMinCutPos=pCutPos2[0];
				int idxMaxCutPos=pCutPos2[0];
				for ( int k=1; k<pPopulation[iOrgan1].nSeqSize; ++k )
				{
					if( idxMaxCutPos<pCutPos2[k] )  idxMaxCutPos=pCutPos2[k];
					if( idxMinCutPos>pCutPos2[k] )  idxMinCutPos=pCutPos2[k];
				}
				//the recombination of child 1
				for ( int k=0; k<pPopulation[iOrgan1].nSeqSize; ++k )
				{
					CSeqData& FatherSeq1 = pPopulation[iOrgan1].pSequence[k].sequence.getSequenceContext();
					CSeqData& FatherSeq2 = pPopulation[iOrgan2].pSequence[k].sequence.getSequenceContext();
					CSeqData& OffSpringSeq = m_pPopulation[i].pSequence[k].sequence.getSequenceContext();
					OffSpringSeq.clear();
					for( int j=0; j<idxCutPos1+1; ++j )
					{
						OffSpringSeq.push_back( FatherSeq1[j] );
					}
					for( int j=0; j<pCutPos2[k]-idxMinCutPos; ++j )
					{
						OffSpringSeq.push_back( StruSeqElem(GENESPACE, -1) );
					}
					for( int j=pCutPos2[k]; j<(int)FatherSeq2.size(); ++j )
					{
						OffSpringSeq.push_back( FatherSeq2[j] );
					}
				}
				//the recombination of child 2
				for ( int k=0; k<pPopulation[iOrgan1].nSeqSize; ++k )
				{
					CSeqData& FatherSeq1 = pPopulation[iOrgan1].pSequence[k].sequence.getSequenceContext();
					CSeqData& FatherSeq2 = pPopulation[iOrgan2].pSequence[k].sequence.getSequenceContext();
					CSeqData& OffSpringSeq = m_pPopulation[i+1].pSequence[k].sequence.getSequenceContext();
					OffSpringSeq.clear();
					for( int j=0; j<pCutPos2[k]; ++j )
					{
						OffSpringSeq.push_back( FatherSeq2[j] );
					}
					for( int j=0; j<idxMaxCutPos-pCutPos2[k]; ++j )
					{
						OffSpringSeq.push_back( StruSeqElem(GENESPACE, -1) );
					}
					for( int j=idxCutPos1+1; j<(int)FatherSeq1.size(); ++j )
					{
						OffSpringSeq.push_back( FatherSeq1[j] );
					}
				}
				SAFE_DELETE_ARRAY(pCutPos2);

				//find the best two as the children
				CCompSPS sps[4];
				ResizeOrgan( m_pPopulation[i] );
				sps[0].sps = SPScore( m_pPopulation[i] );
				sps[0].pOrgan = &(m_pPopulation[i]);

				ResizeOrgan( m_pPopulation[i+1] );
				sps[1].sps = SPScore( m_pPopulation[i+1] );
				sps[1].pOrgan = &(m_pPopulation[i+1]);

				sps[2].sps = pPopulation[iOrgan1].score;
				sps[2].pOrgan = &(pPopulation[iOrgan1]);

				sps[3].sps = pPopulation[iOrgan2].score;
				sps[3].pOrgan = &(pPopulation[iOrgan2]);

				qsort( sps, 4, sizeof(CCompSPS), compScore );
				m_pPopulation[i] = *sps[0].pOrgan;
				m_pPopulation[i+1] = *sps[1].pOrgan;
			}
		}
	}

	SAFE_DELETE_ARRAY(pPopulation);
}

void		CMSAGAAlgorithm::VerticalRecombination()
{
	COrganism  *pPopulation( new COrganism[MSA_GA_POPULATION_SIZE] );
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		pPopulation[i] =m_pPopulation[i];
	}

	for ( int i=1; i<MSA_GA_POPULATION_SIZE-1; i+=2 )
	{
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < MSA_GA_VERTICAL_RECOMB_RATIO )
		{
			//randomly choose the organism to recombine with.
			int iOrgan1 = CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_POPULATION_SIZE);
			int iOrgan2 = CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_POPULATION_SIZE);

			if( iOrgan1==iOrgan2) 
			{
				m_pPopulation[i] = pPopulation[iOrgan1];
				m_pPopulation[i+1] = pPopulation[iOrgan2];
			}
			else
			{
				//randomly choose the sequence position to recombine with.
				int *pCutPos2 = new int[pPopulation[iOrgan1].nSeqSize];
				int idxCutPos1 = CGlobalSpace::m_sUtility.getRandomNumber(pPopulation[iOrgan1].pSequence[0].sequence.getLen());
				for ( int k=0; k<pPopulation[iOrgan1].nSeqSize; ++k )
				{
					//find the cut point. we count the none GENSPACE character number to locate the same cut point of the sequence in another organism
					int	chNum=0;
					for ( int j=0; j<(idxCutPos1+1); ++j )
					{
						if ( pPopulation[iOrgan1].pSequence[k].sequence.getSequenceContext().at(j).m_iCode !=CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
						{
							++chNum;
						}
					}
					//find the cut point of the sequence in another organism
					int chNum2=0;
					int idxCutPos2=0;
					while( chNum2<chNum )
					{
						if ( pPopulation[iOrgan2].pSequence[k].sequence.getSequenceContext().at(idxCutPos2).m_iCode !=CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
						{
							++chNum2;
						}
						++idxCutPos2;
					}
					pCutPos2[k]= idxCutPos2;
				}
				//vertical recombination at the two sequence
				int idxMinCutPos=pCutPos2[0];
				int idxMaxCutPos=pCutPos2[0];
				for ( int k=1; k<pPopulation[iOrgan1].nSeqSize; ++k )
				{
					if( idxMaxCutPos<pCutPos2[k] )  idxMaxCutPos=pCutPos2[k];
					if( idxMinCutPos>pCutPos2[k] )  idxMinCutPos=pCutPos2[k];
				}
				//the recombination of child 1
				for ( int k=0; k<pPopulation[iOrgan1].nSeqSize; ++k )
				{
					CSeqData& FatherSeq1 = pPopulation[iOrgan1].pSequence[k].sequence.getSequenceContext();
					CSeqData& FatherSeq2 = pPopulation[iOrgan2].pSequence[k].sequence.getSequenceContext();
					CSeqData& OffSpringSeq = m_pPopulation[i].pSequence[k].sequence.getSequenceContext();
					OffSpringSeq.clear();
					for( int j=0; j<idxCutPos1+1; ++j )
					{
						OffSpringSeq.push_back( FatherSeq1[j] );
					}
					for( int j=0; j<pCutPos2[k]-idxMinCutPos; ++j )
					{
						OffSpringSeq.push_back( StruSeqElem(GENESPACE, -1) );
					}
					for( int j=pCutPos2[k]; j<(int)FatherSeq2.size(); ++j )
					{
						OffSpringSeq.push_back( FatherSeq2[j] );
					}
				}
				//the recombination of child 2
				for ( int k=0; k<pPopulation[iOrgan1].nSeqSize; ++k )
				{
					CSeqData& FatherSeq1 = pPopulation[iOrgan1].pSequence[k].sequence.getSequenceContext();
					CSeqData& FatherSeq2 = pPopulation[iOrgan2].pSequence[k].sequence.getSequenceContext();
					CSeqData& OffSpringSeq = m_pPopulation[i+1].pSequence[k].sequence.getSequenceContext();
					OffSpringSeq.clear();
					for( int j=0; j<pCutPos2[k]; ++j )
					{
						OffSpringSeq.push_back( FatherSeq2[j] );
					}
					for( int j=0; j<idxMaxCutPos-pCutPos2[k]; ++j )
					{
						OffSpringSeq.push_back( StruSeqElem(GENESPACE, -1) );
					}
					for( int j=idxCutPos1+1; j<(int)FatherSeq1.size(); ++j )
					{
						OffSpringSeq.push_back( FatherSeq1[j] );
					}
				}
				ResizeOrgan(m_pPopulation[i]);
				ResizeOrgan(m_pPopulation[i+1]);
				SAFE_DELETE_ARRAY(pCutPos2);
			}
		}
	}
	
	SAFE_DELETE_ARRAY(pPopulation);
}

void		CMSAGAAlgorithm::Selection2004()
{
	//use "roulette wheel selection" algorithm, we may implement "tournament selection" in the selection procedure later as the paper does.

	//select the best organism as the first organism of the next generation
	//exchange the organism if the best one is NOT the first one
	int bestScore = -INT_MAX;
	int bestOrganismNo = -1;
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( bestScore < m_pPopulation[i].score )  
		{
			bestScore = m_pPopulation[i].score;
			bestOrganismNo = i;
		}
	}

	//use the score field of TmpPath to store the proportion of the score
	COrganism  *pPopulation( new COrganism[MSA_GA_POPULATION_SIZE] );
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		pPopulation[i] =m_pPopulation[i];
	}

	//reserve the best one before the selection
	if( bestOrganismNo!=0 )
	{
		m_pPopulation[0] = pPopulation[bestOrganismNo];
		m_pPopulation[bestOrganismNo] =  pPopulation[0];
	}

	bestScore = -INT_MAX;
	bestOrganismNo = -1;
	for ( int i=1; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( bestScore < m_pPopulation[i].score )  
		{
			bestScore = m_pPopulation[i].score;
			bestOrganismNo = i;
		}
	}

	//reserve the second best one before the selection
	if( bestOrganismNo!=0 )
	{
		m_pPopulation[1] = pPopulation[bestOrganismNo];
		m_pPopulation[bestOrganismNo] =  pPopulation[1];
	}

	int lowestScore=0; 
	int highestScore=0;
	highestScore = lowestScore = pPopulation[0].score;
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( pPopulation[i].score<lowestScore )
		{
			lowestScore = pPopulation[i].score;
		}
		if ( pPopulation[i].score>highestScore )
		{
			highestScore = pPopulation[i].score;
		}
	}

	int scoreRange = highestScore - lowestScore;

	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		pPopulation[i].adjustScore = pPopulation[i].score - lowestScore;//change all score as positive number
		//if the adjust score is zero, there may be ERROR occur, we just add one point to the adjustScore
		++(pPopulation[i].adjustScore);
	}

	int totalScore = 0;
	//int arrBoundary[MSA_GA_POPULATION_SIZE][2];
	int *pArrBoundary = new int(MSA_GA_POPULATION_SIZE*2);
	//use the score field of TmpPath to store the proportion of the score
	for (int i = 0; i<MSA_GA_POPULATION_SIZE; ++i)
	{
		pArrBoundary[i*2+0] = totalScore;
		totalScore += pPopulation[i].adjustScore;
		pArrBoundary[i*2+1] = totalScore - 1;
	}

	if( totalScore==0 ) 		
	{
		SAFE_DELETE_ARRAY( pPopulation );
		return;
	}

	//select the organism of next generation with roulette-wheel method
	for ( int i=2; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		int idx = CGlobalSpace::m_sUtility.getRandomNumber(totalScore);
		int k=0;
		while (idx > pArrBoundary[k*2+1]) ++k;
		//create next generation with all new copy!
		m_pPopulation[i] = pPopulation[k];
	}

	SAFE_DELETE_ARRAY(pArrBoundary);
	SAFE_DELETE_ARRAY( pPopulation );
}

void		CMSAGAAlgorithm::Selection()
{
	//use "roulette wheel selection" algorithm, we may implement "tournament selection" in the selection procedure later as the paper does.

	//select the best organism as the first organism of the next generation
	//exchange the organism if the best one is NOT the first one
	int bestScore = -INT_MAX;
	int bestOrganismNo = -1;
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( bestScore < m_pPopulation[i].score )  
		{
			bestScore = m_pPopulation[i].score;
			bestOrganismNo = i;
		}
	}

	//use the score field of TmpPath to store the proportion of the score
	COrganism  *pPopulation( new COrganism[MSA_GA_POPULATION_SIZE] );
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		pPopulation[i] =m_pPopulation[i];
	}

	//reserve the best one before the selection
	if( bestOrganismNo!=0 )
	{
		m_pPopulation[0] = pPopulation[bestOrganismNo];
	}
		
	int lowestScore=0; 
	int highestScore=0;
	highestScore = lowestScore = pPopulation[0].score;
	for ( int i=1; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		if ( pPopulation[i].score<lowestScore )
		{
			lowestScore = pPopulation[i].score;
		}
		if ( pPopulation[i].score>highestScore )
		{
			highestScore = pPopulation[i].score;
		}
	}

	int scoreRange = highestScore - lowestScore;

	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		pPopulation[i].adjustScore = pPopulation[i].score - lowestScore;//change all score as positive number
		//if the adjust score is zero, there may be ERROR occur, we just add one point to the adjustScore
		++(pPopulation[i].adjustScore);
	}

	int totalScore = 0;
	//int arrBoundary[MSA_GA_POPULATION_SIZE][2];
	int *pArrBoundary = new int[MSA_GA_POPULATION_SIZE*2];
	//use the score field of TmpPath to store the proportion of the score
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		pArrBoundary[i*2+0] = totalScore;
		totalScore += pPopulation[i].adjustScore;
		pArrBoundary[i*2+1] = totalScore - 1;
	}

	if( totalScore==0 ) 		
	{
		SAFE_DELETE_ARRAY( pPopulation );
		return;
	}

	//select the organism of next generation with roulette-wheel method
	for ( int i=1; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		int idx = CGlobalSpace::m_sUtility.getRandomNumber(totalScore);
		int k=0;
		while (idx > pArrBoundary[k*2+1]) ++k;
		//create next generation with all new copy!
		m_pPopulation[i] = pPopulation[k];
	}

	SAFE_DELETE_ARRAY(pArrBoundary);
	SAFE_DELETE_ARRAY( pPopulation );
}

bool    CMSAGAAlgorithm::testOrgLen( COrganism* vpPopulation )
{
	int iShortestAcidLen = INT_MAX;
	int idxShortestSeq = -1;
#ifdef DEBUG
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		for ( int j=0; j<vpPopulation[i].nSeqSize; ++j)
		{
			int iLen = 0;
			for( int k=0; k<vpPopulation[i].pSequence[j].sequence.getLen(); ++k )
			{
				if ( vpPopulation[i].pSequence[j].sequence.getSequenceContext().at(k).m_iCode!=GENESPACECODE )
				{
					++iLen;
				}
			}
			if(iShortestAcidLen>iLen)	iShortestLen=iLen;
		} 
	}
	_ASSERT( iShortestAcidLen==m_idbgShortestInputAcidNum )
	return true;
#else
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		for ( int j=0; j<vpPopulation[i].nSeqSize; ++j)
		{
			int iLen = 0;
			for( int k=0; k<vpPopulation[i].pSequence[j].sequence.getLen(); ++k )
			{
				if ( vpPopulation[i].pSequence[j].sequence.getSequenceContext().at(k).m_iCode!=CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
				{
					++iLen;
				}
			}
			if(iShortestAcidLen>iLen)
			{
				iShortestAcidLen=iLen;
				idxShortestSeq = j;
			}
		} 
	}
	if( iShortestAcidLen!=m_idbgShortestInputAcidNum )
	{
		char ch[100];
		sprintf_s(ch,  "generation=%i, maxscore=%f, shortestseqlen=%i, id=%i, m_idbgShortestInputAcidNum=%i" ,
			          m_iRun, m_fMaxScore, iShortestAcidLen, idxShortestSeq, m_idbgShortestInputAcidNum );
		CGlobalSpace::m_sEventLog.writeEvent(ch);				
		return false;
	}
	return true;
#endif
}

int		CMSAGAAlgorithm::CountSeqCharNum( const CSequence& vSeq )
{
	int iLen = 0;
	for( int k=0; k<vSeq.getLen(); ++k )
	{
		if ( vSeq.getSequenceContext().at(k).m_iCode!=CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
		{
			++iLen;
		}
	}
	return iLen;
}

void		CMSAGAAlgorithm::InitPopulationREV(const std::vector<CSequence>& vSequences, int vCol )
{
	m_pPopulation = new COrganism[MSA_GA_POPULATION_SIZE];
	m_OrganismSize = vSequences.size();
	int *arrCharIdx = new int[vSequences.size()];
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		CSeqData  *tmpSeqData = new CSeqData[vSequences.size()];
		m_pPopulation[i].pSequence = new CMSA_GASeq[vSequences.size()];
		m_pPopulation[i].nSeqSize=vSequences.size();
		m_pPopulation[i].score = 0;
		
		for ( int j=0; j<(int)vSequences.size(); ++j )
		{
			m_pPopulation[i].pSequence[j].charNum = CountSeqCharNum( vSequences[j] );
			arrCharIdx[j] = 0;
		}

		bool		bSeqEnd = false;
		while(!bSeqEnd)
		{
			for (  int j=0; j<(int)vSequences.size(); ++j  )
			{
				if ( arrCharIdx[j]==m_pPopulation[i].pSequence[j].charNum )
				{
					bSeqEnd = true;
					break;
				}
			}
			if(!bSeqEnd)
			{
				int n = CGlobalSpace::m_sUtility.getRandomNumber(vSequences.size())+1;
				int  idxRow=0;
				while( n!=0 )
				{
					if ( n<2 )
					{
						tmpSeqData[idxRow].push_back( StruSeqElem( vSequences[idxRow].getSequenceContext().at(arrCharIdx[idxRow]).m_char, vSequences[idxRow].getSequenceContext().at(arrCharIdx[idxRow]).m_index ) );
						arrCharIdx[idxRow] += 1;
						n=0;  //end the column.
					}
					else
					{
						int idx = n%2;
						if(1==idx)
						{
							tmpSeqData[idxRow].push_back( StruSeqElem( vSequences[idxRow].getSequenceContext().at(arrCharIdx[idxRow]).m_char, vSequences[idxRow].getSequenceContext().at(arrCharIdx[idxRow]).m_index ) );
							arrCharIdx[idxRow] += 1;
						}
						else
						{
							tmpSeqData[idxRow].push_back( StruSeqElem(GENESPACE,CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE)) );
						}
						n/=2;
					}
				}
			}
		}

		for ( int j=0; j<(int)vSequences.size(); ++j )
		{
			while( arrCharIdx[j]<m_pPopulation[i].pSequence[j].charNum )
			{
				tmpSeqData[j].push_back( StruSeqElem( vSequences[j].getSequenceContext().at(arrCharIdx[j]).m_char, vSequences[j].getSequenceContext().at(arrCharIdx[j]).m_index ) );
				arrCharIdx[j] += 1;
			}
			m_pPopulation[i].pSequence[j].sequence = CSequence( tmpSeqData[j], vSequences[j].getName(), vSequences[j].getTitle() );
		}
		if(!CheckSumCharNum( m_pPopulation[i]  ))
		{
			CGlobalSpace::m_sEventLog.writeEvent("Init Population ERROR!");
		}
		SAFE_DELETE_ARRAY(tmpSeqData);
	}
	SAFE_DELETE_ARRAY(arrCharIdx);

	ArrangeSequences();
	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA: Population Initialized with REV method!" );
}

void		CMSAGAAlgorithm::InitPopulation2009( const std::vector<CSequence>& vSequences, int vCol )
{
	m_pPopulation = new COrganism[MSA_GA_POPULATION_SIZE];
	m_OrganismSize = vSequences.size();
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		m_pPopulation[i].pSequence = new CMSA_GASeq[vSequences.size()];
		m_pPopulation[i].nSeqSize=vSequences.size();
		m_pPopulation[i].score = 0;

		for ( int j=0; j<(int)vSequences.size(); ++j )
		{
			m_pPopulation[i].pSequence[j].charNum = CountSeqCharNum( vSequences[j] );
			m_pPopulation[i].pSequence[j].gapOffset =  CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_OFFSET_FACTOR_X*m_iLongestSeqSize);
			CSeqData  tmpSeqData;
			for ( int k=0; k<m_pPopulation[i].pSequence[j].gapOffset; ++k )
			{
				tmpSeqData.push_back( StruSeqElem(GENESPACE,CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE)) );
			}
			for ( int k=0; k<vSequences[j].getLen(); ++k )
			{
				tmpSeqData.push_back( StruSeqElem( vSequences[j].getSequenceContext().at(k).m_char, vSequences[j].getSequenceContext().at(k).m_index ) );
			}
			m_pPopulation[i].pSequence[j].sequence = CSequence( tmpSeqData, vSequences[j].getName(), vSequences[j].getTitle() );
		}
	}
	ArrangeSequences();

	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA: Population Initialized using method in paper 2009!" );
}

void		CMSAGAAlgorithm::InitPopulation2007( const std::vector<CSequence>& vSequences, int vCol )
{
	CSequence  *pNWAlignmentMatrix (new CSequence[vSequences.size()*(vSequences.size()-1)]);
	
	//use NW algorithm to align every two sequences in the vSequences
	std::auto_ptr<CNWAlgorithm>	pNWAlgorithm( new CNWAlgorithm );
	for ( int i=0; i<(int)vSequences.size()-1; ++i )
	{
		for ( int j=i+1; j<(int)vSequences.size(); ++j )
		{
			std::vector<CSequence> VectEachCompOrg, VectEachCompResult;
			VectEachCompOrg.push_back( vSequences[i] );
			VectEachCompOrg.push_back( vSequences[j] );
			pNWAlgorithm->Align( VectEachCompResult, VectEachCompOrg );
			//create the Alignment Matrix
			pNWAlignmentMatrix[ POS(i, j-1, vSequences.size()-1) ] = VectEachCompResult[0];
			pNWAlignmentMatrix[ POS(j, i, vSequences.size()-1) ]  = VectEachCompResult[1];
		}
	}

	m_pPopulation = new COrganism[MSA_GA_POPULATION_SIZE];
	m_OrganismSize = vSequences.size();
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		m_pPopulation[i].pSequence = new CMSA_GASeq[vSequences.size()];
		m_pPopulation[i].nSeqSize=vSequences.size();
		m_pPopulation[i].score = 0;

		for ( int j=0; j<(int)vSequences.size(); ++j )
		{
			m_pPopulation[i].pSequence[j].charNum = CountSeqCharNum( vSequences[j] );
			m_pPopulation[i].pSequence[j].gapOffset =  CGlobalSpace::m_sUtility.getRandomNumber(MSA_GA_OFFSET_FACTOR_X*m_iLongestSeqSize);
			CSeqData  tmpSeqData;
			for ( int k=0; k<m_pPopulation[i].pSequence[j].gapOffset; ++k )
			{
				tmpSeqData.push_back( StruSeqElem(GENESPACE,CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE)) );
			}
			CSequence&  NWAlignSeq = pNWAlignmentMatrix[POS(j, CGlobalSpace::m_sUtility.getRandomNumber(vSequences.size()-1), vSequences.size()-1)];
			for ( int k=0; k<NWAlignSeq.getLen(); ++k )
			{
				tmpSeqData.push_back( StruSeqElem( NWAlignSeq.getSequenceContext()[k].m_char, NWAlignSeq.getSequenceContext()[k].m_index ) );
			}
			m_pPopulation[i].pSequence[j].sequence = CSequence( tmpSeqData, NWAlignSeq.getName(), NWAlignSeq.getTitle() );
		}
	}
	SAFE_DELETE_ARRAY(pNWAlignmentMatrix);

	//delete void GENSPACE of the sequences
	ArrangeSequences();

	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA: Population Initialized using method in paper 2007!" );
}

int		CMSAGAAlgorithm::GetLongestSeqLen( const COrganism& og )
{
	//find the length of the longest sequence
	int iLongestLength = og.pSequence[0].sequence.getLen();
	for ( int j=1; j<og.nSeqSize; ++j )
	{
		if ( iLongestLength<og.pSequence[j].sequence.getLen() )
		{
			iLongestLength = og.pSequence[j].sequence.getLen();
		}
	}
	return iLongestLength;
}

void     CMSAGAAlgorithm::ResizeOrgan( COrganism& vOrgan )
{
	CSeqData::iterator	*itrSeq(new CSeqData::iterator[vOrgan.nSeqSize]);
	for ( int j=0; j<vOrgan.nSeqSize; ++j )
	{
		itrSeq[j] = vOrgan.pSequence[j].sequence.getSequenceContext().begin();
	}

	int iLongestLength = GetLongestSeqLen( vOrgan );
	for ( int j=0; j<iLongestLength; ++j )
	{
		bool		bSpace = true;
		for ( int k=0; k<vOrgan.nSeqSize; ++k )
		{
			if ( itrSeq[k]!=vOrgan.pSequence[k].sequence.getSequenceContext().end() && itrSeq[k]->m_iCode != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
			{
				bSpace = false;
				break;
			}
		}
		if( bSpace )	
		{
			for ( int k=0; k<vOrgan.nSeqSize; ++k )
			{
				if ( itrSeq[k]!=vOrgan.pSequence[k].sequence.getSequenceContext().end() )
				{
					itrSeq[k] = vOrgan.pSequence[k].sequence.getSequenceContext().erase( itrSeq[k] );
				}
			}
		}
		else
		{
			for ( int k=0; k<vOrgan.nSeqSize; ++k )
			{
				if ( itrSeq[k]!=vOrgan.pSequence[k].sequence.getSequenceContext().end() )
				{
					++itrSeq[k];
				}
			}
		}
	}

	SAFE_DELETE_ARRAY(itrSeq);

	//recalculate the length after it has changed.
	iLongestLength = GetLongestSeqLen( vOrgan );

	//adding GENSPACE at the shorter sequence tail.
	for ( int j=0; j<vOrgan.nSeqSize; ++j )
	{
		if ( vOrgan.pSequence[j].sequence.getLen()<iLongestLength )
		{
			int orglen = vOrgan.pSequence[j].sequence.getLen();
			for ( int k=0; k<(iLongestLength-orglen); ++k )
			{
				vOrgan.pSequence[j].sequence.getSequenceContext().push_back( StruSeqElem(GENESPACE, -1) );
			}				
		}
	}
}

void		CMSAGAAlgorithm::ArrangeSequences()
{
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		ResizeOrgan( m_pPopulation[i] );
	}
}

float		CMSAGAAlgorithm::SPScore2009( const COrganism& vOrgan )
{
	float spScore = 0;
	//pair score
	for ( int j=0; j<vOrgan.pSequence[0].sequence.getLen(); ++j )
	{
		for ( int k=0; k<vOrgan.nSeqSize-1; ++k )
		{
			for ( int m=k+1; m<vOrgan.nSeqSize; ++m )
			{
				int a =  vOrgan.pSequence[k].sequence.getSequenceContext().at(j).m_iCode;
				int b =  vOrgan.pSequence[m].sequence.getSequenceContext().at(j).m_iCode;
				if ( a != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) && b != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
				{
#ifdef CLUSTALW_SP_SCORE
					spScore += 100*CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getSubMatrixScore(MSAGA_MAT_TYPE, a, b)*m_SeqWeight[k] * m_SeqWeight[m];
#else
					spScore += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getSubMatrixScore(MSAGA_MAT_TYPE, a, b);
#endif
				}
			}
		}
	}
	spScore = spScore * 2 / ((vOrgan.nSeqSize - 1)*vOrgan.nSeqSize);
	//gap penalty
	float penalty=0;
	for ( int j=0; j<vOrgan.nSeqSize; ++j )
	{
		int gapNum=0;
		int idx=0;
		while( idx<vOrgan.pSequence[j].sequence.getLen() )
		{
			while( idx<vOrgan.pSequence[j].sequence.getLen() && vOrgan.pSequence[j].sequence.getSequenceContext().at(idx).m_iCode != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )	++idx;
			while( idx<vOrgan.pSequence[j].sequence.getLen() && vOrgan.pSequence[j].sequence.getSequenceContext().at(idx).m_iCode == CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )	
			{
				++gapNum;
				++idx;
			}
			if( gapNum>0 )
			{
				penalty += (CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(MSAGA_MAT_TYPE) + (gapNum - 1)*CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapExtendCost(MSAGA_MAT_TYPE));
			}
			gapNum=0;
		}
	}
	penalty = penalty / vOrgan.nSeqSize;
	return (spScore + penalty);
}

int		CMSAGAAlgorithm::SPScore( const COrganism& vOrgan )
{
	int spScore = 0;
	//pair score
	for ( int j=0; j<vOrgan.pSequence[0].sequence.getLen(); ++j )
	{
		for ( int k=0; k<vOrgan.nSeqSize-1; ++k )
		{
			for ( int m=k+1; m<vOrgan.nSeqSize; ++m )
			{
				int a =  vOrgan.pSequence[k].sequence.getSequenceContext().at(j).m_iCode;
				int b =  vOrgan.pSequence[m].sequence.getSequenceContext().at(j).m_iCode;
				if ( a != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) && b != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
				{
					spScore += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getSubMatrixScore( MSAGA_MAT_TYPE, a, b );
				}
			}
		}
	}
	//gap penalty
	for ( int j=0; j<vOrgan.nSeqSize; ++j )
	{
		int gapNum=0;
		int idx=0;
		while( idx<vOrgan.pSequence[j].sequence.getLen() )
		{
			while( idx<vOrgan.pSequence[j].sequence.getLen() && vOrgan.pSequence[j].sequence.getSequenceContext().at(idx).m_iCode != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )	++idx;
			while( idx<vOrgan.pSequence[j].sequence.getLen() && vOrgan.pSequence[j].sequence.getSequenceContext().at(idx).m_iCode == CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )	
			{
				++gapNum;
				++idx;
			}
			if( gapNum>0 )
			{
				spScore += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(MSAGA_MAT_TYPE) + (gapNum-1)*CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapExtendCost(MSAGA_MAT_TYPE);
			}
			gapNum=0;
		}
	}
	return spScore;
}

bool		CMSAGAAlgorithm::Fitness()
{
	if ( m_OrganismSize<0 )
	{
		throw CAppException( DEF_EXCEPTION_INVALID_PARAMETER, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "Bad parameter!" );
	}

	m_fMaxScore = -INT_MAX;
	//use SP score as the fitness
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		m_pPopulation[i].score = SPScore( m_pPopulation[i] );
		if( m_fMaxScore<m_pPopulation[i].score )	m_fMaxScore=m_pPopulation[i].score;
	}

	return false;//always run the max generation of evolution.
}

bool		CMSAGAAlgorithm::Fitness2009()
{
	if ( m_OrganismSize<0 )
	{
		throw CAppException( DEF_EXCEPTION_INVALID_PARAMETER, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "Bad parameter!" );
	}

	m_fMaxScore = -INT_MAX;
	//use SP score as the fitness
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		m_pPopulation[i].score = SPScore2009( m_pPopulation[i] );
		if( m_fMaxScore<m_pPopulation[i].score )	m_fMaxScore=m_pPopulation[i].score;
	}

	return false;//always run the max generation of evolution.
}

void		CMSAGAAlgorithm::testSPScore( const std::vector<CSequence>& vSequences )
{
	COrganism	org;
	org.pSequence = new CMSA_GASeq[vSequences.size()];
	org.nSeqSize = vSequences.size();
	for ( int i=0; i<vSequences.size(); ++i )
	{
		org.pSequence[i].sequence = vSequences[i];
	}
	ResizeOrgan( org );
	
	int spScore = SPScore(org);

	char ch[50];
	sprintf_s( ch, "The sp score of original aln is %i", spScore );
	CGlobalSpace::m_sEventLog.writeEvent( ch );

	SAFE_DELETE_ARRAY(org.pSequence);	
}

};