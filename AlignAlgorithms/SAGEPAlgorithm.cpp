#include<memory>
#include <ctime>
#include"Common.h"
#include"SAGEPAlgorithm.h"
#include"GlobalSpace.h"
#include "GeneDecoder.h"

namespace SeqAnsis
{

void	CSAGEPAlgorithm::Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences )
{
	if( vSequences.size()!=2 )
	{
		CGlobalSpace::m_sEventLog.writeEvent("Wrong sequences size to align by SAGEP algorithm!");
		return;
	}

	int nSeq1 = vSequences[0].getSequenceContext().size();
	int nSeq2 = vSequences[1].getSequenceContext().size();

	m_FuncSet = std::string("A");
	m_HeadSet = std::string("ARBD");
	m_TailSet = std::string("RBD");

	int idxRow=0;
	int idxCol=0;

	while( idxRow<nSeq1 && idxCol<nSeq2 )
	{
		initPopulation( idxRow, idxCol );
		int iEvolutionCount=0;
		m_iMaxScore=INT_MIN;
		int maxScore = INT_MIN;
		int iCount = 0;
		if( !fitness( vSequences ) )
		{
			while( iEvolutionCount<EVOLUTION_GENERATION_NUM && iCount<GEP_NO_ADV_GENERATION_NUM )	
			{
				selection();
				reproduction();
				++iEvolutionCount;
				if(fitness(vSequences))	break;	
				if( m_iMaxScore>maxScore )
				{
					iCount=0;
					maxScore = m_iMaxScore;
				}
				else
				{
					++iCount;
				}
			}
		}

		char ch[50];
		sprintf_s(ch,  "generation=%i, maxscore=%i" , iEvolutionCount, m_iMaxScore );
		CGlobalSpace::m_sEventLog.writeEvent(ch);

		translateChromosome(idxRow, idxCol, vAlignedSequences, vSequences);
	}

}

void	CSAGEPAlgorithm::translateChromosome( int& voRow, int& voCol, std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences )
{
	//find the best chromosome.
	int bestScore = -INT_MAX;
	int bestPathNo = -1;

	for ( int i=0; i<POPULATION_NUM; ++i )
	{
		if ( bestScore < m_Population[i]->score )
		{
			bestScore = m_Population[i]->score;
			bestPathNo = i;
		}
	}
	//translate the GEP expression to aligned sequence
	CGeneDecoder  GeneDecoder;
	int idx0 = m_Population[bestPathNo]->x-1;//idx0 point to seq1
	int idx1 = m_Population[bestPathNo]->y-1;//idx1 point to seq2

	//link all the sub path with 'A' function.
	std::string strPath;
	int idxPath = 0;
	for ( int j=0; j<GENE_NUM_PER_CHROMOSOME; ++j )
	{
		GeneDecoder.Decode( &(m_Population[bestPathNo]->Path[j*GENE_LEN]), GENE_LEN );
		GeneDecoder.GetPath( strPath );
	}

	bool bFinished = false;
	int iPathIdx = 0;
	int nSeq1 = vSequences[0].getSequenceContext().size();
	int nSeq2 = vSequences[1].getSequenceContext().size();
	const CSeqData& Seq1 = vSequences[0].getSequenceContext();
	const CSeqData& Seq2 = vSequences[1].getSequenceContext();

	if ( vAlignedSequences.empty() )
	{
		CSeqData tmp;
		vAlignedSequences.push_back( CSequence( tmp, vSequences[0].getName(), vSequences[0].getTitle() ) );
		vAlignedSequences.push_back( CSequence( tmp, vSequences[1].getName(), vSequences[1].getTitle() ) );
	}

	CSeqData& ReSeq1 = vAlignedSequences[0].getSequenceContext();
	CSeqData& ReSeq2 = vAlignedSequences[1].getSequenceContext();

	while( iPathIdx<(int)(strPath.length()) && !bFinished )
	{
		switch( strPath[iPathIdx] )
		{
		case 'R':
			{
				++idx0;
				if( idx0<nSeq1 )
				{
					//the element of seq1 match space
					ReSeq1.push_back(Seq1[idx0]);
					ReSeq2.push_back(StruSeqElem(GENESPACE,-1));
				}
				else
				{
					bFinished = true;
				}
			}
			break;
		case 'B':
			{
				++idx1;
				if( idx1<nSeq2 )
				{
					ReSeq1.push_back(StruSeqElem(GENESPACE,-1));
					ReSeq2.push_back(Seq2[idx1]);
				}
				else
				{
					bFinished = true;
				}
			}
			break;
		case 'D':
			{
				++idx0;
				++idx1;
				if( idx1<nSeq2 && idx0<nSeq1 )
				{
					ReSeq1.push_back(Seq1[idx0]);
					ReSeq2.push_back(Seq2[idx1]);
				}
				else
				{
					bFinished = true;
				}
			}
			break;
		}
		++iPathIdx;
	}



	voRow = idx0+1;
	voCol = idx1+1;
}

void	CSAGEPAlgorithm::initPopulation(int vX, int vY)
{
	m_Population.clear();
	for ( int i=0; i<POPULATION_NUM; ++i )
	{
		std::auto_ptr<CChromosome>	pChromosome( new CChromosome );
		pChromosome->x = vX;
		pChromosome->y = vY;
		for( int k=0; k<GENE_NUM_PER_CHROMOSOME; ++k )
		{
			for ( int j=0; j<GENE_HEAD_LEN; ++j )
			{
				pChromosome->Path[j+k*GENE_LEN] = m_HeadSet[CGlobalSpace::m_sUtility.getRandomNumber(m_HeadSet.length())];
			}
			for ( int j=GENE_HEAD_LEN; j<GENE_LEN; ++j )
			{
				pChromosome->Path[j+k*GENE_LEN] = m_TailSet[CGlobalSpace::m_sUtility.getRandomNumber(m_TailSet.length())];
			}
		}
		pChromosome->score = pChromosome->adjustScore = -INT_MAX;

		m_Population.push_back( pChromosome );
	} 
}

bool	CSAGEPAlgorithm::fitness(const std::vector<CSequence>& vSequences)
{
	int nSeq1 = vSequences[0].getSequenceContext().size();
	int nSeq2 = vSequences[1].getSequenceContext().size();
	const CSeqData& Seq1 = vSequences[0].getSequenceContext();
	const CSeqData& Seq2 = vSequences[1].getSequenceContext();

	int bestScore = m_Population[0]->score;
	CGeneDecoder GeneDecoder;
	for ( int i=0; i<POPULATION_NUM; ++i )
	{
		int idx0 = m_Population[i]->x-1;//idx0 point to seq1
		int idx1 = m_Population[i]->y-1;//idx1 point to seq2

		//link all the sub path with 'A' function.
		std::string		strPath;
		for ( int j=0; j<GENE_NUM_PER_CHROMOSOME; ++j )
		{
			GeneDecoder.Decode( &(m_Population[i]->Path[j*GENE_LEN]), GENE_LEN );
			GeneDecoder.GetPath( strPath );
		}

		//score
		int score=0;
		bool bFinished = false;
		int iPathIdx = 0;
		char preCH = 'D';
		while( iPathIdx<(int)(strPath.length()) && !bFinished )
		{
			switch( strPath[iPathIdx] )
			{
			case 'R':
				{
					++idx0;
					if( idx0<nSeq1 )
					{
						if( preCH ==  'R' )
						{
							score += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapExtendCost(SAGEP_MAT_TYPE);
						}
						else
						{
							score += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(SAGEP_MAT_TYPE);
						}
						preCH =  'R';
					}
					else
					{
						bFinished = true;
					}
				}
				break;
			case 'B':
				{	
					++idx1;
					if( idx1<nSeq2 )
					{
						if( preCH ==  'B' )
						{
							score += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapExtendCost(SAGEP_MAT_TYPE);
						}
						else
						{
							score += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(SAGEP_MAT_TYPE);
						}
						preCH =  'B';
					}
					else
					{
						bFinished = true;
					}
				}
				break;
			case 'D':
				{
					++idx0;
					++idx1;
					if( idx1<nSeq2 && idx0<nSeq1 )
					{
						score += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getSubMatrixScore(  SAGEP_MAT_TYPE, Seq1[idx0].m_iCode, Seq2[idx1].m_iCode );
						preCH =  'D';
					}
					else
					{
						bFinished = true;
					}
				}
				break;
			}
			++iPathIdx;
		}

		m_Population[i]->score = score;
		if(m_iMaxScore<m_Population[i]->score)  
		{
			m_iMaxScore = m_Population[i]->score;
		}
#if 0
		if ( bestScore < m_Population[i]->score )
		{
			bestScore = m_Population[i]->score;
		}
#endif
	}
	return false;//bestScore >= MIN_SCORE? true:false;
}

void	CSAGEPAlgorithm::selection()
{
	//select the best chromosome as the first chromosome of the next generation
	//exchange the chromosome if the best one is NOT the first one
	int bestScore = -INT_MAX;
	int bestChromosomeNo = -1;
	for ( int i=0; i<POPULATION_NUM; ++i )
	{
		if ( bestScore < m_Population[i]->score )  
		{
			bestScore = m_Population[i]->score;
			bestChromosomeNo = i;
		}
	}

	//exchange chromosome to ensure the best one is always in the next generation
	if ( 0!=bestChromosomeNo )
	{
		std::shared_ptr<CChromosome>  pChrom = m_Population[0];
		m_Population[0] = m_Population[bestChromosomeNo];
		m_Population[bestChromosomeNo] = pChrom;
		bestChromosomeNo = 0;
	}

	//acquire the score range
	int lowestScore=0; 
	int highestScore=0;
	highestScore = lowestScore = m_Population[0]->score;
	for ( int i=0; i<POPULATION_NUM; ++i )
	{
		if ( m_Population[i]->score<lowestScore )
		{
			lowestScore = m_Population[i]->score;
		}
		if ( m_Population[i]->score>highestScore )
		{
			highestScore = m_Population[i]->score;
		}
	}

	int scoreRange = highestScore - lowestScore;
	for ( int i=0; i<POPULATION_NUM; ++i )
	{
		m_Population[i]->adjustScore = m_Population[i]->score + abs(lowestScore);//change all score as positive number
		//if the adjust score is zero, there may be ERROR occur, we just add one point to the adjustScore
		++(m_Population[i]->adjustScore);
	}

	int totalScore = 0;
	std::vector<std::shared_ptr<CChromosome>>  TmpPath;
	int arrBoundary[POPULATION_NUM][2];
	//use the score field of TmpPath to store the proportion of the score
	for ( int i=0; i<POPULATION_NUM; ++i )
	{
		arrBoundary[i][0] = totalScore;
		totalScore += m_Population[i]->adjustScore;
		arrBoundary[i][1] = totalScore - 1;
		//copy the chromosome pointer to tmp path
		TmpPath.push_back(m_Population[i]);
	}

	//select the chromosomes of next generation with roulette-wheel method
	for ( int i=1; i<POPULATION_NUM; ++i )
	{
		int idx = CGlobalSpace::m_sUtility.getRandomNumber(totalScore);
		int k=0;
		while( idx > arrBoundary[k][1] ) ++k;
		//create next generation chromosome with all new copy!
		*(m_Population[i]) = *(TmpPath[k]);
	}
}

void	CSAGEPAlgorithm::reproduction()
{
	mutation();
	transposition();
	recombination();
}

void	CSAGEPAlgorithm::recombination()
{
	recombinationP1();
	recombinationP2();
	if ( GENE_NUM_PER_CHROMOSOME>1 )  recombinationGene();
}

void	CSAGEPAlgorithm::recombinationP1()
{
	if ( CGlobalSpace::m_sUtility.getRandomNumber() < RECOMB_RATE_P1 )
	{
		int p1 = CGlobalSpace::m_sUtility.getRandomNumber(CHROMOSOME_LEN);
		int iRecombLen = CHROMOSOME_LEN - p1;
		//randomly choose the chromosome to recombine with.
		int ichromosome0 = CGlobalSpace::m_sUtility.getRandomNumber(POPULATION_NUM-1) + 1;
		int ichromosome1 = CGlobalSpace::m_sUtility.getRandomNumber(POPULATION_NUM-1) + 1;
		if ( ichromosome0 != ichromosome1 )
		{
			char *pRecomb = new char[iRecombLen];
			memcpy( pRecomb, &(m_Population[ichromosome0]->Path[p1]), iRecombLen );
			memcpy( &(m_Population[ichromosome0]->Path[p1]), &(m_Population[ichromosome1]->Path[p1]), iRecombLen );
			memcpy( &(m_Population[ichromosome1]->Path[p1]), pRecomb, iRecombLen );
			SAFE_DELETE_ARRAY(pRecomb);
		}
	}
}

void	CSAGEPAlgorithm::recombinationP2()
{
	if ( CGlobalSpace::m_sUtility.getRandomNumber() < RECOMB_RATE_P2 )
	{
		int p1 = CGlobalSpace::m_sUtility.getRandomNumber(CHROMOSOME_LEN);
		int p2 = CGlobalSpace::m_sUtility.getRandomNumber(CHROMOSOME_LEN);
		//we define p1<=p2
		if( p1 > p2 )
		{
			int tmp = p1;
			p1 = p2;
			p2 = tmp;
		}

		int iRecombLen = p2 - p1 + 1;

		//randomly choose the chromosome to recombine with.
		int ichromosome0 = CGlobalSpace::m_sUtility.getRandomNumber(POPULATION_NUM-1) + 1;
		int ichromosome1 = CGlobalSpace::m_sUtility.getRandomNumber(POPULATION_NUM-1) + 1;
		if ( ichromosome0 != ichromosome1 )
		{
			char *pRecomb = new char[iRecombLen];
			memcpy( pRecomb, &(m_Population[ichromosome0]->Path[p1]), iRecombLen );
			memcpy( &(m_Population[ichromosome0]->Path[p1]), &(m_Population[ichromosome1]->Path[p1]), iRecombLen );
			memcpy( &(m_Population[ichromosome1]->Path[p1]), pRecomb, iRecombLen );
			SAFE_DELETE_ARRAY(pRecomb);
		}
	}
}

void	CSAGEPAlgorithm::recombinationGene()
{
	if ( CGlobalSpace::m_sUtility.getRandomNumber() < RECOMB_RATE_GENE )
	{
		//randomly choose the chromosome to recombine with.
		int ichromosome0 = CGlobalSpace::m_sUtility.getRandomNumber(POPULATION_NUM-1) + 1;
		int ichromosome1 = CGlobalSpace::m_sUtility.getRandomNumber(POPULATION_NUM-1) + 1;

		if ( ichromosome0 != ichromosome1 )
		{
			//randomly choose the gene to recombine with
			int iGeneNo = CGlobalSpace::m_sUtility.getRandomNumber(GENE_NUM_PER_CHROMOSOME);
			char *pRecombGene = new char[GENE_LEN];
			memcpy( pRecombGene, &(m_Population[ichromosome0]->Path[iGeneNo*GENE_LEN]), GENE_LEN );
			memcpy( &(m_Population[ichromosome0]->Path[iGeneNo*GENE_LEN]), &(m_Population[ichromosome1]->Path[iGeneNo*GENE_LEN]), GENE_LEN );
			memcpy( &(m_Population[ichromosome1]->Path[iGeneNo*GENE_LEN]), pRecombGene, GENE_LEN );
			SAFE_DELETE_ARRAY(pRecombGene);
		}
	}
}

void	CSAGEPAlgorithm::mutation()
{
	for ( int i=1; i<POPULATION_NUM; ++i )
	{
		//the path mutation
		for( int k=0; k<GENE_NUM_PER_CHROMOSOME; ++k )
		{
			for( int j=0; j<GENE_HEAD_LEN; ++j )
			{
				if ( float(rand())/RAND_MAX < MUTATION_RATE_GENE )
				{//the mutation occur
					m_Population[i]->Path[j+k*GENE_LEN] = m_HeadSet[CGlobalSpace::m_sUtility.getRandomNumber(m_HeadSet.length())];
				}
			}
			for( int j=GENE_HEAD_LEN; j<GENE_LEN; ++j )
			{
				if ( float(rand())/RAND_MAX < MUTATION_RATE_GENE )
				{//the mutation occur
					m_Population[i]->Path[j+k*GENE_LEN] = m_TailSet[CGlobalSpace::m_sUtility.getRandomNumber(m_TailSet.length())];
				}
			}
		}
	}
}

void	CSAGEPAlgorithm::transposition()
{
	transposition_IS();
	transposition_RIS();
	if ( GENE_NUM_PER_CHROMOSOME>1 )  transposition_Gene();
}

void	CSAGEPAlgorithm::transposition_IS()
{
	for ( int i=1; i<POPULATION_NUM; ++i )
	{
		//randomly choose the chromosome to transposition
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < TRANS_RATE_IS )
		{
			//IS transposition operation
			for( int k=0; k<GENE_NUM_PER_CHROMOSOME; ++k )
			{
				//randomly choose the transposon start point.
				int iISStartPoint = k*GENE_LEN + CGlobalSpace::m_sUtility.getRandomNumber(GENE_LEN);
				//randomly choose the insert point except for the root position.
				int iISInsertPoint = k*GENE_LEN + CGlobalSpace::m_sUtility.getRandomNumber(GENE_HEAD_LEN-1) + 1;
				//randomly choose the transposon length
				int iTransLen =CGlobalSpace::m_sUtility.getRandomNumber( (k+1)*GENE_LEN-iISStartPoint );

				char *pTransposon = new char[iTransLen]; 
				memcpy( pTransposon, &(m_Population[i]->Path[iISStartPoint]), iTransLen );
				char *pHead = new char[GENE_HEAD_LEN];
				memcpy( pHead, &(m_Population[i]->Path[k*GENE_LEN]), GENE_HEAD_LEN );
				for ( int j=iISInsertPoint; j<(iISInsertPoint+iTransLen); ++j )
				{
					if ( j<(GENE_HEAD_LEN+k*GENE_LEN) )
					{
						m_Population[i]->Path[j] = pTransposon[j-iISInsertPoint];
					}
					else
					{
						break;//only the head part can be inserted.
					}
				}
				//move the gene part to the tail of the inserted part
				if( (k*GENE_LEN+GENE_HEAD_LEN-iISInsertPoint-iTransLen)>0 )
				{
					memcpy( &(m_Population[i]->Path[iISInsertPoint+iTransLen]), pHead+(iISInsertPoint-k*GENE_LEN),  k*GENE_LEN+GENE_HEAD_LEN-iISInsertPoint-iTransLen );
				}
				SAFE_DELETE_ARRAY(pTransposon);
				SAFE_DELETE_ARRAY(pHead);
			}
		}
	}
}

void	CSAGEPAlgorithm::transposition_RIS()
{
	for ( int i=1; i<POPULATION_NUM; ++i )
	{
		//randomly choose the chromosome to transposition
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < TRANS_RATE_RIS )
		{
			//RIS transposition operation
			for( int n=0; n<GENE_NUM_PER_CHROMOSOME; ++n )
			{
				//randomly choose the transposon start point.
				int iRISStartPoint = n*GENE_LEN + CGlobalSpace::m_sUtility.getRandomNumber(GENE_HEAD_LEN);
				//find if the current elem is a function
				bool bFunc = false;
				for ( int k=0; k<(int)(m_FuncSet.length()); ++k )
				{
					if ( m_Population[i]->Path[iRISStartPoint] == m_FuncSet[k] )
					{
						bFunc = true;
					}
				}

				//if the current elem is not a function, we search the head started from this elem
				if ( !bFunc )
				{
					for ( int j=(iRISStartPoint+1); j<(n*GENE_LEN + GENE_HEAD_LEN); ++j )
					{
						for ( int k=0; k<(int)(m_FuncSet.length()); ++k )
						{
							if ( m_Population[i]->Path[j] == m_FuncSet[k] )
							{
								bFunc = true;
								iRISStartPoint = j;
								break;
							}
						}
						if( bFunc )  break;
					}
				}
				//if we can NOT found function elem in the head, we do nothing in this chromosome
				if( !bFunc )  continue;


				//Found the function elem in the start of the transposon.
				//randomly choose the transposon length
				int iTransLen = CGlobalSpace::m_sUtility.getRandomNumber(n*GENE_LEN+GENE_LEN-iRISStartPoint);

				char *pTransposon = new char[iTransLen]; 
				memcpy( pTransposon, &(m_Population[i]->Path[iRISStartPoint]), iTransLen );
				char *pHead = new char[GENE_HEAD_LEN];
				memcpy( pHead, &(m_Population[i]->Path[n*GENE_LEN]), GENE_HEAD_LEN );
				for ( int j=0; j<iTransLen; ++j )
				{
					if ( j<GENE_HEAD_LEN )
					{
						m_Population[i]->Path[j] = pTransposon[j];
					}
					else
					{
						break;//only the head part can be inserted.
					}
				}
				//move the gene part to the tail of the inserted part
				if(  (n*GENE_LEN+GENE_HEAD_LEN-iRISStartPoint-iTransLen)>0 )
				{
					memcpy( &(m_Population[i]->Path[iRISStartPoint+iTransLen]), pHead+(iRISStartPoint-n*GENE_LEN),  n*GENE_LEN+GENE_HEAD_LEN-iRISStartPoint-iTransLen );
				}
				SAFE_DELETE_ARRAY(pTransposon);
				SAFE_DELETE_ARRAY(pHead);
			}
		}
	}
}

void	CSAGEPAlgorithm::transposition_Gene()
{
	for ( int i=1; i<POPULATION_NUM; ++i )
	{
		//randomly choose the chromosome to transposition
		if ( CGlobalSpace::m_sUtility.getRandomNumber() < TRANS_RATE_GENE )
		{
			int iGeneNo = CGlobalSpace::m_sUtility.getRandomNumber(GENE_NUM_PER_CHROMOSOME-1) + 1;
			char *pGene = new char[GENE_LEN];
			memcpy( pGene, &(m_Population[i]->Path[iGeneNo*GENE_LEN]), GENE_LEN );
			for ( int j=(iGeneNo*GENE_LEN-1); j>=0; --j )
			{
				m_Population[i]->Path[j+GENE_LEN] = m_Population[i]->Path[j];
			}
			memcpy( &(m_Population[i]->Path[0]), pGene, GENE_LEN );
			SAFE_DELETE_ARRAY(pGene);
		}
	}
}

}