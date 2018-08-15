#include <ctime>
#include <cuda.h>
#include <cuda_runtime.h>
#include <builtin_types.h>

#include"Common.h"
#include"GlobalSpace.h"
#include "MSA-GA-CUDA.h"

namespace SeqAnsis
{

CMSAGA_CUDA_Algorithm::CMSAGA_CUDA_Algorithm():m_pPopulation(NULL), m_iMaxScore(-INT_MAX), m_iRun(0)
{

}

CMSAGA_CUDA_Algorithm::~CMSAGA_CUDA_Algorithm()
{
	SAFE_DELETE_ARRAY(m_pPopulation);
}

void CMSAGA_CUDA_Algorithm::Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences )
{
	if( vSequences.size()>5 )
	{
		CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA-CUDA: The sequence size been aligned should below or equal 5!" );
		return;
	}

	m_bdbgFirstRun = true;
	//find the longest sequence in the sequence vector
	int iLongestSeqSize = -1;
	m_idbgShortestInputAcidNum = INT_MAX;
	int idxShortestSeq = -1;
	for ( int i=0; i<(int)vSequences.size(); ++i )
	{
		m_iLongestSeqSize = iLongestSeqSize<vSequences[i].getLen()?vSequences[i].getLen():iLongestSeqSize;
		if ( m_idbgShortestInputAcidNum>vSequences[i].getLen() )
		{
			m_idbgShortestInputAcidNum =vSequences[i].getLen();
			idxShortestSeq = i;
		}
	}
	char ch[100];
	sprintf_s(ch,  "shortestSeqAcidNum=%i, the No=%i" , m_idbgShortestInputAcidNum,  idxShortestSeq);
	CGlobalSpace::m_sEventLog.writeEvent(ch);


	//define the column of the sequence matrix
	int nCol = (int)( m_iLongestSeqSize*(MSA_GA_SCALING_FACTOR_K+MSA_GA_OFFSET_FACTOR_X) );

	//population initialization
	InitPopulation( vSequences, nCol );

	//Data transfer from CPU to GPU and evolution on GPU
	Evolution();











	//calculate the fitness of every organism. if not meet the requirement ,then evolution begin.
	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA: Population evolution begin!" );
	m_iRun = 0;
	bool  bStop=false;
	int iCount=0;
	int maxScore=INT_MIN;
	while(!Fitness()  && m_iRun<MSA_GA_GENERATION_NUM && iCount<MSA_GA_NO_ADV_GENERATION_NUM)
	{
		if( m_iMaxScore>maxScore )
		{
			iCount=0;
			maxScore = m_iMaxScore;
		}
		else
		{
			++iCount;
		}
		if ( m_bdbgFirstRun )
		{
			sprintf_s(ch,  "generation=%i, maxscore=%i" , m_iRun, m_iMaxScore );
			CGlobalSpace::m_sEventLog.writeEvent(ch);
			m_bdbgFirstRun = false;
		}
		++m_iRun;
		Evolution();
	}
	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA: Population evolution end!" );

	sprintf_s(ch,  "generation=%i, maxscore=%i" , m_iRun, m_iMaxScore );
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

	if( !CheckSumCharNum( m_pPopulation[idxBestOrganism] ) )
	{
		CGlobalSpace::m_sEventLog.writeEvent("CheckSumCharNum  ERROR!");
	}

	for ( int i=0; i<m_pPopulation[idxBestOrganism].nSeqSize; ++i )
	{
		vAlignedSequences.push_back( m_pPopulation[idxBestOrganism].pSequence[i].sequence );
	}
}

void		CMSAGA_CUDA_Algorithm::InitPopulation( const std::vector<CSequence>& vSequences, int vCol )
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
				tmpSeqData.push_back( StruSeqElem(GENESPACE,-1) );
			}
			for ( int k=0; k<vSequences[j].getLen(); ++k )
			{
				tmpSeqData.push_back( StruSeqElem( vSequences[j].getSequenceContext().at(k).m_char, vSequences[j].getSequenceContext().at(k).m_index ) );
			}
			m_pPopulation[i].pSequence[j].sequence = CSequence( tmpSeqData, vSequences[j].getName(), vSequences[j].getTitle() );
		}
	}
	ArrangeSequences();

	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA-CUDA: Population Initialized!" );
}

int		CMSAGA_CUDA_Algorithm::CountSeqCharNum( const CSequence& vSeq )
{
	int iLen = 0;
	for( int k=0; k<vSeq.getLen(); ++k )
	{
		if ( vSeq.getSequenceContext().at(k).m_iCode!=GENESPACECODE )
		{
			++iLen;
		}
	}
	return iLen;
}

void     CMSAGA_CUDA_Algorithm::ResizeOrgan( COrganism& vOrgan )
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
			if ( itrSeq[k]!=vOrgan.pSequence[k].sequence.getSequenceContext().end() && itrSeq[k]->m_iCode != GENESPACECODE )
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

void		CMSAGA_CUDA_Algorithm::ArrangeSequences()
{
	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
	{
		ResizeOrgan( m_pPopulation[i] );
	}
}

int		CMSAGA_CUDA_Algorithm::GetLongestSeqLen( const COrganism& og )
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

void		CMSAGA_CUDA_Algorithm::Evolution()
{
	//The evolution main code with cuda 

	// Allocate host memory for organs
	unsigned int size_organs = m_pPopulation->nSeqSize*MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	unsigned int mem_size_organs = sizeof(float) * size_organs;
	float *h_organs = (float *)malloc(mem_size_organs);

	//initialize the h_organs with population
	int n=0;
	for( int i=0; i<m_pPopulation->nSeqSize; ++i )
	{
		for ( int j=0; j<MSA_GA_POPULATION_SIZE; ++j )
		{
			const CSeqData& seq = m_pPopulation[j].pSequence[i].sequence.getSequenceContext();
			if( seq.size()>MSA_GA_CUDA_MAX_LEN )
			{
				free(h_organs);
				h_organs  = NULL;
				throw   CAppException( DEF_EXCEPTION_INDEX_OUT_OF_RANGE,DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK
					,__EXCEPTION_SITE__ ,"The max length of the sequence exceed the MSA_GA_CUDA_MAX_LEN" );
			}
			for( int k=0; k<seq.size(); ++k )
			{
				h_organs[n] = float(seq[k].m_iCode);
				++n;
			}
			for( int k=seq.size(); k<MSA_GA_CUDA_MAX_LEN; ++k )
			{
				h_organs[n] = MSA_GA_CUDA_SPACE;
				++n;
			}
		}
	}
	//copy the data to GPU and computing
	float *d_organs;
	cudaError_t error;
	error = cudaMalloc((void **) &d_organs, mem_size_organs);
	if (error != cudaSuccess)
	{
		throw   CAppException( DEF_EXCEPTION_UNEXPECTED,DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK
			,__EXCEPTION_SITE__ ,"cudaMalloc d_organs returned error code!" );
	}

	// copy host memory to device
	error = cudaMemcpy(d_organs, h_organs, mem_size_organs, cudaMemcpyHostToDevice);
	if (error != cudaSuccess)
	{
		throw   CAppException( DEF_EXCEPTION_UNEXPECTED,DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK
			,__EXCEPTION_SITE__ ,"cudaMemcpy (d_organs,h_organs) returned error code!" );
	}

	//fitness
	dim3 thread(m_pPopulation->nSeqSize, MSA_GA_CUDA_MAX_LEN);
	dim3 grid( 1, MSA_GA_POPULATION_SIZE );

	int block_size = 32;
	// Performs warmup operation using matrixMul CUDA kernel
	if (block_size == 16)
	{
		FitnessCUDA<16><<< grid, threads >>>(d_C, d_A, d_B, dimsA.x, dimsB.x);
	}
	else
	{
		FitnessCUDA<32><<< grid, threads >>>(d_C, d_A, d_B, dimsA.x, dimsB.x);
	}


}

};
