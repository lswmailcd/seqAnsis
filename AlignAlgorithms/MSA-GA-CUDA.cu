#include <ctime>
#include <cuda.h>
#include <cuda_runtime.h>
#include <builtin_types.h>
#include <crtdbg.h>

#include"Common.h"
#include"GlobalSpace.h"
#include "Timer.h"
#include"sortingNetworks_common.cuh"
#include "MSA-GA-CUDA.h"

namespace SeqAnsis
{

#define	SubMatPos(r,c)	(r)*(r+1)/2+(c)

CMSAGA_CUDA_Algorithm::CMSAGA_CUDA_Algorithm():m_pPopulation(NULL), m_fMaxScore(-INT_MAX), m_iRun(0), m_pSubMatDevice(NULL)
{

}

CMSAGA_CUDA_Algorithm::~CMSAGA_CUDA_Algorithm()
{
	SAFE_DELETE_ARRAY(m_pPopulation);
}

void	CMSAGA_CUDA_Algorithm::SetAlignParams(int nPopulationNum, int nNoAdvGenerationNum, int nMaxOrgLen, const std::vector<float>& SeqWeight)
{
	int n_seq = (int)SeqWeight.size();
	MSA_GA_CUDA_MAX_LEN = nMaxOrgLen;
	R_MSA_GA_CUDA_MAX_LEN = 1.0f / MSA_GA_CUDA_MAX_LEN;
	MSA_GA_POPULATION_SIZE = nPopulationNum;
	MSA_GA_NO_ADV_GENERATION_NUM = nNoAdvGenerationNum;
	POPULATION_WIDTH = MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	R_POPULATION_WIDTH = 1.0f / POPULATION_WIDTH;
	//Init weight of sequences
	m_SeqWeight.clear();
	m_SeqWeight.resize(SeqWeight.size());

#ifdef WEIGHT_SUM_SP_SCORE
	float sum = 0.0f;
	for (int i = 0; i < (int)SeqWeight.size() - 1; ++i)
	{
		for (int j = i + 1; j < (int)SeqWeight.size(); ++j)
		{
			sum += SeqWeight[i] * SeqWeight[j];
		}
	}
#endif

	float sum = 0.0f;
	for (int i = 0; i < (int)SeqWeight.size(); ++i)
	{
		sum += SeqWeight[i];
	}

	float similarity = 1.0-sum/n_seq;

#ifdef WEIGHT_SUM_SP_SCORE
	for (int i = 0; i < (int)SeqWeight.size(); ++i)
	{
		m_SeqWeight[i] = SeqWeight[i] / sqrt(sum);
	}
#endif

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

#if 0
	sum = 0.0f;
	for (int i = 0; i < (int)SeqWeight.size() - 1; ++i)
	{
		for (int j = i + 1; j < (int)SeqWeight.size(); ++j)
		{
			sum += m_SeqWeight[i] * m_SeqWeight[j];
		}
	}
	sum += 0;
#endif
}

void CMSAGA_CUDA_Algorithm::Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences )
{
	m_OrginSequences = vSequences;
	int n_seq = m_OrginSequences.size();
	m_rnPairs = 2.0f / (n_seq*(n_seq - 1));
	m_bdbgFirstRun = true;
	//find the  SP score of the original sequences
	int score = -1;
	//score = SPScore(vSequences);
	char ch[100];
	sprintf_s(ch, "The input sequences SP score=%i (-1 means not scoring!)", score);
	CGlobalSpace::m_sEventLog.writeEvent(ch);

	//population initialization
	InitPopulation( vSequences );

	//init all data needed for device computing.
	InitDevice();
	
	//calculate the fitness of every organism. if not meet the requirement ,then evolution begin.
	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA-CUDA: Population evolution begin!" );
	CTimer  time;
	double dStartTime = time.getCurrentTime();

////============FOR CPU  VERIFICATION================//
//	short*	score = new short[MSA_GA_POPULATION_SIZE];
//	for ( int i=0; i<MSA_GA_POPULATION_SIZE; ++i )
//	{
//		m_pPopulation[i].score = SPScore( m_pPopulation[i] );
//		score[i] = m_pPopulation[i].score;
//	}
//	std::auto_ptr<SeqAnsis::CFileWriter>		pFileWriter( new SeqAnsis::CFileWriter( "cpuScore.txt" ) );
//	pFileWriter->openFile();
//	pFileWriter->OutputVector( score,  MSA_GA_POPULATION_SIZE);
//	pFileWriter->closeFile();
//	SAFE_DELETE_ARRAY(score);
//	
//	FitnessDeviceUnlimit();
//	dbgWriteFile("m_pOrganScore_Device.txt", m_pOrganScore_Device,MSA_GA_POPULATION_SIZE,MSA_GA_POPULATION_SIZE);
////============FOR CPU  VERIFICATION================//

	m_iRun = 0;
	int iCount = 0;
	float maxScore = INT_MIN;
	while(!FitnessDeviceUnlimit() && m_iRun<MSA_GA_GENERATION_NUM )//&& iCount<MSA_GA_NO_ADV_GENERATION_NUM)
	{
		Evolution();
		if (m_bestScore>maxScore)
		{
			iCount = 0;
			maxScore = m_bestScore;
		}
		else
		{
			++iCount;
		}
		++m_iRun;		
		if (m_bdbgFirstRun)
		{
			sprintf_s(ch, "generation=%i, maxscore=%f", 1, maxScore);
			CGlobalSpace::m_sEventLog.writeEvent(ch);
			m_bdbgFirstRun = false;
		}
	}
	double dInterval = time.getCurrentTime() - dStartTime;
	sprintf_s(ch, "the process time of CUDA-MSA-GA algorithm is %f", dInterval);
	SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent(ch);

	//dbgCheckSequence();
	sprintf_s(ch, "generation=%i, maxscore=%f, iCount=%i", m_iRun, maxScore, iCount);
	CGlobalSpace::m_sEventLog.writeEvent(ch);
	dbgWriteSeqsFromDevice2File( "The final align sequences.txt" );
	//dbgWriteFile("OrganScore_Sorted_Device.txt", m_pOrganScore_Align_Sorted_Device, MSA_GA_POPULATION_SIZE, MSA_GA_POPULATION_SIZE);
	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA-CUDA: Population evolution end!" );
	ReadAlignedSeqsFromDevice( vAlignedSequences );

	ClearDevice();
}

void     CMSAGA_CUDA_Algorithm::ReadAlignedSeqsFromDevice( std::vector<CSequence>& vAlignedSequences )
{
	unsigned int size_organs = m_pPopulation->nSeqSize*MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	unsigned int mem_size_organs = sizeof(short) * size_organs;
	short *h_organs = (short *)malloc(mem_size_organs);
	//read sequences from device to host
	CUDA_DEV2HOST(   h_organs, m_pPopulationDevice[m_curPopulationIndex], mem_size_organs );

	//fill the alignedSequences with best aligned one
	int nWidth = MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	for( int i=0; i<m_pPopulation->nSeqSize; ++i )
	{
		int idx = 0;
		struct StruSeqElem  seqElem;
		CSeqData  seqData;
		for ( int j=0; j<MSA_GA_CUDA_MAX_LEN; ++j )
		{
			if( h_organs[ i*nWidth+j ] != MSA_GA_CUDA_SPACE )
			{
				seqElem.m_iCode = h_organs[ i*nWidth+j ];
				seqElem.m_char = CGlobalSpace::m_sAlignParams.getAminoAcidInt2CharCode( seqElem.m_iCode );
				if( seqElem.m_char != GENESPACE || seqElem.m_char!=NONEGENE )
				{
					seqElem.m_index = idx;
					++idx;
				}
				else
				{
					seqElem.m_index = -1;
				}
				seqData.push_back( seqElem );
			}
			else
			{
				break;
			}
		}
		CSequence seq( seqData,  m_pPopulation[0].pSequence[i].sequence.getName(), m_pPopulation[0].pSequence[i].sequence.getTitle() );
		vAlignedSequences.push_back( seq );
	}
	free(h_organs);
	h_organs=NULL;
}

__device__ inline void devMOD(int M, int N, float revN, int& quotient, int& remainder)
{
	quotient = floor(M*revN);
	remainder = M - quotient*N;
}

template <int BLOCK_SIZE> __global__ void
MSAGA_CALC_SEQ_LEN( bool bRecomb, unsigned int *d_OrganLenNextGen, unsigned int *d_OrganLen, short *d_population, float *d_rand_organ_idx, float fChance, int nOrganWidth, int nPopulationNum)
{	
	// Thread index
	int tx = threadIdx.x;

	if ( tx==0 )
	{
		d_OrganLenNextGen[tx] = d_OrganLen[tx];
		return;
	}

	int  curIdx = -1;
	if (bRecomb)
	{
		if ((nPopulationNum - 1) % 2 != 0 && tx == (nPopulationNum - 1))
		{//总个体为奇数，最后一个个体直接复制
			d_OrganLenNextGen[tx] = d_OrganLen[tx];
			return;
		}
		else
		{
			int  curOrder = tx;
			if (tx % 2 == 0)//如果当前个体是偶数序号个体
			{
				curOrder = tx - 1;//使用相邻奇数个体序号进行计算
			}
			curIdx = curOrder / 2;//half of the population size
		}
	}
	else
	{
		curIdx = tx;
	}

	if (d_rand_organ_idx[curIdx]<fChance)
	{	
		int len=0;
		int i=0;
		while( d_population[tx*nOrganWidth+i]!=MSA_GA_CUDA_SPACE && i<nOrganWidth )
		{
			++len;
			++i;
		}
		d_OrganLenNextGen[tx]=len;
	}
	else
	{
		d_OrganLenNextGen[tx]=d_OrganLen[tx];
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_CALC_SEQ_LEN_UNLIMIT(bool bRecomb, unsigned int *d_OrganLenNextGen, unsigned int *d_OrganLen, short *d_population, float *d_rand_organ_idx, float fChance, int nOrganWidth, int nWidth, int nPopulationNum)
{	
	// Block index
	int bx = blockIdx.x;
	//    int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int orgIdx = bx*blockDim.x + ty*nWidth + tx;
	if (orgIdx >= nPopulationNum)  return;

	if (orgIdx == 0)
	{
		d_OrganLenNextGen[orgIdx] = d_OrganLen[orgIdx];
		return;
	}

	int  curIdx = -1;
	if (bRecomb)
	{
		//result=(MSA_GA_POPULATION_SIZE - 1) % 2
		int j = (nPopulationNum - 1) * 0.5f;
		int result = (nPopulationNum - 1) - j * 2;
		if (result != 0 && orgIdx == (nPopulationNum - 1))
		{//总个体为奇数，最后一个个体直接复制
			d_OrganLenNextGen[orgIdx] = d_OrganLen[orgIdx];
			return;
		}
		else
		{
			int  curOrder = orgIdx;
			//result = orgIdx % 2
			int j = orgIdx * 0.5f;
			int result = orgIdx - j * 2;
			if (result == 0)//如果当前个体是偶数序号个体
			{
				curOrder = orgIdx - 1;//使用相邻奇数个体序号进行计算
			}
			curIdx = curOrder / 2;//half of the population size
		}
	}
	else
	{
		curIdx = orgIdx;
	}

	if (d_rand_organ_idx[curIdx]<fChance)
	{	
		int len=0;
		int i=0;
		while (d_population[orgIdx*nOrganWidth + i] != MSA_GA_CUDA_SPACE && i<nOrganWidth)
		{
			++len;
			++i;
		}
		d_OrganLenNextGen[orgIdx] = len;
	}
	else
	{
		d_OrganLenNextGen[orgIdx]=d_OrganLen[orgIdx];
	}
}

#if 0
__device__ inline void ComparatorShort(
    short &keyA,
    short &valA,
    short &keyB,
    short &valB,
    short dir
)
{
    short t;

    if ((keyA > keyB) == dir)
    {
        t = keyA;
        keyA = keyB;
        keyB = t;
        t = valA;
        valA = valB;
        valB = t;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Monolithic bitonic sort kernel for short arrays fitting into shared memory
////////////////////////////////////////////////////////////////////////////////
__global__ void bitonicSortShared(
    short *d_DstKey,
    short *d_DstVal,
    short *d_SrcKey,
    short *d_SrcVal,
    short arrayLength,
    short dir
)
{
    //Shared memory storage for one or more short vectors
    __shared__ short s_key[SHARED_SIZE_LIMIT];
    __shared__ short s_val[SHARED_SIZE_LIMIT];

    //Offset to the beginning of subbatch and load data
    d_SrcKey += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
    d_SrcVal += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
    d_DstKey += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
    d_DstVal += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
    s_key[threadIdx.x +                       0] = d_SrcKey[                      0];
    s_val[threadIdx.x +                       0] = d_SrcVal[                      0];
    s_key[threadIdx.x + (SHARED_SIZE_LIMIT / 2)] = d_SrcKey[(SHARED_SIZE_LIMIT / 2)];
    s_val[threadIdx.x + (SHARED_SIZE_LIMIT / 2)] = d_SrcVal[(SHARED_SIZE_LIMIT / 2)];
	
    for (uint size = 2; size < arrayLength; size <<= 1)
    {
        //Bitonic merge
        uint ddd = dir ^ ((threadIdx.x & (size / 2)) != 0);

        for (uint stride = size / 2; stride > 0; stride >>= 1)
        {
            __syncthreads();
            uint pos = 2 * threadIdx.x - (threadIdx.x & (stride - 1));
            ComparatorShort(
                s_key[pos +      0], s_val[pos +      0],
                s_key[pos + stride], s_val[pos + stride],
                ddd
            );
        }
    }
	
    //ddd == dir for the last bitonic merge step
    {
        for (uint stride = arrayLength / 2; stride > 0; stride >>= 1)
        {
            __syncthreads();
            uint pos = 2 * threadIdx.x - (threadIdx.x & (stride - 1));
            ComparatorShort(
                s_key[pos +      0], s_val[pos +      0],
                s_key[pos + stride], s_val[pos + stride],
                dir
            );
        }
    }

    __syncthreads();
    d_DstKey[                      0] = s_key[threadIdx.x +                       0];
    d_DstVal[                      0] = s_val[threadIdx.x +                       0];
    d_DstKey[(SHARED_SIZE_LIMIT / 2)] = s_key[threadIdx.x + (SHARED_SIZE_LIMIT / 2)];
    d_DstVal[(SHARED_SIZE_LIMIT / 2)] =  s_val[threadIdx.x + (SHARED_SIZE_LIMIT / 2)];
}
#endif

void		CMSAGA_CUDA_Algorithm::InitPopulation( const std::vector<CSequence>& vSequences )
{
	int iLongestSeqSize = -INT_MAX;
	for (int i = 0; i<(int)vSequences.size(); ++i)
	{
		m_iLongestSeqSize = iLongestSeqSize<vSequences[i].getLen() ? vSequences[i].getLen() : iLongestSeqSize;

	}
	if (m_pPopulation)	SAFE_DELETE_ARRAY(m_pPopulation);
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
	ArrangeSequences();//在CPU上整理序列

	//load the sequence alignment from HOST to DEVICE
	unsigned int size_organs = m_pPopulation->nSeqSize*MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	unsigned int mem_size_organs = sizeof(short) * size_organs;
	short *h_organs = (short *)malloc(mem_size_organs);
	unsigned int mem_size_organLen = sizeof(unsigned int)*MSA_GA_POPULATION_SIZE;
	unsigned int *h_organLen = (unsigned int *)malloc(mem_size_organLen);
	int n=0;
	for ( int j=0; j<MSA_GA_POPULATION_SIZE; ++j )
	{
		const CSeqData& seq = m_pPopulation[j].pSequence[0].sequence.getSequenceContext();
		h_organLen[j] = seq.size();
	}

	for( int i=0; i<m_pPopulation->nSeqSize; ++i )
	{
		for ( int j=0; j<MSA_GA_POPULATION_SIZE; ++j )
		{
			const CSeqData& seq = m_pPopulation[j].pSequence[i].sequence.getSequenceContext();
			if( seq.size()>MSA_GA_CUDA_MAX_LEN )
			{
				free(h_organs);
				h_organs  = NULL;
				free(h_organLen);
				h_organLen=NULL;
				throw   CAppException( DEF_EXCEPTION_INDEX_OUT_OF_RANGE,DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK
					,__EXCEPTION_SITE__ ,"The max length of the sequence exceed the MSA_GA_CUDA_MAX_LEN" );
			}
			for( int k=0; k<seq.size(); ++k )
			{
				h_organs[n] = short(seq[k].m_iCode);
				++n;
			}
			for( int k=seq.size(); k<MSA_GA_CUDA_MAX_LEN; ++k )
			{
				h_organs[n] = MSA_GA_CUDA_SPACE;
				++n;
			}
		}
	}
	//copy the data to Device
	m_curPopulationIndex = 0;
	CUDA_MALLOC( (void**)&m_pPopulationDevice[m_curPopulationIndex],  h_organs, mem_size_organs );
	CUDA_MALLOC( (void**)&m_pPopulationDevice[1-m_curPopulationIndex],  mem_size_organs );
	m_curOrganLenIndex=0;
	CUDA_MALLOC( (void**)&m_pOrganLenDevice[m_curOrganLenIndex], h_organLen,  mem_size_organLen );
	CUDA_MALLOC( (void**)&m_pOrganLenDevice[1-m_curPopulationIndex],  mem_size_organLen );
	free(h_organs);
	h_organs=NULL;
	free(h_organLen);
	h_organLen=NULL;

	//在网格中的有效计算单元的个数
	m_nCellNum = POPULATION_WIDTH*m_pPopulation->nSeqSize;

	//dbgWriteSeqsFromDevice2File("The init align sequences.txt");

	CGlobalSpace::m_sEventLog.writeEvent( "MSA-GA-CUDA: Population Initialized and load to Device!" );
}

int		CMSAGA_CUDA_Algorithm::CountSeqCharNum( const CSequence& vSeq )
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

__global__ void  MSAGA_SPS_UNIT_CUDA_UNLIMIT(float *d_score, short *d_population, unsigned int *d_populationLen, short *d_subMat, float *d_seqWeight, int nSeq, int nWidth, int nPopulationWidth, float rnPopulationWidth, int nOrgWidth, float rnOrgWidth, int iGeneSpace, float rnPairs, float gapOpenCost, float gapExtendCost)
{
	// Block index
	int bx = blockIdx.x;
	//    int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = ty*nWidth + bx*blockDim.x + tx;
#if 0
	int orgIdx = curPos%nPopulationWidth / nOrgWidth;
	if (orgIdx >= MSA_GA_POPULATION_SIZE)  return;
	int seqIdx = curPos / nPopulationWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;
#endif
	int seqIdx, r;
	devMOD(curPos, nPopulationWidth, rnPopulationWidth, seqIdx, r);
	int orgIdx = r*rnOrgWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;

	//if the unit is invalid, we simply set the score of the unit to zero and return.
	if (d_population[curPos] == MSA_GA_CUDA_SPACE)
	{
		d_score[curPos] = 0;
		return;
	}

	//gap panalty,must calculate before the pair score!!!!!
	float panaltyScore = 0;
	//gapOpenCost = gapOpenCost + log10(d_populationLen[orgIdx]);
	if (d_population[curPos] == iGeneSpace)
	{
#if 0
		int nGapPair = 1;
		if (orgIdx < nSeq - 1)
		{
			for (int i = seqIdx + 1; i < nSeq; ++i)
			{
				if (d_population[i*nPopulationWidth + orgIdx*nOrgWidth + posIdx] == iGeneSpace)  ++nGapPair;
			}
		}

		bool bLeftGap=false, bRightGap=false;
		int curPosIdx = posIdx;
		for (int i = 1; i < 8; ++i)
		{
			--curPosIdx;
			if (curPosIdx <= 0) break;
			if (d_population[seqIdx*nPopulationWidth + orgIdx*nOrgWidth + curPosIdx - 1] == iGeneSpace && d_population[seqIdx*nPopulationWidth + orgIdx*nOrgWidth + curPosIdx] != iGeneSpace)
			{
				bLeftGap = true;
				break;
			}
		}
		curPosIdx = posIdx;
		for (int i = 1; i < 8; ++i)
		{
			++curPosIdx;
			if (curPosIdx >= d_populationLen[orgIdx]) break;
			if (d_population[seqIdx*nPopulationWidth + orgIdx*nOrgWidth + curPosIdx + 1] == iGeneSpace && d_population[seqIdx*nPopulationWidth + orgIdx*nOrgWidth + curPosIdx] != iGeneSpace)
			{
				bRightGap = true;
				break;
			}
		}
		if (bLeftGap || bRightGap)  gapExtendCost = gapExtendCost*2.0f;
#endif
		int nGapPair = 1;
		if (posIdx == 0)
		{//the start of the seq is a GAP
			panaltyScore = gapOpenCost;// / nGapPair;
		}
		else
		{
			if (d_population[curPos - 1] == iGeneSpace)
			{//gap extension
				panaltyScore = gapExtendCost;
			}
			else
			{//gap open
				panaltyScore = gapOpenCost;// / nGapPair;
			}
		}
	}
	d_score[curPos] = panaltyScore / nSeq;// *nSeq*(nSeq - 1)*0.5f;

	if (seqIdx == nSeq - 1)
	{//如果是最后一条序列，由于序列对得分(pair score)已经在前面序列中进行了计算，不需要计算最后一个序列的序列对得分，
		return;
	}

	//pair score
	float colScore = 0;
	short	aCode = d_population[curPos];
	for (int i = seqIdx + 1; i<nSeq; ++i)
	{
		int idx = i*nPopulationWidth + orgIdx*nOrgWidth + posIdx;
		short	bCode = d_population[idx];
		if (aCode != iGeneSpace && bCode != iGeneSpace)
		{
			if (aCode<bCode)
			{
				colScore += d_subMat[SubMatPos(bCode, aCode)] * rnPairs;// *d_seqWeight[curPos] * d_seqWeight[idx] * rnPairs;//
			}
			else
			{
				colScore += d_subMat[SubMatPos(aCode, bCode)] * rnPairs;// *d_seqWeight[curPos] * d_seqWeight[idx] * rnPairs;//
			}
		}
	}

	d_score[curPos] += colScore;
}

template <int BLOCK_SIZE> __global__ void
MSAGA_SPS_UNIT_CUDA(short *d_score, short *d_population, short *d_subMat, int nHeight, int nWidth, int blockWidth, int GAPCODE, int gapOpenCost, int gapExtendCost )
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int curPos = ty*nWidth+bx*blockWidth+tx;
	
	//if the unit is invalid, we simply set the score of the unit to zero and return.
	if (d_population[curPos] == MSA_GA_CUDA_SPACE)
	{
		d_score[curPos] = 0;
		return;
	}

	//gap panalty,must calculate before the pair score!!!!!
	short panaltyScore = 0;
	if( d_population[curPos]==GAPCODE )
	{
		if(tx==0)
		{//the start of the seq is a GAP
			panaltyScore = gapOpenCost;
		}
		else
		{
			if( d_population[curPos-1]==GAPCODE )
			{//gap extension
				panaltyScore = gapExtendCost;
			}
			else
			{//gap open
				panaltyScore = gapOpenCost;
			}
		}
	}
	d_score[curPos] = panaltyScore;

	 if(ty==nHeight-1)  
	 {//如果是最后一条序列，由于序列对得分(pair score)已经在前面序列中进行了计算，不需要计算最后一个序列的序列对得分，
		 return;
	 }
	
	//pair score
	short colScore=0;
	short	aCode = d_population[curPos];
	for( int i=ty+1; i<nHeight; ++i )
	{
		short	bCode =  d_population[i*nWidth+bx*blockWidth+tx];
		if( aCode!=GAPCODE && bCode!=GAPCODE )
		{
			if ( aCode<bCode )
			{
				colScore += d_subMat[SubMatPos(bCode, aCode)];
			}
			else
			{
				colScore += d_subMat[SubMatPos(aCode, bCode)];
			}
		}
	}

	d_score[curPos] += colScore;
}

__global__ void MSAGA_SPS_COLSUM_CUDA_UNLIMIT(float *d_score_colsum, float *d_score, int nSeq, int nWidth, int nPopuationNum, int nPopulationWidth, int nOrgWidth, float rnOrgWidth)
{
	// Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = ty*nWidth + bx*blockDim.x + tx;
	int orgIdx = curPos*rnOrgWidth;
	if (orgIdx >= nPopuationNum)  return;
	int posIdx = curPos - orgIdx*nOrgWidth;

	float score_col_sum = 0;

	for (int i = 0; i<nSeq; ++i)
	{
		score_col_sum += d_score[i*nPopulationWidth + orgIdx*nOrgWidth + posIdx];
	}

	d_score_colsum[curPos] = score_col_sum;
}

template <int BLOCK_SIZE> __global__ void
MSAGA_SPS_COLSUM_CUDA(short *d_score_colsum, short *d_score, int nHeight, int nWidth, int blockWidth)
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	short score_col_sum = 0;

	for( int i=0; i<nHeight; ++i )
	{
		score_col_sum += d_score[i*nWidth+bx*blockWidth+tx];
	}

	//__syncthreads();

	d_score_colsum[ty*nWidth+bx*blockWidth+tx] = score_col_sum;
}

__global__ void MSAGA_SPS_ORGAN_CUDA_UNLIMIT(float *d_score_organ, float *d_score_colsum, int nWidth, int nOrgWidth, int nPopulationNum)
{
	// Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = ty*nWidth + bx*blockDim.x + tx;
	if (curPos >= nPopulationNum)  return;

	float score_organ = 0;
	for (int i = 0; i<nOrgWidth; ++i)
	{
		score_organ += d_score_colsum[curPos*nOrgWidth + i];
	}

	d_score_organ[curPos] = score_organ;
}

template <int BLOCK_SIZE> __global__ void
MSAGA_SPS_ORGAN_CUDA(short *d_score_organ, short *d_score_colsum, int nSeqMaxLen, int blockWidth)
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
//    int ty = threadIdx.y;

	short score_organ = 0;
	for( int i=0; i<nSeqMaxLen; ++i )
	{
		score_organ += d_score_colsum[(bx*blockWidth+tx)*nSeqMaxLen+i];
	}

	//__syncthreads();

	d_score_organ[bx*blockWidth+tx] = score_organ;
}

template <int BLOCK_SIZE> __global__ void
MSAGA_ALIGN_CUDA_UNLIMIT(float *d_score_organ_align, uint *d_organ_index_align, float *d_score_organ, int len, int nWidth)
{
	// Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int orgIdx = ty*nWidth + bx*blockDim.x + tx;
	if (orgIdx<len)
	{
		d_score_organ_align[orgIdx] = d_score_organ[orgIdx];
		d_organ_index_align[orgIdx] = orgIdx;
	}
	else
	{
		d_score_organ_align[orgIdx] = -1e4;
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_ALIGN_CUDA(short *d_score_organ_align, uint *d_organ_index_align, short *d_score_organ, int len)
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
//    int ty = threadIdx.y;

	int curPos = bx*blockDim.x+tx;
	if( curPos<len )
	{
		d_score_organ_align[curPos] = d_score_organ[curPos];
		d_organ_index_align[curPos] = curPos;
	}
	else
	{
		d_score_organ_align[curPos] = -1e4;
	}
}

//the arrangement only occur at the changed organ,
//so the gap recorder only record the changed organ.
template <int BLOCK_SIZE> __global__ void
MSAGA_GAP_RECORDER(bool bRecomb, short *d_gapRecorder, short *d_population, float* d_rand_organ_idx, float fChance, int GAPCODE, int nWidth, int nSeq, int nPopulationNum)
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	if ( bx==0 )
	{
		return;
	}

	int  curIdx = -1;
	if (bRecomb)
	{
		//r=(MSA_GA_POPULATION_SIZE - 1)%2
		int q, r;
		devMOD(nPopulationNum - 1, 2, 0.5f, q, r);
		if (r = 0 && (bx == nPopulationNum - 1))
		{//总个体为奇数，最后一个个体不处理
			return;
		}
		else
		{
			int  curOrder = bx;
			//r=bx%2
			int q, r;
			devMOD(bx, 2, 0.5f, q, r);
			if (r == 0)//如果当前个体是偶数序号个体
			{
				curOrder = bx - 1;//使用相邻奇数个体序号进行计算
			}
			curIdx = curOrder / 2;//half of the population size
		}
	}
	else
	{
		curIdx = bx;
	}

	if (d_rand_organ_idx[curIdx]<fChance)//只处理变动了的个体
	{
		bool   bGapColumn = true;
		for( int i=0; i<nSeq; ++i )
		{//对当前列的所有残基进行遍历，检查是否全列都为GAP
			if(d_population[bx*blockDim.x+i*nWidth+tx] != GAPCODE&&d_population[bx*blockDim.x+i*nWidth+tx] != MSA_GA_CUDA_SPACE)
			{
				bGapColumn = false;
			}
		}
		if( bGapColumn )
		{//如果当前列全列都为GAP或是未用空间，则标记此列为0
			d_gapRecorder[bx*blockDim.x+tx] = 0;
		}
		else
		{//如果当前列不全为GAP，则标记此列为1
			d_gapRecorder[bx*blockDim.x+tx] = 1;
		}
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_GAP_RANGE( short *d_gapRange, short *d_GapRecorder, float* d_rand_organ_idx, float fChance  )
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	if( d_rand_organ_idx[bx]<fChance )
	{
		int countGap=0;
		for( int i=0; i<tx; ++i )
		{
			int prePos=bx*blockDim.x + i;
			int curPos=prePos+1;
			if( d_GapRecorder[prePos]==1 && d_GapRecorder[curPos]==0 )
			{
				++countGap;
			}
		}
		d_gapRange[bx*blockDim.x+tx] = countGap;
	}

}

template <int BLOCK_SIZE> __global__ void
MSAGA_GAP_LEN( short *d_GapLenRecorder, short *d_GapRecorder, float* d_rand_organ_idx, float fChance, int nWidth )
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	if( d_rand_organ_idx[bx]<fChance  )
	{
		int gapNum=tx+1;
		int count=0;
		int prePos=0;
		int curPos=prePos+1;
		int gapStartPos=0;
		int gapLen=0;
		while( count<gapNum && curPos<nWidth && d_GapRecorder[curPos] !=MSA_GA_CUDA_SPACE )
		{
			if( d_GapRecorder[prePos]==1 && d_GapRecorder[curPos]==0 )
			{
				gapStartPos = curPos;
			}
			if( d_GapRecorder[prePos]==0 && d_GapRecorder[curPos]==1 )
			{
				gapLen += curPos-gapStartPos;
				++count;
			}
			++prePos;
			++curPos;
		}
		d_GapLenRecorder[bx*blockDim.x+tx]=gapLen;
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_CLEAR_POPULATON_NEXT_GEN( short *d_populationNextGen, float* d_rand_organ_idx, float fChance, int nWidth )
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	if(d_rand_organ_idx[bx]<fChance)//处理变动了的个体
	{
		d_populationNextGen[bx*blockDim.x+ty*nWidth+tx]=MSA_GA_CUDA_SPACE;
	}
}

__global__	void	MSAGA_DELETE_SPACE_AT_TAIL_UNLIMIT(bool bRecomb, short *d_populationNextGen, short *d_population, float* d_rand_organ_idx, float fChance, int nWidth, int nPopulationWidth, float rnPopulationWidth, int nOrgWidth, float rnOrgWidth, int nSeq, int iGeneSpace, int nCellNum, int nPopulationNum)
{
	// Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = bx*blockDim.x + ty*nWidth + tx;
	if(curPos>=nCellNum)  return;
#if 0
	int orgIdx = int(curPos%nPopulationWidth / nOrgWidth);
	if (orgIdx >= MSA_GA_POPULATION_SIZE)   return;
	int seqIdx = int(curPos / nPopulationWidth);
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;
#endif
	int seqIdx = floor(curPos*rnPopulationWidth);
	int i = curPos - seqIdx*nPopulationWidth;
	int orgIdx = i*rnOrgWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;

	if (orgIdx == 0)
	{
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}

	int  curIdx = -1;
	if (bRecomb)
	{
		int q, r;
		devMOD(nPopulationNum - 1, 2, 0.5f, q, r);
		if (r!= 0 && orgIdx == (nPopulationNum - 1))
		{//总个体为奇数，最后一个个体直接复制
			d_populationNextGen[curPos] = d_population[curPos];
			return;
		}
		else
		{
			int  curOrder = orgIdx;
			int q, r;
			devMOD(orgIdx, 2, 0.5f, q, r);
			if (r == 0)//如果当前个体是偶数序号个体
			{
				curOrder = orgIdx - 1;//使用相邻奇数个体序号进行计算
			}
			curIdx = curOrder * 0.5f;//curIdx = curOrder / 2;//half of the population size
		}
	}
	else
	{
		curIdx = orgIdx;
	}

	if (d_rand_organ_idx[curIdx]<fChance)//处理变动了的个体
	{
		for (int i = 0; i < nSeq; i++)
		{
			if( d_population[orgIdx*nOrgWidth + i*nPopulationWidth + posIdx] != iGeneSpace)
			{//如果本列有残基或是MSA_GA_CUDA_SPACE，则直接复制
				d_populationNextGen[curPos] = d_population[curPos];
				return;
			}
		}
		d_populationNextGen[curPos] = MSA_GA_CUDA_SPACE;
	}
	else//没有变动的个体直接复制
	{
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

__global__	void	MSAGA_DELETE_SPACE_UNLIMIT( bool bRecomb, short *d_populationNextGen, short *d_population, float* d_rand_organ_idx, float fChance, short *d_gapRecorder, int nWidth, int nPopulationNum, int nPopulationWidth, float rnPopulationWidth, int nOrgWidth, float rnOrgWidth, int nCellNum )
{
	// Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = bx*blockDim.x+ty*nWidth+tx;
	if (curPos >= nCellNum)  return;
#if 0
	int orgIdx = int(curPos%nPopulationWidth/nOrgWidth);
	if (orgIdx >= MSA_GA_POPULATION_SIZE)   return;
	int seqIdx = int(curPos/nPopulationWidth);
	int posIdx = curPos-seqIdx*nPopulationWidth-orgIdx*nOrgWidth;
#endif
	int seqIdx = floor(curPos*rnPopulationWidth);
	int i = curPos - seqIdx*nPopulationWidth;
	int orgIdx = i*rnOrgWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;

	if ( orgIdx==0 )
	{
		d_populationNextGen[curPos]=d_population[curPos];
		return;
	}

	int  curIdx = -1;
	if (bRecomb)
	{
		int q, r;
		devMOD(nPopulationNum - 1, 2, 0.5f, q, r);
		if (r != 0 && orgIdx == (nPopulationNum - 1))
		{//总个体为奇数，最后一个个体直接复制
			d_populationNextGen[curPos] = d_population[curPos];
			return;
		}
		else
		{
			int  curOrder = orgIdx;
			devMOD(orgIdx, 2, 0.5f, q, r);
			if (r == 0)//如果当前个体是偶数序号个体
			{
				curOrder = orgIdx - 1;//使用相邻奇数个体序号进行计算
			}
			curIdx = curOrder / 2;//half of the population size
		}
	}
	else
	{
		curIdx = orgIdx;
	}

	if (d_rand_organ_idx[curIdx]<fChance)//处理变动了的个体
	{
		if (d_gapRecorder[orgIdx*nOrgWidth + posIdx] != 0)
		{//当前列是有残基的列
			int gapLen=0;
			for( int i=0; i<posIdx; ++i )
			{
				if( d_gapRecorder[orgIdx*nOrgWidth+i]==0 )
				{
					++gapLen;
				}
			}
			//将当前位置的数据移动到空格之前	
			d_populationNextGen[curPos-gapLen]=d_population[curPos];
		}
	}
	else//没有变动的个体直接复制
	{
		d_populationNextGen[curPos]=d_population[curPos];
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_DELETE_SPACE(  short *d_populationNextGen, short *d_population, float* d_rand_organ_idx, float fChance, short *d_gapRecorder, int nWidth )
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int pos=bx*blockDim.x+ty*nWidth+tx;
	if(d_rand_organ_idx[bx]<fChance)//处理变动了的个体
	{
		if( d_gapRecorder[bx*blockDim.x+tx]==1 )//当前列不是空格且不是未用空间
		{
			int gapLen=0;
			for( int i=0; i<tx; ++i )
			{
				if( d_gapRecorder[bx*blockDim.x+i]==0 )
				{
					++gapLen;
				}
			}
			//将当前位置的数据移动到空格之前	
			d_populationNextGen[pos-gapLen]=d_population[pos];
		}
	}
	else//没有变动的个体直接复制
	{
		d_populationNextGen[pos]=d_population[pos];
	}
}

__global__	void		MSAGA_FILL_SPACE_AT_TAIL_UNLIMIT( bool bRecomb, short *d_populationNextGen, short *d_population, float* d_rand_organ_idx, float fChance, int nWidth, short SPACE,  int nOrgWidth, float rnOrgWidth, int nPopulationWidth, float rnPopulationWidth, int nSeq, int nCellNum, int nPopulationNum )
{
	// Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = bx*blockDim.x+ty*nWidth+tx;
	if(curPos>=nCellNum)  return;
#if 0
	int orgIdx = int(curPos%nPopulationWidth/nOrgWidth);
	if (orgIdx >= MSA_GA_POPULATION_SIZE)  return;
	int seqIdx = int(curPos/nPopulationWidth);
	int posIdx = curPos-seqIdx*nPopulationWidth-orgIdx*nOrgWidth;
#endif
	int seqIdx = floor(curPos*rnPopulationWidth);
	int i = curPos - seqIdx*nPopulationWidth;
	int orgIdx = i*rnOrgWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;


	if ( orgIdx==0 )
	{
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}

	short curIdx = -1;
	if (bRecomb)
	{
		int d, r;
		devMOD(nPopulationNum - 1, 2, 0.5f, d, r);
		if (r != 0 && orgIdx == (nPopulationNum - 1))
		{//总个体为奇数，最后一个个体直接复制
			d_populationNextGen[curPos] = d_population[curPos];
			return;
		}
		else
		{
			int  curOrder = orgIdx;
			devMOD(orgIdx, 2, 0.5f, d, r);
			if (r == 0)//如果当前个体是偶数序号个体
			{
				curOrder = orgIdx - 1;//使用相邻奇数个体序号进行计算
			}
			//curIdx = curOrder / 2;//half of the population size
			curIdx = curOrder *0.5f;
		}
	}
	else
	{
		curIdx = orgIdx;
	}

	if (d_rand_organ_idx[curIdx]<fChance)//处理变动了的个体
	{
		if( d_population[curPos]==MSA_GA_CUDA_SPACE )
		{
			bool  bAddSpace = false;
			for( int i=0; i<nSeq; ++i )
			{
				if( d_population[orgIdx*nOrgWidth+i*nPopulationWidth+posIdx]!=MSA_GA_CUDA_SPACE )
				{
					bAddSpace = true;
					break;
				}
			}

			if(bAddSpace)
			{
				d_populationNextGen[curPos] = SPACE;
			}
			else
			{
				d_populationNextGen[curPos] = MSA_GA_CUDA_SPACE;
			}
		}
		else
		{
			d_populationNextGen[curPos] = d_population[curPos];
		}
	}
	else
	{
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_FILL_SPACE_AT_TAIL( short *d_populationNextGen, short *d_population, float* d_rand_organ_idx, float fChance, int nWidth, short SPACE )
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int curPos = bx*blockDim.x+ty*nWidth+tx;
	if(d_rand_organ_idx[bx]<fChance)//处理变动了的个体
	{
		if( d_population[curPos]==MSA_GA_CUDA_SPACE )
		{
			bool  bAddSpace = false;
			for( int i=0; i<blockDim.y; ++i )
			{
				if( d_population[bx*blockDim.x+i*nWidth+tx]!=MSA_GA_CUDA_SPACE )
				{
					bAddSpace = true;
					break;
				}
			}

			if(bAddSpace)
			{
				d_populationNextGen[curPos] = SPACE;
			}
			else
			{
				d_populationNextGen[curPos] = d_population[curPos];
			}
		}
		else
		{
			d_populationNextGen[curPos] = d_population[curPos];
		}
	}
	else
	{
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_DELETE_SPACE_COLUMN( short *d_populationNextGen, short *d_population, short *d_gapRecorder, float *d_rand_organ_idx, float fChance, int nWidth )
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int curPos = bx*blockDim.x+ty*nWidth+tx;
	if(d_rand_organ_idx[bx]<fChance)//处理变动了的个体
	{
		int nGap=0;
		for( int i=0; i<=tx; ++i )
		{//检测[0，tx]区间内全是空位的数目nGap，Tx处由Tx+nGap位置的数据填充.
			if( d_gapRecorder[bx*blockDim.x+i]==0 )  ++nGap; 
		}
		if( tx+nGap<blockDim.x )
		{
			d_populationNextGen[curPos] = d_population[curPos+nGap];
		}
		else
		{//如果tx+nGap所指位置超过了整个序列长度，则用未使用空间占位符填充
			d_populationNextGen[curPos] = MSA_GA_CUDA_SPACE;
		}
	}
	else
	{
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_LOCATE_GAP( short* d_Gap_Tx, short*  d_mark, float* d_rand_Tx_number, float *d_rand_organ_idx,  int seqMaxLen, int nWidth, float fChance, int GAPCODE )
{
	// Block index
//    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
//    int ty = threadIdx.y;

	//only handle the randomly chosen organs.
	int n=0;
	d_Gap_Tx[tx]=-1;
	if(d_rand_organ_idx[tx]<fChance)
	{
		for( int i=0; i<seqMaxLen; ++i )
		{
			if( d_mark[tx*seqMaxLen+i]==1 )
			{
				++n;
			}
		}
		int idx=d_rand_Tx_number[tx]*n;//根据总的连续空位个数，随机选择一个连续空位作为操作对象
		n=-1;
		for( int i=0; i<seqMaxLen; ++i )
		{
			if( d_mark[tx*seqMaxLen+i]==1 )
			{
				++n;
			}
			if( n==idx )//找到应该插入GAP的连续空位的位置
			{
				d_Gap_Tx[tx] = i;
				break;
			}
		}
		
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_LOCATE_GAP_UNLIMIT(short* d_Gap_Tx, short*  d_mark, float* d_rand_Tx_number, float *d_rand_organ_idx, int seqMaxLen, int nWidth, float fChance, int GAPCODE, int nCellNum)
{
	// Block index
	    int bx = blockIdx.x;
	//    int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int orgIdx = bx*blockDim.x + ty*nWidth + tx;
	if (orgIdx >= nCellNum)  return;

	//only handle the randomly chosen organs.
	int n = 0;
	d_Gap_Tx[orgIdx] = -1;
	if (d_rand_organ_idx[orgIdx]<fChance)
	{
		for (int i = 0; i<seqMaxLen; ++i)
		{
			if (d_mark[orgIdx*seqMaxLen + i] == 1)
			{
				++n;
			}
		}
		int idx = d_rand_Tx_number[orgIdx] * n;//根据总的连续空位个数，随机选择一个连续空位作为操作对象
		n = -1;
		for (int i = 0; i<seqMaxLen; ++i)
		{
			if (d_mark[orgIdx*seqMaxLen + i] == 1)
			{
				++n;
			}
			if (n == idx)//找到应该插入GAP的连续空位的位置
			{
				d_Gap_Tx[orgIdx] = i;
				break;
			}
		}

	}
}

__global__	void	MSAGA_MARK_GAP_START_UNLIMIT( short*  d_mark, short* d_population, float *d_rand_organ_idx,  float* d_rand_Ty, int nSeqNum, int nWidth, float fChance, int GAPCODE,  int nOrgWidth, float rnOrgWidth, int nPopulationWidth, float rnPopulationWidth, int nCellNum )
{
	// Block index
    int bx = blockIdx.x;
	
	// Thread index
    int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = bx*blockDim.x+ty*nWidth+tx;
	if( curPos>=nCellNum )  return;
#if 0
	int orgIdx = curPos%nPopulationWidth/nOrgWidth;
	if (orgIdx >= MSA_GA_POPULATION_SIZE)  return;
	int seqIdx = curPos/nPopulationWidth;
	int posIdx = curPos-seqIdx*nPopulationWidth-orgIdx*nOrgWidth;
#endif
	int seqIdx = floor(curPos*rnPopulationWidth);
	int i = curPos - seqIdx*nPopulationWidth;
	int orgIdx = i*rnOrgWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;

	//only handle the randomly chosen organs.
	if(d_rand_organ_idx[orgIdx]<fChance)
	{
		int idxTy = d_rand_Ty[orgIdx]*(nSeqNum-1);//确定该个体中哪个序列进行标记
		int curOrganPos=orgIdx*nOrgWidth+idxTy*nPopulationWidth+posIdx;
		if(posIdx==0)//当前位置在序列头部
		{
			if( d_population[curOrganPos]==GAPCODE )
			{
				d_mark[orgIdx*nOrgWidth + posIdx] = 1;
			}
			else
			{
				d_mark[orgIdx*nOrgWidth + posIdx] = 0;
			}
		}
		else//当前位置不在序列头部时
		{
			if(  d_population[curOrganPos]==GAPCODE && d_population[curOrganPos-1]!=GAPCODE )//当前位置为GAP，而前一位置不是GAP，即确定当前位置为空位开始位置
			{
				d_mark[orgIdx*nOrgWidth + posIdx] = 1;
			}
			else
			{
				d_mark[orgIdx*nOrgWidth + posIdx] = 0;
			}
		}
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_MARK_GAP_START( short*  d_mark, short* d_population, float *d_rand_organ_idx,  float* d_rand_Ty, int nSeqNum, int nWidth, float fChance, int GAPCODE )
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
//    int ty = threadIdx.y;

	//only handle the randomly chosen organs.
	if(d_rand_organ_idx[bx]<fChance)
	{
		int curPos = bx*blockDim.x+tx;
		int idxTy = d_rand_Ty[bx]*(nSeqNum-1);//确定该个体中哪个序列进行标记
		int curOrganPos=bx*blockDim.x+idxTy*nWidth+tx;
		if(tx==0)//当前位置在序列头部
		{
			if( d_population[curOrganPos]==GAPCODE )
			{
				d_mark[curPos] = 1;
			}
			else
			{
				d_mark[curPos] = 0;
			}
		}
		else//当前位置不在序列头部时
		{
			if(  d_population[curOrganPos]==GAPCODE && d_population[curOrganPos-1]!=GAPCODE )//当前位置为GAP，而前一位置不是GAP，即确定当前位置为空位开始位置
			{
				d_mark[curPos] = 1;
			}
			else
			{
				d_mark[curPos] = 0;
			}
		}
	}
}

__global__ void  MSAGA_RECOMB_VERTICAL_UNLIMIT(  short *d_populationNextGen, short *d_population, short* d_pOrgIdx0, short* d_pOrgIdx1, short* d_pRecombPos0, short* d_pRecombPos1, float *d_rand_organ_idx, float  fRecombVertical, short *d_Pos1_min, short *d_Pos1_max, int geneSpace, int nOrgWidth, float rnOrgWidth, int nPopulationWidth, float rnPopulationWidth, int nPopulation, int nWidth, int nCellNum  )
{
	// Block index    
	int bx = blockIdx.x;
	//    int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = ty*nWidth+bx*blockDim.x+tx;
	if(curPos>=nCellNum)  return;
	int seqIdx = floor(curPos*rnPopulationWidth);
	int i = curPos - seqIdx*nPopulationWidth;
	int orgIdx = i*rnOrgWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;
#if 0
	int orgIdx = curPos%nPopulationWidth/nOrgWidth;
	if( orgIdx>=nPopulation )  return;
	int seqIdx = curPos/nPopulationWidth;
	int posIdx = curPos-seqIdx*nPopulationWidth-orgIdx*nOrgWidth;
#endif

	if(orgIdx==0) 
	{//最优个体直接进入下一代
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}	

	int r, d;
	devMOD(nPopulation-1, 2, 0.5f, d, r);
	if (r != 0 && orgIdx == nPopulation - 1)//如果待处理个体总数为奇数，则最后一个个体直接进入下一代
	{
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}

	int  curOrder=orgIdx;
	devMOD(orgIdx, 2, 0.5f, d, r);
	if( r==0 )//如果当前个体是偶数序号个体
	{
		curOrder=orgIdx-1;//使用相邻奇数个体序号进行计算
	}

	int curIdx = curOrder/2;//half of the population size=nRecomb
	//d_pBx[curIdx]=0;
	if(d_rand_organ_idx[curIdx]<fRecombVertical)
	{
		//d_pBx[curIdx]=1;
		int nRecomb = (nPopulation-1)/2;
		if( orgIdx%2==0 )//如果当前个体是偶数序号个体
		{
			//d_populationNextGen[curPos] = d_population[curPos];
			if( posIdx<=d_pRecombPos1[curIdx+seqIdx*nRecomb] )
			{//copy directly from organ1,使用P1中相应位置数据
				d_populationNextGen[curPos] = d_population[d_pOrgIdx1[curIdx]*nOrgWidth+seqIdx*nPopulationWidth+posIdx];
			}
			else if( posIdx>d_Pos1_max[curIdx] )
			{
				int posIdx1=d_pRecombPos0[curIdx]+posIdx-d_Pos1_max[curIdx];
				int copyPos = d_pOrgIdx0[curIdx]*nOrgWidth+seqIdx*nPopulationWidth+posIdx1;
				if( posIdx1<nOrgWidth )
				{
					d_populationNextGen[curPos] = d_population[copyPos];
				}
				else
				{//如果P0所有数据均已经复制，则后面的空间全都应该标记为未使用
					d_populationNextGen[curPos] = MSA_GA_CUDA_SPACE;
				}
			}
			else
			{//insert gap,填入SPACE1个空格
				d_populationNextGen[curPos] = geneSpace;
			}
		}
		else//如果当前个体是奇数序号个体
		{
			//d_populationNextGen[curPos] = d_population[curPos];
			if( posIdx<=d_pRecombPos0[curIdx] )
			{//copy directly from organ0,直接复制个体P0相应位置的数据
				d_populationNextGen[curPos] = d_population[d_pOrgIdx0[curIdx]*nOrgWidth+seqIdx*nPopulationWidth+posIdx];
			}
			else if( posIdx>(d_pRecombPos0[curIdx]+d_pRecombPos1[curIdx+seqIdx*nRecomb]-d_Pos1_min[curIdx]) )
			{
				int posIdx1=d_pRecombPos1[curIdx+seqIdx*nRecomb]+posIdx-(d_pRecombPos0[curIdx]+d_pRecombPos1[curIdx+seqIdx*nRecomb]-d_Pos1_min[curIdx]);
				int copyPos = d_pOrgIdx1[curIdx]*nOrgWidth+seqIdx*nPopulationWidth+posIdx1;
				if( posIdx1<nOrgWidth )
				{	
					d_populationNextGen[curPos] =d_population[copyPos];
				}
				else
				{
					d_populationNextGen[curPos] = MSA_GA_CUDA_SPACE;
				}
			}
			else
			{//insert gap，填入recombPos1-min(recombPos1)个空格
				d_populationNextGen[curPos] = geneSpace;
			}
		}
	}
	else//不满足垂直杂交条件的个体直接复制
	{
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_RECOMB_VERTICAL(  short *d_populationNextGen, short *d_population, short* d_pOrgIdx0, short* d_pOrgIdx1, short* d_pRecombPos0, short* d_pRecombPos1, float *d_rand_organ_idx, float  fRecombVertical, int nWidth,  int nPopulation, int organNum, short *d_Pos1_min, short *d_Pos1_max, int GENESPACE, short NONE )
{
	// Block index    
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int curPos = ty*nWidth+bx*blockDim.x+tx;
	if(bx==0) 
	{//最优个体直接进入下一代
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}	

	//!!!!!!!!!!!!!!!!!!!!!!!!!需要区分奇偶个体，分别处理！！！！！//	
	if( (nPopulation-1)%2!=0 && bx==nPopulation-1 )//如果待处理个体总数为奇数，则最后一个个体直接进入下一代
	{
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}

	int  curOrder=bx;
	if( bx%2==0 )//如果当前个体是偶数序号个体
	{
		curOrder=bx-1;//使用相邻奇数个体序号进行计算
	}

	int curIdx = curOrder/2;//half of the population size=nRecomb
	//d_pBx[curIdx]=0;
	if(d_rand_organ_idx[curIdx]<fRecombVertical)
	{
		//d_pBx[curIdx]=1;
		int nRecomb = (nPopulation-1)/2;
		if( bx%2==0 )//如果当前个体是偶数序号个体
		{
			//d_populationNextGen[curPos] = d_population[curPos];
			if( tx<=d_pRecombPos1[curIdx+ty*nRecomb] )
			{//copy directly from organ1,使用P1中相应位置数据
				d_populationNextGen[curPos] = d_population[d_pOrgIdx1[curIdx]*blockDim.x+ty*nWidth+tx];
			}
			else if( tx>d_Pos1_max[curIdx] )
			{
				int tx1=d_pRecombPos0[curIdx]+tx-d_Pos1_max[curIdx];
				int copyPos = d_pOrgIdx0[curIdx]*blockDim.x+ty*nWidth+tx1;
				if( tx1<blockDim.x )
				{
						d_populationNextGen[curPos] = d_population[copyPos];
				}
				else
				{//如果P0所有数据均已经复制，则后面的空间全都应该标记为未使用
						d_populationNextGen[curPos] = GENESPACE;
				}
			}
			else
			{//insert gap,填入SPACE1个空格
				d_populationNextGen[curPos] = GENESPACE;
			}
		}
		else//如果当前个体是奇数序号个体
		{
			//d_populationNextGen[curPos] = d_population[curPos];
			if( tx<=d_pRecombPos0[curIdx] )
			{//copy directly from organ0,直接复制个体P0相应位置的数据
				d_populationNextGen[curPos] = d_population[d_pOrgIdx0[curIdx]*blockDim.x+ty*nWidth+tx];
			}
			else if( tx>(d_pRecombPos0[curIdx]+d_pRecombPos1[curIdx+ty*nRecomb]-d_Pos1_min[curIdx]) )
			{
				int tx1=d_pRecombPos1[curIdx+ty*nRecomb]+tx-(d_pRecombPos0[curIdx]+d_pRecombPos1[curIdx+ty*nRecomb]-d_Pos1_min[curIdx]);
				int copyPos = d_pOrgIdx1[curIdx]*blockDim.x+ty*nWidth+tx1;
				if( tx1<blockDim.x )
				{	
					d_populationNextGen[curPos] =d_population[copyPos];
				}
				else
				{
					d_populationNextGen[curPos] = GENESPACE;
				}
			}
			else
			{//insert gap，填入recombPos1-min(recombPos1)个空格
				d_populationNextGen[curPos] = GENESPACE;
			}
		}
	}
	else//不满足垂直杂交条件的个体直接复制
	{
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_RECOMB_VERTICAL_SPACENUM( short *d_SpaceOrg0, short *d_SpaceOrg1, short *d_pOffset_max, short *d_pOffset_min, short *d_pOffset, int nSeqNum, int nWidth )
{
	// Block index
//    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	d_SpaceOrg0[ty*nWidth+tx] = abs(d_pOffset_min[tx])+d_pOffset[ty*nWidth+tx];
	d_SpaceOrg1[ty*nWidth+tx] = d_pOffset_max[tx]-d_pOffset[ty*nWidth+tx];
}

template <int BLOCK_SIZE> __global__ void
MSAGA_RECOMB_VERTICAL_POS1_MIN_MAX(short *d_pPos1_min, short *d_pPos1_max, short *d_pRecombPos1, int nSeqNum, float fRecombVertical, float *d_rand_organ_idx, int nWidth)
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
//    int ty = threadIdx.y;

	d_pPos1_min[tx] = -10000;
	d_pPos1_max[tx] = 10000;

	if(d_rand_organ_idx[tx]<fRecombVertical )
	{
		short max = -10000;
		short min = 10000;

		for( int i=0; i<nSeqNum; ++i )
		{
			if( d_pRecombPos1[i*nWidth+tx]>max ) 
			{
				max =  d_pRecombPos1[i*nWidth+tx];
			}
			if( d_pRecombPos1[i*nWidth+tx]<min )
			{
				min =  d_pRecombPos1[i*nWidth+tx];
			}
		}

		d_pPos1_min[tx] = min;
		d_pPos1_max[tx] = max;
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_RECOMB_VERTICAL_POS1_MIN_MAX_UNLIMIT(short *d_pPos1_min, short *d_pPos1_max, short *d_pRecombPos1, int nSeqNum, float fRecombVertical, float *d_rand_organ_idx, int nWidth, int nRecomb)
{
	// Block index
	int bx = blockIdx.x;
	//    int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = bx*blockDim.x + ty*nWidth + tx;
	if (curPos >= nRecomb)   return;

	d_pPos1_min[curPos] = -10000;
	d_pPos1_max[curPos] = 10000;

	if (d_rand_organ_idx[curPos]<fRecombVertical)
	{
		short max = -10000;
		short min = 10000;

		for (int i = 0; i<nSeqNum; ++i)
		{
			if (d_pRecombPos1[i*nRecomb + curPos]>max)
			{
				max = d_pRecombPos1[i*nRecomb + curPos];
			}
			if (d_pRecombPos1[i*nRecomb + curPos]<min)
			{
				min = d_pRecombPos1[i*nRecomb + curPos];
			}
		}

		d_pPos1_min[curPos] = min;
		d_pPos1_max[curPos] = max;
	}
}


__global__ void	MSAGA_RECOMB_VERTICAL_RECOMB_POS1_UNLIMIT(short* d_pRecombPos0, short* d_pRecombPos1, short* d_pOrgIdx0, short* d_pOrgIdx1, float *d_rand_organ_idx, float fRecombVertical, int nWidth, short* d_population, int nPopulationWidth, int nOrgWidth, int nRecomb, float rnRecomb, short iGeneSpace, int nSeq)
{
	//block index
	int bx = blockIdx.x;
	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = bx*blockDim.x + ty*nWidth+tx;
	if (curPos >= nRecomb*nSeq)  return;

	int curSeq = curPos*rnRecomb;//当前位置所对应的个体中的当前序列
	int curOrg = curPos-curSeq*nRecomb;//当前位置所对应的个体编号
	d_pRecombPos1[curPos] = -2;
	int n = 0;
	if (d_rand_organ_idx[curOrg] < fRecombVertical)
	{
		//统计个体每个序列从开头到杂交点之前的残基个数
		for (int i = 0; i <= d_pRecombPos0[curOrg]; i++)
		{
			if (d_population[d_pOrgIdx0[curOrg] * nOrgWidth + curSeq*nPopulationWidth + i] != iGeneSpace)
			{
				++n;
			}
		}

		//查找第二个个体相应杂交点位置
		short j = 0;
		short pos1 = 0;
		if (d_population[d_pOrgIdx0[curOrg] * nOrgWidth + curSeq*nPopulationWidth + 0] == iGeneSpace && n == 0)
		{//从第一个个体到杂交点处没有残基，全为空格，则posOrg1为第二个个体中第一个残基前一位置
			while (d_population[d_pOrgIdx1[curOrg] * nOrgWidth + curSeq*nPopulationWidth + pos1] == iGeneSpace)
			{
				++pos1;
			}
		}
		else
		{
			while (j < n)
			{
				if (d_population[d_pOrgIdx1[curOrg] * nOrgWidth + curSeq*nPopulationWidth + pos1] != iGeneSpace)
				{
					++j;
				}
				++pos1;
			}
		}
		d_pRecombPos1[curPos] = pos1 - 1;
	}
}

__global__ void	MSAGA_RECOMB_VERTICAL_RECOMB_POS0_UNLIMIT(short* d_pRecombPos0, short* d_pOrgIdx0, short* d_pOrgIdx1, float *d_rand_organ_idx, float* d_rand_Tx, float* d_rand_org0, float* d_rand_org1, float fRecombVertical, unsigned int *d_organLen, int nWidth, int nRecomb, int nPopulationNum)
{
	//Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = bx*blockDim.x + ty*nWidth + tx;
	if (curPos >= nRecomb)	return;

	d_pRecombPos0[curPos] = -1;
	if (d_rand_organ_idx[curPos] < fRecombVertical)//每个curPos对应一对奇偶序号的个体
	{
		//随机生成一对进行垂直杂交的个体
		d_pOrgIdx0[curPos] = short(d_rand_org0[curPos] * (nPopulationNum - 1));
		d_pOrgIdx1[curPos] = short(d_rand_org1[curPos] * (nPopulationNum - 1));

		//随机生成一个杂交点
		d_pRecombPos0[curPos] = d_rand_Tx[curPos] * (d_organLen[d_pOrgIdx0[curPos]] - 1);
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_RECOMB_VERTICAL_RECOMB_POS(   short* d_pRecombPos0, short* d_pRecombPos1, short* d_pOrgIdx0, short* d_pOrgIdx1, short *d_population, float *d_rand_organ_idx, float* d_rand_Tx,   float* d_rand_org0, float* d_rand_org1, float fRecombVertical, int nPopulation, int nBlockWidth, int iGeneSpace, unsigned int *d_organLen, int nWidth, int nRecomb)
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int curPos = bx*blockDim.x+ty*nWidth+tx;
	int curPopulation = curPos%nRecomb;//当前位置所对应的个体编号
	int curSeqNo = curPos/nRecomb;//当前位置所对应的个体中的当前序列

	d_pOrgIdx0[curPopulation]=-1;
	d_pOrgIdx1[curPopulation]=-1;
	d_pRecombPos0[curSeqNo*nRecomb+curPopulation]=-100;
	d_pRecombPos1[curSeqNo*nRecomb+curPopulation]=-100;
	
	if(d_rand_organ_idx[curPopulation]<fRecombVertical)//每个tx对应一对奇偶序号的个体
	{
		d_pOrgIdx0[curPopulation] = d_rand_org0[curPopulation]*(nPopulation-1);//随机生成一对进行垂直杂交的个体
		d_pOrgIdx1[curPopulation] = d_rand_org1[curPopulation]*(nPopulation-1);
		int posOrg0 = d_rand_Tx[curPopulation]*(d_organLen[curPopulation]-1);//随机生成一个杂交点
		//统计从第一个个体开头到杂交点之前的残基个数
		int n=0;
		for( int i=0; i<=posOrg0; ++i )
		{
			if( d_population[d_pOrgIdx0[curPopulation]*nBlockWidth+curSeqNo*(nPopulation*nBlockWidth)+i]!=iGeneSpace )
			{
				++n;
			}
		}

		//查找第二个个体相应杂交点位置
		int j=0;
		int posOrg1=0;
		if( d_population[d_pOrgIdx0[curPopulation]*nBlockWidth+curSeqNo*(nPopulation*nBlockWidth)+0]==iGeneSpace && n==0 )//从第一个个体到杂交点处没有残基，全为空格，则posOrg1为第二个个体中第一个残基前一位置
		{
			while( d_population[d_pOrgIdx1[curPopulation]*nBlockWidth+curSeqNo*(nPopulation*nBlockWidth)+j]==iGeneSpace  )
			{
				++j;
			}
			posOrg1=j-1;
		}
		else
		{
			while( j<n )
			{
				if( d_population[d_pOrgIdx1[curPopulation]*nBlockWidth+curSeqNo*(nPopulation*nBlockWidth)+posOrg1]!=iGeneSpace )
				{
					++j;
				}
				++posOrg1;
			}
			if(posOrg1!=0)
			{
				--posOrg1;
			}
		}
		
		d_pRecombPos0[curSeqNo*nRecomb+curPopulation] = posOrg0;	
		d_pRecombPos1[curSeqNo*nRecomb+curPopulation] = posOrg1;
	}
}

__global__	void	MSAGA_RECOMB_HORIZEN_UNLIMIT(short *d_populationNextGen, short *d_population, float *d_rand_organ_idx, float* d_rand_Ty, float* d_rand_org0, float* d_rand_org1, float fRecombHorizen, int nOrgWidth, float rnOrgWidth, int nPopulationWidth, float rnPopulationWidth, int nWidth, int nPopulation, int nOrgNum, int nCellNum)
{
	// Block index
	int bx = blockIdx.x;
	//    int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = bx*blockDim.x+ty*nWidth+tx;
	if (curPos >= nCellNum)   return;
#if 0
	int orgIdx = (curPos%nPopulationWidth)/nOrgWidth;
	int seqIdx = curPos/nPopulationWidth;
	int posIdx = curPos-seqIdx*nPopulationWidth-orgIdx*nOrgWidth;
#endif
	int seqIdx = floor(curPos*rnPopulationWidth);
	int i = curPos - seqIdx*nPopulationWidth;
	int orgIdx = i*rnOrgWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;

	if(orgIdx==0) 
	{//最优个体直接进入下一代
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}		

	//result = 

	if( (nPopulation-1)%2!=0 && orgIdx==nPopulation-1 )//如果待处理个体总数为奇数，则最后一个个体直接复制
	{
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}

	int  curOrder=orgIdx;
	if( orgIdx%2==0 )//如果当前个体是偶数序号个体
	{
		curOrder=orgIdx-1;//使用相邻奇数个体序号进行计算
	}

	int curIdx = curOrder/2;//half of the population size
	if(d_rand_organ_idx[curIdx]<fRecombHorizen)//对满足水平杂交条件的个体
	{
		int bxOrg0 = d_rand_org0[curIdx]*(nPopulation-1);
		int bxOrg1 = d_rand_org1[curIdx]*(nPopulation-1);
		int selectSeqIdx = d_rand_Ty[curIdx]*(nOrgNum-1);

		if( orgIdx%2==0 )//如果当前个体是偶数序号个体
		{
			if(seqIdx==selectSeqIdx)
			{//exchange the cur Seq
				d_populationNextGen[curPos] = d_population[bxOrg0*nOrgWidth+selectSeqIdx*nPopulationWidth+posIdx];
			}
			else
			{
				d_populationNextGen[curPos] = d_population[bxOrg1*nOrgWidth+seqIdx*nPopulationWidth+posIdx];
			}
		}
		else//如果当前个体是奇数序号个体
		{
			if(seqIdx==selectSeqIdx)
			{//exchange the cur Seq
				d_populationNextGen[curPos] = d_population[bxOrg1*nOrgWidth+selectSeqIdx*nPopulationWidth+posIdx];
			}
			else
			{
				d_populationNextGen[curPos] = d_population[bxOrg0*nOrgWidth+seqIdx*nPopulationWidth+posIdx];
			}
		}
	}
	else
	{//不满足杂交条件的个体直接拷贝
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_RECOMB_HORIZEN(   short *d_populationNextGen, short *d_population, float *d_rand_organ_idx, float* d_rand_Ty, float* d_rand_org0, float* d_rand_org1, float fRecombHorizen, int organNum, int nPopulation )
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int nBlockWidth = blockDim.x;
	int curPos=bx*nBlockWidth+ty*(nPopulation*nBlockWidth)+tx;
	if(bx==0) 
	{//最优个体直接进入下一代
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}	

	if( (nPopulation-1)%2!=0 && bx==nPopulation-1 )//如果待处理个体总数为奇数，则最后一个个体直接复制
	{
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}

	int  curOrder=bx;
	if( bx%2==0 )//如果当前个体是偶数序号个体
	{
		curOrder=bx-1;//使用相邻奇数个体序号进行计算
	}

	int curIdx = curOrder/2;//half of the population size
	if(d_rand_organ_idx[curIdx]<fRecombHorizen)//对满足水平杂交条件的个体
	{
		int bxOrg0 = d_rand_org0[curIdx]*(nPopulation-1);
		int bxOrg1 = d_rand_org1[curIdx]*(nPopulation-1);
		int selectTy = d_rand_Ty[curIdx]*(organNum-1);

		if( bx%2==0 )//如果当前个体是偶数序号个体
		{
			if(ty==selectTy)
			{//exchange the cur Seq
				d_populationNextGen[curPos] = d_population[bxOrg0*nBlockWidth+selectTy*(nPopulation*nBlockWidth)+tx];
			}
			else
			{
				d_populationNextGen[curPos] = d_population[bxOrg1*nBlockWidth+ty*(nPopulation*nBlockWidth)+tx];
			}
		}
		else//如果当前个体是奇数序号个体
		{
			if(ty==selectTy)
			{//exchange the cur Seq
				d_populationNextGen[curPos] = d_population[bxOrg1*nBlockWidth+selectTy*(nPopulation*nBlockWidth)+tx];
			}
			else
			{
				d_populationNextGen[curPos] = d_population[bxOrg0*nBlockWidth+ty*(nPopulation*nBlockWidth)+tx];
			}
		}
	}
	else
	{//不满足杂交条件的个体直接拷贝
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

__global__ void MSAGA_GAP_REDUCE_UNLIMIT(  short *d_populationNextGen, short *d_population, float *d_rand_organ_idx,  float* d_rand_Ty, short* d_gapPos_Tx, int nSeqNum, int nWidth, float fReduceGapChance, int GAPCODE, int nOrgWidth, float rnOrgWidth, int nPopulationWidth, float rnPopulationWidth, int nCellNum )
{
	// Block index
	int bx = blockIdx.x;
	//    int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = bx*blockDim.x+ty*nWidth+tx;
	if(curPos>=nCellNum)  return;
#if 0
	int nPopulationWidth = nOrgWidth*nPopulation;
	int orgIdx = curPos%nPopulationWidth/nOrgWidth;
	if (curPos >= nPopulation)  return;
	int seqIdx = curPos/nPopulationWidth;
	int posIdx = curPos-seqIdx*nPopulationWidth-orgIdx*nOrgWidth;
#endif
	int seqIdx = floor(curPos*rnPopulationWidth);
	int i = curPos - seqIdx*nPopulationWidth;
	int orgIdx = i*rnOrgWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;

	if(orgIdx==0) 
	{//最优个体直接进入下一代
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}	

	if(d_rand_organ_idx[orgIdx]<fReduceGapChance)
	{
		int reduceTy = d_rand_Ty[orgIdx]*(nSeqNum-1);
		if( seqIdx==reduceTy )
		{
			int reduceTx = d_gapPos_Tx[orgIdx];
			if( reduceTx!=-1 )//当有空格可用于减少时
			{
				if( posIdx<reduceTx )
				{//direct copy
					d_populationNextGen[curPos] = d_population[curPos];
				}
				else if( posIdx<(nOrgWidth-1) )
				{
					d_populationNextGen[curPos] = d_population[curPos+1];
				}
				else//最后一个数据点直接填为未使用的空间记号
				{
					d_populationNextGen[curPos] = MSA_GA_CUDA_SPACE;
				}
			}
			else
			{
				d_populationNextGen[curPos] = d_population[curPos];
			}
		}
		else
		{//just copy to next generation
			d_populationNextGen[curPos] = d_population[curPos];
		}
	}
	else
	{//just copy to next generation
		d_populationNextGen[curPos] = d_population[curPos];
	}

}

template <int BLOCK_SIZE> __global__ void
MSAGA_GAP_REDUCE(  short *d_populationNextGen, short *d_population, float *d_rand_organ_idx,  float* d_rand_Ty, short* d_gapPos_Tx, int nSeqNum, int nWidth, float fReduceGapChance, int GAPCODE )
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int curPos=bx*blockDim.x+ty*nWidth+tx;
	if(bx==0) 
	{//最优个体直接进入下一代
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}	

	if(d_rand_organ_idx[bx]<fReduceGapChance)
	{
		int reduceTy = d_rand_Ty[bx]*(nSeqNum-1);
		if( ty==reduceTy )
		{
			int reduceTx = d_gapPos_Tx[bx];
			if( reduceTx!=-1 )//当有空格可用于减少时
			{
				if( tx<reduceTx )
				{//direct copy
					d_populationNextGen[curPos] = d_population[curPos];
				}
				else if( tx<(blockDim.x-1) )
				{
					d_populationNextGen[curPos] = d_population[curPos+1];
				}
				else//最后一个数据点直接填为未使用的空间记号
				{
					d_populationNextGen[curPos] = MSA_GA_CUDA_SPACE;
				}
			}
			else
			{
				d_populationNextGen[curPos] = d_population[curPos];
			}
		}
		else
		{//just copy to next generation
			d_populationNextGen[curPos] = d_population[curPos];
		}
	}
	else
	{//just copy to next generation
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

__global__ void MSAGA_GAP_INSERT_UNLIMIT(short *d_populationNextGen, short *d_population, float *d_rand_organ_idx, float* d_rand_Ty, float* d_rand_Tx, int populationNum, int nSeq, int nBlockWidth, float fOpenGapChance, int GAPCODE, int nWidth, int nOrgWidth, float rnOrgWidth, int nPopulationWidth, float rnPopulationWidth, unsigned int* d_OrgLen, int nCellNum)
{
	// Block index
	int bx = blockIdx.x;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curPos = bx*blockDim.x + ty*nWidth + tx;
	if (curPos >= nCellNum)  return;
#if 0
	int orgIdx = curPos%nPopulationWidth/nOrgWidth;
	int seqIdx = curPos/nPopulationWidth;
	int posIdx = curPos-seqIdx*nPopulationWidth-orgIdx*nOrgWidth;
	int curPopulationPos = orgIdx*nOrgWidth+seqIdx*nPopulationWidth+posIdx;
#endif
	int seqIdx = floor(curPos*rnPopulationWidth);
	int i = curPos - seqIdx*nPopulationWidth;
	int orgIdx = i*rnOrgWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;

	if(orgIdx==0) 
	{//最优个体直接进入下一代
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}	

	if(d_rand_organ_idx[orgIdx]<fOpenGapChance)
	{//open gap
		int insertSeqIdx = d_rand_Ty[orgIdx]*(nSeq-1);
		if( seqIdx==insertSeqIdx )//当前序列为需要插入GAP的序列
		{
			int insertPosIdx = d_rand_Tx[orgIdx]*(d_OrgLen[orgIdx]-1);
			if(  insertPosIdx==0 )//当插入位置在序列头部
			{
				if( posIdx==insertPosIdx )
				{//直接在下一代此序列头部插入GAP
					d_populationNextGen[curPos] = GAPCODE;
				}
				else//在下一代中，原序列其它残基依次向后一位
				{
					d_populationNextGen[curPos] = d_population[curPos - 1];
				}
			}
			else//当插入位置不在序列头部
			{
				if( posIdx<insertPosIdx )
				{//当前位置位于插入点前，则下一代直接复制上一代相应残基。
					d_populationNextGen[curPos] = d_population[curPos];
				}
				else if( posIdx==insertPosIdx )
				{//当前位置为插入点，则插入GAP
					d_populationNextGen[curPos] = GAPCODE;
				}
				else
				{//当前位置在插入点后，则在下一代中，原序列相应残基依次向后移动一位。
					d_populationNextGen[curPos] = d_population[curPos - 1];
				}
			}
		}
		else//当前序列不是要插入GAP的序列，则直接复制.
		{
			d_populationNextGen[curPos] = d_population[curPos];
		}
	}
	else//当前进化个体不需要GapInsertMutation，则直接复制。
	{
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_GAP_INSERT(  short *d_populationNextGen, short *d_population, float *d_rand_organ_idx,  float* d_rand_Ty, float* d_rand_Tx, int populationNum, int nSeq, int nBlockWidth, float fOpenGapChance, int GAPCODE )
{
	// Block index
    int bx = blockIdx.x;
// int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int curPos=bx*nBlockWidth+ty*(populationNum*nBlockWidth)+tx;
	if(bx==0) 
	{//最优个体直接进入下一代
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}	

	if(d_rand_organ_idx[bx]<fOpenGapChance)
	{//open gap
		int insertTy = d_rand_Ty[bx]*(nSeq-1);
		if( ty==insertTy )//当前序列为需要插入GAP的序列
		{
			int insertTx = d_rand_Tx[bx]*(nBlockWidth-1);
			if(d_population[bx*nBlockWidth+ty*(populationNum*nBlockWidth)+insertTx]==MSA_GA_CUDA_SPACE)
			{//如果插入位置超出有效数据所在区域，则insertTx=(nBlockWidth-1)-insertTx为新的插入位置（因为留出的空间只占全部空间的20%左右，因此第二次Tx必然指向有数据的区域）
				insertTx = nBlockWidth-1-insertTx;
			}
			if(  insertTx==0 )//当插入位置在序列头部
			{
				if( tx==insertTx )
				{//直接在下一代此序列头部插入GAP
					d_populationNextGen[curPos] = GAPCODE;
				}
				else//在下一代中，原序列其它残基依次向后一位
				{
					d_populationNextGen[curPos] = d_population[curPos-1];
				}
			}
			else//当插入位置不在序列头部
			{
				if( tx<insertTx )
				{//当前位置位于插入点前，则下一代直接复制上一代相应残基。
					d_populationNextGen[curPos] = d_population[curPos];
				}
				else if( tx==insertTx )
				{//当前位置为插入点，则插入GAP
					d_populationNextGen[curPos] = GAPCODE;
				}
				else
				{//当前位置在插入点后，则在下一代中，原序列相应残基依次向后移动一位。
					d_populationNextGen[curPos] = d_population[curPos-1];
				}
			}
		}
		else//当前序列不是要插入GAP的序列，则直接复制.
		{
			d_populationNextGen[curPos] = d_population[curPos];
		}
	}
	else//当前进化个体不需要GapInsertMutation，则直接复制。
	{
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

__global__	void	MSAGA_GAP_EXTENSION_UNLIMIT(short *d_populationNextGen, short *d_population, float *d_rand_organ_idx,  float* d_rand_Ty, short* d_gapExtendPos_Tx, int populationNum, int nSeq, float fOpenGapChance, int GAPCODE, int nWidth, int nOrgWidth, float rnOrgWidth, int nPopulationWidth, float rnPopulationWidth, int nCellNum )
{
	// Block index
    int bx = blockIdx.x;
// int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int curPos = bx*blockDim.x+ty*nWidth+tx;
	if(curPos>=nCellNum)  return;
#if 0
	int orgIdx = curPos%nPopulationWidth/nOrgWidth;
	if (orgIdx >= populationNum)  return;
	int seqIdx = curPos/nPopulationWidth;
	int posIdx = curPos-seqIdx*nPopulationWidth-orgIdx*nOrgWidth;
#endif
	int seqIdx = floor(curPos*rnPopulationWidth);
	int i = curPos - seqIdx*nPopulationWidth;
	int orgIdx = i*rnOrgWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - orgIdx*nOrgWidth;

	if(orgIdx==0) 
	{//最优个体直接进入下一代
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}	

	if(d_rand_organ_idx[orgIdx]<fOpenGapChance)
	{//open gap
		int insertTy = d_rand_Ty[orgIdx]*(nSeq-1);
		if( seqIdx==insertTy )//当前序列为需要插入GAP的序列
		{
			int insertTx = d_gapExtendPos_Tx[orgIdx];
			if(  insertTx==0 )//当插入位置在序列头部
			{
				if( posIdx==insertTx )
				{//直接在下一代此序列头部插入GAP
					d_populationNextGen[curPos] = GAPCODE;
				}
				else//在下一代中，原序列其它残基依次向后一位
				{
					d_populationNextGen[curPos] = d_population[curPos-1];
				}
			}
			else//当插入位置不在序列头部
			{

				if( posIdx<insertTx )
				{//当前位置位于插入点前，则下一代直接复制上一代相应残基。
					d_populationNextGen[curPos] = d_population[curPos];
				}
				else if( posIdx==insertTx )
				{//当前位置为插入点，则插入GAP
					d_populationNextGen[curPos] = GAPCODE;
				}
				else
				{//当前位置在插入点后，则在下一代中，原序列相应残基依次向后移动一位。
					d_populationNextGen[curPos] = d_population[curPos-1];
				}
			}
		}
		else//当前序列不是要插入GAP的序列，则直接复制.
		{
			d_populationNextGen[curPos] = d_population[curPos];
		}
	}
	else//当前进化个体不需要GapInsertMutation，则直接复制。
	{
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

template <int BLOCK_SIZE> __global__ void
MSAGA_GAP_EXTENSION(  short *d_populationNextGen, short *d_population, float *d_rand_organ_idx,  float* d_rand_Ty, short* d_gapExtendPos_Tx, int populationNum, int nSeq, int nBlockWidth, float fOpenGapChance, int GAPCODE )
{
	// Block index
    int bx = blockIdx.x;
// int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int curPos=bx*nBlockWidth+ty*(populationNum*nBlockWidth)+tx;
	if(bx==0) 
	{//最优个体直接进入下一代
		d_populationNextGen[curPos] = d_population[curPos];
		return;
	}	

	if(d_rand_organ_idx[bx]<fOpenGapChance)
	{//open gap
		int insertTy = d_rand_Ty[bx]*(nSeq-1);
		if( ty==insertTy )//当前序列为需要插入GAP的序列
		{
			int insertTx = d_gapExtendPos_Tx[bx];
			//if(d_population[bx*nBlockWidth+ty*(populationNum*nBlockWidth)+insertTx]==MSA_GA_CUDA_SPACE)
			//{//如果插入位置超出有效数据所在区域，则insertTx=(nBlockWidth-1)-insertTx为新的插入位置（因为留出的空间只占全部空间的20%左右，因此第二次Tx必然指向有数据的区域）
			//	insertTx = nBlockWidth-1-insertTx;
			//}
			if(  insertTx==0 )//当插入位置在序列头部
			{
				if( tx==insertTx )
				{//直接在下一代此序列头部插入GAP
					d_populationNextGen[curPos] = GAPCODE;
				}
				else//在下一代中，原序列其它残基依次向后一位
				{
					d_populationNextGen[curPos] = d_population[curPos-1];
				}
			}
			else//当插入位置不在序列头部
			{

				if( tx<insertTx )
				{//当前位置位于插入点前，则下一代直接复制上一代相应残基。
					d_populationNextGen[curPos] = d_population[curPos];
				}
				else if( tx==insertTx )
				{//当前位置为插入点，则插入GAP
					d_populationNextGen[curPos] = GAPCODE;
				}
				else
				{//当前位置在插入点后，则在下一代中，原序列相应残基依次向后移动一位。
					d_populationNextGen[curPos] = d_population[curPos-1];
				}
			}
		}
		else//当前序列不是要插入GAP的序列，则直接复制.
		{
			d_populationNextGen[curPos] = d_population[curPos];
		}
	}
	else//当前进化个体不需要GapInsertMutation，则直接复制。
	{
		d_populationNextGen[curPos] = d_population[curPos];
	}
}

__global__  void  MSAGA_SELECT_NEXT_GEN_UNLIMIT(short *d_populationNextGen, short *d_population, uint *d_organ_index, int nOrgWidth, float rnOrgWidth, int nPopulationWidth, float rnPopulationWidth, int nWidth, int nCellNum)
{
	// Block index
    int bx = blockIdx.x;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int curPos = bx*blockDim.x+ty*nWidth+tx;
	if (curPos >= nCellNum)  return;
	//curPos = nPopulationWidth*j+i
//#if 0
	int j = floor(curPos*rnPopulationWidth);
	int i = curPos - j*nPopulationWidth;
	int destOrgIdx = i*rnOrgWidth;
//#endif

//	int destOrgIdx = curPos%nPopulationWidth / nOrgWidth;

	int srcOrgIdx = d_organ_index[destOrgIdx];
#if 0
	int seqIdx = curPos/nPopulationWidth;
	int posIdx = curPos - seqIdx*nPopulationWidth - destOrgIdx*nOrgWidth;
	d_populationNextGen[curPos] = d_population[srcOrgIdx*nOrgWidth + seqIdx*nPopulationWidth + posIdx];
#endif
	d_populationNextGen[curPos] = d_population[srcOrgIdx*nOrgWidth + curPos - destOrgIdx*nOrgWidth];
}


template <int BLOCK_SIZE> __global__ void
MSAGA_SELECT_NEXT_GEN( short *d_populationNextGen, short *d_population, uint *d_organ_index,  int nWidth)
{
	// Block index
    int bx = blockIdx.x;
	
	// Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

	int offsetBlock = d_organ_index[bx];
	d_populationNextGen[bx*blockDim.x+ty*nWidth+tx] = d_population[offsetBlock*blockDim.x+ty*nWidth+tx];
}

template <int BLOCK_SIZE> __global__ void
MSAGA_SELECT_MAX_CUDA(uint *d_selected_organ_idx, uint* d_organ_idx, short* d_organ_score, float* d_random_organ_idx, short organ_num, short candidateLen)
{
	// Block index
    int bx = blockIdx.x;
//    int by = blockIdx.y;
	
	// Thread index
    int tx = threadIdx.x;
//    int ty = threadIdx.y;

	int curPos = bx*blockDim.x+tx;

	if( curPos==0 )
	{
		d_selected_organ_idx[curPos] = d_organ_idx[0];//the best one directly selected into next generation.
		return;
	}

	int curIdxPos = curPos*candidateLen;
	short IdxCandidate= d_organ_idx[(int)(d_random_organ_idx[curIdxPos]*(organ_num-1))];
	short maxScore = d_organ_score[IdxCandidate];
	for( int i=1; i<candidateLen; ++i )
	{
		short idx = d_organ_idx[(int)d_random_organ_idx[curIdxPos+i]*(organ_num-1)];
		if(maxScore<d_organ_score[idx])
		{
			maxScore = d_organ_score[idx];
			IdxCandidate = idx;
		}
	}
	d_selected_organ_idx[curPos] = IdxCandidate;
}

__global__ void MSAGA_SELECT_MAX_CUDA_UNLIMIT(uint *d_selected_organ_idx, uint* d_organ_idx, float* d_organ_score, float* d_random_organ_idx, int organ_num, short candidateLen, int nWidth)
{
	// Block index
	int bx = blockIdx.x;
	//    int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curOrgIdx = bx*blockDim.x + ty*nWidth + tx;
	if (curOrgIdx >= organ_num)  return;
	//int orgIdx = curPos%nPopulationWidth / nOrgWidth;
	if (curOrgIdx == 0)
	{
		d_selected_organ_idx[curOrgIdx] = d_organ_idx[0];//the best one directly selected into next generation.
		return;
	}

	int curIdxPos = curOrgIdx*candidateLen;
	short IdxCandidate = d_organ_idx[(int)(d_random_organ_idx[curIdxPos] * (organ_num - 1))];
	short maxScore = d_organ_score[IdxCandidate];
	for (int i = 1; i<candidateLen; ++i)
	{
		short idx = d_organ_idx[(int)d_random_organ_idx[curIdxPos + i] * (organ_num - 1)];
		if (maxScore<d_organ_score[idx])
		{
			maxScore = d_organ_score[idx];
			IdxCandidate = idx;
		}
	}
	d_selected_organ_idx[curOrgIdx] = IdxCandidate;
}

__global__ void MSAGA_RECORD_NEXT_GEN_LEN_UNLIMIT(unsigned int* d_orgLenNextGen, unsigned int* d_orgLen, uint* d_SelectedOrgIdx, short nWidth, int nCellNum)
{
	// Block index
	int bx = blockIdx.x;
	//    int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int curOrgIdx = bx*blockDim.x + ty*nWidth + tx;
	if (curOrgIdx >= nCellNum)  return;
	d_orgLenNextGen[curOrgIdx] = d_orgLen[d_SelectedOrgIdx[curOrgIdx]];
}

void		CMSAGA_CUDA_Algorithm::CUDA_DEV2HOST( void *h_mem, const void * const d_mem,  unsigned int size_mem )
{
	 // copy device memory to host
	cudaError_t error;
    error = cudaMemcpy(h_mem, d_mem, size_mem, cudaMemcpyDeviceToHost);

    if (error != cudaSuccess)
    {
		throw   CAppException( DEF_EXCEPTION_UNEXPECTED,DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK
			,__EXCEPTION_SITE__ ,"cudaMemcpy  Device2Host Failured !" );
    }
}
 
void		CMSAGA_CUDA_Algorithm::CUDA_MALLOC( void **d_mem, unsigned int size_mem )
{
	cudaError_t error;
	error = cudaMalloc(d_mem, size_mem);
	if (error != cudaSuccess)
	{
		throw   CAppException( DEF_EXCEPTION_UNEXPECTED,DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK
			,__EXCEPTION_SITE__ ,"cudaMalloc Failured !" );
	}
}

void		CMSAGA_CUDA_Algorithm::CUDA_MALLOC( void **d_mem, void *h_mem, unsigned int size_mem )
{
	cudaError_t error;
	error = cudaMalloc(d_mem, size_mem);
	if (error != cudaSuccess)
	{
		throw   CAppException( DEF_EXCEPTION_UNEXPECTED,DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK
			,__EXCEPTION_SITE__ ,"cudaMalloc Failured !" );
	}
	error = cudaMemcpy(*d_mem, h_mem, size_mem, cudaMemcpyHostToDevice);
	if (error != cudaSuccess)
	{
		throw   CAppException( DEF_EXCEPTION_UNEXPECTED,DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK
			,__EXCEPTION_SITE__ ,"cudaMemcpy  Host2Device Failured!" );
	}
}

void		CMSAGA_CUDA_Algorithm::CUDA_FREE( void *d_mem )
{
	cudaError_t error;
	error = cudaFree(d_mem);
	if (error != cudaSuccess)
	{
		throw   CAppException( DEF_EXCEPTION_UNEXPECTED,DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK
			,__EXCEPTION_SITE__ ,"cudaFree Failured !" );
	}
}

void	CMSAGA_CUDA_Algorithm::SetSubstitutionMat2GPU(  SubMatrixType type  )
{
	if( NULL!=m_pSubMatDevice )
	{
		if( type != m_curSubMatType )
		{
			cudaFree( m_pSubMatDevice );
			CUDA_MALLOC( (void**)&m_pSubMatDevice, CGlobalSpace::m_sAlignParams.getSubMatrix(type), sizeof(short)*CGlobalSpace::m_sAlignParams.getSubMatrixSize(type) );
		}
	}
	else
	{
		m_curSubMatType = type;
		CUDA_MALLOC( (void**)&m_pSubMatDevice, CGlobalSpace::m_sAlignParams.getSubMatrix(type), sizeof(short)*CGlobalSpace::m_sAlignParams.getSubMatrixSize(type) );
	}
}

void		CMSAGA_CUDA_Algorithm::MutationGapReductionDeviceUnlimit()
{
	RandomNumberGeneratorDevice( m_pMutation_RandOrgan_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTy_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTx_Device, MSA_GA_POPULATION_SIZE );

	//mark the gap start position with 1, or else 0.
	dim3 thread( MSA_GA_CUDA_BLOCK_SIZE, MSA_GA_CUDA_BLOCK_SIZE );
	dim3 grid( ceil( (float)m_nCellNum/(thread.x*thread.y) ), 1 );
	MSAGA_MARK_GAP_START_UNLIMIT<<< grid, thread>>>( m_pGapStartRecorder_Device, m_pPopulationDevice[m_curPopulationIndex], 
													m_pMutation_RandOrgan_Device, m_pMutation_RandTy_Device, 
													m_pPopulation->nSeqSize, thread.x*grid.x, MSA_GA_EXTENSION_MUTATION, 
													CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), 
													MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN, POPULATION_WIDTH, R_POPULATION_WIDTH, m_nCellNum);
	//dbgWriteFile( "gapStartRecorder.txt", m_pGapStartRecorder_Device, MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN, MSA_GA_CUDA_MAX_LEN  );
	grid.x = ceil((float)MSA_GA_POPULATION_SIZE / (thread.x*thread.y));
	MSAGA_LOCATE_GAP_UNLIMIT<32><<< grid, thread>>>( m_pGapPos_Tx_Device, m_pGapStartRecorder_Device, m_pMutation_RandTx_Device, m_pMutation_RandOrgan_Device,
													 MSA_GA_CUDA_MAX_LEN, grid.x*thread.x, MSA_GA_REDUCTION_MUTATION,
													 CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), MSA_GA_POPULATION_SIZE);

	//delete one gap at Tx
	grid.x = ceil((float)m_nCellNum / (thread.x*thread.y));
	//dbgWriteFile( "gapPos_Tx.txt", m_pGapPos_Tx_Device, MSA_GA_POPULATION_SIZE, MSA_GA_POPULATION_SIZE  );
	MSAGA_GAP_REDUCE_UNLIMIT<<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex], 
												m_pMutation_RandOrgan_Device, m_pMutation_RandTy_Device, m_pGapPos_Tx_Device,
												m_pPopulation->nSeqSize, thread.x*grid.x, MSA_GA_REDUCTION_MUTATION,
												CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE),  
												MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN, POPULATION_WIDTH, R_POPULATION_WIDTH, m_nCellNum);
	m_curPopulationIndex = 1-m_curPopulationIndex;

	//arrange the organs.
	ArrangeSeqDeviceUnlimit( m_pMutation_RandOrgan_Device, MSA_GA_REDUCTION_MUTATION );
}

void		CMSAGA_CUDA_Algorithm::MutationGapReductionDevice()
{
	RandomNumberGeneratorDevice( m_pMutation_RandOrgan_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTy_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTx_Device, MSA_GA_POPULATION_SIZE );
	
	//mark the gap start position with 1, or else 0.
	dim3 thread(MSA_GA_CUDA_MAX_LEN, m_pPopulation->nSeqSize );
	dim3 grid( MSA_GA_POPULATION_SIZE, 1 );
	MSAGA_MARK_GAP_START<32><<< grid, thread>>>( m_pGapStartRecorder_Device, m_pPopulationDevice[m_curPopulationIndex], 
		                                                                                        m_pMutation_RandOrgan_Device, m_pMutation_RandTy_Device, 
																								m_pPopulation->nSeqSize, MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN, MSA_GA_REDUCTION_MUTATION, 
																								CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE));

	thread.x = MSA_GA_POPULATION_SIZE;
	thread.y = 1;
	grid.x = 1;
	MSAGA_LOCATE_GAP<32><<< grid, thread>>>( m_pGapPos_Tx_Device, m_pGapStartRecorder_Device, m_pMutation_RandTx_Device, m_pMutation_RandOrgan_Device,
		                                                                               MSA_GA_CUDA_MAX_LEN, MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN, MSA_GA_REDUCTION_MUTATION,
																					   CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE));

	//delete one gap at Tx
	thread.x=MSA_GA_CUDA_MAX_LEN;
	thread.y=m_pPopulation->nSeqSize;
	grid.x=MSA_GA_POPULATION_SIZE;
	grid.y=1;
	MSAGA_GAP_REDUCE<32><<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex], 
		                                                                               m_pMutation_RandOrgan_Device, m_pMutation_RandTy_Device, m_pGapPos_Tx_Device,
																					   m_pPopulation->nSeqSize, MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN, MSA_GA_REDUCTION_MUTATION,
																					   CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE));
	m_curPopulationIndex = 1-m_curPopulationIndex;

	//arrange the organs.
	ArrangeSeqDevice( m_pMutation_RandOrgan_Device, MSA_GA_REDUCTION_MUTATION );
}

void		CMSAGA_CUDA_Algorithm::MutationGapExtensionDeviceUnlimit()
{
	RandomNumberGeneratorDevice( m_pMutation_RandOrgan_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTy_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTx_Device, MSA_GA_POPULATION_SIZE );
	
	//mark the gap start position with 1, or else 0.
	dim3 thread( MSA_GA_CUDA_BLOCK_SIZE, MSA_GA_CUDA_BLOCK_SIZE );
	dim3 grid( ceil( (float)m_nCellNum/(thread.x*thread.y) ), 1 );
	MSAGA_MARK_GAP_START_UNLIMIT<<< grid, thread>>>( m_pGapStartRecorder_Device, m_pPopulationDevice[m_curPopulationIndex], 
													 m_pMutation_RandOrgan_Device, m_pMutation_RandTy_Device, 
													 m_pPopulation->nSeqSize, thread.x*grid.x, MSA_GA_EXTENSION_MUTATION, 
													 CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), 
													 MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN, POPULATION_WIDTH, R_POPULATION_WIDTH, m_nCellNum);
	
	grid.x = ceil((float)MSA_GA_POPULATION_SIZE / (thread.x*thread.y));
	MSAGA_LOCATE_GAP_UNLIMIT<32><<< grid, thread>>>( m_pGapPos_Tx_Device, m_pGapStartRecorder_Device, m_pMutation_RandTx_Device, m_pMutation_RandOrgan_Device,
		                                             MSA_GA_CUDA_MAX_LEN, thread.x*grid.x, MSA_GA_EXTENSION_MUTATION,
													 CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), MSA_GA_POPULATION_SIZE);
	
	//extend the gap at Tx
#if 0
	thread.x=MSA_GA_CUDA_BLOCK_SIZE;
	thread.y=thread.x;
	grid.x=ceil( (float)MSA_GA_CUDA_MAX_LEN*m_pPopulation->nSeqSize*MSA_GA_POPULATION_SIZE/(thread.x*thread.y) );
	grid.y=1;
#endif
	MSAGA_GAP_EXTENSION_UNLIMIT<<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex], 
													m_pMutation_RandOrgan_Device, m_pMutation_RandTy_Device, m_pGapPos_Tx_Device,												
													MSA_GA_POPULATION_SIZE, m_OrganismSize, MSA_GA_EXTENSION_MUTATION, 												
													CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), 												
													thread.x*grid.x, MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN,
													POPULATION_WIDTH, R_POPULATION_WIDTH, m_nCellNum);
	m_curPopulationIndex = 1-m_curPopulationIndex;
	
	//arrange the organs.
	ArrangeSeqDeviceUnlimit( m_pMutation_RandOrgan_Device, MSA_GA_EXTENSION_MUTATION );

}

void		CMSAGA_CUDA_Algorithm::MutationGapExtensionDevice()
{
	RandomNumberGeneratorDevice( m_pMutation_RandOrgan_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTy_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTx_Device, MSA_GA_POPULATION_SIZE );

	//mark the gap start position with 1, or else 0.
	dim3 thread(MSA_GA_CUDA_MAX_LEN, m_pPopulation->nSeqSize );
	dim3 grid( MSA_GA_POPULATION_SIZE, 1 );
	MSAGA_MARK_GAP_START<32><<< grid, thread>>>( m_pGapStartRecorder_Device, m_pPopulationDevice[m_curPopulationIndex], 
		                                                                                        m_pMutation_RandOrgan_Device, m_pMutation_RandTy_Device, 
																								m_pPopulation->nSeqSize, MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN, MSA_GA_EXTENSION_MUTATION, 
																								CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE));

	thread.x = MSA_GA_POPULATION_SIZE;
	thread.y = 1;
	grid.x = 1;
	MSAGA_LOCATE_GAP<32><<< grid, thread>>>( m_pGapPos_Tx_Device, m_pGapStartRecorder_Device, m_pMutation_RandTx_Device, m_pMutation_RandOrgan_Device,
		                                                                               MSA_GA_CUDA_MAX_LEN, MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN, MSA_GA_EXTENSION_MUTATION,
																					   CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE));

	//extend the gap at Tx
	thread.x=MSA_GA_CUDA_MAX_LEN;
	thread.y=m_pPopulation->nSeqSize;
	grid.x=MSA_GA_POPULATION_SIZE;
	grid.y=1;
	MSAGA_GAP_EXTENSION<32><<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex], 
		                                                                                     m_pMutation_RandOrgan_Device, m_pMutation_RandTy_Device, m_pGapPos_Tx_Device,
																					         MSA_GA_POPULATION_SIZE, m_OrganismSize, MSA_GA_CUDA_MAX_LEN, MSA_GA_EXTENSION_MUTATION,
																					         CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE));
	m_curPopulationIndex = 1-m_curPopulationIndex;

	//arrange the organs.
	ArrangeSeqDevice( m_pMutation_RandOrgan_Device, MSA_GA_EXTENSION_MUTATION );
}

void		CMSAGA_CUDA_Algorithm::MutationGapInsertDeviceUnlimit()
{
	RandomNumberGeneratorDevice( m_pMutation_RandOrgan_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTx_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTy_Device, MSA_GA_POPULATION_SIZE );

	//locate the position of gap insert seq and position.
	dim3 thread( MSA_GA_CUDA_BLOCK_SIZE, MSA_GA_CUDA_BLOCK_SIZE );
	dim3 grid(ceil((float)m_nCellNum/(thread.x*thread.y)), 1);
	MSAGA_GAP_INSERT_UNLIMIT<<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex], 
												 m_pMutation_RandOrgan_Device, m_pMutation_RandTy_Device, m_pMutation_RandTx_Device,
												 MSA_GA_POPULATION_SIZE, m_OrganismSize, MSA_GA_CUDA_MAX_LEN, MSA_GA_BLOCK_MUTATION,
												 CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE),
												 thread.x*grid.x, MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN, POPULATION_WIDTH, R_POPULATION_WIDTH,
												 m_pOrganLenDevice[m_curOrganLenIndex], m_nCellNum);
	m_curPopulationIndex = 1-m_curPopulationIndex;
	//dbgCheckSequence();
	//arrange the organs.
	ArrangeSeqDeviceUnlimit( m_pMutation_RandOrgan_Device, MSA_GA_BLOCK_MUTATION );
	//dbgCheckSequence();
}

void		CMSAGA_CUDA_Algorithm::MutationGapInsertDevice()
{
	//gap insert mutation
	RandomNumberGeneratorDevice( m_pMutation_RandOrgan_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTx_Device, MSA_GA_POPULATION_SIZE );
	RandomNumberGeneratorDevice( m_pMutation_RandTy_Device, MSA_GA_POPULATION_SIZE );

	//locate the position of gap insert seq and position.
	dim3 thread(MSA_GA_CUDA_MAX_LEN, m_pPopulation->nSeqSize );
	dim3 grid( MSA_GA_POPULATION_SIZE, 1 );
	MSAGA_GAP_INSERT<32><<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex], 
		                                                                              m_pMutation_RandOrgan_Device, m_pMutation_RandTy_Device, m_pMutation_RandTx_Device,
																					  MSA_GA_POPULATION_SIZE, m_OrganismSize, MSA_GA_CUDA_MAX_LEN, MSA_GA_BLOCK_MUTATION,
																					  CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE));
	m_curPopulationIndex = 1-m_curPopulationIndex;


	//arrange the organs.
	ArrangeSeqDevice( m_pMutation_RandOrgan_Device, MSA_GA_BLOCK_MUTATION );
}

bool		CMSAGA_CUDA_Algorithm::FitnessDeviceUnlimit()
{
	dim3 thread(MSA_GA_CUDA_BLOCK_SIZE, MSA_GA_CUDA_BLOCK_SIZE);
	dim3 grid(ceil((float)MSA_GA_CUDA_MAX_LEN*m_pPopulation->nSeqSize*MSA_GA_POPULATION_SIZE / (thread.x*thread.y)), 1);
//	dbgWriteFile("score_before.txt", m_pOrganScore_Device, MSA_GA_POPULATION_SIZE, MSA_GA_POPULATION_SIZE);
//	dbgWriteSeqsFromDevice2File("before scoring.txt");
	//sp unit score
	MSAGA_SPS_UNIT_CUDA_UNLIMIT << < grid, thread >> >(m_pUnitScore_Device, m_pPopulationDevice[m_curPopulationIndex], m_pOrganLenDevice[m_curPopulationIndex], m_pSubMatDevice, m_pSeqWeight_Device,
													  m_pPopulation->nSeqSize, thread.x*grid.x, POPULATION_WIDTH, R_POPULATION_WIDTH, MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN,
													  CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), m_rnPairs,
													  CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(MSAGA_MAT_TYPE), 
													  CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapExtendCost(MSAGA_MAT_TYPE));

//	dbgWriteFile(std::string("UnitScore.txt"), m_pUnitScore_Device, MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE*m_pPopulation->nSeqSize, MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE*m_pPopulation->nSeqSize);

	//column score
	grid.x = ceil((float)MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE / (thread.x*thread.y));
	MSAGA_SPS_COLSUM_CUDA_UNLIMIT << < grid, thread >> >(m_pColumnScore_Device, m_pUnitScore_Device, m_pPopulation->nSeqSize,
		thread.x*grid.x, MSA_GA_POPULATION_SIZE, POPULATION_WIDTH, MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN);

//	dbgWriteFile(std::string("ColScore.txt"), m_pUnitScore_Device, MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE*m_pPopulation->nSeqSize, MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE*m_pPopulation->nSeqSize);

	grid.x = ceil((float)MSA_GA_POPULATION_SIZE / (thread.x*thread.y));
	MSAGA_SPS_ORGAN_CUDA_UNLIMIT << < grid, thread >> >(m_pOrganScore_Device, m_pColumnScore_Device, thread.x*grid.x, MSA_GA_CUDA_MAX_LEN, MSA_GA_POPULATION_SIZE);

//	dbgWriteFile("score_after.txt", m_pOrganScore_Device, MSA_GA_POPULATION_SIZE, MSA_GA_POPULATION_SIZE);

	return false;
}

void		CMSAGA_CUDA_Algorithm::InitArrangement()
{
	CUDA_MALLOC( (void**)&m_pGapRecorder_Device,  sizeof(short)*MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN );
	CUDA_MALLOC( (void**)&m_pRangeRecorder_Device,  sizeof(short)*MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN );
	CUDA_MALLOC( (void**)&m_pGapLenRecorder_Device,  sizeof(short)*MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN );
}

void		CMSAGA_CUDA_Algorithm::InitSelection()
{
	float size=log10((float)MSA_GA_POPULATION_SIZE)/log10((float)2);
	m_nOrganSize_Align = pow(2.0f, ceil(size));

	unsigned int size_organ_score_align = m_nOrganSize_Align;//SHARED_SIZE_LIMIT;
	unsigned int mem_size_organ_score_align =  sizeof(float) * size_organ_score_align;
	unsigned int mem_size_organ_index_align = sizeof(uint) * size_organ_score_align;;

	CUDA_MALLOC( (void**)&m_pOrganScore_Align_Device,  mem_size_organ_score_align );
	CUDA_MALLOC( (void**)&m_pOrganIndex_Align_Device,  mem_size_organ_index_align );
	CUDA_MALLOC( (void**)&m_pOrganScore_Align_Sorted_Device,  mem_size_organ_score_align );
	CUDA_MALLOC( (void**)&m_pOrganIndex_Align_Sorted_Device,  mem_size_organ_index_align );

	//using tournament selection method to select the organ
	m_nRandom_organ_idx = sizeof(float)*MSA_GA_CUDA_SELECTION_CANDIDATE*MSA_GA_POPULATION_SIZE;
	CUDA_MALLOC( (void**)&m_pRandomOrganIdx_Device,  m_nRandom_organ_idx );

	CUDA_MALLOC( (void**)&m_pSelectedOrganIdx_Device,  sizeof(uint)*MSA_GA_POPULATION_SIZE );
}

void		CMSAGA_CUDA_Algorithm::InitFitness()
{
	//unit score
	unsigned int size_score_unit = m_pPopulation->nSeqSize*MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	CUDA_MALLOC((void**)&m_pUnitScore_Device, sizeof(float)*size_score_unit);
	//load seq weight
	unsigned int size_seq_weight = size_score_unit;
	unsigned int mem_size_seq_weight = sizeof(float)*size_seq_weight;
	CUDA_MALLOC((void**)&m_pSeqWeight_Device, mem_size_seq_weight);
	float *pSeqWeight = new float[size_seq_weight];
	for (int i = 0; i < m_pPopulation->nSeqSize; ++i)
	{
		for (int j = 0; j < POPULATION_WIDTH; ++j)
		{
			pSeqWeight[i*POPULATION_WIDTH + j] = m_SeqWeight[i];
		}
	}

	CUDA_MALLOC((void**)&m_pSeqWeight_Device, pSeqWeight, mem_size_seq_weight);
	SAFE_DELETE_ARRAY(pSeqWeight);

	//column score
	unsigned int size_score_colsum = MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	CUDA_MALLOC((void**)&m_pColumnScore_Device, sizeof(float)*size_score_colsum);
	//organ score
	unsigned int size_score_population = MSA_GA_POPULATION_SIZE;
	CUDA_MALLOC((void**)&m_pOrganScore_Device, sizeof(float)*size_score_population);
}

void		CMSAGA_CUDA_Algorithm::InitRecomb()
{
	int	nRecomb = (MSA_GA_POPULATION_SIZE-1)/2;//half of the population number

	CUDA_MALLOC( (void**)&m_pRecomb_RandOrgan_Device,  sizeof(float)*nRecomb );
	CUDA_MALLOC( (void**)&m_pRecomb_RandTy_Device,  sizeof(float)*nRecomb );

	//cudaMemset( (void*)m_pRecomb_RandOrgan_Device, 1, sizeof(float)*nRecomb );
	//cudaMemset( (void*)m_pRecomb_RandTy_Device, 2, sizeof(float)*nRecomb );
	//dbgWriteFile( "m_pRecomb_RandOrgan_Device.txt", m_pRecomb_RandOrgan_Device, nRecomb, nRecomb );
	//dbgWriteFile( "m_pRecomb_RandTy_Device.txt", m_pRecomb_RandTy_Device, nRecomb, nRecomb );

	CUDA_MALLOC( (void**)&m_pRecomb_RandOrganIdx0_Device,  sizeof(float)*nRecomb );
	CUDA_MALLOC( (void**)&m_pRecomb_RandOrganIdx1_Device,  sizeof(float)*nRecomb );

	//for vertical recombination
	CUDA_MALLOC( (void**)&m_pPos1_Min_Device,  sizeof(short)*nRecomb );
	CUDA_MALLOC( (void**)&m_pPos1_Max_Device,  sizeof(short)*nRecomb );

	CUDA_MALLOC( (void**)&m_pRecomb_RandTx_Device,  sizeof(float)*nRecomb );

	CUDA_MALLOC( (void**)&m_pOrganIndex0_Device,  sizeof(short)*nRecomb );
	CUDA_MALLOC( (void**)&m_pOrganIndex1_Device,  sizeof(short)*nRecomb );

	CUDA_MALLOC( (void**)&m_pVerticalRecomb_pos0_Device,  sizeof(short)*nRecomb );
	CUDA_MALLOC( (void**)&m_pVerticalRecomb_pos1_Device,  sizeof(short)*nRecomb*m_pPopulation->nSeqSize );
}

void		CMSAGA_CUDA_Algorithm::InitMutation()
{
	CUDA_MALLOC( (void**)&m_pMutation_RandOrgan_Device,  sizeof(float)*MSA_GA_POPULATION_SIZE );
	CUDA_MALLOC( (void**)&m_pMutation_RandTx_Device,  sizeof(float)*MSA_GA_POPULATION_SIZE );
	CUDA_MALLOC( (void**)&m_pMutation_RandTy_Device,  sizeof(float)*MSA_GA_POPULATION_SIZE );
	CUDA_MALLOC( (void**)&m_pGapStartRecorder_Device,  sizeof(short)*MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN );
	CUDA_MALLOC( (void**)&m_pGapPos_Tx_Device,  sizeof(short)*MSA_GA_POPULATION_SIZE );
}

void     CMSAGA_CUDA_Algorithm::ClearDevice()
{
	CUDA_FREE((void*)m_pPopulationDevice[m_curPopulationIndex]);
	CUDA_FREE((void*)m_pPopulationDevice[1-m_curPopulationIndex]);
	CUDA_FREE((void*)m_pOrganLenDevice[m_curOrganLenIndex]);
	CUDA_FREE((void*)m_pOrganLenDevice[1-m_curOrganLenIndex]);

	curandDestroyGenerator( m_RandomNumberGen );

	CUDA_FREE((void*)m_pSubMatDevice);

	//FITNESS
	//unit score
	CUDA_FREE((void*)m_pUnitScore_Device);
	//column score
	CUDA_FREE((void*)m_pColumnScore_Device);
	//organ score
	CUDA_FREE((void*)m_pOrganScore_Device);
	//weight
	CUDA_FREE((void*)m_pSeqWeight_Device);

	//SELECTION 
	CUDA_FREE((void*)m_pOrganScore_Align_Device );
	CUDA_FREE((void*)m_pOrganIndex_Align_Device );
	CUDA_FREE((void*)m_pOrganScore_Align_Sorted_Device );
	CUDA_FREE((void*)m_pOrganIndex_Align_Sorted_Device );
	CUDA_FREE((void*)m_pRandomOrganIdx_Device );
	CUDA_FREE((void*)m_pSelectedOrganIdx_Device );

	//MUTATION
	CUDA_FREE((void*)m_pMutation_RandOrgan_Device );
	CUDA_FREE((void*)m_pMutation_RandTx_Device );
	CUDA_FREE((void*)m_pMutation_RandTy_Device );
	CUDA_FREE((void*)m_pGapStartRecorder_Device );
	CUDA_FREE((void*)m_pGapPos_Tx_Device );

	//RECOMB
	CUDA_FREE((void*)m_pRecomb_RandOrgan_Device );
	CUDA_FREE((void*)m_pRecomb_RandTy_Device );

	CUDA_FREE((void*)m_pRecomb_RandOrganIdx0_Device );
	CUDA_FREE((void*)m_pRecomb_RandOrganIdx1_Device );

	//for vertical recombination
	CUDA_FREE((void*)m_pPos1_Min_Device );
	CUDA_FREE((void*)m_pPos1_Max_Device );

	CUDA_FREE((void*)m_pRecomb_RandTx_Device );

	CUDA_FREE((void*)m_pOrganIndex0_Device );
	CUDA_FREE((void*)m_pOrganIndex1_Device );
	CUDA_FREE((void*)m_pVerticalRecomb_pos0_Device );
	CUDA_FREE((void*)m_pVerticalRecomb_pos1_Device );

	//ALIGNMENT
	CUDA_FREE((void*)m_pGapRecorder_Device );
	CUDA_FREE((void*)m_pRangeRecorder_Device );
	CUDA_FREE((void*)m_pGapLenRecorder_Device );
}

void		CMSAGA_CUDA_Algorithm::InitDevice()
{
	//init the random number generator
	curandCreateGenerator( &m_RandomNumberGen,  CURAND_RNG_PSEUDO_XORWOW );

	//set the substitution matrix.
	SetSubstitutionMat2GPU( MSAGA_MAT_TYPE );
	InitFitness();
	InitSelection();
	InitMutation();
	InitRecomb();
	InitArrangement();

	//计算个体的初始长度
	dim3 thread(MSA_GA_CUDA_BLOCK_SIZE, MSA_GA_CUDA_BLOCK_SIZE);
	dim3 grid(ceil((float)MSA_GA_POPULATION_SIZE / (thread.x*thread.y)), 1);
	MSAGA_CALC_SEQ_LEN_UNLIMIT<32> << < grid, thread >> >(false, m_pOrganLenDevice[1 - m_curOrganLenIndex], m_pOrganLenDevice[m_curOrganLenIndex],
														  m_pPopulationDevice[m_curPopulationIndex], m_pRandomOrganIdx_Device, 1,
														  MSA_GA_CUDA_MAX_LEN, grid.x*thread.x, MSA_GA_POPULATION_SIZE);
	m_curOrganLenIndex = 1-m_curOrganLenIndex;
}

void		CMSAGA_CUDA_Algorithm::RandomNumberGeneratorDevice( float *pRandNum, int nRand )
{
	curandStatus_t  curandResult = curandGenerateUniform( m_RandomNumberGen, pRandNum, nRand );
    if (curandResult != CURAND_STATUS_SUCCESS)
    {
		throw   CAppException( DEF_EXCEPTION_UNEXPECTED,DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK
			,__EXCEPTION_SITE__ ,"Could not generate random numbers! ");
    }
}

void		CMSAGA_CUDA_Algorithm::SelectionDeviceUnlimit()
{
	//align the organ score and index with  power-of-two array lengths
	//the non used score initialized with -1e4;
	dim3 thread(m_nOrganSize_Align, 1);
	dim3 grid(1, 1);
	if (m_nOrganSize_Align >= SHARED_SIZE_LIMIT)
	{
		thread.x = MSA_GA_CUDA_BLOCK_SIZE;
		thread.y = MSA_GA_CUDA_BLOCK_SIZE;
		grid.x = ceil((float)m_nOrganSize_Align / (MSA_GA_CUDA_BLOCK_SIZE*MSA_GA_CUDA_BLOCK_SIZE));
	}
	//dbgWriteFile("OrganScore_Org.txt", m_pOrganScore_Device, MSA_GA_POPULATION_SIZE, MSA_GA_POPULATION_SIZE);
	MSAGA_ALIGN_CUDA_UNLIMIT<32> << < grid, thread >> >(m_pOrganScore_Align_Device, m_pOrganIndex_Align_Device, m_pOrganScore_Device, MSA_GA_POPULATION_SIZE, grid.x*thread.x);
	
	//sort the score using bitonic sort method
#if 0
	thread.x = SHARED_SIZE_LIMIT/2;
	thread.y = 1;
	grid.x =1;
	grid.y = 1;
    bitonicSortShared<<< grid, thread >>>(m_pOrganScore_Align_Sorted_Device, m_pOrganIndex_Align_Sorted_Device, 
		                                  m_pOrganScore_Align_Device, m_pOrganIndex_Align_Device, 
										  SHARED_SIZE_LIMIT, 0 );
#endif
	//dbgWriteFile("score_align.txt", m_pOrganScore_Align_Device, m_nOrganSize_Align, m_nOrganSize_Align);
	if (m_nOrganSize_Align < SHARED_SIZE_LIMIT)
	{
		bitonicSortShared(m_pOrganScore_Align_Sorted_Device, m_pOrganIndex_Align_Sorted_Device,
						  m_pOrganScore_Align_Device, m_pOrganIndex_Align_Device, 1, m_nOrganSize_Align, 0);
	}
	else
	{
		bitonicSort(m_pOrganScore_Align_Sorted_Device, m_pOrganIndex_Align_Sorted_Device,
			m_pOrganScore_Align_Device, m_pOrganIndex_Align_Device, 1, m_nOrganSize_Align, 0);
	}
	//dbgWriteFile("aligned_sorted_score.txt", m_pOrganScore_Align_Sorted_Device, m_nOrganSize_Align, m_nOrganSize_Align);

	CUDA_DEV2HOST(&m_bestScore, m_pOrganScore_Align_Sorted_Device, sizeof(float));

	//using tournament selection method to select the organ
	// Generate random numbers
	RandomNumberGeneratorDevice( m_pRandomOrganIdx_Device, MSA_GA_CUDA_SELECTION_CANDIDATE*MSA_GA_POPULATION_SIZE );

	//select the max one of the candidate
#if 0
	thread.x = MSA_GA_POPULATION_SIZE;
	thread.y = 1;
	grid.x = 1;
	grid.y = 1;
	MSAGA_SELECT_MAX_CUDA<32><<< grid, thread>>>(m_pSelectedOrganIdx_Device, m_pOrganIndex_Align_Sorted_Device, m_pOrganScore_Align_Sorted_Device, 
		                                         m_pRandomOrganIdx_Device, MSA_GA_POPULATION_SIZE, MSA_GA_CUDA_SELECTION_CANDIDATE );
#endif
	thread.x = MSA_GA_CUDA_BLOCK_SIZE;
	thread.y = thread.x;
	grid.x = ceil((float)MSA_GA_POPULATION_SIZE / (thread.x*thread.y));
	grid.y = 1;
	MSAGA_SELECT_MAX_CUDA_UNLIMIT <<<grid, thread >>>(m_pSelectedOrganIdx_Device, m_pOrganIndex_Align_Sorted_Device, m_pOrganScore_Align_Sorted_Device,
														m_pRandomOrganIdx_Device, MSA_GA_POPULATION_SIZE, MSA_GA_CUDA_SELECTION_CANDIDATE, grid.x*thread.x);
	//dbgWriteFile("m_pSelectedOrganIdx_Device.txt", m_pSelectedOrganIdx_Device, MSA_GA_POPULATION_SIZE, MSA_GA_POPULATION_SIZE);
	//dbgWriteSeqsFromDevice2File("before next gen.txt");
	//produce the next generation of organs.
	thread.x=MSA_GA_CUDA_BLOCK_SIZE;
	thread.y=thread.x;
	grid.x = ceil((float)m_nCellNum / (thread.x*thread.y));
	grid.y = 1;
	MSAGA_SELECT_NEXT_GEN_UNLIMIT<<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex], 
													 m_pSelectedOrganIdx_Device, MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN,
													 POPULATION_WIDTH, R_POPULATION_WIDTH, grid.x*thread.x, m_nCellNum);
	m_curPopulationIndex = 1-m_curPopulationIndex;
	
//	dbgWriteSeqsFromDevice2File("after next gen.txt");

	//dbgCheckSequence();

	//record the next generation length
	thread.x = MSA_GA_CUDA_BLOCK_SIZE;
	thread.y = thread.x;
	grid.x = ceil((float)MSA_GA_POPULATION_SIZE / (thread.x*thread.y));
	grid.y = 1;
	MSAGA_RECORD_NEXT_GEN_LEN_UNLIMIT << <grid, thread >> >(m_pOrganLenDevice[1 - m_curOrganLenIndex], m_pOrganLenDevice[m_curOrganLenIndex], m_pSelectedOrganIdx_Device, grid.x*thread.x, MSA_GA_POPULATION_SIZE);
	m_curOrganLenIndex = 1 - m_curOrganLenIndex;
	//dbgWriteFile<unsigned int>("m_pOrganLen.txt", m_pOrganLenDevice[m_curOrganLenIndex], MSA_GA_POPULATION_SIZE, MSA_GA_POPULATION_SIZE);
	//dbgWriteSeqsFromDevice2File("checkPass.txt");
}

void		CMSAGA_CUDA_Algorithm::RecombinationDeviceUnlimit()
{
	RecombinationHorizentalDeviceUnlimit();
	//dbgCheckSequence();
	RecombinationVerticalDeviceUnlimit();
}

void		CMSAGA_CUDA_Algorithm::RecombinationDevice()
{
	//dbgWriteSeqsFromDevice2File( "seqResult0.txt" );
	RecombinationHorizentalDevice();
	//dbgWriteSeqsFromDevice2File( "seqResult0.txt" );
    RecombinationVerticalDevice();
	//dbgWriteSeqsFromDevice2File( "seqResult1.txt" );
	//static int count=0;
	//++count;
	//if(count==2)
	//{
	//	int i=0;
	//}
}

void		CMSAGA_CUDA_Algorithm::RecombinationVerticalDeviceUnlimit()
{
	int	nRecomb = (MSA_GA_POPULATION_SIZE-1)/2;
	RandomNumberGeneratorDevice( m_pRecomb_RandOrgan_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandTx_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandOrganIdx0_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandOrganIdx1_Device, nRecomb );

	//find the recomb pos of org0 and org1.
#if 0
	dim3 thread(nRecomb, 1);
	dim3 grid(1, 1 );
	MSAGA_RECOMB_VERTICAL_RECOMB_POS0_UNLIMIT<<< grid, thread>>>( m_pVerticalRecomb_pos0_Device, m_pOrganIndex0_Device, m_pOrganIndex1_Device,
																  m_pRecomb_RandOrgan_Device, m_pRecomb_RandTx_Device, 
																  m_pRecomb_RandOrganIdx0_Device, m_pRecomb_RandOrganIdx1_Device,
																  MSA_GA_VERTICAL_RECOMB_RATIO, m_pOrganLenDevice[m_curOrganLenIndex], thread.x*grid.x);
#endif 
	dim3 thread(MSA_GA_CUDA_BLOCK_SIZE, MSA_GA_CUDA_BLOCK_SIZE);
	dim3 grid(ceil((float)nRecomb/(thread.x*thread.y)), 1);
	MSAGA_RECOMB_VERTICAL_RECOMB_POS0_UNLIMIT << < grid, thread >> >(m_pVerticalRecomb_pos0_Device, m_pOrganIndex0_Device, m_pOrganIndex1_Device,
																	 m_pRecomb_RandOrgan_Device, m_pRecomb_RandTx_Device,
																	 m_pRecomb_RandOrganIdx0_Device, m_pRecomb_RandOrganIdx1_Device,
																	 MSA_GA_VERTICAL_RECOMB_RATIO, m_pOrganLenDevice[m_curOrganLenIndex], thread.x*grid.x, nRecomb, MSA_GA_POPULATION_SIZE);

	grid.x = ceil((float)nRecomb*m_pPopulation->nSeqSize / (thread.x*thread.y));
	MSAGA_RECOMB_VERTICAL_RECOMB_POS1_UNLIMIT<<< grid, thread >>>(m_pVerticalRecomb_pos0_Device, m_pVerticalRecomb_pos1_Device, m_pOrganIndex0_Device, m_pOrganIndex1_Device,
																	m_pRecomb_RandOrgan_Device, MSA_GA_VERTICAL_RECOMB_RATIO, thread.x*grid.x, 
																	m_pPopulationDevice[m_curPopulationIndex], POPULATION_WIDTH, MSA_GA_CUDA_MAX_LEN, nRecomb, 1.0f / nRecomb,
																	CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), m_pPopulation->nSeqSize);

	//dbgWriteFile("m_pVerticalRecomb_pos0_Device.txt", m_pVerticalRecomb_pos0_Device, nRecomb, nRecomb);
	//dbgWriteFile("m_pVerticalRecomb_pos1_Device.txt", m_pVerticalRecomb_pos1_Device, nRecomb*m_pPopulation->nSeqSize, nRecomb);
	//dbgWriteFile("m_pOrganIndex0_Device.txt", m_pOrganIndex0_Device, nRecomb, nRecomb);
	//dbgWriteFile("m_pOrganIndex1_Device.txt", m_pOrganIndex1_Device, nRecomb, nRecomb);

	grid.x = ceil((float)nRecomb/(thread.x*thread.y));
	//find the min and max value of org1 recomb position 
	MSAGA_RECOMB_VERTICAL_POS1_MIN_MAX_UNLIMIT<32><<< grid, thread>>>( m_pPos1_Min_Device, m_pPos1_Max_Device, m_pVerticalRecomb_pos1_Device, m_pPopulation->nSeqSize,
															           MSA_GA_VERTICAL_RECOMB_RATIO, m_pRecomb_RandOrgan_Device, grid.x*thread.x, nRecomb );

	cudaMemset( (void*)m_pPopulationDevice[1-m_curPopulationIndex], MSA_GA_CUDA_SPACE, m_nCellNum*sizeof(short) );
	grid.x = ceil((float)m_nCellNum/(thread.x*thread.y));
	MSAGA_RECOMB_VERTICAL_UNLIMIT<<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex],  
													  m_pOrganIndex0_Device, m_pOrganIndex1_Device, m_pVerticalRecomb_pos0_Device, m_pVerticalRecomb_pos1_Device,													
													  m_pRecomb_RandOrgan_Device, MSA_GA_VERTICAL_RECOMB_RATIO, m_pPos1_Min_Device, m_pPos1_Max_Device,														
													  CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN, 
													  POPULATION_WIDTH, R_POPULATION_WIDTH, MSA_GA_POPULATION_SIZE, thread.x*grid.x, m_nCellNum);
	m_curPopulationIndex = 1-m_curPopulationIndex;
	//dbgCheckSequence();
	ArrangeSeqDeviceUnlimit( m_pRecomb_RandOrgan_Device, MSA_GA_VERTICAL_RECOMB_RATIO, true );
}

void		CMSAGA_CUDA_Algorithm::RecombinationVerticalDevice()
{
	int	nRecomb = (MSA_GA_POPULATION_SIZE-1)/2;
	RandomNumberGeneratorDevice( m_pRecomb_RandOrgan_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandTx_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandOrganIdx0_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandOrganIdx1_Device, nRecomb );

	dim3 thread;
	dim3 grid;
	if( nRecomb*m_pPopulation->nSeqSize>1024 )
	{
		thread.x = 1024/m_pPopulation->nSeqSize;
		thread.y = m_pPopulation->nSeqSize;
		grid.x = ceil((float)nRecomb*m_pPopulation->nSeqSize/1024);
		grid.y = 1;
	}
	else
	{
		thread.x = nRecomb;
		thread.y = m_pPopulation->nSeqSize;
		grid.x = 1;
		grid.y = 1;
	}
	//find the recomb pos of org0 and org1.
	MSAGA_RECOMB_VERTICAL_RECOMB_POS<32><<< grid, thread>>>( m_pVerticalRecomb_pos0_Device, m_pVerticalRecomb_pos1_Device, m_pOrganIndex0_Device, m_pOrganIndex1_Device,
																														  m_pPopulationDevice[m_curPopulationIndex], m_pRecomb_RandOrgan_Device, m_pRecomb_RandTx_Device,
																														  m_pRecomb_RandOrganIdx0_Device, m_pRecomb_RandOrganIdx1_Device,
																														  MSA_GA_VERTICAL_RECOMB_RATIO, MSA_GA_POPULATION_SIZE, MSA_GA_CUDA_MAX_LEN,
																														  CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), m_pOrganLenDevice[m_curOrganLenIndex], grid.x*thread.x, nRecomb );

	//dbgWriteFile( std::string("Recomb_Organ0.txt"), m_pOrganIndex0_Device, nRecomb, nRecomb );
	//dbgWriteFile( std::string("Recomb_Organ1.txt"), m_pOrganIndex1_Device, nRecomb, nRecomb );
	//dbgWriteFile( std::string("Recomb_pos0.txt"), m_pVerticalRecomb_pos0_Device, nRecomb*m_pPopulation->nSeqSize, nRecomb );
	//dbgWriteFile( std::string("Recomb_pos1.txt"), m_pVerticalRecomb_pos1_Device, nRecomb*m_pPopulation->nSeqSize, nRecomb );
	
	thread.x = nRecomb;
	thread.y = 1;
	grid.x = 1;
	grid.y = 1;
	//find the min and max value of org1 recomb position 
	MSAGA_RECOMB_VERTICAL_POS1_MIN_MAX<32><<< grid, thread>>>( m_pPos1_Min_Device, m_pPos1_Max_Device, m_pVerticalRecomb_pos1_Device, m_pPopulation->nSeqSize,
		                                                                                                                          MSA_GA_VERTICAL_RECOMB_RATIO, m_pRecomb_RandOrgan_Device,  nRecomb );

	//writeDebugFile( std::string("m_pPos1_Min_Device.txt"), m_pPos1_Min_Device, nRecomb );
	//writeDebugFile( std::string("m_pPos1_Max_Device.txt"), m_pPos1_Max_Device, nRecomb );

	//dbgWriteSeqsFromDevice2File( "seqResult0.txt" );
	thread.x = MSA_GA_CUDA_MAX_LEN;
	thread.y = m_pPopulation->nSeqSize;
	grid.x = MSA_GA_POPULATION_SIZE;
	grid.y = 1;
	MSAGA_RECOMB_VERTICAL<32><<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex],  
		                                                                                         m_pOrganIndex0_Device, m_pOrganIndex1_Device, m_pVerticalRecomb_pos0_Device, m_pVerticalRecomb_pos1_Device,
		                                                                                         m_pRecomb_RandOrgan_Device, MSA_GA_VERTICAL_RECOMB_RATIO, 
																								 MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN, MSA_GA_POPULATION_SIZE, m_pPopulation->nSeqSize, 
																								 m_pPos1_Min_Device, m_pPos1_Max_Device,
																								 CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), MSA_GA_CUDA_SPACE );
	m_curPopulationIndex = 1-m_curPopulationIndex;
	
	//dbgWriteSeqsFromDevice2File( "seqResult1.txt" );
	//writeDebugFile( std::string("BX_VERTICAL.txt"), m_pVerticalRecomb_pos1_Device, nRecomb );
	//dbgWriteSeqsFromDevice2File( "seqResult1.txt" );

	//arrange the organs.
	ArrangeSeqDevice( m_pRecomb_RandOrgan_Device, 1 );
}


void		CMSAGA_CUDA_Algorithm::RecombinationHorizentalDeviceUnlimit()
{
	int	nRecomb = (MSA_GA_POPULATION_SIZE-1)/2;
	RandomNumberGeneratorDevice( m_pRecomb_RandOrgan_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandTy_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandOrganIdx0_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandOrganIdx1_Device, nRecomb );

	dim3 thread(MSA_GA_CUDA_BLOCK_SIZE, MSA_GA_CUDA_BLOCK_SIZE );
	dim3 grid(ceil((float)m_nCellNum/(thread.x*thread.y)), 1);
	MSAGA_RECOMB_HORIZEN_UNLIMIT<<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex],
													 m_pRecomb_RandOrgan_Device,m_pRecomb_RandTy_Device, m_pRecomb_RandOrganIdx0_Device, 
													 m_pRecomb_RandOrganIdx1_Device, MSA_GA_HORIZENTAL_RECOMB_RATIO, 
													 MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN, POPULATION_WIDTH, R_POPULATION_WIDTH,
													 thread.x*grid.x, MSA_GA_POPULATION_SIZE, m_OrganismSize, m_nCellNum);
	m_curPopulationIndex = 1-m_curPopulationIndex;
	//dbgCheckSequence();
	ArrangeSeqDeviceUnlimit( m_pRecomb_RandOrgan_Device, MSA_GA_HORIZENTAL_RECOMB_RATIO, true );
}

void		CMSAGA_CUDA_Algorithm::RecombinationHorizentalDevice()
{
	int	nRecomb = (MSA_GA_POPULATION_SIZE-1)/2;
	RandomNumberGeneratorDevice( m_pRecomb_RandOrgan_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandTy_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandOrganIdx0_Device, nRecomb );
	RandomNumberGeneratorDevice( m_pRecomb_RandOrganIdx1_Device, nRecomb );

	dim3 thread(MSA_GA_CUDA_MAX_LEN, m_pPopulation->nSeqSize );
	dim3 grid(MSA_GA_POPULATION_SIZE, 1 );
	MSAGA_RECOMB_HORIZEN<32><<< grid, thread>>>( m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex],
	                                                                                             m_pRecomb_RandOrgan_Device,m_pRecomb_RandTy_Device, m_pRecomb_RandOrganIdx0_Device, m_pRecomb_RandOrganIdx1_Device,
																								 MSA_GA_HORIZENTAL_RECOMB_RATIO, m_pPopulation->nSeqSize,  MSA_GA_POPULATION_SIZE );
	m_curPopulationIndex = 1-m_curPopulationIndex;

	//arrange the organs.
	ArrangeSeqDevice( m_pRecomb_RandOrgan_Device, MSA_GA_HORIZENTAL_RECOMB_RATIO );
}

void		CMSAGA_CUDA_Algorithm::MutationDeviceUnlimit()
{
	MutationGapInsertDeviceUnlimit();
	MutationGapExtensionDeviceUnlimit();
	MutationGapReductionDeviceUnlimit();
}

void		CMSAGA_CUDA_Algorithm::MutationDevice()
{
   //dbgWriteSeqsFromDevice2File( "seqResult0.txt" );
	MutationGapInsertDevice();
	MutationGapExtensionDevice();
	MutationGapReductionDevice();
	//dbgWriteSeqsFromDevice2File( "seqResult1.txt" );
}

void		CMSAGA_CUDA_Algorithm::Evolution()
{
	//The evolution main code with cuda
	//dbgCheckSequence();
	SelectionDeviceUnlimit();
	//dbgCheckSequence();
	MutationDeviceUnlimit();
	//dbgCheckSequence();
	RecombinationDeviceUnlimit();
	//dbgCheckSequence();
}

void     CMSAGA_CUDA_Algorithm::testSPScoreCPU( const COrganism& vOrgan )
{
	//calculate col score
	// Allocate host memory for organs
	unsigned int size_colscore = m_pPopulation->nSeqSize*MSA_GA_CUDA_MAX_LEN;
	unsigned int mem_size_colscore = sizeof(short) * size_colscore;
	short *colscore = (short *)malloc(mem_size_colscore);
	memset(colscore, 0, mem_size_colscore);

	for( int i=0; i<vOrgan.pSequence[0].sequence.getLen(); ++i )
	{
		for( int j=0; j<vOrgan.nSeqSize-1; ++j )
		{
			short score=0;
			for( int k=j+1; k<vOrgan.nSeqSize; ++k )
			{
				int a =  vOrgan.pSequence[j].sequence.getSequenceContext().at(i).m_iCode;
				int b =  vOrgan.pSequence[k].sequence.getSequenceContext().at(i).m_iCode;
				if ( a != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) && b != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
				{
					score += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getSubMatrixScore( MSAGA_MAT_TYPE, a, b );
				}
			}
			colscore[ j*MSA_GA_CUDA_MAX_LEN+i ] = score;
		}
	}
	std::auto_ptr<SeqAnsis::CFileWriter>		pFileWriter( new SeqAnsis::CFileWriter( std::string("seqColScore.txt") ) );
	pFileWriter->OutputVector( colscore,  size_colscore);
}

void		CMSAGA_CUDA_Algorithm::testColSumScoreCPU( const COrganism& vOrgan )
{
	//calculate col score
	// Allocate host memory for organs
	unsigned int size_colscore = m_pPopulation->nSeqSize*MSA_GA_CUDA_MAX_LEN;
	unsigned int mem_size_colscore = sizeof(short) * size_colscore;
	short *colscore = (short *)malloc(mem_size_colscore);
	memset(colscore, 0, mem_size_colscore);

	unsigned int size_colsumscore = MSA_GA_CUDA_MAX_LEN;
	unsigned int mem_size_colsumscore = sizeof(short) * size_colsumscore;
	short *colsumscore = (short *)malloc(mem_size_colsumscore);
	memset(colsumscore, 0, mem_size_colsumscore);

	for( int i=0; i<vOrgan.pSequence[0].sequence.getLen(); ++i )
	{
		for( int j=0; j<vOrgan.nSeqSize-1; ++j )
		{
			short score=0;
			for( int k=j+1; k<vOrgan.nSeqSize; ++k )
			{
				int a =  vOrgan.pSequence[j].sequence.getSequenceContext().at(i).m_iCode;
				int b =  vOrgan.pSequence[k].sequence.getSequenceContext().at(i).m_iCode;
				if ( a != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) && b != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) )
				{
					score += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getSubMatrixScore( MSAGA_MAT_TYPE, a, b );
				}
			}
			colscore[ j*MSA_GA_CUDA_MAX_LEN+i ] = score;
		}
	}

	for( int i=0; i<vOrgan.pSequence[0].sequence.getLen(); ++i )
	{
		for( int j=0; j<vOrgan.nSeqSize; ++j )
		{
			colsumscore[i] += colscore[ j*MSA_GA_CUDA_MAX_LEN+i ];
		}
	}

	std::auto_ptr<SeqAnsis::CFileWriter>		pFileWriter( new SeqAnsis::CFileWriter( std::string("seqColSumScore.txt") ) );
	pFileWriter->OutputVector( colsumscore,  size_colsumscore);
}

int		CMSAGA_CUDA_Algorithm::SPScore(const std::vector<CSequence>& vInputSequences)
{
	int spScore = 0;
	//pair score
	for (int j = 0; j<vInputSequences[0].getLen(); ++j)
	{
		for (int k = 0; k<vInputSequences.size() - 1; ++k)
		{
			for (int m = k + 1; m<vInputSequences.size(); ++m)
			{
				int a = vInputSequences[k].getSequenceContext().at(j).m_iCode;
				int b = vInputSequences[m].getSequenceContext().at(j).m_iCode;
				if (a != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE) && b != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE))
				{
					spScore += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getSubMatrixScore(MSAGA_MAT_TYPE, a, b);
				}
			}
		}
	}
	//gap penalty
	for (int j = 0; j<vInputSequences.size(); ++j)
	{
		int gapNum = 0;
		int idx = 0;
		while (idx<vInputSequences[j].getLen())
		{
			while (idx<vInputSequences[j].getLen() && vInputSequences[j].getSequenceContext().at(idx).m_iCode != CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE))	++idx;
			while (idx<vInputSequences[j].getLen() && vInputSequences[j].getSequenceContext().at(idx).m_iCode == CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE))
			{
				++gapNum;
				++idx;
			}
			if (gapNum>0)
			{
				spScore += CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(MSAGA_MAT_TYPE) + (gapNum - 1)*CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapExtendCost(MSAGA_MAT_TYPE);
			}
			gapNum = 0;
		}
	}
	return spScore;
}

int		CMSAGA_CUDA_Algorithm::SPScore( const COrganism& vOrgan )
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

void		CMSAGA_CUDA_Algorithm::ArrangeSeqDevice(float*	pRandOrgan_Device, float randLimit)
{

}

void		CMSAGA_CUDA_Algorithm::ArrangeSeqDeviceUnlimit( float*	pRandOrgan_Device, float randLimit,	bool bRecomb )
{
	//对要进行空格去除的下一代个体存储空间全部置为未使用
	cudaMemset( (void*)m_pPopulationDevice[1-m_curPopulationIndex], MSA_GA_CUDA_SPACE, m_nCellNum*sizeof(short) );

	//在较短序列尾部填充空格，补齐长度
	dim3 thread(MSA_GA_CUDA_BLOCK_SIZE, MSA_GA_CUDA_BLOCK_SIZE);
	dim3 grid( ceil( (float)m_nCellNum/(thread.x*thread.y) ), 1 );
	MSAGA_FILL_SPACE_AT_TAIL_UNLIMIT<<< grid, thread>>>( bRecomb, m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex], 
														 pRandOrgan_Device, randLimit, thread.x*grid.x, 													
														 CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), 													
														 MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN,
														 POPULATION_WIDTH, R_POPULATION_WIDTH,
														 m_pPopulation->nSeqSize, m_nCellNum, MSA_GA_POPULATION_SIZE);
	m_curPopulationIndex = 1-m_curPopulationIndex;
	//dbgCheckSequence();
	//记录全为空格的列	
	thread.x=MSA_GA_CUDA_MAX_LEN;	
	thread.y=1;	
	grid.x=MSA_GA_POPULATION_SIZE;	
	MSAGA_GAP_RECORDER<32><<< grid, thread>>>( bRecomb, m_pGapRecorder_Device, m_pPopulationDevice[m_curPopulationIndex], pRandOrgan_Device, randLimit, 	
											   CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE),												 	
											   MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN, m_OrganismSize, MSA_GA_POPULATION_SIZE );													

	//删除全列为空格的列(不处理尾部空格)
	thread.x = MSA_GA_CUDA_BLOCK_SIZE;	
	thread.y = thread.x;	
	grid.x = ceil( (float)m_nCellNum/(thread.x*thread.y) );	
	grid.y = 1;		
	cudaMemset( (void*)m_pPopulationDevice[1-m_curPopulationIndex], MSA_GA_CUDA_SPACE, m_nCellNum*sizeof(short) );
	MSAGA_DELETE_SPACE_UNLIMIT<<< grid, thread>>>( bRecomb, m_pPopulationDevice[1-m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex],	
												   pRandOrgan_Device, randLimit, m_pGapRecorder_Device, thread.x*grid.x,
												   MSA_GA_POPULATION_SIZE, POPULATION_WIDTH, R_POPULATION_WIDTH, MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN, m_nCellNum);
	m_curPopulationIndex = 1-m_curPopulationIndex;	
	//dbgCheckSequence();
	//单独删除尾部空格
#if 0
	thread.x = MSA_GA_CUDA_BLOCK_SIZE;
	thread.y = thread.x;
	grid.x = ceil((float)MSA_GA_CUDA_MAX_LEN*m_pPopulation->nSeqSize*MSA_GA_POPULATION_SIZE / (thread.x*thread.y));
	grid.y = 1;

	cudaMemset((void*)m_pPopulationDevice[1 - m_curPopulationIndex], MSA_GA_CUDA_SPACE, m_nCellNum*sizeof(short));
	MSAGA_DELETE_SPACE_AT_TAIL_UNLIMIT <<< grid, thread >>>(bRecomb, m_pPopulationDevice[1 - m_curPopulationIndex], m_pPopulationDevice[m_curPopulationIndex],
															pRandOrgan_Device, randLimit, thread.x*grid.x, POPULATION_WIDTH, R_POPULATION_WIDTH,
															MSA_GA_CUDA_MAX_LEN, R_MSA_GA_CUDA_MAX_LEN, m_pPopulation->nSeqSize, 
															CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode(GENESPACE), m_nCellNum, MSA_GA_POPULATION_SIZE);
	m_curPopulationIndex = 1 - m_curPopulationIndex;
#endif
	//dbgCheckSequence();
	//重新计算整理后的个体长度
#if 0
	thread.x = MSA_GA_CUDA_BLOCK_SIZE;
	thread.y = MSA_GA_CUDA_BLOCK_SIZE;
#endif
	grid.x = ceil((float)MSA_GA_POPULATION_SIZE / (thread.x*thread.y));
	grid.y=1;
	MSAGA_CALC_SEQ_LEN_UNLIMIT<32><<< grid, thread>>>( bRecomb, m_pOrganLenDevice[1-m_curOrganLenIndex], m_pOrganLenDevice[m_curOrganLenIndex], 
													   m_pPopulationDevice[m_curPopulationIndex], pRandOrgan_Device, randLimit, MSA_GA_CUDA_MAX_LEN, 
													   grid.x*thread.x, MSA_GA_POPULATION_SIZE);
	m_curOrganLenIndex = 1-m_curOrganLenIndex;
}

//==================FOR CPU VERIFICATION PURPOSE=====================//
void		CMSAGA_CUDA_Algorithm::Mutation()
{
	GapInsertMutation();
	//GapExtensionMutation();
	//GapReductionMutation();
}

void		CMSAGA_CUDA_Algorithm::GapInsertMutation()
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

int        CMSAGA_CUDA_Algorithm::getSeqLen( const CSequence& Seq )
{
	const CSeqData& SeqData = Seq.getSequenceContext();
	int seqLen = SeqData.size();
	int count=0;
	for( int i=0; i<seqLen; ++i )
	{
		if( SeqData[i].m_char != GENESPACE&&SeqData[i].m_char!=NONEGENE )
		{
			++count;
		}
	}
	return count;
}

void     CMSAGA_CUDA_Algorithm::dbgCheckSequence()
{
	unsigned int size_organs = m_pPopulation->nSeqSize*MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	unsigned int mem_size_organs = sizeof(short) * size_organs;
	short *h_organs = (short *)malloc(mem_size_organs);
	//read sequences from device to host
	CUDA_DEV2HOST(   h_organs, m_pPopulationDevice[m_curPopulationIndex], mem_size_organs );

	//每一个个体中的每一个序列进行长度校验
	int nWidth = MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	for( int k=0; k<MSA_GA_POPULATION_SIZE; ++k )
	{
		for( int i=0; i<m_pPopulation->nSeqSize; ++i )
		{
			int nSeqLen=0;
			for ( int j=0; j<MSA_GA_CUDA_MAX_LEN; ++j )
			{
				if( h_organs[ k*MSA_GA_CUDA_MAX_LEN+i*nWidth+j ] != MSA_GA_CUDA_SPACE )
				{
					char acid = CGlobalSpace::m_sAlignParams.getAminoAcidInt2CharCode( h_organs[k*MSA_GA_CUDA_MAX_LEN+i*nWidth+j ] );
					if( acid!='-' )
					{
						++nSeqLen;
					}
				}
				else
				{
					while( j<MSA_GA_CUDA_MAX_LEN && h_organs[ k*MSA_GA_CUDA_MAX_LEN+i*nWidth+j ] == MSA_GA_CUDA_SPACE )
					{
						++j;
					}
					if( j<MSA_GA_CUDA_MAX_LEN && h_organs[ k*MSA_GA_CUDA_MAX_LEN+i*nWidth+j ] != MSA_GA_CUDA_SPACE )
					{
						char ch[256];
						sprintf_s(ch,  "The pos %i of seq %i of the population %i is checksum error!" , j, i,  k);
						CGlobalSpace::m_sEventLog.writeEvent(ch);
						dbgWriteSeqsFromDevice2File( "checkError.txt" );
					}
				}
			}
			if( nSeqLen!=m_OrginSequences[i].getLen() )
			{
				char ch[256];
				sprintf_s(ch,  "The seq %i of the population %i is wrong, seq len in device=%i, org seq len=%i  The evolution run is %i" , i,  k, nSeqLen, m_OrginSequences[i].getLen(), m_iRun);
				CGlobalSpace::m_sEventLog.writeEvent(ch);
				dbgWriteSeqsFromDevice2File( "checkError.txt" );
			}
		}
	}

	free(h_organs);
	h_organs=NULL;

	//dbgWriteSeqsFromDevice2File( "checkPass.txt" );
}

void     CMSAGA_CUDA_Algorithm::dbgWriteSeqsFromDevice2File( const std::string& filename )
{
	std::fstream   fileOut;
	if( !fileOut.is_open() )
	{
		try
		{
			fileOut.open( filename.c_str(), std::ios::out|std::ios::trunc );
		}
		catch (CAppException const& e)
		{
			throw	e;
		}
	}
	unsigned int size_organs = m_pPopulation->nSeqSize*MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	unsigned int mem_size_organs = sizeof(short) * size_organs;
	short *h_organs = (short *)malloc(mem_size_organs);
	//read sequences from device to host
	CUDA_DEV2HOST(   h_organs, m_pPopulationDevice[m_curPopulationIndex], mem_size_organs );

	//将每一个个体打印到文件
	int nWidth = MSA_GA_CUDA_MAX_LEN*MSA_GA_POPULATION_SIZE;
	for( int k=0; k<MSA_GA_POPULATION_SIZE; ++k )
	{
		fileOut<<"==========Population  "<<k<<"=========="<<std::endl;
		for( int i=0; i<m_pPopulation->nSeqSize; ++i )
		{
			for ( int j=0; j<MSA_GA_CUDA_MAX_LEN; ++j )
			{
				if( h_organs[ k*MSA_GA_CUDA_MAX_LEN+i*nWidth+j ] != MSA_GA_CUDA_SPACE )
				{
					char acid = CGlobalSpace::m_sAlignParams.getAminoAcidInt2CharCode( h_organs[k*MSA_GA_CUDA_MAX_LEN+i*nWidth+j ] );
					fileOut<<acid;
				}
				else
				{
					break;
				}
			}
			fileOut<<std::endl;//换行输出同一个个体的下一个序列
		}
		fileOut<<std::endl<<std::endl;//空一行再输出下一个个体
	}
	fileOut.close();
	free(h_organs);
	h_organs=NULL;
}

};
