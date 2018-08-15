//This is the implementation of MSA-GA algorithm based on the paper of 
//2007=A simple genetic algorithm for multiple sequence alignment=Genetics and Molecular Research
#pragma once

//#define		MSA_GA_REVERSION
#include <curand.h>

#include "FileWriter.h"
#include"AlignAlgorithmBase.h"
#include"Sequence.h"
#include"sortingNetworks_common.h"
#include "MSA-GA.h"

namespace SeqAnsis
{
#if 0
	const int  MSA_GA_CUDA_MAX_LEN=256;
	const float R_MSA_GA_CUDA_MAX_LEN = 1.0f / MSA_GA_CUDA_MAX_LEN;
	const int POPULATION_WIDTH = MSA_GA_POPULATION_SIZE*MSA_GA_CUDA_MAX_LEN;
	const float R_POPULATION_WIDTH = 1.0f / POPULATION_WIDTH;
#endif
	const short  MSA_GA_CUDA_SPACE = -1;//sign the position as not used one.
	const short MSA_GA_CUDA_BLOCK_SIZE = 32;
	const short  MSA_GA_CUDA_SELECTION_CANDIDATE = 10;
	const short WEIGHT_SCALE = 10;

	class CMSAGA_CUDA_Algorithm:public CAlignAlgorithmBase
	{
	public:
		CMSAGA_CUDA_Algorithm();
		virtual void    Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences );
		void	SetAlignParams(int nPopulationNum, int nNoAdvGenerationNum, int nMaxOrgLen, const std::vector<float>& SeqWeight );
		~CMSAGA_CUDA_Algorithm();
	private:
		void		InitPopulation( const std::vector<CSequence>& vSequences );
		void		InitDevice();
		void		InitFitness();
		void		InitSelection();
		void     InitMutation();
		void		InitRecomb();
		void     InitArrangement();

		void     ClearDevice();

		bool		FitnessDevice();
		bool		FitnessDeviceUnlimit();
		
		void		MutationDevice();
		void		MutationDeviceUnlimit();
		void     MutationGapInsertDevice();
		void     MutationGapInsertDeviceUnlimit();
		void		MutationGapExtensionDevice();
		void		MutationGapExtensionDeviceUnlimit();
		void     MutationGapReductionDevice();
		void     MutationGapReductionDeviceUnlimit();
		
		void		RecombinationDevice();
		void		RecombinationDeviceUnlimit();
		void		RecombinationHorizentalDevice();
		void		RecombinationHorizentalDeviceUnlimit();
		void		RecombinationVerticalDevice();
		void		RecombinationVerticalDeviceUnlimit();

		void     ArrangeSeqDevice( float*	pRandOrgan_Device, float randLimit );
		void     ArrangeSeqDeviceUnlimit(float*	pRandOrgan_Device, float randLimit, bool bRecomb=false);

		void		ReadAlignedSeqsFromDevice( std::vector<CSequence>& vAlignedSequences );

		int		SPScore( const COrganism& vOrgan );
		int		SPScore2009( const COrganism& vOrgan );
		int		SPScore(const std::vector<CSequence>& vInputSequences);
		void		Evolution();
		void		SelectionDevice();
		void		SelectionDeviceUnlimit();
		void		Selection2004();
		void		Selection2009();
		void		Recombination();
		void		Recombination2004();
		void		HorizentalRecombination();
		void		VerticalRecombination();
		void		Recombination2009();
		void		HorizentalRecombination2009();
		void		VerticalRecombination2009();
		void		Mutation2004();
		void		GapInsertMutation2004();
		void		GapExtensionMutation2004();
		void		GapReductionMutation2004();
		void		Mutation();
		void		GapInsertMutation();
		void		GapExtensionMutation();
		void		GapReductionMutation();
		int        SelectRandomGap(CSeqData& vSeqSelected);
		void     ResizeOrgan( COrganism& vOrgan );//arrange the  sequences into equal length sequences by adding GENSPACE character at the tail of the shorter sequences.
		void		ArrangeSequences();
		int		GetLongestSeqLen( const COrganism& og );

		bool    testOrgLen( COrganism* vpPopulation );
		int		CountSeqCharNum( const CSequence& vSeq );
		bool		CheckSumCharNum( const COrganism& vOG );

		void		SetSubstitutionMat2GPU( SubMatrixType t );

		void		CUDA_MALLOC( void **d_mem, void *h_mem, unsigned int size_mem );
		void		CUDA_MALLOC( void **d_mem, unsigned int size_mem );

		void		CUDA_FREE( void *d_mem );

		void		CUDA_DEV2HOST( void *h_mem, const void * const d_mem,  unsigned int size_mem  );

		void		testSPScoreCPU( const COrganism& vOrgan );
		void		testColSumScoreCPU( const COrganism& vOrgan );

		template<typename T> void	dbgWriteFile(const std::string& fileName, const T * const d_data, int dataSize, int nDataPerLine)
		{
			int nLine = (dataSize + nDataPerLine - 1) / nDataPerLine;
			unsigned int dataMemSize = sizeof(T)*dataSize;
			T *h_data = (T*)malloc(dataMemSize);
			CUDA_DEV2HOST(h_data, d_data, dataMemSize);
			std::auto_ptr<CFileWriter>		pFileWriter(new CFileWriter(fileName));
			pFileWriter->openFile();
			for (int i = 0; i<(nLine - 1); ++i)
			{
				pFileWriter->OutputVector(h_data + i*nDataPerLine, nDataPerLine);
			}

			//最后一行
			pFileWriter->OutputVector(h_data + (nLine - 1)*nDataPerLine, dataSize - nDataPerLine*(nLine - 1));
			pFileWriter->closeFile();
			free(h_data);
			h_data = NULL;
		}

		void     RandomNumberGeneratorDevice( float *pRandNum, int nRand );
		int        getSeqLen( const CSequence& Seq );
		void     dbgCheckSequence();
		void     dbgWriteSeqsFromDevice2File( const std::string& filename );

		std::vector<CSequence>		m_OrginSequences;
		COrganism   *m_pPopulation;
		int                m_OrganismSize;//record the sequence number in each organism.
		int				m_iLongestSeqSize;
		int				m_iRun;
		int                m_fMaxScore;

		int                m_idbgHcombLen;
		int                m_idbgShortestInputAcidNum;
		bool             m_bdbgFirstRun;

		curandGenerator_t		m_RandomNumberGen;
		short *m_pSubMatDevice;
		short *m_pWeightSeq;
		short *m_pPopulationDevice[2];
		int	   m_curPopulationIndex;
		SubMatrixType		m_curSubMatType;
		float*		m_pOrganScore_Align_Device;
		uint*		m_pOrganIndex_Align_Device;
		float*		m_pOrganScore_Align_Sorted_Device;
		uint*		m_pOrganIndex_Align_Sorted_Device;

		float		m_bestScore;

		float*		m_pUnitScore_Device;
		float*		m_pColumnScore_Device;
		float*		m_pOrganScore_Device;

		int			m_nOrganSize_Align;
		int         m_nCellNum;//实际使用的计算单元个数，=POPULATION_WIDTH*序列个数

		unsigned int	*m_pOrganLenDevice[2];
		int			m_curOrganLenIndex;

		float *		m_pRandomOrganIdx_Device;
		unsigned int		m_nRandom_organ_idx;
		uint*		m_pSelectedOrganIdx_Device;

		float*		m_pMutation_RandOrgan_Device;
		float*		m_pMutation_RandTx_Device;
		float*		m_pMutation_RandTy_Device;

		float*		m_pRecomb_RandOrgan_Device;
		float*		m_pRecomb_RandTx_Device;
		float*		m_pRecomb_RandTy_Device;
		float*		m_pRecomb_RandOrganIdx0_Device;
		float*		m_pRecomb_RandOrganIdx1_Device;

		short*		m_pGapRecorder_Device;
		short*		m_pRangeRecorder_Device;
		short*      m_pGapLenRecorder_Device;

		short*      m_pGapStartRecorder_Device;
		short*      m_pGapPos_Tx_Device;

		short*		m_pPos1_Min_Device;
		short*		m_pPos1_Max_Device;
		short*      m_pOrganIndex0_Device;
		short*      m_pOrganIndex1_Device;
		short*		m_pVerticalRecomb_pos0_Device;
		short*		m_pVerticalRecomb_pos1_Device;

		float*		m_pSeqWeight_Device;

		/////////////////////
		int m_noAdvGenerationNum;
		std::vector<float>	m_SeqWeight;
		float m_rnPairs;
		//////////////////////

		////////////////////
		int MSA_GA_POPULATION_SIZE;
		int  MSA_GA_CUDA_MAX_LEN;
		float R_MSA_GA_CUDA_MAX_LEN;
		int POPULATION_WIDTH;
		float R_POPULATION_WIDTH;
		int MSA_GA_NO_ADV_GENERATION_NUM;

		const int MSA_GA_GENERATION_NUM = 1024;
		const float	MSA_GA_SCALING_FACTOR_K = 1.2f;
		const float MSA_GA_OFFSET_FACTOR_X = 0.2f;
		const float MSA_GA_HORIZENTAL_RECOMB_RATIO = 0.3f;
		const float MSA_GA_VERTICAL_RECOMB_RATIO = 0.5f;
		const float MSA_GA_BLOCK_MUTATION = 0.1f;
		const float MSA_GA_EXTENSION_MUTATION = 0.05f;
		const float MSA_GA_REDUCTION_MUTATION = 0.05f;
		const int   MSA_GA_MAX_INSERT_GAP = 3;
		const float MSA_GA_RESERVED_RATIO = 0.1f;

	};

}