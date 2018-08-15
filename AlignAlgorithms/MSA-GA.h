//This is the implementation of MSA-GA algorithm based on the paper of 
//2007=A simple genetic algorithm for multiple sequence alignment=Genetics and Molecular Research
#pragma once

//#define		MSA_GA_REVERSION

#include"AlignAlgorithmBase.h"
#include"Sequence.h"

namespace SeqAnsis
{
#define  MSAGA_MAT_TYPE   BLOSUM62//PAM250//

//#define CLUSTALW_SP_SCORE
//#define WEIGHT_SUM_SP_SCORE

	const int  MSA_GA_POPULATION_SIZE_1 = 2048;//must bigger than 1024

	class CMSA_GASeq
	{
	public:
		int      gapOffset;
		CSequence	sequence;
		int      charNum;//the amino acid count number
	};

	class COrganism
	{
	public:
		CMSA_GASeq	  *pSequence;
		int       nSeqSize;			//the number of the Sequences in one organism
		float      score;				//the score of this organism
		int      adjustScore;   //the score for selection process to avoid negative score

		~COrganism()
		{
			SAFE_DELETE_ARRAY(pSequence);
		}

		COrganism()
		{
			score=0;
			adjustScore = 0;
			nSeqSize=-INT_MAX;
			pSequence = NULL;
		}

		COrganism(const COrganism& vIn)
		{
			nSeqSize = vIn.nSeqSize;
			score = vIn.score;
			adjustScore = vIn.adjustScore;

			if(NULL!=pSequence)  SAFE_DELETE_ARRAY(pSequence);

			pSequence = new CMSA_GASeq[ nSeqSize ];
			for ( int i=0; i<nSeqSize; ++i )
			{
				pSequence[i] = vIn.pSequence[i];
			}
		}

		COrganism& operator=( const COrganism& vRight )
		{
			if( this == &vRight )	return *this;

			nSeqSize = vRight.nSeqSize;
			score = vRight.score;
			adjustScore = vRight.adjustScore;

			if(NULL!=pSequence)  SAFE_DELETE_ARRAY(pSequence);

			pSequence = new CMSA_GASeq[ nSeqSize ];
			for ( int i=0; i<nSeqSize; ++i )
			{
				pSequence[i] = vRight.pSequence[i];
			}

			return *this;
		}
	};

	class		CCompSPS
	{
	public:
		int sps;
		COrganism *pOrgan;
		friend int operator-(const CCompSPS& lhs, const CCompSPS& rhs)
		{
			return	lhs.sps - rhs.sps;
		}
	};

class CMSAGAAlgorithm:public CAlignAlgorithmBase
{
public:
	CMSAGAAlgorithm();
	virtual void    Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences );
	void	SetAlignParams(int nPopulationNum, const std::vector<float>& SeqWeight);
	~CMSAGAAlgorithm();
	void		testSPScore( const std::vector<CSequence>& vSequences );
private:
	void		InitPopulation2007( const std::vector<CSequence>& vSequences, int vCol );
	void		InitPopulation2009( const std::vector<CSequence>& vSequences, int vCol );
	void		InitPopulationREV(const std::vector<CSequence>& vSequences, int vCol );
	bool    Fitness2009();
	bool    Fitness();
	int		SPScore( const COrganism& vOrgan );
	float		SPScore2009( const COrganism& vOrgan );
	void		Evolution();
	void		Selection();
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
	void		dbgCheckSequences(const std::vector<CSequence>& vSequences);

	COrganism   *m_pPopulation;
	int                m_OrganismSize;//record the sequence number in each organism.
	int				m_iLongestSeqSize;
	int				m_iRun;
	float                m_fMaxScore;

	int                m_idbgHcombLen;
	int                m_idbgShortestInputAcidNum;
	bool             m_bdbgFirstRun;

	int MSA_GA_POPULATION_SIZE;
	std::vector<float>	m_SeqWeight;

	const int MSA_GA_GENERATION_NUM = 4096;
	const int MSA_GA_NO_ADV_GENERATION_NUM = 1000;	
	const float	MSA_GA_SCALING_FACTOR_K = 1.2f;
	const float MSA_GA_OFFSET_FACTOR_X = 0.2f;
	const float MSA_GA_HORIZENTAL_RECOMB_RATIO = 0.2f;
	const float MSA_GA_VERTICAL_RECOMB_RATIO = 0.5f;
	const float MSA_GA_BLOCK_MUTATION = 0.2f;
	const float MSA_GA_EXTENSION_MUTATION = 0.05f;
	const float MSA_GA_REDUCTION_MUTATION = 0.05f;
	const int   MSA_GA_MAX_INSERT_GAP = 3;
	const float MSA_GA_RESERVED_RATIO = 0.1f;
};

}