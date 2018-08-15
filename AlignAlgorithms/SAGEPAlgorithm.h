//This is the implementation of needleman_wunsch algorithm
#pragma once

#include"AlignAlgorithmBase.h"
#include"Sequence.h"

namespace SeqAnsis
{

#define  SAGEP_MAT_TYPE   BLOSUM62

const int GEP_NO_ADV_GENERATION_NUM = 50;
const int POPULATION_NUM =200;	
const int EVOLUTION_CIRCLE_NUM = 1;	

//mutation rate	
const float MUTATION_RATE_COORD1 = 0.2f;	
const float MUTATION_RATE_COORD2 = 0.02f;	
const float MUTATION_RATE_GENE = 0.05f;	

//transposition rate	
const float TRANS_RATE_IS = 0.1f;	
const float TRANS_RATE_RIS = 0.1f;	
const float TRANS_RATE_GENE = 0.1f;	

//recombination rate	
const float RECOMB_RATE_P1 = 0.3f;	
const float RECOMB_RATE_P2 = 0.3f;	
const float RECOMB_RATE_GENE = 0.1f;	

const int EVOLUTION_GENERATION_NUM = 2000;	

const int MIN_SCORE = 10;	

const int GENE_HEAD_LEN = 15;	
const int MAX_FUNC_OP_NUM = 2;	
const int GENE_TAIL_LEN = GENE_HEAD_LEN*(MAX_FUNC_OP_NUM-1) + 1;	
const int GENE_LEN = GENE_HEAD_LEN + GENE_TAIL_LEN;	
const int GENE_NUM_PER_CHROMOSOME =10;
const int CHROMOSOME_LEN = GENE_NUM_PER_CHROMOSOME*GENE_LEN; 	

class CChromosome
{
public:
	int      x,y;//the start position of the path.
	char     Path[CHROMOSOME_LEN];
	int      score;     //the chromosome original score after fitness function.
	int      adjustScore;  //the chromosome score adjusted for selection operation.

	CChromosome()
	{
		score=-INT_MAX;
		adjustScore=-INT_MAX;
	}

	CChromosome(const CChromosome& vIn)
	{
		x = vIn.x;
		y = vIn.y;
		score = vIn.score;
		adjustScore = vIn.adjustScore;
		for ( int i=0; i<CHROMOSOME_LEN; ++i )
		{
			Path[i] = vIn.Path[i];
		}
	}

	CChromosome& operator=( const CChromosome& vRight )
	{
		if( this == &vRight )	return *this;

		x = vRight.x;
		y = vRight.y;
		score = vRight.score;
		adjustScore = vRight.adjustScore;
		for ( int i=0; i<CHROMOSOME_LEN; ++i )
		{
			Path[i] = vRight.Path[i];
		}
		return *this;
	}
};

class CSAGEPAlgorithm:public CAlignAlgorithmBase
{
public:
	virtual void    Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences );
private:
	void	initPopulation(int vX, int vY);
	bool	fitness( const std::vector<CSequence>& vSequences );
	void	selection();
	void	reproduction();
	void    mutation();
	void	transposition();
	void	transposition_IS();
	void	transposition_RIS();
	void	transposition_Gene();
	void	recombination();
	void	recombinationP1();
	void	recombinationP2();	
	void	recombinationGene();

	void    translateChromosome( int& voRow, int& voCol, std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences );

	std::vector<std::shared_ptr<CChromosome>>	m_Population;
	std::string		m_FuncSet;
	std::string		m_HeadSet;
	std::string		m_TailSet;

	int				m_iMaxScore;

};

}