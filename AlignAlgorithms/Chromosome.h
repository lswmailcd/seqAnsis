#pragma once

#include <memory>

namespace	SeqAnsis
{

const char	FuncSet[FUNCSET_LEN] = {'A'};	
const char	HeadSet[HEADSET_LEN] = {'A','R','B','D'};	
const char	TailSet[TAILSET_LEN] = {'R','B','D'};	

class CChromosome
{
public:
	CChromosome();
	~CChromosome();
	CChromosome& operator=( const CChromosome& vLeft );

	void	init( const std::string& vFuncSet, const std::string& vTerminalSet, int vGeneNum, int vGeneHeadLen );
	
	int		score() const {return m_score;}
	void	score( int s ){ m_score=s; }
	int		getChromosomeLength() const { return m_chromosomeLength; }
	bool	getPath( char *pPath, int iPathLen ) const;

private:
	std::auto_ptr<char>     m_pData;//the character expression of chromosome.
	int		m_geneHeadLen;
	int     m_geneNum;



	int     m_headLength;

	int      m_score;     //the chromosome original score after fitness function.
	
	int      m_tailLength;
	int      m_geneNum;
	int	     m_chromosomeLength;
};

}