#include <string>
#include "Chromosome.h"
#include "AppException.h"

namespace	SeqAnsis
{

CChromosome::CChromosome():m_geneHeadLen(0)
{
}

CChromosome::~CChromosome()
{
	SAFE_DELETE_ARRAY(m_pPath);
}

void	CChromosome::init( const std::string& vFuncSet, const std::string& vTerminalSet, int vGeneNum, int vGeneHeadLen )
{
	if(vGeneHeadLen<0 || vGeneNum<=0 || vTerminalSet.empty())
	{
		throw	CAppException( DEF_EXCEPTION_INVALID_PARAMETER, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "Chromosome init param is illegal in SAGEP algorithm!" );
	}

	m_geneHeadLen = vGeneHeadLen;
	m_geneNum = vGeneNum;

	//calculate gene length
	m_geneLen = m_geneHeadLen*((int)(vTerminalSet.length())-1)+1;



	for ( int j=0; j<GENE_HEAD_LEN; ++j )
		{
			g_Path[i].Path[j+k*GENE_LEN] = g_HeadSet[rand()%HEADSET_LEN];
		}
		for ( int j=GENE_HEAD_LEN; j<GENE_LEN; ++j )
		{
			g_Path[i].Path[j+k*GENE_LEN] = g_TailSet[rand()%TAILSET_LEN];
		}
	}
	g_Path[i].score = g_Path[i].adjustScore = -INT_MAX;	
}

CChromosome& CChromosome::operator=( const CChromosome& left )
{
	if( this == &left )	return *this;

	m_x = left.x();
	m_y = left.y();
	m_score = left.score();
	m_chromosomeLength = left.getChromosomeLength();
	left.getPath( m_pPath, m_chromosomeLength );
	return *this;
}

bool CChromosome::getPath( char *pPath, int iPathLen ) const
{
	if( NULL==pPath || iPathLen<m_chromosomeLength )  return false;

	memcpy( pPath, m_pPath, m_chromosomeLength );
	return true;
}

}