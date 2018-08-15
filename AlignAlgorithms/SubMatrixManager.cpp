#include "SubMatrices.h"
#include "SubMatrixManager.h"

namespace	SeqAnsis
{

#define	SubMatPos(r,c)	(r)*(r+1)/2+(c)

short    SubMatrixManager::getSubMatrixSize(SubMatrixType vType)
{
	short  size=0;
	switch( vType )
	{
	case DNASUBMAT:
		size=8;
		break;
	case PAM20:
	case PAM60:
	case PAM250:
	case PAM350:
	case BLOSUM30:
	case BLOSUM45:
	case BLOSUM62:
	case BLOSUM80:
		size=276;
		break;
	}
	return size;
}

short *     SubMatrixManager::getSubMatrix(SubMatrixType vType)
{
	short *pMat = NULL;
	switch( vType )
	{
	case DNASUBMAT:
		pMat=dnaSubMat;
		break;
	case PAM20:
		pMat = pam20mt;
		break;
	case PAM60:
		pMat = pam60mt;
		break;
	case PAM250:
		pMat=pam250mt;
		break;
	case PAM350:
		pMat = pam350mt;
		break;
	case BLOSUM30:
		pMat=blosum30mt;
		break;
	case BLOSUM45:
		pMat=blosum45mt;
		break;
	case BLOSUM62:
		pMat=blosum62mt2;
		break;
	case BLOSUM80:
		pMat = blosum80mt;
		break;
	}
	return pMat;
}

float		SubMatrixManager::getGapOpenCost(SubMatrixType vType)
{
	float	cost = 0;
	switch( vType )
	{
	case DNASUBMAT:
		cost = -1.0f;
		break;
	case PAM20:
		cost = -10.0f;
		break;
	case PAM60:
		cost = -20.0f;
		break;
	case PAM250:
		cost = -10.0f;
		break;
	case PAM350:
		cost = -10.0f;
		break;
	case BLOSUM30:
		cost = -10.0f;
		break;
	case BLOSUM45:
		cost = -10.0f;
		break;
	case BLOSUM62:
		cost = -10.0f;
		break;
	case BLOSUM80:
		cost = -10.0f;
		break;
	}
	return cost;
}

float		SubMatrixManager::getGapExtendCost(SubMatrixType vType)
{
	float	cost = 0;
	switch( vType )
	{
	case DNASUBMAT:
		cost = -1.0f;
		break;
	case PAM20:
		cost = -0.2f;
		break;
	case PAM60:
		cost = -0.5f;
		break;
	case PAM250:
		cost = -0.2f;
		break;
	case PAM350:
		cost = -0.2f;
		break;
	case BLOSUM30:
		cost=-0.2f;
		break;
	case BLOSUM45:
		cost=-0.2f;
		break;
	case BLOSUM62:
		cost = -0.5f;
		break;
	case BLOSUM80:
		cost = -0.2f;
		break;
	}
	return cost;
}

float		SubMatrixManager::getSubMatrixScore(SubMatrixType vType, int aCode, int bCode)
{
	switch( vType )
	{
	case DNASUBMAT:
		m_pCurSubMat = dnaSubMat;
		break;
	case PAM20:
		m_pCurSubMat = pam20mt;
		break;
	case PAM60:
		m_pCurSubMat = pam60mt;
		break;
	case PAM250:
		m_pCurSubMat = pam250mt;
		break;
	case PAM350:
		m_pCurSubMat = pam350mt;
		break;
	case BLOSUM30:
		m_pCurSubMat = blosum30mt;
		break;
	case BLOSUM45:
		m_pCurSubMat = blosum45mt;
		break;
	case BLOSUM62:
		m_pCurSubMat = blosum62mt2;
		break;
	case BLOSUM80:
		m_pCurSubMat = blosum80mt;
		break;
	}

	float	score;
	if ( aCode<bCode )
	{
		score = (float)m_pCurSubMat[SubMatPos(bCode, aCode)];
	}
	else
	{
		score = (float)m_pCurSubMat[SubMatPos(aCode, bCode)];
	}
	return	score;
}

}