#pragma once

#include <vector>
#include <string>

namespace SeqAnsis
{

typedef std::vector<short> CrossRef;
typedef std::vector<short> Matrix;

typedef enum
{
	DNASUBMAT,
	PAM20,
	PAM60,
	PAM250,
	PAM350,
	BLOSUM30,
	BLOSUM45,
	BLOSUM62,
	BLOSUM80
}SubMatrixType;

class SubMatrixManager
{
	public:
        /* Functions */
		SubMatrixManager(){};
		~SubMatrixManager(){};
        
		float	getSubMatrixScore(SubMatrixType vType, int aCode, int bCode);
		float    getGapOpenCost(SubMatrixType vType);
		float    getGapExtendCost(SubMatrixType vType);
		short   *getSubMatrix(SubMatrixType vType);
		short    getSubMatrixSize(SubMatrixType vType);
	private:
		short*	m_pCurSubMat;
};

}

