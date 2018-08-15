#pragma once

#include <string>
#include <hash_map>
#include "SubMatrixManager.h"

namespace SeqAnsis
{

const char NONEGENE='*';

class CAlignParams
{
	public:
		CAlignParams();
		const std::string& getAminoAcidCodes(){ return m_AminoAcidCodes; }
		int	getAminoAcidCodesNum(){return (int)m_AminoAcidCodes.length()-2;}
		short	getAminoAcidChar2IntCode( char a );
		char					getAminoAcidInt2CharCode( short  iCode );

		void setMaxAllowedSeqLength(int num){m_maxAllowedSeqLength = num;}
        int getMaxAllowedSeqLength(){return m_maxAllowedSeqLength;}
		unsigned long	getUniqueSequenceIdentifier(){ return m_curAvailableIdentifier++;}

	   bool getDNAFlag(){return m_bDNAFlag;}
       void setDNAFlag(bool value);
	   void setDNAParams();
       void setProtParams();

	   short   *getSubMatrix(SubMatrixType vType);
	   short    getSubMatrixSize(SubMatrixType vType);

	   SubMatrixManager	   m_SubstitutionMatMgr;

	private:
		std::string m_AminoAcidCodes;
		std::hash_map<char, short>		m_AminoAcidChar2IntCodeMap;
		int		m_maxAllowedSeqLength;
		static unsigned long  m_curAvailableIdentifier;
		bool  m_bDNAFlag;
};

}