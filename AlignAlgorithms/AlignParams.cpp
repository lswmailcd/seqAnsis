#include"Common.h"
#include"AlignParams.h"

namespace SeqAnsis
{

unsigned long	CAlignParams::m_curAvailableIdentifier=1;

CAlignParams::CAlignParams():m_AminoAcidCodes("ABCDEFGHIKLMNPQRSTVWXYZ-"), m_maxAllowedSeqLength(INT_MAX)
{
	for ( int i=0; i<(int)m_AminoAcidCodes.size(); ++i )
	{
		m_AminoAcidChar2IntCodeMap.insert( std::make_pair( m_AminoAcidCodes.at(i), (short)i ) );
	}
	
}

short	CAlignParams::getAminoAcidChar2IntCode( char a )
{
	return	m_AminoAcidChar2IntCodeMap[a];
}

char		CAlignParams::getAminoAcidInt2CharCode( short  iCode )
{
	if ( iCode>=0 && iCode<(short)m_AminoAcidCodes.size() )
	{
		return	m_AminoAcidCodes[iCode];
	}
	return NONEGENE;
}

short*  CAlignParams::getSubMatrix(SubMatrixType vType)
{
	return		m_SubstitutionMatMgr.getSubMatrix( vType );
}

short    CAlignParams::getSubMatrixSize(SubMatrixType vType)
{
	return		m_SubstitutionMatMgr.getSubMatrixSize(vType);
}

void CAlignParams::setDNAFlag(bool value)
{
#if 0
	if(value == true)
	{
		setDNAParams();
	}
	else
	{
		setProtParams();
	}    
	DNAFlag = value;
#endif
}

void CAlignParams::setDNAParams()
{
#if 0
	gapOpen       = DNAGapOpen;
	gapExtend     = DNAGapExtend;
	PWGapOpen  = DNAPWGapOpen;
	PWGapExtend  = DNAPWGapExtend;
	ktup           = DNAKtup;
	window         = DNAWindow;
	signif         = DNASignif;
	windowGap       = DNAWindowGap;
#endif
}

void CAlignParams::setProtParams()
{
#if 0
	gapOpen       = AAGapOpen;
	gapExtend     = AAGapExtend;
	PWGapOpen  = AAPWGapOpen;
	PWGapExtend  = AAPWGapExtend;
	ktup           = AAKtup;
	window         = AAWindow;
	signif         = AASignif;
	windowGap       = AAWindowGap;
#endif
}

}