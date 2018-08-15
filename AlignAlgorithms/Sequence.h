//This is the definition of DNA or PROTEIN sequence used in SeqAnsis
#pragma once

#include<string>
#include "Common.h"
#include "GlobalSpace.h"

namespace SeqAnsis
{

const char GENESPACE='-';
const int   SEQUENCETOOBIG = -900;

struct StruSeqElem
{
	StruSeqElem(){}
	StruSeqElem( char vCh, int vIndex )
	{ 
		m_char=vCh; 
		m_index=vIndex; 
		m_iCode=CGlobalSpace::m_sAlignParams.getAminoAcidChar2IntCode( m_char ); 		
	}
	char	 m_char;//�ַ���ʾ�Ĳл�
	int	 m_iCode;//���ֱ�ʾ�Ĳл�
	int     m_index;//�л��������е���ţ���0Ϊ��ʼ���
};

typedef	std::vector<StruSeqElem>  CSeqData;	

class CSequence
{
public:
	CSequence(){}
	 CSequence(const CSeqData& vSeq, const std::string& vName, const std::string& vTitle);
	 bool   isEmpty();
	 const CSeqData& getSequenceContext() const  {return m_Sequence;} 
	 CSeqData& getSequenceContext()  {return m_Sequence;} 
	 const std::string&	getName() const; 
	 const std::string&	getTitle() const;
	 int	getLen() const { return m_Sequence.size(); }
private:
	CSeqData				m_Sequence;
	std::string				m_SeqName;
	std::string				m_SeqTitle;
	unsigned long	    m_seqIdentifier;
};

}
