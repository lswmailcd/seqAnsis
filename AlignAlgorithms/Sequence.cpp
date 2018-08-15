#include"Sequence.h"
#include"GlobalSpace.h"

namespace SeqAnsis
{

CSequence::CSequence(const CSeqData& vSeq, const std::string& vName, const std::string& vTitle):m_Sequence(vSeq), m_SeqName(vName), m_SeqTitle(vTitle)
{
	m_seqIdentifier=CGlobalSpace::m_sAlignParams.getUniqueSequenceIdentifier();
}

bool CSequence::isEmpty()
{
	return (m_Sequence.size() == 0);   
}

const std::string& CSequence::getName() const
{
    return m_SeqName;
}

const std::string& CSequence::getTitle() const
{
    return m_SeqTitle;
}

}