#pragma once

#include"Utility.h"
#include"EventLog.h"
#include"AlignParams.h"

namespace SeqAnsis
{

class CGlobalSpace
{
public:
	enum ParserType
	{
		PARSER_TYPE_TEXT = 2000,
		PARSER_TYPE_XML  = 2001,
	};

public:
	static CUtility    m_sUtility;
	static CEventLog   m_sEventLog;
	static CAlignParams		m_sAlignParams;

protected:
	~CGlobalSpace(void);

private:
	CGlobalSpace(void);
};

}
