// Timer.h: interface for the CTimer class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TIMER_H__0859A443_B838_11D3_9AC6_0000B4B5B268__INCLUDED_)
#define AFX_TIMER_H__0859A443_B838_11D3_9AC6_0000B4B5B268__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

namespace	SeqAnsis
{

class CTimer  
{
public:
	CTimer();
	virtual ~CTimer();
   
	bool   begin();
	double end();
	double getCurrentTime();

	inline bool isActive() {return m_Active;}

private :
    int m_Initialized;
	bool m_Active;
    __int64 m_Frequency;
    __int64 m_BeginTime;
};

}

#endif // !defined(AFX_TIMER_H__0859A443_B838_11D3_9AC6_0000B4B5B268__INCLUDED_)
