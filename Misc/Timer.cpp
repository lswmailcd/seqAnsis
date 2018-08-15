#include "Windows.h"
#include "Timer.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

namespace	SeqAnsis
{

CTimer::CTimer() : m_Active(false), m_BeginTime(0)
{
    m_Initialized = QueryPerformanceFrequency((LARGE_INTEGER *)&m_Frequency);
}

bool CTimer::begin()
{
	if (!m_Initialized )
    	return false; // error - couldn't get frequency

	m_Active = true;
	// get the starting counter value
	if (QueryPerformanceCounter((LARGE_INTEGER *)&m_BeginTime) == 0)
		return false;
	else
		return true;
}

//FUNCTION: compute the time elapsed and stop timer
double CTimer::end()
{
	if (!m_Initialized )
	    return 0.0; // error - couldn't get frequency

	__int64 iEndTime,iElapsed;

	// get the ending counter value
	QueryPerformanceCounter((LARGE_INTEGER *)&iEndTime);

	// determine the elapsed counts
	iElapsed = iEndTime - m_BeginTime;
	m_Active = false;

	// convert counts to time in seconds and return it
	return (double)iElapsed / (double)m_Frequency;
}

//FUNCTION: compute the time elapsed without stopping timer
double CTimer::getCurrentTime()
{
	if (!m_Initialized )
		return 0.0; // error - couldn't get frequency

	__int64 iEndTime,iElapsed;

	// get the ending counter value
	QueryPerformanceCounter((LARGE_INTEGER *)&iEndTime);

	// determine the elapsed counts
	iElapsed = iEndTime - m_BeginTime;

	// convert counts to time in seconds and return it
	return (double)iElapsed / (double)m_Frequency;
}

CTimer::~CTimer()
{

}

}