#include "CursorControler.h"

namespace  SeqAnsis
{

void	CCursor::pushPos()
{
	StruCursorPos  PrevPos;
	PrevPos.m_fX = m_fX;
	PrevPos.m_fY = m_fY;
	m_PrevPosStack.push( PrevPos );
}

void	CCursor::popPos()
{
	if ( !m_PrevPosStack.empty() )
	{
		StruCursorPos  PrevPos = m_PrevPosStack.top();
		m_fX = PrevPos.m_fX;
		m_fY = PrevPos.m_fY;
		m_PrevPosStack.pop();
	}

}

}