#include"2DRect.h"

namespace SeqAnsis
{

C2DRect::C2DRect()
{

}

C2DRect::~C2DRect() 
{

}

void	C2DRect::draw( CDC* vpDC )
{
	if( !m_bDraw )  return;

	CBrush backBrush1(RGB(getRed(),getGreen(),getBlue())); 
	CBrush *pOldBrush=vpDC->SelectObject(&backBrush1);
	vpDC->Rectangle( m_iXcoord, m_iYcoord, m_iXcoord+m_iWidth, m_iYcoord+m_iHeight );
	vpDC->SelectObject(pOldBrush);
}

}