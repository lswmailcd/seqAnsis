#include"2DText.h"

namespace SeqAnsis
{

C2DText::C2DText()
{

}

C2DText::~C2DText()
{

}

void	C2DText::draw( CDC* vpDC )
{
	if( !m_bDraw )  return;

	CFont font;
	VERIFY(font.CreatePointFont(120, _T("����"), vpDC));//�������������ʽ  100Ϊ�ָ� 
	CFont*def_font = vpDC->SelectObject(&font);  //ѡ����������PDC 
	vpDC->SetBkMode(TRANSPARENT);//�������屳��Ϊ͸��

	CRect rect( m_iXcoord, m_iYcoord, m_iXcoord+m_iWidth, m_iYcoord+m_iHeight );
	vpDC->SetTextColor(RGB(getRed(),getGreen(),getBlue()));
	vpDC->DrawText( m_char, &rect, DT_SINGLELINE|DT_CENTER|DT_VCENTER );

	vpDC->SelectObject(def_font);//�ָ�PDC��ȱʡ����  
	// Done with the font.Delete the font object. 
	font.DeleteObject();//�ͷ�font����  
}

}