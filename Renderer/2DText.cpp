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
	VERIFY(font.CreatePointFont(120, _T("宋体"), vpDC));//创建字体宋体格式  100为字高 
	CFont*def_font = vpDC->SelectObject(&font);  //选择该字体进入PDC 
	vpDC->SetBkMode(TRANSPARENT);//设置字体背景为透明

	CRect rect( m_iXcoord, m_iYcoord, m_iXcoord+m_iWidth, m_iYcoord+m_iHeight );
	vpDC->SetTextColor(RGB(getRed(),getGreen(),getBlue()));
	vpDC->DrawText( m_char, &rect, DT_SINGLELINE|DT_CENTER|DT_VCENTER );

	vpDC->SelectObject(def_font);//恢复PDC的缺省字体  
	// Done with the font.Delete the font object. 
	font.DeleteObject();//释放font对象  
}

}