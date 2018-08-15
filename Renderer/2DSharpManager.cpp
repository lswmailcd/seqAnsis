#pragma once

#include"Common.h"
#include"2DSharpManager.h"

namespace SeqAnsis
{

CSharpManager::CSharpManager():m_iMAX_RANGE_X(1920), m_iMAX_RANGE_Y(1080), m_iSTART_COORD_X(0), m_iSTART_COORD_Y(0)
, m_fCurUnifiedStartPosY(0.0f), m_bDividScreen(false),m_fPageUpDownOffset(UNFIIED_PAGEUPDOWN_OFFSET)
{
	m_Page2Disp.setPage(0.0f, 1.0f);
}

void    CSharpManager::release()
{
	for( std::vector<I2DSharp*>::iterator itr=m_2DSharpList.begin();
		 itr!=m_2DSharpList.end();
		 ++itr)
	{
		SAFE_DELETE((*itr));
	}
}

CSharpManager::~CSharpManager()
{
	release();
}

void	CSharpManager::draw(CDC* vpDC)
{
	coordMapping();

	for( std::vector<I2DSharp*>::iterator itr=m_2DSharpList.begin();
		 itr!=m_2DSharpList.end();
		 ++itr)
	{
		(*itr)->draw(vpDC);
	}
	//if ( m_bDividScreen )
	//{
	//	::MessageBox(NULL, "Press button to continue!", "ÌבÊ¾",MB_OK );
	//	m_bDividScreen = false;
	//	CBrush backBrush(RGB(255,255,255));
	//	RECT rc;
	//	rc.left=m_iSTART_COORD_X;
	//	rc.top=m_iSTART_COORD_Y;
	//	rc.right=m_iSTART_COORD_X+m_iMAX_RANGE_X; 
	//	rc.bottom=m_iSTART_COORD_Y+m_iMAX_RANGE_Y;
	//	vpDC->FillRect(&rc,&backBrush); 
	//	draw( vpDC );
	//}
}

void	CSharpManager::coordMapping()
{
	for( std::vector<I2DSharp*>::iterator itr=m_2DSharpList.begin();
		itr!=m_2DSharpList.end();
		++itr)
	{
		(*itr)->m_bDraw = true;
		(*itr)->m_fUnifiedYcoordDisp = (*itr)->m_fUnifiedYcoord - m_Page2Disp.getTop();
		//we suppose the sharp are all in one page in one screen draw.
		if( ((*itr)->m_fUnifiedYcoordDisp+(*itr)->m_fUnifiedHeight) < 0.0f || (*itr)->m_fUnifiedYcoordDisp > UNIFIEDCOORD_BOTTOM )
		{
			(*itr)->m_bDraw = false;
		}
	}
	mapUnifiedCoord2Display();
}

void    CSharpManager::coordMappingDisp()
{
	//find the max unified y coord.
	float fUnifiedMaxBottom=findMaxUnifiedBottomY();

	m_bDividScreen = (fUnifiedMaxBottom-m_fCurUnifiedStartPosY)>1.0f;

	for( std::vector<I2DSharp*>::iterator itr=m_2DSharpList.begin();
		itr!=m_2DSharpList.end();
		++itr)
	{
		(*itr)->m_bDraw = true;
		(*itr)->m_fUnifiedYcoordDisp = (*itr)->m_fUnifiedYcoord-m_fCurUnifiedStartPosY;
		//we suppose the sharp are all in one page in one screen draw.
		if( (*itr)->m_fUnifiedYcoordDisp < 0.0f || (*itr)->m_fUnifiedYcoordDisp > UNIFIEDCOORD_BOTTOM )
		{
			(*itr)->m_bDraw = false;
		}
	}

	if( !m_2DSharpList.empty() && m_bDividScreen )	
	{
		m_fCurUnifiedStartPosY += UNIFIEDCOORD_BOTTOM;
		m_fCurUnifiedStartPosY += UNIFIEDCOORD_OFFSET;
	}

	mapUnifiedCoord2Display();
}

float	CSharpManager::findMaxUnifiedBottomY()
{
	//find the max unified y coord.
	float fUnifiedMaxBottom=0.0f;
	for( std::vector<I2DSharp*>::iterator itr=m_2DSharpList.begin();
		itr!=m_2DSharpList.end();
		++itr)
	{
		(*itr)->m_fUnifiedYcoordDisp = (*itr)->m_fUnifiedYcoord;//init the disp coord as the original one
		float fBottom = (*itr)->m_fUnifiedYcoord + (*itr)->m_fUnifiedHeight;
		fUnifiedMaxBottom = (fBottom - fUnifiedMaxBottom)>0.0f ? fBottom : fUnifiedMaxBottom;
	}

	return fUnifiedMaxBottom;
}

void    CSharpManager::coordMappingList()
{
	//find the max unified y coord.
	float fUnifiedMaxBottom=findMaxUnifiedBottomY();

	if( fUnifiedMaxBottom>m_Page2Disp.getBottom() )
	{
		m_Page2Disp.setPage( fUnifiedMaxBottom-1.0f, fUnifiedMaxBottom );
	}
}

void	CSharpManager::mapUnifiedCoord2Display()
{
	//map the logic coordination to screen coordination.
	for( std::vector<I2DSharp*>::iterator itr=m_2DSharpList.begin();
		itr!=m_2DSharpList.end();
		++itr)
	{
		//we suppose the sharp are all in one page in one screen draw.
		if( (*itr)->m_fUnifiedYcoordDisp > -1.0f )
		{
			(*itr)->m_iXcoord = int(m_iMAX_RANGE_X*(*itr)->m_fUnifiedXcoord);
			(*itr)->m_iYcoord = int(m_iMAX_RANGE_Y*(*itr)->m_fUnifiedYcoordDisp);
			(*itr)->m_iWidth = int(m_iMAX_RANGE_X*(*itr)->m_fUnifiedWidth);
			(*itr)->m_iHeight = int(m_iMAX_RANGE_Y*(*itr)->m_fUnifiedHeight);
		}
	}
}

void	CSharpManager::addSharp( I2DSharp* vpSharp )
{
	if(NULL==vpSharp)	return;

	m_2DSharpList.push_back(vpSharp);
}

void	CSharpManager::pageUp()
{
	m_Page2Disp.setPage( m_Page2Disp.getTop()-UNFIIED_PAGEUPDOWN_OFFSET, m_Page2Disp.getBottom()-UNFIIED_PAGEUPDOWN_OFFSET );
}

void	CSharpManager::pageDown()
{
	m_Page2Disp.setPage( m_Page2Disp.getTop()+UNFIIED_PAGEUPDOWN_OFFSET, m_Page2Disp.getBottom()+UNFIIED_PAGEUPDOWN_OFFSET );	
}

void	CSharpManager::scrollUp()
{
	m_Page2Disp.setPage( m_Page2Disp.getTop()-UNFIIED_SCROLL_OFFSET, m_Page2Disp.getBottom()-UNFIIED_SCROLL_OFFSET );
}

void	CSharpManager::scrollDown()
{
	m_Page2Disp.setPage( m_Page2Disp.getTop()+UNFIIED_SCROLL_OFFSET, m_Page2Disp.getBottom()+UNFIIED_SCROLL_OFFSET );
}

void	CSharpManager::beginDraw()
{

}

void	CSharpManager::endDraw()
{
	coordMappingList();
}

}