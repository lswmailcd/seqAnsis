#pragma once
#include"2DSharp.h"
#include"2DText.h"
#include"2DRect.h"
#include "PageControler.h"

namespace SeqAnsis
{

const float UNIFIEDCOORD_OFFSET = -0.1f;
const float UNIFIEDCOORD_BOTTOM = 1.0f;
const float UNFIIED_PAGEUPDOWN_OFFSET = 0.9f;
const float UNFIIED_SCROLL_OFFSET = 0.1f;

class CSharpManager
{
public:
	CSharpManager();
	~CSharpManager();
	void	draw(CDC* vpDC);
	void	addSharp( I2DSharp* vpSharp );
	void    setRange( int vX, int vY );
	void    release();
	void    setRenderArea(int vX, int vY, int vW, int vH){ m_iSTART_COORD_X=vX; m_iSTART_COORD_Y=vY; m_iMAX_RANGE_X=vW; m_iMAX_RANGE_Y=vH;}
	void	pageUp();
	void    pageDown();
	void	scrollUp();
	void    scrollDown();
	void    beginDraw();
	void    endDraw();
protected:
	std::vector<I2DSharp*>	m_2DSharpList;
private:
	void    coordMappingDisp();
	void    coordMappingList();
	void    coordMapping();
	float	findMaxUnifiedBottomY();
	void	mapUnifiedCoord2Display();
	int		m_iMAX_RANGE_X;
	int		m_iMAX_RANGE_Y; 
	int		m_iSTART_COORD_X;
	int		m_iSTART_COORD_Y;
	int     m_iOverflowAreaTop;
	int     m_iOverflowAreaBottom;
	float   m_fCurUnifiedStartPosY;
	bool    m_bDividScreen;
	float   m_fPageUpDownOffset;
	CPage	m_Page2Disp;
};

}