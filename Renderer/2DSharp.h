#pragma once

#ifndef _WIN32_WINNT
#define _WIN32_WINNT  0x0601
#endif

#include<afxwin.h>
#include<vector>

namespace SeqAnsis
{

class I2DSharp
{
public:
	I2DSharp();
	virtual	~I2DSharp(){}
	virtual	void	draw( CDC* pDC ) = 0;
	virtual void    reLocate(int vOffsetX, int vOffsetY);
	virtual void	setColor( int r, int g, int b ){ m_r=r; m_g=g; m_b=b; }
	virtual int     getRed(){return m_r;}
	virtual int     getGreen(){return m_g;}
	virtual int     getBlue(){return m_b;}
public:
	float	m_fUnifiedWidth;
	float	m_fUnifiedHeight;
	float	m_fUnifiedXcoord;
	float	m_fUnifiedYcoord;

	float	m_fUnifiedYcoordDisp;

	int		m_iXcoord;
	int		m_iYcoord;
	int		m_iWidth;
	int		m_iHeight;
	bool    m_bDraw;
private:	
	int     m_r;
	int     m_g;
	int     m_b;
};

}