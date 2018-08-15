#include"2DSharp.h"

namespace SeqAnsis
{

I2DSharp::I2DSharp():m_r(0), m_g(0), m_b(0), m_bDraw(true)
{
}

void	I2DSharp::reLocate( int vOffsetX, int vOffsetY )
{
	m_iXcoord += vOffsetX;
	m_iYcoord += vOffsetY;
}

}