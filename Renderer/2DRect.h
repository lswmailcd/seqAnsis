#pragma once

#include"2DSharp.h"

namespace SeqAnsis
{

class C2DRect:public I2DSharp	
{
public:
	C2DRect();
	virtual ~C2DRect();
	virtual	void	draw(CDC* pDC);

	std::string		m_Text;
};

}