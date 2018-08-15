#pragma once

#include<string>
#include"2DSharp.h"

namespace SeqAnsis
{

class C2DText:public I2DSharp	
{
public:
	C2DText();
	virtual ~C2DText();
	virtual	void	draw( CDC* vpDC );
	
	char	m_char;
};

}