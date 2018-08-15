#pragma once

namespace SeqAnsis
{

class CPage	
{	
public:	
	float	getTop(){ return m_fTop; }	
	float	getBottom(){ return m_fBottom; }	
	void    setPage( float vTop, float vBottom );	
private:	
	float               m_fTop;	
	float               m_fBottom;
};	

}