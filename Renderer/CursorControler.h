#pragma once

#include <stack>

namespace SeqAnsis
{

struct StruCursorPos 
{
	float	m_fX;
	float	m_fY;
};

class CCursor
{
public:
	float	getXPos(){ return m_fX; }
	float	getYPos(){ return m_fY; }
	void    setXPos( float vX ){m_fX=vX;}
	void    setYPos( float vY ){m_fY=vY;}
	void	setPos(float vX, float vY){ m_fX=vX; m_fY=vY; }
	float   getStepX(){return m_fStepX;}
	float   getStepY(){return m_fStepY;}
	void	setStep(float vX, float vY){ m_fStepX=vX; m_fStepY=vY; }
	void    moveX( int vX ){ m_fX+=(m_fStepX*vX); }
	void    moveY( int vY ){ m_fY+=(m_fStepY*vY); }
	void    move(int vX, int vY){ m_fX+=(m_fStepX*vX); m_fY+=(m_fStepY*vY); }
	void    pushPos();
	void    popPos();
private:
	float               m_fX;
	float               m_fY;
	float               m_fStepX;
	float               m_fStepY;
	std::stack<StruCursorPos>	m_PrevPosStack;
	float               m_fPrevX;
	float               m_fPrevY;
};

}