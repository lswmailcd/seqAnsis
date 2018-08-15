#pragma once
#include"2DSharpManager.h"
#include"Sequence.h"
#include"CursorControler.h"

namespace SeqAnsis
{

const float START_POS_X = 0.01f;
const float START_POS_Y = 0.0f;
const float END_POS_X = 0.95f;
const float END_POS_Y = 1.0f;
const float TITLE_WIDTH = 0.1f;
const int   TITLE_CONTENT_INTERVAL = 2;//two char space interval
const float PURE_TEXT_WIDTH=0.005f;
const float PURE_TEXT_HEIGHT=0.02f;
const float RECT_TEXT_WIDTH=0.01f;
const float RECT_TEXT_HEIGHT=PURE_TEXT_HEIGHT;
const int   RECT_TEXT_PER_LINE = 80;
const int   PRINT_LINE_INTERVAL = 6;

class CSeqRenderer
{
public:
	CSeqRenderer();
	~CSeqRenderer(){}
	void	clear();
	void	draw(CDC *vpDC);
	void    setRenderArea(int vX, int vY, int vW, int vH);

	void    printInputSeqs(const std::vector<CSequence>& vSeqs, const std::string& fileName);
	void		printNWAlignedSeqs( const std::vector<CSequence>& vSeqs );
	void		printSWAlignedSeqs( const std::vector<CSequence>& vSeqs );
	void		printBlastAlignedSeqs( const std::vector<CSequence>& vSeqs );
	void		printSAGEPAlignedSeqs( const std::vector<CSequence>& vSeqs );
	void		printMSAGAlignedSeqs( const std::vector<CSequence>& vSeqs, const std::string& strDrawText  );
	void		pageUp();
	void    pageDown();
	void		OnScrollUp();
	void		OnScrollDown();
private:
	void    countSeqPair( const std::vector<CSequence>& vSeqs, int& voNEq, int& voNGap );
	void    printSeqPair( const std::vector<CSequence>& vSeqs  );
	void    printSeqPairIndex( const std::vector<CSequence>& vSeqs, int vCurPos, bool vDoubleNum=true, bool vUpper=true);
	void    printSeqPairContent( const std::vector<CSequence>& vSeqs, int vCurPos, int vTxtRed, int vTxtGreen, int vTxtBlue, int vBkRed, int vBkGreen, int vBkBlue);
	
	void		printMSASeqs( const std::vector<CSequence>& vSeqs );

	void    drawTextWithRect( const char& vChar, int vTxtRed, int vTxtGreen, int vTxtBlue, int vRectRed, int vRectGreen, int vRectBlue );
	void  	drawText( const std::string& vString, int vRed, int vGreen, int vBlue, bool vVertical=false );
	CSharpManager		m_SharpMgr;
	CCursor				m_Cursor;
};

}