#include"Common.h"
#include"SeqRenderer.h"
#include"GlobalSpace.h"

namespace SeqAnsis
{

CSeqRenderer::CSeqRenderer()
{
	m_Cursor.setPos(START_POS_X,START_POS_Y);
	m_Cursor.setStep(PURE_TEXT_WIDTH,PURE_TEXT_HEIGHT );
}

void    CSeqRenderer::setRenderArea(int vX, int vY, int vW, int vH)
{
	m_SharpMgr.setRenderArea( vX, vY, vW, vH );
}

void    CSeqRenderer::drawTextWithRect( const char& vChar, int vTxtRed, int vTxtGreen, int vTxtBlue, int vRectRed, int vRectGreen, int vRectBlue )
{
	SeqAnsis::C2DRect *pRect = new SeqAnsis::C2DRect;			
	pRect->setColor(vRectRed,vRectGreen,vRectBlue);					
	pRect->m_fUnifiedXcoord=m_Cursor.getXPos();					
	pRect->m_fUnifiedYcoord=m_Cursor.getYPos();	
	pRect->m_fUnifiedWidth = m_Cursor.getStepX();
	pRect->m_fUnifiedHeight = m_Cursor.getStepY();
	m_SharpMgr.addSharp( pRect );

	std::string str;
	str.push_back(vChar);
	drawText( str, vTxtRed, vTxtGreen, vTxtBlue );
}

void    CSeqRenderer::drawText( const std::string& vString, int vRed, int vGreen, int vBlue, bool vVertical /*=false*/  )
{
	if ( vVertical )
	{
		for( int i=0; i<(int)(vString.size()); ++i )
		{
			SeqAnsis::C2DText *pText = new SeqAnsis::C2DText;			
			pText->m_char = vString.at(i);							
			pText->m_fUnifiedXcoord=m_Cursor.getXPos();					
			pText->m_fUnifiedYcoord=m_Cursor.getYPos();	
			pText->m_fUnifiedWidth = m_Cursor.getStepX();
			pText->m_fUnifiedHeight = m_Cursor.getStepY();
			pText->setColor( vRed, vGreen, vBlue );
			m_SharpMgr.addSharp( pText );
			m_Cursor.moveY(1);
		}
	}
	else
	{
		for( int i=0; i<(int)(vString.size()); ++i )
		{
			SeqAnsis::C2DText *pText = new SeqAnsis::C2DText;			
			pText->m_char = vString.at(i);							
			pText->m_fUnifiedXcoord=m_Cursor.getXPos();					
			pText->m_fUnifiedYcoord=m_Cursor.getYPos();	
			pText->m_fUnifiedWidth = m_Cursor.getStepX();
			pText->m_fUnifiedHeight = m_Cursor.getStepY();
			pText->setColor( vRed, vGreen, vBlue );
			m_SharpMgr.addSharp( pText );
			m_Cursor.moveX(1);
		}
	}
}

void	CSeqRenderer::printSeqPairContent( const std::vector<CSequence>& vSeqs, int vCurPos, int vTxtRed, int vTxtGreen, int vTxtBlue, int vBkRed, int vBkGreen, int vBkBlue)
{
	const CSeqData& Seq1=vSeqs[0].getSequenceContext();
	const CSeqData& Seq2=vSeqs[1].getSequenceContext();

	drawTextWithRect( Seq1[vCurPos].m_char, vTxtRed, vTxtGreen, vTxtBlue, vBkRed, vBkGreen, vBkBlue );

	m_Cursor.move(-1, 1);
	drawTextWithRect( Seq2[vCurPos].m_char, vTxtRed, vTxtGreen, vTxtBlue, vBkRed, vBkGreen, vBkBlue );
}

void    CSeqRenderer::printSAGEPAlignedSeqs( const std::vector<CSequence>& vSeqs )
{
	if( vSeqs.size()==2 )
	{
		const CSeqData& Seq1=vSeqs[0].getSequenceContext();
		const CSeqData& Seq2=vSeqs[1].getSequenceContext();
		if( Seq1.size()==Seq2.size() )
		{
			//if the current line is not first line, step another new line to display
			m_Cursor.setStep( PURE_TEXT_WIDTH, PURE_TEXT_HEIGHT );
			if( abs(m_Cursor.getYPos())>EPS )
			{
				m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+PRINT_LINE_INTERVAL*m_Cursor.getStepY() );
			}
			else
			{
				m_Cursor.setPos( START_POS_X, START_POS_Y );
			}

			m_SharpMgr.beginDraw();
			drawText( "SAGEP alignment result:", 0,0,0 );

			int nEq, nGap;
			countSeqPair( vSeqs, nEq, nGap );
			char ch[100];
			sprintf_s( ch, "identity  %d,  gap  %d,  seq length  %d", nEq, nGap, Seq1.size() );
			drawText( ch, 0,0,0 );

			m_Cursor.moveY( 2 );

			printSeqPair( vSeqs );
			m_SharpMgr.endDraw();

		}
		else
		{
			CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "two sequences size is not equal!" );
		}
	}
	else
	{
		CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "Wrong sequences number" );
	}
}

void		CSeqRenderer::printMSASeqs( const std::vector<CSequence>& vSeqs )
{
	if ( (int)(vSeqs.size())<2 )
	{
		CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "sequences number is <2!" );
		return;
	}

	int titleCharNum = (int)(TITLE_WIDTH/m_Cursor.getStepX());
	for ( int i=0;  i<(int)vSeqs.size(); ++i )
	{
		if( titleCharNum<(int)vSeqs[i].getName().size() )
		{
			CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "sequences name is too long to print!" );
			break;
		}
	}

	//find longest sequence length
	int iLongestSeqLen = 0;
	for ( int i=0;  i<(int)vSeqs.size(); ++i )
	{
		if( iLongestSeqLen<vSeqs[i].getLen() )
		{
			iLongestSeqLen = vSeqs[i].getLen();
		}
	}

	m_Cursor.setStep( RECT_TEXT_WIDTH, RECT_TEXT_HEIGHT );
	int iCharNum = (int)(((END_POS_X-START_POS_X)-TITLE_WIDTH-TITLE_CONTENT_INTERVAL*m_Cursor.getStepX())/m_Cursor.getStepX());
	int iLineNum = (int)ceil((float)iLongestSeqLen/(float)iCharNum);

	bool  bCharEqual = false;

	int iCurPos = 0;
	bool  bFirstRow=true;
	bool  bFirstNoPrinted = true;
	while( iCurPos<iLongestSeqLen )
	{
		//print the seq name
		for(int i=0; i<(int)vSeqs.size(); ++i)
		{
			//new line, Seq name
			m_Cursor.setStep( PURE_TEXT_WIDTH, PURE_TEXT_HEIGHT );	
			if( bFirstRow )	
			{	
				if(bFirstNoPrinted)
				{
					m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+3*m_Cursor.getStepY() );
					bFirstNoPrinted = false;
				}
				else
				{
					m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+(vSeqs.size()+1)*m_Cursor.getStepY() );
				}
				bFirstRow=false;	
			}	
			else	
			{	
				m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+m_Cursor.getStepY() );	
			}	
			drawText( vSeqs[i].getName(), 0,0,0 );	
		}

		//step back n-1 lines (n=seqSize),  prepare for content printing of one new line	
		m_Cursor.setStep( RECT_TEXT_WIDTH, RECT_TEXT_HEIGHT );	
		m_Cursor.setPos( m_Cursor.getXPos()+TITLE_CONTENT_INTERVAL*m_Cursor.getStepX()	
			, m_Cursor.getYPos()-((int)vSeqs.size()-1)*m_Cursor.getStepY() );	
		bFirstRow = true;

		while ( !((m_Cursor.getXPos()+m_Cursor.getStepX())>END_POS_X) && iCurPos<iLongestSeqLen )
		{
			//find the equal character at current position
			bool bEqual = true;
			if( iCurPos<vSeqs[0].getLen() )
			{
				char ch = vSeqs[0].getSequenceContext()[iCurPos].m_char;
				for ( int i=1; i<(int)vSeqs.size(); ++i )
				{
					if( iCurPos<(int)vSeqs[i].getSequenceContext().size() &&  ch !=vSeqs[i].getSequenceContext()[iCurPos].m_char )
					{
						bEqual = false;
						break;
					}
				}
			}
			else
			{
				bEqual = false;
			}

			//draw current character 
			for ( int i=0; i<(int)vSeqs.size(); ++i )
			{
				if( iCurPos<vSeqs[i].getLen() )
				{
					if( bEqual )
					{
						drawTextWithRect( vSeqs[i].getSequenceContext()[iCurPos].m_char, 0, 0, 0, 0, 255, 0 );
					}
					else
					{
						drawTextWithRect( vSeqs[i].getSequenceContext()[iCurPos].m_char, 0, 0, 0, 255, 0, 0 );
					}
					//the drawTextWithRect will automatic step one at X direction, so we rollback one step at X.
					//step one at Y
					m_Cursor.move(-1,1);
				}
				else
				{//if current sequence is end , step one at Y direction and draw next sequence.
					m_Cursor.moveY(1);
				}
			}
			m_Cursor.move(1, -(int)vSeqs.size() );
			++iCurPos;
		}
	}
}

void    CSeqRenderer::printSeqPair( const std::vector<CSequence>& vSeqs  )
{
	if ( (int)(vSeqs.size())!=2 )
	{
		CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "sequences number is NOT two!" );
		return;
	}
	const CSeqData& Seq1=vSeqs[0].getSequenceContext();
	const CSeqData& Seq2=vSeqs[1].getSequenceContext();

	int titleCharNum = (int)(TITLE_WIDTH/m_Cursor.getStepX());

	if( titleCharNum<(int)(vSeqs[0].getName().size()) || titleCharNum<(int)(vSeqs[1].getName().size()) )
	{
		CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "sequences name is to long to print!" );
	}

	m_Cursor.setStep( RECT_TEXT_WIDTH, RECT_TEXT_HEIGHT );
	int iCharNum = (int)(((END_POS_X-START_POS_X)-TITLE_WIDTH-TITLE_CONTENT_INTERVAL*m_Cursor.getStepX())/m_Cursor.getStepX());
	int iLineNum = (int)ceil((float)Seq1.size()/(float)iCharNum);

	bool  bCharEqual = false;
	
	int iCurPos = 0;
	bool  bFirstRow=true;
	while( iCurPos<(int)(Seq1.size()) )
	{
		bool  bFirstNoPrinted = false;
		//new line, Seq1 name
		m_Cursor.setStep( PURE_TEXT_WIDTH, PURE_TEXT_HEIGHT );	
		if( bFirstRow )	
		{	
			m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+m_Cursor.getStepY() );	
			bFirstRow=false;	
		}	
		else	
		{	
			m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+10*m_Cursor.getStepY() );	
		}	
		drawText( vSeqs[0].getName(), 0,0,0 );	

		//new line, Seq2 name	
		m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+m_Cursor.getStepY() );	
		drawText( vSeqs[1].getName(), 0,0,0 );	

		//step back one line, prepare for content printing of one new line	
		m_Cursor.setStep( RECT_TEXT_WIDTH, RECT_TEXT_HEIGHT );	
		m_Cursor.setPos( m_Cursor.getXPos()+TITLE_CONTENT_INTERVAL*m_Cursor.getStepX()	
			, m_Cursor.getYPos()-m_Cursor.getStepY() );	

		while ( !((m_Cursor.getXPos()+m_Cursor.getStepX())>END_POS_X) && iCurPos<(int)(Seq1.size()) )
		{
			if( Seq1[iCurPos].m_char==Seq2[iCurPos].m_char )
			{
				if( !bFirstNoPrinted || !bCharEqual )
				{
					m_Cursor.pushPos();
					printSeqPairIndex( vSeqs, iCurPos, bFirstNoPrinted );
					m_Cursor.popPos();

					printSeqPairContent(vSeqs, iCurPos, 0,0,0,0, 255, 0 );

					m_Cursor.pushPos();
					m_Cursor.moveX(-1);
					printSeqPairIndex( vSeqs, iCurPos, bFirstNoPrinted, false );
					m_Cursor.popPos();

					bCharEqual = true;
					bFirstNoPrinted = true;
				}
				else
				{
					printSeqPairContent(vSeqs, iCurPos, 0,0,0,0, 255, 0 );
				}
			}
			else
			{
				if( !bFirstNoPrinted || bCharEqual )
				{
					m_Cursor.pushPos();
					printSeqPairIndex( vSeqs, iCurPos, bFirstNoPrinted );
					m_Cursor.popPos();

					printSeqPairContent(vSeqs, iCurPos, 0,0,0,255,0,0 );

					m_Cursor.pushPos();
					m_Cursor.moveX(-1);
					printSeqPairIndex( vSeqs, iCurPos, bFirstNoPrinted, false );
					m_Cursor.popPos();

					bCharEqual = false;
					bFirstNoPrinted = true;
				}
				else
				{
					printSeqPairContent(vSeqs, iCurPos, 0,0,0,255,0,0 );
				}
			}
			m_Cursor.moveY( -1 );	
			++iCurPos;
		}
		//print the sequence pair line tail index number
		--iCurPos;
		m_Cursor.pushPos();	
		m_Cursor.moveX(-1);	
		printSeqPairIndex( vSeqs, iCurPos, false );	
		m_Cursor.popPos();	

		m_Cursor.pushPos();	
		m_Cursor.move(-1,1);	
		printSeqPairIndex( vSeqs, iCurPos, false, false );	
		m_Cursor.popPos();	
		++iCurPos;
	}
}

void    CSeqRenderer::printSeqPairIndex( const std::vector<CSequence>& vSeqs, int vCurPos, bool vDoubleNum/*=true*/, bool vUpper/*=true*/)
{
	const CSeqData& Seq1=vSeqs[0].getSequenceContext();
	const CSeqData& Seq2=vSeqs[1].getSequenceContext();
	if(vUpper)
	{//index number is above the sequence pair
		if(vDoubleNum && vCurPos>0)
		//draw the index of prior seq tail
		{
			char ch[10];
			sprintf_s( ch, "%d", Seq1[vCurPos-1].m_index );
			std::string str1( ch );
			m_Cursor.move( -1, -(int)(str1.length()) );
			drawText( ch, 0,0,0, true );
			m_Cursor.moveX(1);
		}

		char ch[10];
		sprintf_s( ch, "%d", Seq1[vCurPos].m_index );
		std::string str1( ch );
		m_Cursor.moveY( -(int)(str1.length()) );
		drawText( ch, 0,0,0, true );
	}
	else
	{//index number is under the sequence pair
		//draw the index of prior seq tail
		if(vDoubleNum && vCurPos>0)
		{
			char ch[10];
			sprintf_s( ch, "%d", Seq2[vCurPos-1].m_index );
			m_Cursor.pushPos();
			m_Cursor.move(-1,1);
			drawText( ch, 0,0,0, true );
			m_Cursor.popPos();
		}

		char ch[10];
		sprintf_s( ch, "%d", Seq2[vCurPos].m_index );
		std::string str1( ch );
		m_Cursor.moveY(1);
		drawText( ch, 0,0,0, true );
	}
}

void	CSeqRenderer::printNWAlignedSeqs( const std::vector<CSequence>& vSeqs )
{
	if( vSeqs.size()==2 )
	{
		const CSeqData& Seq1=vSeqs[0].getSequenceContext();
		const CSeqData& Seq2=vSeqs[1].getSequenceContext();
		if( Seq1.size()==Seq2.size() )
		{
			//if the current line is not first line, step another new line to display
			m_Cursor.setStep( PURE_TEXT_WIDTH, PURE_TEXT_HEIGHT );
			if( abs(m_Cursor.getYPos())>EPS )
			{
				m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+PRINT_LINE_INTERVAL*m_Cursor.getStepY() );
			}
			else
			{
				m_Cursor.setPos( START_POS_X, START_POS_Y );
			}

			int nEq, nGap;
			countSeqPair( vSeqs, nEq, nGap );

			m_SharpMgr.beginDraw();
			drawText( "Needleman-Wunsch alignment result:", 0,0,0 );
			char ch[100];
			sprintf_s( ch, "identity  %d,  gap  %d,  seq length  %d", nEq, nGap, Seq1.size() );
			drawText( ch, 0,0,0 );
			m_Cursor.moveY( 2 );
			printSeqPair( vSeqs );
			m_SharpMgr.endDraw();

		}
		else
		{
			CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "two sequences size is not equal!" );
		}
	}
	else
	{
		CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "Wrong sequences number" );
	}
}

void	CSeqRenderer::printSWAlignedSeqs( const std::vector<CSequence>& vSeqs )
{
	if( vSeqs.size()==2 )
	{
		const CSeqData& Seq1=vSeqs[0].getSequenceContext();
		const CSeqData& Seq2=vSeqs[1].getSequenceContext();
		if( Seq1.size()==Seq2.size() )
		{
			//if the current line is not first line, step another new line to display
			m_Cursor.setStep( PURE_TEXT_WIDTH, PURE_TEXT_HEIGHT );
			if( abs(m_Cursor.getYPos())>EPS )
			{
				m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+PRINT_LINE_INTERVAL*m_Cursor.getStepY() );
			}
			else
			{
				m_Cursor.setPos( START_POS_X, START_POS_Y );
			}

			m_SharpMgr.beginDraw();
			drawText( "Smith-Waterman alignment result:", 0,0,0 );

			m_Cursor.moveY( 2 );

			printSeqPair( vSeqs );
			m_SharpMgr.endDraw();

		}
		else
		{
			CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "two sequences size is not equal!" );
		}
	}
	else
	{
		CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "Wrong sequences number" );
	}
}

void		CSeqRenderer::printMSAGAlignedSeqs( const std::vector<CSequence>& vSeqs, const std::string& strDrawText )
{
	if( vSeqs.size()>1 )
	{
			//if the current line is not first line, step another new line to display
			m_Cursor.setStep( PURE_TEXT_WIDTH, PURE_TEXT_HEIGHT );
			if( abs(m_Cursor.getYPos())>EPS )
			{
				m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+PRINT_LINE_INTERVAL*m_Cursor.getStepY() );
			}
			else
			{
				m_Cursor.setPos( START_POS_X, START_POS_Y );
			}

			m_SharpMgr.beginDraw();
			drawText( strDrawText, 0,0,0 );
			printMSASeqs( vSeqs );
			m_SharpMgr.endDraw();
	}
	else
	{
		CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "Wrong sequences number" );
	}
}

void		CSeqRenderer::printBlastAlignedSeqs( const std::vector<CSequence>& vSeqs )
{
	if( vSeqs.size()%2==0 )
	{
		//if the current line is not first line, step another new line to display
		m_Cursor.setStep( PURE_TEXT_WIDTH, PURE_TEXT_HEIGHT );
		if( abs(m_Cursor.getYPos())>EPS )	
		{	
			m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+PRINT_LINE_INTERVAL*m_Cursor.getStepY() );	
		}	
		else	
		{	
			m_Cursor.setPos( START_POS_X, START_POS_Y );	
		}	

		m_SharpMgr.beginDraw();

		drawText( "BLAST alignment result:", 0,0,0 );	

		for ( int i=0; i<(int)(vSeqs.size()); i+=2 )
		{
			m_Cursor.setStep( PURE_TEXT_WIDTH, PURE_TEXT_HEIGHT );
			m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+PRINT_LINE_INTERVAL*m_Cursor.getStepY() );
			char strSeqNoInfo[100];
			sprintf_s( strSeqNoInfo, "Sequence pair no: %d ", i/2+1 );
			drawText( strSeqNoInfo, 0,0,0 );
			std::vector<CSequence> SeqPair;
			SeqPair.push_back( vSeqs[i] );
			SeqPair.push_back( vSeqs[i+1] );
			printSeqPair( SeqPair );				
		}

		m_SharpMgr.endDraw();

	}
	else
	{
		CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "Wrong sequences number" );
	}	
}

void	CSeqRenderer::printInputSeqs(const std::vector<CSequence>& vSeqs, const std::string& fileName)
{
	if ( (int)(vSeqs.size())==0 )	return;
	
	//if the current line is not first line, step another new line to display
	m_Cursor.setStep( PURE_TEXT_WIDTH, PURE_TEXT_HEIGHT );
	if( abs(m_Cursor.getYPos())>EPS )	
	{	
		m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+PRINT_LINE_INTERVAL*m_Cursor.getStepY() );	
	}	
	else	
	{	
		m_Cursor.setPos( START_POS_X, START_POS_Y );	
	}	

	m_SharpMgr.beginDraw();

	drawText("Input sequences is listed below:" + fileName, 0, 0, 0);

	for ( int i=0; i<(int)(vSeqs.size()); ++i )
	{
		int titleCharNum = (int)(TITLE_WIDTH/PURE_TEXT_WIDTH);

		if( titleCharNum<(int)(vSeqs[i].getName().size()) )
		{
			CGlobalSpace::m_sEventLog.writeWarning( __EXCEPTION_SITE__, "sequences name is to long to print!" );
		}

		m_Cursor.setStep( RECT_TEXT_WIDTH, RECT_TEXT_HEIGHT );
		int iLineNum = (int)ceil((float)vSeqs[i].getSequenceContext().size()/(float)RECT_TEXT_PER_LINE);

		//new line
		m_Cursor.setStep( PURE_TEXT_WIDTH, PURE_TEXT_HEIGHT );
		m_Cursor.setPos( START_POS_X, m_Cursor.getYPos()+3*m_Cursor.getStepY() );
		drawText( vSeqs[i].getName(), 0,0,0 );

		float fCurNameEndPosX = m_Cursor.getXPos();

		m_Cursor.setXPos( fCurNameEndPosX+TITLE_CONTENT_INTERVAL*m_Cursor.getStepX());
		int iNo=0;
		int iCurPos=0;
		while( iCurPos<int(vSeqs[i].getSequenceContext().size()) )
		{
			if( (m_Cursor.getXPos()+m_Cursor.getStepX()) <END_POS_X)
			{
				m_Cursor.setYPos( m_Cursor.getYPos()-m_Cursor.getStepY());
				char ch[10];
				sprintf_s( ch,"%d",iNo );
				drawText( ch, 0,0,0 );
				++iNo;
				m_Cursor.setPos( m_Cursor.getXPos()-m_Cursor.getStepX(), m_Cursor.getYPos()+m_Cursor.getStepY());
				std::string str;
				str.push_back( vSeqs[i].getSequenceContext().at(iCurPos).m_char );
				drawText( str, 0,0,0 );
				++iCurPos;
			}
			else
			{
				m_Cursor.setPos( fCurNameEndPosX+TITLE_CONTENT_INTERVAL*m_Cursor.getStepX(), m_Cursor.getYPos()+2*m_Cursor.getStepY());
			}
		}
	}
	m_SharpMgr.endDraw();
}

void	CSeqRenderer::countSeqPair( const std::vector<CSequence>& vSeqs, int& voNEq, int& voNGap )
{
	const CSeqData& Seq1=vSeqs[0].getSequenceContext();
	const CSeqData& Seq2=vSeqs[1].getSequenceContext();
	int nSeq = Seq1.size();
	voNEq=0;
	voNGap=0;
	for ( int i=0; i<nSeq; ++i )
	{
		voNEq += (Seq1[i].m_char == Seq2[i].m_char)?1:0;
		voNGap += (Seq1[i].m_char == GENESPACE )?1:0 + (Seq2[i].m_char == GENESPACE )?1:0;		
	}
}

void	CSeqRenderer::draw(CDC *vpDC)
{
	m_SharpMgr.draw(vpDC);
}

void	CSeqRenderer::pageUp()
{
	m_SharpMgr.pageUp();
}

void	CSeqRenderer::pageDown()
{
	m_SharpMgr.pageDown();
}

void	CSeqRenderer::OnScrollUp()
{
	m_SharpMgr.scrollUp();
}

void	CSeqRenderer::OnScrollDown()
{
	m_SharpMgr.scrollDown();
}

}