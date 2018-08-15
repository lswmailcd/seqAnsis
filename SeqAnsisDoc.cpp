
// SeqAnsisDoc.cpp : CSeqAnsisDoc 类的实现
//

#include "stdafx.h"
// SHARED_HANDLERS 可以在实现预览、缩略图和搜索筛选器句柄的
// ATL 项目中进行定义，并允许与该项目共享文档代码。
#ifndef SHARED_HANDLERS
#include "SeqAnsis.h"
#endif

#include "SeqAnsisDoc.h"

#include <propkey.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#include<fstream>
#include<memory>
#include"NWAlgorithm.h"
#include"SWAlgorithm.h"
#include "BlastAlgorithm.h"
#include "SAGEPAlgorithm.h"
#include "MSA-GA.h"
#include "MSA-GA-CUDA.h"
#include"GlobalSpace.h"
#include"FileReader.h"
#include"PhylogeneticTree.h"

// CSeqAnsisDoc

IMPLEMENT_DYNCREATE(CSeqAnsisDoc, CDocument)

BEGIN_MESSAGE_MAP(CSeqAnsisDoc, CDocument)
	ON_COMMAND(ID_FILE_OPEN, &CSeqAnsisDoc::OnFileOpen)
	ON_COMMAND(ID_SA_NW, &CSeqAnsisDoc::OnSaNw)
	ON_COMMAND(ID_SA_SW, &CSeqAnsisDoc::OnSaSw)
	ON_COMMAND(ID_SA_BLAST, &CSeqAnsisDoc::OnSaBlast)
	ON_COMMAND(ID_32783, &CSeqAnsisDoc::OnSAGEP)
	ON_COMMAND(ID_MSA_GA, &CSeqAnsisDoc::OnMsaGa)
	ON_COMMAND(ID_CUDA_MSAGA, &CSeqAnsisDoc::OnCudaMsaga)
END_MESSAGE_MAP()


// CSeqAnsisDoc 构造/析构

CSeqAnsisDoc::CSeqAnsisDoc()
{
	// TODO: 在此添加一次性构造代码
}

CSeqAnsisDoc::~CSeqAnsisDoc()
{
}

BOOL CSeqAnsisDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: 在此添加重新初始化代码
	// (SDI 文档将重用该文档)

	return TRUE;
}




// CSeqAnsisDoc 序列化

void CSeqAnsisDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: 在此添加存储代码
	}
	else
	{
		// TODO: 在此添加加载代码
	}
}

#ifdef SHARED_HANDLERS

// 缩略图的支持
void CSeqAnsisDoc::OnDrawThumbnail(CDC& dc, LPRECT lprcBounds)
{
	// 修改此代码以绘制文档数据
	dc.FillSolidRect(lprcBounds, RGB(255, 255, 255));

	CString strText = _T("TODO: implement thumbnail drawing here");
	LOGFONT lf;

	CFont* pDefaultGUIFont = CFont::FromHandle((HFONT) GetStockObject(DEFAULT_GUI_FONT));
	pDefaultGUIFont->GetLogFont(&lf);
	lf.lfHeight = 36;

	CFont fontDraw;
	fontDraw.CreateFontIndirect(&lf);

	CFont* pOldFont = dc.SelectObject(&fontDraw);
	dc.DrawText(strText, lprcBounds, DT_CENTER | DT_WORDBREAK);
	dc.SelectObject(pOldFont);
}

// 搜索处理程序的支持
void CSeqAnsisDoc::InitializeSearchContent()
{
	CString strSearchContent;
	// 从文档数据设置搜索内容。
	// 内容部分应由“;”分隔

	// 例如:  strSearchContent = _T("point;rectangle;circle;ole object;")；
	SetSearchContent(strSearchContent);
}

void CSeqAnsisDoc::SetSearchContent(const CString& value)
{
	if (value.IsEmpty())
	{
		RemoveChunk(PKEY_Search_Contents.fmtid, PKEY_Search_Contents.pid);
	}
	else
	{
		CMFCFilterChunkValueImpl *pChunk = NULL;
		ATLTRY(pChunk = new CMFCFilterChunkValueImpl);
		if (pChunk != NULL)
		{
			pChunk->SetTextValue(PKEY_Search_Contents, value, CHUNK_TEXT);
			SetChunkValue(pChunk);
		}
	}
}

#endif // SHARED_HANDLERS

// CSeqAnsisDoc 诊断

#ifdef _DEBUG
void CSeqAnsisDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CSeqAnsisDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// CSeqAnsisDoc 命令


void CSeqAnsisDoc::OnFileOpen()
{
	// TODO: 在此添加命令处理程序代码
	CFileDialog dlg(TRUE, _T("*"), NULL, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
		_T("所有文件(*.*)|*.*||"), NULL);

	if(dlg.DoModal() == IDOK)
	{
		m_filePathName = std::string(CT2CA(dlg.GetPathName())); // 文件路径
		m_fileName = std::string(CT2CA(dlg.GetFileTitle()));
		try
		{
			std::auto_ptr<SeqAnsis::CFileReader> pFileReader(new SeqAnsis::CFileReader(m_filePathName));
			int iFirstSeq=0;
			std::string OffendingSeq;
			m_Sequences.clear();
			pFileReader->readSeqs( m_Sequences, iFirstSeq, OffendingSeq );
			m_SeqRenderer.printInputSeqs(m_Sequences, m_filePathName);

			POSITION pos = GetFirstViewPosition();
			GetNextView(pos)->Invalidate();
			
		}
		catch(SeqAnsis::CAppException const& e)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeException(e);
			return;
		}
		catch(...)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("CSeqAnsisDoc::OnFileOpen() failure!");
		}
	} 
}


void CSeqAnsisDoc::OnSaNw()
{
	// TODO: 在此添加命令处理程序代码
	if( !m_Sequences.empty() )
	{
		try
		{
			std::auto_ptr<SeqAnsis::CNWAlgorithm> pNWAlgorithm( new SeqAnsis::CNWAlgorithm );
			if( m_Sequences.size()==2 )
			{
				int len1 = m_Sequences[0].getSequenceContext().size();
				int len2 = m_Sequences[1].getSequenceContext().size();
				m_NWAlignedSequences.clear();
				pNWAlgorithm->Align( m_NWAlignedSequences, m_Sequences );

				CRect rc;
				POSITION pos = GetFirstViewPosition();
				GetNextView(pos)->GetClientRect(&rc);
				m_SeqRenderer.setRenderArea( rc.left, rc.top, rc.Width(), rc.Height() );

				m_SeqRenderer.printNWAlignedSeqs( m_NWAlignedSequences );
				std::auto_ptr<SeqAnsis::CFileWriter>		pSeqFileWriter( new SeqAnsis::CFileWriter( std::string("seqTest_NW.txt") ) );
				pSeqFileWriter->seqOutput( m_NWAlignedSequences );

				pos = GetFirstViewPosition();
				GetNextView(pos)->Invalidate();
			}
			else
			{
				SeqAnsis::CGlobalSpace::m_sEventLog.writeWarning(__EXCEPTION_SITE__, "sequence number is not right to be handled with Needman-Wunsch algorithm.");
			}
		}
		catch(SeqAnsis::CAppException const& e)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeException(e);
			return;
		}
		catch(...)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("CSeqAnsisDoc::OnSaNw() failure!");
		}
	}
	else
	{
		SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("Sequence is empty to handled by Needleman-wunsch method!");
	}
}


void CSeqAnsisDoc::OnSaSw()
{
	// TODO: 在此添加命令处理程序代码
	if( !m_Sequences.empty() )
	{
		try
		{
			std::auto_ptr<SeqAnsis::CSWAlgorithm> pSWAlgorithm( new SeqAnsis::CSWAlgorithm );
			if( m_Sequences.size()==2 )
			{
				int len1 = m_Sequences[0].getSequenceContext().size();
				int len2 = m_Sequences[1].getSequenceContext().size();
				m_SWAlignedSequences.clear();

				pSWAlgorithm->Align( m_SWAlignedSequences, m_Sequences );

				CRect rc;
				POSITION pos = GetFirstViewPosition();
				GetNextView(pos)->GetClientRect(&rc);
				m_SeqRenderer.setRenderArea( rc.left, rc.top, rc.Width(), rc.Height() );

				m_SeqRenderer.printSWAlignedSeqs( m_SWAlignedSequences );

				pos = GetFirstViewPosition();
				GetNextView(pos)->Invalidate();
			}
			else
			{
				SeqAnsis::CGlobalSpace::m_sEventLog.writeWarning(__EXCEPTION_SITE__, "sequence number is not right to be handled with Needman-Wunsch algorithm.");
			}
		}
		catch(SeqAnsis::CAppException const& e)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeException(e);
			return;
		}
		catch(...)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("CSeqAnsisDoc::OnSaNw() failure!");
		}
	}
	else
	{
		SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("Sequence is empty to handled by Needleman-wunsch method!");
	}
}


void CSeqAnsisDoc::OnSaBlast()
{
#if 0
	// TODO: 在此添加命令处理程序代码
	if( !m_Sequences.empty() )
	{
		try
		{
			std::auto_ptr<SeqAnsis::CBlastAlgorithm> pBlastAlgorithm( new SeqAnsis::CBlastAlgorithm );
			if( m_Sequences.size()==2 )
			{
				int len1 = m_Sequences[0].getSequenceContext().size();
				int len2 = m_Sequences[1].getSequenceContext().size();
				m_BlastAlignedSequences.clear();

				SeqAnsis::CTimer  time;
				double dStartTime = time.getCurrentTime();
				pBlastAlgorithm->Align( m_BlastAlignedSequences, m_Sequences );
				double dInterval = time.getCurrentTime()-dStartTime;
				char ch[100];
				sprintf_s( ch, "the process time of blast algorithm is %f", dInterval );
				SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent( ch );

				CRect rc;
				POSITION pos = GetFirstViewPosition();
				GetNextView(pos)->GetClientRect(&rc);
				m_SeqRenderer.setRenderArea( rc.left, rc.top, rc.Width(), rc.Height() );

				m_SeqRenderer.printBlastAlignedSeqs( m_BlastAlignedSequences );

				pos = GetFirstViewPosition();
				GetNextView(pos)->Invalidate();
			}
			else
			{
				SeqAnsis::CGlobalSpace::m_sEventLog.writeWarning(__EXCEPTION_SITE__, "sequence number is not right to be handled with BLAST algorithm.");
			}
		}
		catch(SeqAnsis::CAppException const& e)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeException(e);
			return;
		}
		catch(...)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("CSeqAnsisDoc::OnSaBlast() failure!");
		}
	}
	else
	{
		SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("Sequence is empty to handled by BLAST method!");
	}
#endif
}

void CSeqAnsisDoc::OnPageDown()
{
	m_SeqRenderer.pageDown();
}

void CSeqAnsisDoc::OnPageUp()
{
	m_SeqRenderer.pageUp();
}

void CSeqAnsisDoc::OnScrollUp()
{
	m_SeqRenderer.OnScrollUp();
}

void CSeqAnsisDoc::OnScrollDown()
{
	m_SeqRenderer.OnScrollDown();
}


void CSeqAnsisDoc::OnSAGEP()
{
#if 0
	// TODO: 在此添加命令处理程序代码
	if( !m_Sequences.empty() )
	{
		try
		{
			std::auto_ptr<SeqAnsis::CSAGEPAlgorithm> pSAGEPAlgorithm( new SeqAnsis::CSAGEPAlgorithm );
			if( m_Sequences.size()==2 )
			{
				int len1 = m_Sequences[0].getSequenceContext().size();
				int len2 = m_Sequences[1].getSequenceContext().size();
				m_SAGEPAlignedSequences.clear();

				SeqAnsis::CTimer  time;
				double dStartTime = time.getCurrentTime();
				pSAGEPAlgorithm->Align( m_SAGEPAlignedSequences, m_Sequences );
				double dInterval = time.getCurrentTime()-dStartTime;
				char ch[100];
				sprintf_s( ch, "the process time of SAGEP algorithm is %f", dInterval );
				SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent( ch );

				CRect rc;
				POSITION pos = GetFirstViewPosition();
				GetNextView(pos)->GetClientRect(&rc);
				m_SeqRenderer.setRenderArea( rc.left, rc.top, rc.Width(), rc.Height() );

				m_SeqRenderer.printSAGEPAlignedSeqs( m_SAGEPAlignedSequences );
				std::auto_ptr<SeqAnsis::CFileWriter>		pSeqFileWriter( new SeqAnsis::CFileWriter( std::string("seqTest_SAGEP.txt") ) );
				pSeqFileWriter->seqOutput( m_SAGEPAlignedSequences );

				pos = GetFirstViewPosition();
				GetNextView(pos)->Invalidate();
			}
			else
			{
				SeqAnsis::CGlobalSpace::m_sEventLog.writeWarning(__EXCEPTION_SITE__, "sequence number is not right to be handled with SAGEP algorithm.");
			}
		}
		catch(SeqAnsis::CAppException const& e)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeException(e);
			return;
		}
		catch(...)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("CSeqAnsisDoc::OnSaGEP() failure!");
		}
	}
	else
	{
		SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("Sequence is empty to handled by SAGEP method!");
	}
#endif
}


void CSeqAnsisDoc::OnMsaGa()
{
	// TODO: 在此添加命令处理程序代码
	if( !m_Sequences.empty() )
	{
		try
		{
			if (m_DlgCudaMsaga.DoModal() == IDOK)
			{
				std::auto_ptr<SeqAnsis::CMSAGAAlgorithm> pMSAGAAlgorithm(new SeqAnsis::CMSAGAAlgorithm);
				if (m_Sequences.size() > 1)
				{
					SeqAnsis::CPhylogeneticTree TreeWeight;
					std::vector<float>	SeqWeight;
					std::string strWeightTreeFile = m_filePathName.substr(0, m_filePathName.length() - 3) + "ph";
					if (TreeWeight.readTree(strWeightTreeFile, m_Sequences))
					{
						char ch[300];
						sprintf_s(ch, "found phlipy format weight tree file of %s!", strWeightTreeFile.c_str());
						SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent(ch);
						TreeWeight.calcSeqWeights(m_Sequences.size(), SeqWeight);
					}
					else
					{
						char ch[300];
						sprintf_s(ch, "NOT found phlipy format weight tree file of %s!", strWeightTreeFile.c_str());
						SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent(ch);
					}

					pMSAGAAlgorithm->testSPScore(m_Sequences);
					m_MSAGAlignedSequences.clear();
					pMSAGAAlgorithm->SetAlignParams(_wtoi(m_DlgCudaMsaga.m_strPopulationNum.GetBuffer()), SeqWeight);
					pMSAGAAlgorithm->Align(m_MSAGAlignedSequences, m_Sequences);
					CRect rc;
					POSITION pos = GetFirstViewPosition();
					GetNextView(pos)->GetClientRect(&rc);
					m_SeqRenderer.setRenderArea(rc.left, rc.top, rc.Width(), rc.Height());

					m_SeqRenderer.printMSAGAlignedSeqs(m_MSAGAlignedSequences, "Print MSA-GA Seq aln:" + m_filePathName);
					std::auto_ptr<SeqAnsis::CFileWriter>		pSeqFileWriter(new SeqAnsis::CFileWriter(std::string("seqTest_MSAGA.txt")));
					pSeqFileWriter->openFile();
					pSeqFileWriter->seqOutput(m_MSAGAlignedSequences);
					pSeqFileWriter->closeFile();

					pos = GetFirstViewPosition();
					GetNextView(pos)->Invalidate();
				}
				else
				{
					SeqAnsis::CGlobalSpace::m_sEventLog.writeWarning(__EXCEPTION_SITE__, "sequence number is not right to be handled with MSA-GA algorithm.");
				}
			}
		}
		catch(SeqAnsis::CAppException const& e)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeException(e);
			return;
		}
		catch(...)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("CSeqAnsisDoc::OnMsaGa() failure!");
		}
	}
	else
	{
		SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("Sequence is empty to handled by BL MSA-GA method!");
	}
}


void CSeqAnsisDoc::OnCudaMsaga()
{
	// TODO: 在此添加命令处理程序代码
	if( !m_Sequences.empty() )
	{
		try
		{
			if (m_DlgCudaMsaga.DoModal() == IDOK)
			{
				int nPopulationNum = _wtoi(m_DlgCudaMsaga.m_strPopulationNum.GetBuffer());
				int nNoAdvGenerationNum = _wtoi(m_DlgCudaMsaga.m_strNoAdvGenerationNum.GetBuffer());
				int nMaxOrganLen = _wtoi(m_DlgCudaMsaga.m_strMaxOrganLen.GetBuffer());		

				std::auto_ptr<SeqAnsis::CMSAGA_CUDA_Algorithm> pCUDA_MSAGAAlgorithm(new SeqAnsis::CMSAGA_CUDA_Algorithm);
				if (m_Sequences.size() > 1)
				{
					SeqAnsis::CPhylogeneticTree TreeWeight;
					std::vector<float>	SeqWeight;
					std::string strWeightTreeFile = m_filePathName.substr(0, m_filePathName.length() - 3) + "ph";
					if (TreeWeight.readTree(strWeightTreeFile, m_Sequences))
					{
						char ch[300];
						sprintf_s(ch, "found phlipy format weight tree file of %s!", strWeightTreeFile.c_str());
						SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent(ch);
						TreeWeight.calcSeqWeights(m_Sequences.size(), SeqWeight);
					}
					else
					{
						char ch[300];
						sprintf_s(ch, "NOT found phlipy format weight tree file of %s!", strWeightTreeFile.c_str());
						SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent(ch);
					}
					pCUDA_MSAGAAlgorithm->SetAlignParams(nPopulationNum, nNoAdvGenerationNum, nMaxOrganLen, SeqWeight);
					m_CUDA_MSAGAlignedSequences.clear();					
					pCUDA_MSAGAAlgorithm->Align(m_CUDA_MSAGAlignedSequences, m_Sequences);
					CRect rc;
					POSITION pos = GetFirstViewPosition();
					GetNextView(pos)->GetClientRect(&rc);
					m_SeqRenderer.setRenderArea(rc.left, rc.top, rc.Width(), rc.Height());

					m_SeqRenderer.printMSAGAlignedSeqs(m_CUDA_MSAGAlignedSequences, "Print CUDA-MSA-GA Seq aln: " + m_filePathName);
					std::auto_ptr<SeqAnsis::CFileWriter>		pSeqFileWriter(new SeqAnsis::CFileWriter(std::string("seqTest_CUDA_MSAGA.txt")));
					pSeqFileWriter->openFile();
					pSeqFileWriter->seqOutput(m_CUDA_MSAGAlignedSequences);
					pSeqFileWriter->closeFile();

					pos = GetFirstViewPosition();
					GetNextView(pos)->Invalidate();
				}
				else
				{
					SeqAnsis::CGlobalSpace::m_sEventLog.writeWarning(__EXCEPTION_SITE__, "sequence number is not right to be handled with MSA-GA algorithm.");
				}
			}
		}
		catch(SeqAnsis::CAppException const& e)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeException(e);
			return;
		}
		catch(...)
		{
			SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("CSeqAnsisDoc::OnMsaGa() failure!");
		}
	}
	else
	{
		SeqAnsis::CGlobalSpace::m_sEventLog.writeEvent("Sequence is empty to handled by BL MSA-GA method!");
	}
}
