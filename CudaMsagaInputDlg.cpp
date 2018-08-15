// CudaMsagaInputDlg.cpp : implementation file
//

#include "stdafx.h"
#include "SeqAnsis.h"
#include "CudaMsagaInputDlg.h"
#include "afxdialogex.h"


// CCudaMsagaInputDlg dialog

IMPLEMENT_DYNAMIC(CCudaMsagaInputDlg, CDialog)

CCudaMsagaInputDlg::CCudaMsagaInputDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CCudaMsagaInputDlg::IDD, pParent)
{
	m_strPopulationNum = "1024"; 
	m_strNoAdvGenerationNum = "50";
	m_strMaxOrganLen = "256";
}

CCudaMsagaInputDlg::~CCudaMsagaInputDlg()
{
}

void CCudaMsagaInputDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text( pDX, IDC_EDIT_POPULATION_NUM, m_strPopulationNum );
	DDX_Text(pDX, IDC_EDIT_NOADV_GENUM, m_strNoAdvGenerationNum);
	DDX_Text(pDX, IDC_EDIT_ORGAN_MAX_LEN, m_strMaxOrganLen);
}


BEGIN_MESSAGE_MAP(CCudaMsagaInputDlg, CDialog)
END_MESSAGE_MAP()


// CCudaMsagaInputDlg message handlers
