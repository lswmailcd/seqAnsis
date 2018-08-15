#pragma once


// CCudaMsagaInputDlg dialog

class CCudaMsagaInputDlg : public CDialog
{
	DECLARE_DYNAMIC(CCudaMsagaInputDlg)

public:
	CCudaMsagaInputDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~CCudaMsagaInputDlg();

// Dialog Data
	enum { IDD = IDD_DIALOG_CUDA_MSAGA_PARAM };

	CString m_strPopulationNum;
	CString m_strNoAdvGenerationNum;
	CString m_strMaxOrganLen;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};
