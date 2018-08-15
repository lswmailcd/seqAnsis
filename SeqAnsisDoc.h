
// SeqAnsisDoc.h : CSeqAnsisDoc ��Ľӿ�
//


#pragma once
#include<string>
#include<vector>
#include"Sequence.h"
#include"SeqRenderer.h"
#include "FileWriter.h"
#include"CudaMsagaInputDlg.h"

class CSeqAnsisDoc : public CDocument
{
protected: // �������л�����
	CSeqAnsisDoc();
	DECLARE_DYNCREATE(CSeqAnsisDoc)

// ����
public:

// ����
public:

// ��д
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
#ifdef SHARED_HANDLERS
	virtual void InitializeSearchContent();
	virtual void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds);
#endif // SHARED_HANDLERS

// ʵ��
public:
	virtual ~CSeqAnsisDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// ���ɵ���Ϣӳ�亯��
protected:
	DECLARE_MESSAGE_MAP()

#ifdef SHARED_HANDLERS
	// ����Ϊ����������������������ݵ� Helper ����
	void SetSearchContent(const CString& value);
#endif // SHARED_HANDLERS
public:
	afx_msg void OnFileOpen();

private:
	std::string m_filePathName;
	std::string m_fileName;
	std::vector<SeqAnsis::CSequence> m_Sequences;
	std::vector<SeqAnsis::CSequence> m_NWAlignedSequences;
	std::vector<SeqAnsis::CSequence> m_SWAlignedSequences;
	std::vector<SeqAnsis::CSequence> m_BlastAlignedSequences;
	std::vector<SeqAnsis::CSequence> m_SAGEPAlignedSequences;
	std::vector<SeqAnsis::CSequence> m_MSAGAlignedSequences;
	std::vector<SeqAnsis::CSequence> m_CUDA_MSAGAlignedSequences;

	CCudaMsagaInputDlg  m_DlgCudaMsaga;
public:
	SeqAnsis::CSeqRenderer		m_SeqRenderer;
	void	OnPageDown();
	void	OnPageUp();
	void    OnScrollUp();
	void    OnScrollDown();
public:
	afx_msg void OnSaNw();
	afx_msg void OnSaSw();
	afx_msg void OnSaBlast();
	afx_msg void OnSAGEP();
	afx_msg void OnMsaGa();
	afx_msg void OnCudaMsaga();
};
