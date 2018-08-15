
// SeqAnsisDoc.h : CSeqAnsisDoc 类的接口
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
protected: // 仅从序列化创建
	CSeqAnsisDoc();
	DECLARE_DYNCREATE(CSeqAnsisDoc)

// 特性
public:

// 操作
public:

// 重写
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
#ifdef SHARED_HANDLERS
	virtual void InitializeSearchContent();
	virtual void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds);
#endif // SHARED_HANDLERS

// 实现
public:
	virtual ~CSeqAnsisDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// 生成的消息映射函数
protected:
	DECLARE_MESSAGE_MAP()

#ifdef SHARED_HANDLERS
	// 用于为搜索处理程序设置搜索内容的 Helper 函数
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
