
// SeqAnsisView.cpp : CSeqAnsisView ���ʵ��
//

#include "stdafx.h"
// SHARED_HANDLERS ������ʵ��Ԥ��������ͼ������ɸѡ�������
// ATL ��Ŀ�н��ж��壬�����������Ŀ�����ĵ����롣
#ifndef SHARED_HANDLERS
#include "SeqAnsis.h"
#endif

#include "SeqAnsisDoc.h"
#include "SeqAnsisView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CSeqAnsisView

IMPLEMENT_DYNCREATE(CSeqAnsisView, CView)

BEGIN_MESSAGE_MAP(CSeqAnsisView, CView)
	// ��׼��ӡ����
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CSeqAnsisView::OnFilePrintPreview)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_WM_KEYDOWN()
END_MESSAGE_MAP()

// CSeqAnsisView ����/����

CSeqAnsisView::CSeqAnsisView()
{
	// TODO: �ڴ˴���ӹ������

}

CSeqAnsisView::~CSeqAnsisView()
{
}

BOOL CSeqAnsisView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: �ڴ˴�ͨ���޸�
	//  CREATESTRUCT cs ���޸Ĵ��������ʽ

	return CView::PreCreateWindow(cs);
}

// CSeqAnsisView ����

void CSeqAnsisView::OnDraw(CDC* pDC)
{
	CSeqAnsisDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: �ڴ˴�Ϊ����������ӻ��ƴ���
	pDoc->m_SeqRenderer.draw(pDC);
}

void CSeqAnsisView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	switch( nChar )
	{
	case VK_PRIOR:
		{
			GetDocument()->OnPageUp();
		}
		break;
	case VK_NEXT:
		{
			GetDocument()->OnPageDown();
		}
		break;
	case VK_UP:
		{
			GetDocument()->OnScrollUp();
		}
		break;
	case VK_DOWN:
		{
			GetDocument()->OnScrollDown();
		}
	}
	Invalidate();
}


// CSeqAnsisView ��ӡ


void CSeqAnsisView::OnFilePrintPreview()
{
#ifndef SHARED_HANDLERS
	AFXPrintPreview(this);
#endif
}

BOOL CSeqAnsisView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// Ĭ��׼��
	return DoPreparePrinting(pInfo);
}

void CSeqAnsisView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: ��Ӷ���Ĵ�ӡǰ���еĳ�ʼ������
}

void CSeqAnsisView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: ��Ӵ�ӡ����е��������
}

void CSeqAnsisView::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CSeqAnsisView::OnContextMenu(CWnd* /* pWnd */, CPoint point)
{
#ifndef SHARED_HANDLERS
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
}


// CSeqAnsisView ���

#ifdef _DEBUG
void CSeqAnsisView::AssertValid() const
{
	CView::AssertValid();
}

void CSeqAnsisView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CSeqAnsisDoc* CSeqAnsisView::GetDocument() const // �ǵ��԰汾��������
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CSeqAnsisDoc)));
	return (CSeqAnsisDoc*)m_pDocument;
}
#endif //_DEBUG


// CSeqAnsisView ��Ϣ�������
