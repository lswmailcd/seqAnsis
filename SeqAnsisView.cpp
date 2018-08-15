
// SeqAnsisView.cpp : CSeqAnsisView 类的实现
//

#include "stdafx.h"
// SHARED_HANDLERS 可以在实现预览、缩略图和搜索筛选器句柄的
// ATL 项目中进行定义，并允许与该项目共享文档代码。
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
	// 标准打印命令
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CSeqAnsisView::OnFilePrintPreview)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_WM_KEYDOWN()
END_MESSAGE_MAP()

// CSeqAnsisView 构造/析构

CSeqAnsisView::CSeqAnsisView()
{
	// TODO: 在此处添加构造代码

}

CSeqAnsisView::~CSeqAnsisView()
{
}

BOOL CSeqAnsisView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: 在此处通过修改
	//  CREATESTRUCT cs 来修改窗口类或样式

	return CView::PreCreateWindow(cs);
}

// CSeqAnsisView 绘制

void CSeqAnsisView::OnDraw(CDC* pDC)
{
	CSeqAnsisDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: 在此处为本机数据添加绘制代码
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


// CSeqAnsisView 打印


void CSeqAnsisView::OnFilePrintPreview()
{
#ifndef SHARED_HANDLERS
	AFXPrintPreview(this);
#endif
}

BOOL CSeqAnsisView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// 默认准备
	return DoPreparePrinting(pInfo);
}

void CSeqAnsisView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: 添加额外的打印前进行的初始化过程
}

void CSeqAnsisView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: 添加打印后进行的清理过程
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


// CSeqAnsisView 诊断

#ifdef _DEBUG
void CSeqAnsisView::AssertValid() const
{
	CView::AssertValid();
}

void CSeqAnsisView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CSeqAnsisDoc* CSeqAnsisView::GetDocument() const // 非调试版本是内联的
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CSeqAnsisDoc)));
	return (CSeqAnsisDoc*)m_pDocument;
}
#endif //_DEBUG


// CSeqAnsisView 消息处理程序
