#include<memory>
#include"Common.h"
#include"SWAlgorithm.h"
#include"GlobalSpace.h"


namespace SeqAnsis
{

const float fMismatchScore = CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost( SW_MAT_TYPE );
#define W(x)  (1.0f+fMismatchScore*(x))

void CSWAlgorithm::Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences )
{
	if( vSequences.size()!=2 )// || vAlignedSequences.size()!=2 )
	{
		CGlobalSpace::m_sEventLog.writeEvent("Wrong sequences size to align by smith-waterman algorithm!");
		return;
	}

	int nSeq1 = vSequences[0].getSequenceContext().size();
	int nSeq2 = vSequences[1].getSequenceContext().size();
	const CSeqData& Seq1 = vSequences[0].getSequenceContext();
	const CSeqData& Seq2 = vSequences[1].getSequenceContext();

	float *arrScore = new float[(nSeq1+1)*(nSeq2+1)];
	memset( arrScore, 0, sizeof(float)*(nSeq1+1)*(nSeq2+1) );

	int iPathStart=-1, jPathStart=-1;
	float fMaxScore=-1e10f;

		//we define seq1 as the row sequence and the seq2 as the column.
	for ( int row=(nSeq2-1); row>=0; --row )
	{
		for ( int col=(nSeq1-1); col>=0; --col )
		{
			//if the pair is identical
			arrScore[POS(row,col,nSeq1+1)] = (Seq1[col].m_iCode == Seq2[row].m_iCode)? 
				CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getSubMatrixScore( SW_MAT_TYPE, Seq1[col].m_iCode, Seq2[row] .m_iCode):CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(SW_MAT_TYPE);

			////get the H(i,j) with the ending ai&bj
			float maxScore_aibj = arrScore[POS(row,col,nSeq1+1)] + arrScore[POS(row+1,col+1,nSeq1+1)];//(maxSubRowScore>maxSubColScore)?maxSubRowScore:maxSubColScore;

			//get the H(i,j) without ai.(find the max score in the row j)
			float maxSubRowScore = arrScore[POS(row,col+1,nSeq1+1)] - W(1);
			for ( int subcol = col+2; subcol<(nSeq1+1); ++subcol )
			{
				maxSubRowScore = (maxSubRowScore<(arrScore[POS(row,subcol,nSeq1+1)]-W(subcol-col)))? (arrScore[POS(row,subcol,nSeq1+1)]-W(subcol-col)):maxSubRowScore;
			}

			//get the H(i,j) without bj.(find the max score in the column i )
			float maxSubColScore = arrScore[POS(row+1,col,nSeq1+1)] - W(1);
			for ( int subrow = row+2; subrow<(nSeq2+1); ++subrow )
			{
				maxSubColScore = (maxSubColScore<(arrScore[POS(subrow,col,nSeq1+1)]-W(subrow-row)))? (arrScore[POS(subrow,col,nSeq1+1)]-W(subrow-row)):maxSubColScore;
			}			
					
			//get the final result of arrScore[row,col] 
			//which is the maximum of (maxScore_aibj, maxSubRowScore, maxSubColScore)
			arrScore[POS(row,col,nSeq1+1)] = (maxSubRowScore>maxSubColScore)?maxSubRowScore:maxSubColScore;
			arrScore[POS(row,col,nSeq1+1)] = (arrScore[POS(row,col,nSeq1+1)]>maxScore_aibj)?arrScore[POS(row,col,nSeq1+1)]:maxScore_aibj;
			arrScore[POS(row,col,nSeq1+1)] = (arrScore[POS(row,col,nSeq1+1)]<0)?0:arrScore[POS(row,col,nSeq1+1)];
			if ( arrScore[POS(row,col,nSeq1+1)]>fMaxScore )
			{
				fMaxScore = arrScore[POS(row,col,nSeq1+1)];
				iPathStart = row;
				jPathStart = col;
			}
		}
	}

	//find one path in the arrScore 
	//(startRow,startCol)为当前值最大的矩阵元素的行列号。
	//	1\如果(startRow,startCol)的值为零
	//	　　　则结束
	//	   否则
	//	   　　从startRow开始直到行尾，查找最大值所在的列号
	//		 　从startCol开始直到列尾，查找最大值所在的行号

	//		   比较行列最大值，找到最大的值，记录位置i,j
	//			   如果i==startRow且j==startCol,
	//			      则p1[idxSeq1++]=seq1[startCol],p2[idxSeq2++]=seq2[startRow];
	//                startRow=i+1,startCol=j+1,返回1执行
	//	           否则
	//	             如果i>startRow,
	//					从startRow到i-1的长度为len，将p1[idxseq1...idxseq1+len]='-'
	//					p2[idxseq2...idxseq2+len]=seq2[startRow...startRow+len]
	//					p1[idxSeq1++]=seq1[j],p2[idxSeq2++]=seq2[i]
	//				 否则
	//					从startCol到j-1的长度为len，将p2[idxseq1...idxseq1+len]='-'
	//					p1[idxseq1...idxseq1+len]=seq2[startCol...startCol+len]
	//					p1[idxseq1+len+1]=seq1[j],p2[idxSeq2++]=seq2[i]

	//startRow=i+1,startCol=j+1,返回1执行

	int startRow=iPathStart,startCol=jPathStart;//the start point index of the sub matrix.
	int idxSeq1=0,idxSeq2=0;//the index of the rearranged sequences.
	int idxRow=0,idxCol=0;//current max score point in the sub matrix.
	float maxColScore =-1e10f;
	float maxRowScore = maxColScore;
	CSeqData ReSeq1;// = vAlignedSequences[0].getSequenceContext();
	CSeqData ReSeq2;// = vAlignedSequences[1].getSequenceContext();

	while (true)
	{
		//如果(startRow,startCol)的值为零,则结束查找路径
		if ( abs(arrScore[POS(startRow,startCol,nSeq1+1)]) < EPS )
		{
			break;
		}
		else
		{
			//the max score in the same row
			maxColScore = arrScore[POS(startRow,startCol,nSeq1+1)];
			//the max score in the same column
			maxRowScore = maxColScore;
			idxCol = startCol;
			idxRow = startRow;
			//从（startRow，startCol)开始逐列查找直到行尾，查找最大值所在的列号
			for ( int col=(startCol+1); col<nSeq1; ++col )
			{
				if ( maxColScore < arrScore[POS(startRow,col,nSeq1+1)] )
				{
					maxColScore = arrScore[POS(startRow,col,nSeq1+1)];
					idxCol = col;
				}
			}
			//从（startRow，startCol)开始逐行查找直到列尾，查找最大值所在的行号
			for ( int row=(startRow+1); row<nSeq2; ++row )
			{
				if ( maxRowScore<arrScore[POS(row,startCol,nSeq1+1)] )
				{
					maxRowScore = arrScore[POS(row,startCol,nSeq1+1)];
					idxRow = row;
				}
			}
		}
		//比较行列最大值，找到最大的值(当行列最大值相同时，优先考虑行最大值maxColScore)，记录位置idxRow,idxCol
		if ( maxColScore < maxRowScore )
		{
			idxCol = startCol;
		}
		else
		{
			idxRow = startRow;
		}
		//	如果i==startRow且j==startCol,
		//	则p1[idxSeq1++]=seq1[startCol],p2[idxSeq2++]=seq2[startRow];
		//startRow=i+1,startCol=j+1,返回1执行
		if ( idxRow == startRow && idxCol == startCol )
		{
			ReSeq1.push_back(Seq1[startCol]);
			ReSeq2.push_back(Seq2[startRow]);
		}
		else
		{
			//	如果i>startRow,
			//	从startRow到i-1的长度为len，将p1[idxseq1...idxseq1+len]='-'
			//	p2[idxseq2...idxseq2+len]=seq2[startRow...startRow+len]
			//p1[idxSeq1++]=seq1[j],p2[idxSeq2++]=seq2[i]
			if ( idxRow > startRow )
			{
				for ( int i=startRow; i<idxRow; ++i )
				{
					ReSeq1.push_back(StruSeqElem(GENESPACE,-1));
					ReSeq2.push_back(Seq2[i]);
				}
				ReSeq1.push_back(Seq1[idxCol]);
				ReSeq2.push_back(Seq2[idxRow]);
			}
			//否则
			//	从startCol到j-1的长度为len，将p2[idxseq1...idxseq1+len]='-'
			//	p1[idxseq1...idxseq1+len]=seq1[startCol...startCol+len]
			//p1[idxseq1+len+1]=seq1[j],p2[idxSeq2++]=seq2[i]
			else
			{
				for ( int i=startCol; i<idxCol; ++i )
				{
					ReSeq2.push_back(StruSeqElem(GENESPACE, -1));
					ReSeq1.push_back(Seq1[i]);
				}
				ReSeq1.push_back(Seq1[idxCol]);
				ReSeq2.push_back(Seq2[idxRow]);
			}
		}
		//startRow=i+1,startCol=j+1,返回1执行
		startRow = idxRow+1;
		startCol = idxCol+1;
	}

	CSequence  RS1( ReSeq1, vSequences[0].getName(), vSequences[0].getTitle() );
	CSequence  RS2( ReSeq2, vSequences[1].getName(), vSequences[1].getTitle() );
	vAlignedSequences.push_back(RS1);
	vAlignedSequences.push_back(RS2);

	SAFE_DELETE(arrScore);
}

}