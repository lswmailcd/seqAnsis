#include<memory>
#include"Common.h"
#include"NWAlgorithm.h"
#include"GlobalSpace.h"


namespace SeqAnsis
{

void CNWAlgorithm::Align( std::vector<CSequence>& vAlignedSequences, const std::vector<CSequence>& vSequences )
{
	if( vSequences.size()!=2 )// || vAlignedSequences.size()!=2 )
	{
		CGlobalSpace::m_sEventLog.writeEvent("Wrong sequences size to align by needleman-wunsch algorithm!");
		return;
	}

	int nSeq1 = vSequences[0].getSequenceContext().size();
	int nSeq2 = vSequences[1].getSequenceContext().size();
	const CSeqData& Seq1 = vSequences[0].getSequenceContext();
	const CSeqData& Seq2 = vSequences[1].getSequenceContext();

	int *arrScore = new int[(nSeq1+1)*(nSeq2+1)];
	memset( arrScore, 0, sizeof(int)*(nSeq1+1)*(nSeq2+1) );
	
	//initialize the boundary score of score matrix
	for ( int row=nSeq2; row>=0; --row )
	{
		arrScore[POS(row,nSeq1,nSeq1)] = CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(NW_MAT_TYPE)*(nSeq2-row);
	}
	for ( int col=nSeq1; col>=0; --col )
	{
		arrScore[POS(nSeq2,col,nSeq1)] = CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(NW_MAT_TYPE)*(nSeq1-col);
	}

	//we define seq1 as the row sequence and the seq2 as the column.
	for ( int row=(nSeq2-1); row>=0; --row )
	{
		for ( int col=(nSeq1-1); col>=0; --col )
		{
			//if the pair is identical, use substitution matrix to score
			int matchScore = arrScore[POS(row+1,col+1,nSeq1)]+CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getSubMatrixScore(NW_MAT_TYPE, Seq1[col].m_iCode, Seq2[row].m_iCode);
			int deleteScore = arrScore[POS(row+1,col,nSeq1)]+CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(NW_MAT_TYPE);
			int insertScore = arrScore[POS(row,col+1,nSeq1)]+CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(NW_MAT_TYPE);

			arrScore[POS(row,col,nSeq1)] = CGlobalSpace::m_sUtility.max3( matchScore, deleteScore, insertScore );
		}
	}

	//find one path in the arrScore 
	int startRow=0,startCol=0;//the start point index of the sub matrix.
	int idxSeq1=0,idxSeq2=0;//the index of the rearranged sequences.
	int idxRow=0,idxCol=0;//current max score point in the sub matrix.
	int maxColScore = -1;
	int maxRowScore = maxColScore;
	CSeqData ReSeq1;// = vAlignedSequences[0].getSequenceContext();
	CSeqData ReSeq2;// = vAlignedSequences[1].getSequenceContext();
	while( startRow<nSeq2 || startCol<nSeq1 )
	{
		if (startRow<nSeq2 && startCol<nSeq1
			&& arrScore[POS(startRow,startCol,nSeq1)] ==  (arrScore[POS(startRow+1,startCol+1,nSeq1)] + CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getSubMatrixScore(NW_MAT_TYPE, Seq2[startRow].m_iCode, Seq1[startCol].m_iCode)))
		{
			ReSeq1.push_back(Seq1[startCol]);
			ReSeq2.push_back(Seq2[startRow]);
			++startRow;
			++startCol;
		}
		else if(startCol<nSeq1 && (arrScore[POS(startRow,startCol,nSeq1)] == arrScore[POS(startRow,startCol+1,nSeq1)] + CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(NW_MAT_TYPE)))
		{
			ReSeq1.push_back(Seq1[startCol]);
			ReSeq2.push_back(StruSeqElem(GENESPACE, -1));
			++startCol;
		}
		else if(startRow<nSeq2 && (arrScore[POS(startRow,startCol,nSeq1)] == arrScore[POS(startRow+1,startCol,nSeq1)] + CGlobalSpace::m_sAlignParams.m_SubstitutionMatMgr.getGapOpenCost(NW_MAT_TYPE)))
		{
			ReSeq1.push_back(StruSeqElem(GENESPACE, -1));
			ReSeq2.push_back(Seq2[startRow]);
			++startRow;
		}
	}

	CSequence  RS1( ReSeq1, vSequences[0].getName(), vSequences[0].getTitle() );
	CSequence  RS2( ReSeq2, vSequences[1].getName(), vSequences[1].getTitle() );
	vAlignedSequences.push_back(RS1);
	vAlignedSequences.push_back(RS2);

	SAFE_DELETE(arrScore);
}

}