#include "FileWriter.h"
#include "AppException.h"
#include "GlobalSpace.h"

namespace SeqAnsis
{

CFileWriter::CFileWriter( const std::string& vFilePathName ):m_sequenceFileName(vFilePathName)
{

}

void    CFileWriter::openFile() 
{
	if( !m_FileOut.is_open() )
	{
		try
		{
			m_FileOut.open( m_sequenceFileName.c_str(), std::ios::out|std::ios::trunc );
		}
		catch (CAppException const& e)
		{
			CGlobalSpace::m_sEventLog.writeException(e);
		}
	}
}

void    CFileWriter::closeFile()
{
	m_FileOut.close();
}

bool		CFileWriter::seqOutput(const std::vector<CSequence>& vSequences)
{
	if( !m_FileOut.is_open() )
	{
		char ch[128];
		sprintf_s( ch,  "file %s is not open!", m_sequenceFileName.c_str()  );
		throw	CAppException( DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, std::string(ch) );
	}

	for ( int i=0; i<(int)vSequences.size(); ++i )
	{
		m_FileOut<<vSequences[i].getName();
		char ch[SEQ_HEAD_LEN];
		for ( int j=0; j<SEQ_HEAD_LEN; ++j)
		{
			sprintf_s( ch,  " ");
		}
		m_FileOut<<ch;
		for ( int j=0; j<vSequences[i].getLen(); ++j)
		{
			m_FileOut<<vSequences[i].getSequenceContext()[j].m_char;
		}
		m_FileOut<<std::endl;
	}
	m_FileOut<<std::endl;//add one blank line as ending for the balibase score program
	return true;
}

}