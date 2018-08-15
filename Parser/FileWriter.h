#pragma once

#include <fstream>
#include"Sequence.h"

namespace SeqAnsis
{
	const int SEQ_HEAD_LEN =  10;

	class CFileWriter
	{
	public:
		CFileWriter( const std::string& vFilePathName );
		~CFileWriter(){}
		void    openFile();
		void    closeFile();
		bool		seqOutput(const std::vector<CSequence>& vSequences);

		template<typename T> void OutputVector(T* pVect, int size)
		{
			if (!m_FileOut.is_open())
			{
				char ch[128];
				sprintf_s(ch, "file %s is not open!", m_sequenceFileName.c_str());
				throw	CAppException(DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, std::string(ch));
			}

			for (int i = 0; i<size; ++i)
			{
				m_FileOut << pVect[i];
				m_FileOut << " ";
			}
			m_FileOut << std::endl;
		}
	private:
		std::string		    m_sequenceFileName;
		std::fstream	    m_FileOut;
	};

}
