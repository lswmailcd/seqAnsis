#pragma once

#include<ctype.h>
#include<iostream>
#include<vector>
#include"InFileStream.h"

namespace SeqAnsis
{

	class CSeqWeightFileParser
	{
	public:
		CSeqWeightFileParser(const std::string& vFilePath){ m_FileName = vFilePath; }
		~CSeqWeightFileParser(){}

		void	getSeqWeightMatrix(const std::string& vFileName);
		int		getSeqNum(const std::string& vFileName);
	protected:
		char getDelimiter(const std::string& vFileName);
	private:
		std::string m_FileName;
		CInFileStream*  m_pFileIn;
	};

}