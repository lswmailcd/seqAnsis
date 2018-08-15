#include"Common.h"
#include"GlobalSpace.h"
#include"SeqWeightFileParser.h"

const char LF = 0x0a;  //linefeed
const char CR = 0x0d;  //carriage return

namespace SeqAnsis
{

char CSeqWeightFileParser::getDelimiter(const std::string& vFileName)
{
	std::ifstream inFile;
	int type = 0;
	char delim;

	inFile.open(vFileName.c_str(), std::ios::in);
	inFile.seekg(0, std::ios::beg);

	//look for CR or LF or CRLF (or LFCR)
	if (inFile.is_open()) {
		char c;
		while (inFile.get(c)) {
			if (c == CR)
				type |= 1;
			else if (c == LF)
				type |= 2;
			else if (type)
				break;
		}
	}
	inFile.close();

	switch (type) {
	case 1:
		//cout << "file is Mac System 9" << endl;
		delim = '\r';
		break;
	case 2:
		//cout << "file is UNIX" << endl;
		delim = '\n';
		break;
	case 3:
		//cout << "file is DOS" << endl;
		delim = '\n';
		break;
	default: //short or empty file
		//cout << "file is UNIX (default)" << endl;
		delim = '\n';
	}
	return delim;
}

void	CSeqWeightFileParser::getSeqWeightMatrix(const std::string& vFileName)
{
	char delim = getDelimiter(m_FileName);
	std::ifstream fileIn;
	fileIn.open(m_FileName.c_str(), std::ios::in);

	std::string line = "";

	char DIST[] = "DIST";

	do
	{
		std::getline(fileIn, line, delim);
		if (line.substr(16, 4) == DIST)
		{
			float fWeight = atof(line.substr(21, 6));
			CG
		}
	} while (currentSeqNum <vseqNum);
}

int		CSeqWeightFileParser::getSeqNum(const std::string& vFileName)
{

}

}