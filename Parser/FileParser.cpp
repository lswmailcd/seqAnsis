#include"Common.h"
#include"GlobalSpace.h"
#include"FileParser.h"

const char LF = 0x0a;  //linefeed
const char CR = 0x0d;  //carriage return

namespace SeqAnsis
{

CFileParser::CFileParser()
{
    //parseExitCode = OK;
}

CFileParser::~CFileParser()
{
    // Dont do anything!!!!
}


void CFileParser::fillCharTab(void)
{
    register int i;
    register char c;

    for (i = 0; i < 128; chartab[i++] = 0)
        ;
	for (i = 0; i <= CGlobalSpace::m_sAlignParams.getAminoAcidCodesNum() + 1; i++)
    {
		c = CGlobalSpace::m_sAlignParams.getAminoAcidCodes().at(i);
        chartab[(int)c] = chartab[tolower(c)] = c;
    }
}

void CFileParser::freeFileResources(CInFileStream* vpFilePtr)
{
    if(vpFilePtr != 0)
    {
        vpFilePtr->close();
        SAFE_DELETE(vpFilePtr); 
    }
}


char CFileParser::getDelimiter(std::string vFilename)
{
    std::ifstream inFile;
    int type = 0;
    char delim;

    inFile.open(vFilename.c_str(), std::ios::in);
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

}
