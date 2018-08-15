#include <string>
#include <fstream>
#include <iostream>
#include "InFileStream.h"
using namespace std;

const char LF = 0x0a;  //linefeed
const char CR = 0x0d;  //carriage return

namespace SeqAnsis
{

CInFileStream::CInFileStream() : ifstream()
{
    m_Delim = '\n'; // default
    //cout << "CInFileStream() constructor 1" << endl;
}

CInFileStream::CInFileStream(const std::string& vFileName) : ifstream(vFileName, ios::in), m_FileName(vFileName)
{
    m_Delim = findDelimiter();
}

//- copy-constructor: can't copy superclass private members
//- CInFileStream::CInFileStream(const CInFileStream &copy) :
//-     ifstream(static_cast<const ifstream&>(copy))
//- {
//-     cout << "CInFileStream() constructor 3" << endl;
//-     delim = copy.delim;
//- }

void CInFileStream::open(const char *vFilename) 
{

    this->m_FileName = vFilename;
    ifstream::open(vFilename, ios::in);
    if  (ifstream::fail())
        return;
    m_Delim = findDelimiter();
}

//not necessary, but for symmetry to open()
void CInFileStream::close() 
{
    ifstream::close();   
}

//getline with stored delimiter
std::istream& CInFileStream::getline(char *vStr, std::streamsize vNum) 
{
    return ifstream::getline(vStr, vNum, m_Delim);
}

//getline with caller supplied delimiter
std::istream& CInFileStream::getline(char *vStr, std::streamsize vNum, char vDelim) 
{
    return ifstream::getline(vStr, vNum, vDelim);
}


/**
 * Mark 24-1-2007. I added the function findDelimiter to determine if '\r' or
 * '\n' will be used as the line delimiter when parsing the file.
 *
 * 25-01-07,Nigel Brown(EMBL): changed body of loop to check successive chars
 * in case of DOS/Windows
 *
 * 09-02-07,Nigel Brown(EMBL): moved member into new InFileStream subclassed
 * from std::ifstream, so this is called automatically for any file reader
 * that uses InFileStream in place of std::ifstream. Replaced if/then/else
 * with switch.
 */
char CInFileStream::findDelimiter()
{
    ifstream in;
    int type = 0;
    
    in.open(m_FileName.c_str(), ios::in);
    if (in.fail())
        return m_Delim;
    
    in.seekg(0, ios::beg);

    //look for CR or LF or CRLF (or LFCR)
    if (in.is_open()) {
        char c;
        while (in.get(c)) {
            if (c == CR)
                type |= 1;
            else if (c == LF)
                type |= 2;
            else if (type)
                break;
        }
    }
    in.close();

    switch (type) {
	case 1:
	    //cout << "file is Mac System 9" << endl;
	    m_Delim = '\r';
	    break;
	case 2:
	    //cout << "file is UNIX" << endl;
	    m_Delim = '\n';
	    break;
	case 3:
	    //cout << "file is DOS" << endl;
	    m_Delim = '\n';
	    break;
	default: //short or empty file
	    //cout << "file is UNIX (default)" << endl;
	    m_Delim = '\n';
    }
    return m_Delim;
}

}
