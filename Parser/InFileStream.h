#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <memory>

namespace SeqAnsis
{

class CInFileStream : public std::ifstream
{
  public:
    CInFileStream();
    CInFileStream(const std::string& vFileName);

    void open(const char *vFileName);
    void close();

    std::istream& getline(char *vStr, std::streamsize vNum);
    std::istream& getline(char *vStr, std::streamsize vNum, char vDelim);

  protected:
    char findDelimiter();

  private:
    //disable copy-constructor
    CInFileStream(const CInFileStream &vCopy);
    std::string m_FileName;
    char m_Delim;
};

}



