#pragma once

#include<string>
#include<vector>
#include"FileParser.h"
#include"Sequence.h"

namespace SeqAnsis
{

class CFastaFileParser : public CFileParser
{
    public:
        CFastaFileParser(const std::string& vFilePath);
        virtual std::vector<CSequence> getSeqRange(int vfirstSeq, int vnSeqToRead, std::string& vOffendingSeq);
        virtual CSequence getSeq(int vseqNum, std::string *vpOffendingSeq=NULL);
        virtual int countSeqs();
        virtual void getSecStructure(std::vector<char>& vGapPenaltyMask, 
                                     std::vector<char>& vSecStructMask, std::string& vSecStructName, 
                                     int &vStructPenalties, int vLength); 

    private:
        std::string m_FileName;
};

}

