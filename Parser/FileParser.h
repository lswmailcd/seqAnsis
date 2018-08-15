#pragma once

#include<ctype.h>
#include<iostream>
#include<vector>
#include"Sequence.h"
#include"InFileStream.h"

namespace SeqAnsis
{
 
class CFileParser
{
    public:
        /* Functions */
        CFileParser();
        virtual ~CFileParser();
        virtual std::vector<CSequence> getSeqRange(int vFirstSeq, int vNum, std::string& vOffendingSeq) = 0;
        virtual CSequence getSeq(int vSeqNum, std::string *vOffendingSeq=NULL) = 0;
        virtual int countSeqs() = 0; // VIRTUAL 
        virtual void getSecStructure(std::vector<char>& vGapPenaltyMask, 
                                     std::vector<char>& vSecStructMask,
                                     std::string& vSecStructName, int &vStructPenalties, int vLength) = 0;
        void fillCharTab(void);
        char getDelimiter(std::string vFileName);
        /* Attributes */
        char chartab[128];
        //int getParseExitCode() { return parseExitCode; };
        
    protected:
        void freeFileResources(CInFileStream* vpFilePtr);
        CInFileStream*  m_pFileIn;
        int parseExitCode; // reason for returning empty sequence
                           // vector; same as used in FileReader
};

}


