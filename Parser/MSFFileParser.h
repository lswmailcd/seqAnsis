/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef MSFFILEPARSER_H
#define MSFFILEPARSER_H

#include <string>
#include "FileParser.h"

namespace SeqAnsis
{

class CMSFFileParser : public CFileParser
{
    public:
        /* Functions */
        CMSFFileParser(const std::string& filePath);
		virtual std::vector<CSequence> getSeqRange(int vfirstSeq, int vnSeqToRead, std::string& vOffendingSeq);
		virtual CSequence getSeq(int vseqNum, std::string *vpOffendingSeq=NULL);
		virtual int countSeqs();
		virtual void getSecStructure(std::vector<char>& vGapPenaltyMask, 
														std::vector<char>& vSecStructMask, std::string& vSecStructName, 
														int &vStructPenalties, int vLength); 

        /* Attributes */

    private:
        /* Functions */

        /* Attributes */
        std::string fileName;
};

}
#endif


