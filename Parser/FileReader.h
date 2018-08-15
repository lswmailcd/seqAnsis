#pragma once

#include<vector>
#include<string>
#include<memory>
#include"FileParser.h"
#include"FastaFileParser.h"
#include "MSFFileParser.h"
#include"Sequence.h"

namespace SeqAnsis
{

class CFileReader
{
    public:
        CFileReader( const std::string& vFilePathName );
        ~CFileReader();
        void	seqInput(std::vector<CSequence>& vSequences, bool vAppend, std::string& vOffendingSeq);
        void	readSeqs(std::vector<CSequence>& vSequences, int vFirstSeq, std::string& vOffendingSeq);        
		void	profileInput(std::vector<CSequence>& vSequences);
    private:
        void	checkInfile(int& vSeqNum, std::auto_ptr<CFileParser>& vFileParser);
        bool	noEmptySequence(std::vector<CSequence> vSeqVector, std::string& vOffendingSeq);
            
		std::string		m_sequenceFileName;
        CInFileStream*	m_pFileIn;
        //int				structPenalties;
        //std::string		secStructName;
        //std::vector<char>	secStructMask; // Will need to be cleared out before every reading!
        //std::vector<char>	gapPenaltyMask;
        //std::vector<std::string>	formatNames; 
};

}
