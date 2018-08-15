#include"Common.h"
#include"GlobalSpace.h"
#include"FastaFileParser.h"
#include"AppException.h"

namespace SeqAnsis
{

CFastaFileParser::CFastaFileParser(const std::string& vFilePath)
{
    m_FileName = vFilePath;
    fillCharTab(); 
}


/**
 * reads fasta/pearson file in one go instead of calling getSeq for
 * each single sequence.
 *
 * FIXME AW: only CFastaFileParser::getSeqRange is special, rest is the
 * same. should be defined in FileParser and then overloaded in special
 * cases like here
 */
std::vector<CSequence> CFastaFileParser::getSeqRange(int vfirstSeq, int vnSeqToRead, std::string& vOffendingSeq)
{
    std::vector<StruSeqElem> characterSeq;
    std::string name = "";
    std::string title = "";
    std::string blank = "";
    std::string greater = ">";
    //_line[0] = EOS;
    std::vector<CSequence> seqRangeVector;
    
    std::string line;
    
    
    //int i, j;
    int nSeqsRead = 0;
    unsigned char c;
    char delim;
    int currentSeqNum = 0; // Not at any sequence yet!
    
    try
    {
       delim=CFileParser::getDelimiter(m_FileName);
       //cout << "delim = " << delim << endl;
       std::ifstream fileIn;
       fileIn.open(m_FileName.c_str(),std::ios::in);

        // Read in lines until we get to the begining of sequence firstSeq.
        std::string line="";

        do 
		{
			std::getline(fileIn,line,delim);
			if(line.substr(0,1) == greater)
			{
				currentSeqNum++;
			 }
        }while(currentSeqNum <vfirstSeq);
        
        
        while (nSeqsRead < vnSeqToRead)
        {
            // get sequence name from current line (excluded '>' and read up to first ' ' or MAXNAMES
            // remove the first char i.e. '>'
            name=line.substr(1,MAXNAMES);
            //if(name.find(">") != string::npos){
            //  andreas wilm: exit if angle bracket within header?
            //}
            
            while(name.substr(0,1)==" "){
                name=name.substr(1,MAXNAMES);
            }
            //int i;
            //i = name.find(" ");
            if(name.find(" ") != std::string::npos){
                name=name.substr(0,name.find(" "));
            }
            CGlobalSpace::m_sUtility.rTrim(&name); // also replaces linef

            name=CGlobalSpace::m_sUtility.blankToUnderscore(name); // replace blanks with '_'
            
            
            // Read in lines until we get to the begining of next sequence.
            
            title = ""; // No title information
            int iNo=0;
            while(std::getline(fileIn,line,delim) )
			{    
               std::string::const_iterator iterator1 = line.begin();
                while(iterator1 != line.end())
				{
                    // Andreas Wilm (UCD): exit if angle brackets within sequence
                    if(*iterator1=='>' && iterator1!=line.begin()) 
					{                        
                        fileIn.close();
                        seqRangeVector.clear();
						throw CAppException( DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "Bad file format!" );
                        //return seqRangeVector;
                    }
                       
                    if(*iterator1 =='\n' || *iterator1 =='\r' || *iterator1 == '\0' || *iterator1 =='>')
					{
                        break;
                    }
                    c = *iterator1;

                    c = chartab[c];
                    if(c)
					{
						StruSeqElem	elem(c, iNo);
                        characterSeq.push_back( elem );
						++iNo;
                    }
                    iterator1++;
                }
                if(iterator1 != line.end() && *iterator1 == '>')
				{
                    break;
                } 
            }
            
            // check sequence
			if ((int)(characterSeq.size()) > CGlobalSpace::m_sAlignParams.getMaxAllowedSeqLength())
            {
                /* error output*/
                //parseExitCode=SEQUENCETOOBIG;
				if (vOffendingSeq.empty())
                    vOffendingSeq.assign(name);
                fileIn.close();
                seqRangeVector.clear();
				throw CAppException( DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "Sequence length exceed max length!" );
                //return seqRangeVector;
            }
            else if ((int)(characterSeq.size()) == 0)
            {
                //parseExitCode=EMPTYSEQUENCE;
                if (vOffendingSeq.empty())
                    vOffendingSeq.assign(name);
                fileIn.close();
                seqRangeVector.clear();
				throw CAppException( DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "Sequence is empty!" );
                //return seqRangeVector;
            }

            seqRangeVector.push_back(CSequence(characterSeq, name, title ));
            characterSeq.clear();
            nSeqsRead++;
        } // while (nSeqsRead < nSeqsToRead)

        fileIn.close();

        return seqRangeVector;
    }

    catch(...)
    {
        throw CAppException( DEF_EXCEPTION_UNEXPECTED, DEF_EXCEPTION_LEVEL_EXIT_APP, __EXCEPTION_SITE__, "Unexcepted ERROR, must end program!");
    }
}



/**
 * The function getSeq is used to get the sequence 'seqNum' in the file. It returns a
 * sequence object containing the sequence.
 * Deprecated: where possible use faster getSeqRange which reads
 * sequences in one go
 * @param seqNum The number of the sequence to get.
 * @return 
 */
CSequence CFastaFileParser::getSeq(int vseqNum, std::string *vpOffendingSeq)
{
    std::vector<StruSeqElem>  characterSeq;
    std::string name = "";
    std::string title = "";
    std::string blank = "";
    std::string greater = ">";
     
	CGlobalSpace::m_sEventLog.writeEvent("Use of CFastaFileParser::getSeq is deprecated!");
    //int i, j;
    unsigned char c;
    char delim;
    int currentSeqNum = 0; // Not at any sequence yet!
    
    delim=CFileParser::getDelimiter(m_FileName);
    std::ifstream fileIn;
    fileIn.open(m_FileName.c_str(),std::ios::in);

    //////////////////////////////////////////////////    
    //PMcG replace char array with string processing    
    //////////////////////////////////////////////////    

    // Read in lines until we get to the begining of sequence seqNum.    
    std::string line="";

    do 
	{
		std::getline(fileIn,line,delim);
        if(line.substr(0,1) == greater)
		{
			currentSeqNum++;
        }
    }while(currentSeqNum <vseqNum);
        
        
        // get sequence name from current line (excluded '>' and read up to first ' ' or MAXNAMES
        // remove the first char i.e. '>'
     name=line.substr(1,MAXNAMES);
           
     while(name.substr(0,1)==" ")
	 {
          name=name.substr(1,MAXNAMES);
     }

     if(name.find(" ") != std::string::npos)
	 {
          name=name.substr(0,name.find(" "));
     } 
	 name=CGlobalSpace::m_sUtility.blankToUnderscore(name); // replace blanks with '_'

        
     title = ""; // No title information
	 std::string seqLine = "";
	 int iNo=0;
     while(std::getline(fileIn,seqLine,delim) )
	 {
		 std::string::const_iterator iterator1 = seqLine.begin();
         while(iterator1 != seqLine.end())
		 {
			 if(*iterator1 =='\n' || *iterator1 =='\r' || *iterator1 == '\0' || *iterator1 =='>')
			 {
				  break;
			 }
		     c = *iterator1;
             c = chartab[c];
             if(c)
			 {
				 StruSeqElem elem(c, iNo);
				 characterSeq.push_back(elem);
				 ++iNo;
			 }
             iterator1++;
          }
          if(*iterator1 == '>')
		  {
			  break;
          }
     }

     fileIn.close();

	 if ( (int)(characterSeq.size()) > CGlobalSpace::m_sAlignParams.getMaxAllowedSeqLength() )
     {
		 throw CAppException( DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "Sequence length exceed max length!" );
            // return empty seq
            //return CSequence(blank, blank, blank);
     }
     else if ((int)(characterSeq.size()) == 0)
     {
            //parseExitCode=EMPTYSEQUENCE;
            // return empty seq
            //return CSequence(blank, blank, blank);
			throw CAppException( DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "Sequence is empty!" );
      }
       
      return CSequence(characterSeq, name, title);
}

/**
 * The function countSeqs, counts the number of sequences in a file.
 * @return The number of sequences in the file.
 */
int CFastaFileParser::countSeqs()
{
    //char line[1000 + 1];
    int nseqs = 0;
    std::string line2;
    char delim;

    delim=CFileParser::getDelimiter(m_FileName);    
    std::ifstream fileIn;    
	fileIn.open(m_FileName.c_str(),std::ios::in);	
    
    if(!fileIn.is_open())    
    {    
        return 0; // No sequences found!    
    }    
    
    while (std::getline(fileIn,line2,delim))     
	{                 	
        if (line2[0] == '>')    
        {    
            nseqs++;    
        }    
    }    
    fileIn.close();    
    return nseqs;    
}

/**
 * There is no secondary structure information in the Pearson file. This is here to 
 * set the structPenalties to NONE.
 * @param gapPenaltyMask 
 * @param secStructMask 
 * @param secStructName 
 * @param structPenalties 
 * @param length 
 */
void CFastaFileParser::getSecStructure(std::vector<char>& gapPenaltyMask, 
                         std::vector<char>& secStructMask, std::string& secStructName, 
                          int &structPenalties, int length)
{
    structPenalties = NONE;
}

}


