#include "Common.h"
#include "MSFFileParser.h"

namespace SeqAnsis
{

/**
 * MSFFileParser contructor sets up the chartab array.
 * @param filePath 
 * @return 
 */
CMSFFileParser::CMSFFileParser(const std::string& filePath)
{
    fileName = filePath; 
    fillCharTab();
}


    
 std::vector<CSequence> CMSFFileParser::getSeqRange(int vfirstSeq, int vnSeqToRead, std::string& vOffendingSeq)
{
    std::vector<CSequence> seqRangeVector;
    int i;

    for (i=0; i<vnSeqToRead; i++)
    { 
        CSequence tempSeq = getSeq(vfirstSeq + i);
 #if 0
       if (parseExitCode!=OK) {
            seqRangeVector.clear();
            return seqRangeVector;
        }
#endif
        seqRangeVector.push_back(tempSeq);
    }
    return seqRangeVector;
}



/**
 * The function getSeq finds the sequence seqNum in the file and returns it.
 * @param seqNum The number of the sequence in the file to get.
 * @return A sequence object containing the seqNum'th sequence from the file.
 */
CSequence CMSFFileParser::getSeq(int seqNum, std::string *offendingSeq)
{
    char _line[MAXLINE + 1];
    char _sname[MAXNAMES + 1];
    CSeqData  characterSeq;
    std::string name = "";
    std::string title = "";
    std::string blank = "";
    
    _line[0] = EOS;
    int i, j, k;
    unsigned char c;

    try
    {
        m_pFileIn = new CInFileStream;  //nige
        m_pFileIn->open(fileName.c_str());  //nige
        m_pFileIn->seekg(0, std::ios::beg);       
        
        for (i = 0;; i++)
        {
            if (!m_pFileIn->getline(_line, MAXLINE + 1))
            {
                m_pFileIn->close();
				CSeqData seq;
                return CSequence(seq, blank, blank);
            }
            // read the title
            if (CGlobalSpace::m_sUtility.lineType(_line, "//"))
            {
                break;
            }
            // lines...ignore
        }
		int idx = 0;
        while (m_pFileIn->getline(_line, MAXLINE + 1))
        {
            if (!CGlobalSpace::m_sUtility.blankLine(_line))
            {
                for (i = 1; i < seqNum; i++)
                {
                    m_pFileIn->getline(_line, MAXLINE + 1);
                }
                for (j = 0; j <= (int)strlen(_line); j++)
                {
                    if (_line[j] != ' ')
                    {
                        break;
                    }
                }
                for (k = j; k <= (int)strlen(_line); k++)
                {
                    if (_line[k] == ' ')
                    {
                        break;
                    }
                }
                
                // Get the name of the sequence
                strncpy(_sname, _line + j, CGlobalSpace::m_sUtility.MIN(MAXNAMES, k - j));
                _sname[CGlobalSpace::m_sUtility.MIN(MAXNAMES, k - j)] = EOS;
                CGlobalSpace::m_sUtility.rTrim(_sname);
                CGlobalSpace::m_sUtility.blankToUnderscore(_sname);
                name = std::string(_sname);

                for (i = k; i <= MAXLINE; i++)
                {
                    c = _line[i];
                    if (c == '.' || c == '~')
                    {
                        c = '-';
                    }
                    if (c == '*')
                    {
                        c = 'X';
                    }
                    if (c == '\n' || c == EOS)
                    {
                        break;
                    }
                    // EOL 
                    c = chartab[c];
                    if (c)
                    {
						StruSeqElem   elem( c, idx );
                        characterSeq.push_back(elem);
						++idx;
                    }
                }

                for (i = 0;; i++)
                {
                    if (!m_pFileIn->getline(_line, MAXLINE + 1))
                    {
                        m_pFileIn->close();
                        return CSequence(characterSeq, name, title);
                    }
                    if (CGlobalSpace::m_sUtility.blankLine(_line))
                    {
                        break;
                    }
                }
            }
        }
        m_pFileIn->close();
        
        if ((int)characterSeq.size() > CGlobalSpace::m_sAlignParams.getMaxAllowedSeqLength())
        {
            parseExitCode=SEQUENCETOOBIG;
            if (offendingSeq!=NULL)
                offendingSeq->assign(name);
            // return empty seq
			CSeqData seq;
            return CSequence(seq, blank, blank);
        }
        return CSequence(characterSeq, name, title);
    }
    catch(...)
    {
        m_pFileIn->close();
		CGlobalSpace::m_sEventLog.writeEvent( "An exception has occured in the function MSFFileParser::getSeq(), Program needs to terminate." );
        exit(1);
    }
}

/**
 * The function countSeqs counts the number of sequences in the file.
 * @return The number of sequences in the file.
 */
int		CMSFFileParser::countSeqs()
{
    char _line[MAXLINE + 1];
    int _numSeqs;
    
    try
    {
        m_pFileIn = new CInFileStream;  //nige
        m_pFileIn->open(fileName.c_str());  //nige
    
        if(!m_pFileIn->is_open())
        {
            return 0; // No sequences found!
        }
    
        while (m_pFileIn->getline(_line, MAXLINE + 1))
        {
            if (CGlobalSpace::m_sUtility.lineType(_line, "//"))
            {
                break;
            }
        }

        while (m_pFileIn->getline(_line, MAXLINE + 1))
        {
            if (!CGlobalSpace::m_sUtility.blankLine(_line))
            {
                break;
            }
            // Look for next non- blank line
        } 
        _numSeqs = 1;

        while (m_pFileIn->getline(_line, MAXLINE + 1))
        {
            if (CGlobalSpace::m_sUtility.blankLineNumericLabel(_line))
            {
                m_pFileIn->close();
                return _numSeqs;
            }
            _numSeqs++;
        }

        return 0; // if you got to here-funny format/no seqs.
    }
    catch(...)
    {
        m_pFileIn->close();
		CGlobalSpace::m_sEventLog.writeEvent( "An exception has occured in the function MSFFileParser::countSeqs(), Program needs to terminate." );
        exit(1);    
    }
}

/**
 * There is no secondary structure information in MSF files. Set structPenalties to NONE.
 * @param gapPenaltyMask 
 * @param secStructMask 
 * @param secStructName 
 * @param structPenalties 
 * @param length 
 */
void		CMSFFileParser::getSecStructure(std::vector<char>& gapPenaltyMask, std::vector<char>& secStructMask,
															      std::string& secStructName, int &structPenalties, int length)
{
    structPenalties = NONE;
}

}

