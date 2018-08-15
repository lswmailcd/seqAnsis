#include"Common.h"
#include"GlobalSpace.h"
#include"FileReader.h"
#include"AppException.h"

namespace SeqAnsis
{

CFileReader::CFileReader( const std::string& vFilePathName ):m_sequenceFileName(vFilePathName)
{
    m_pFileIn = new CInFileStream();
    //structPenalties = 0; // I think this should be ok.
}

CFileReader::~CFileReader()
{
	SAFE_DELETE(m_pFileIn);
}


/* check if we've read sequences without any information, i.e. header
 *  but no sequence at all
 */
bool CFileReader::noEmptySequence(std::vector<CSequence> vSeqVector, std::string& vOffendingSeq)
{
    std::vector<CSequence>::iterator si;
    for (si = vSeqVector.begin(); si != vSeqVector.end(); si++) 
	{
        if (si->isEmpty())
		{
            vOffendingSeq.assign(si->getName());
            return false;
        }
    }
    return true;
}



/*
 * The function seqInput is called from the interactive menu only. This is because it
 * allows us to append seqs. But this would not happen on command line.
 * It is called twice in the interactive menu, both times with append as false. It calls
 * the readseqs function to do the work.
 */
 void CFileReader::seqInput(std::vector<CSequence>& vSequences, bool vAppend, std::string& vOffendingSeq)
{
	if (vAppend)
	{
		readSeqs(vSequences, vSequences.size()+1, vOffendingSeq);
	}
	else
	{
		readSeqs(vSequences, 1, vOffendingSeq);  //  1 is the first seq to be read 
	}
}

/*
 * The function readSeqs will create the FileParser that is needed and
 * it will read in all of the sequences. Depending on the filetype, it
 * will also check for secondary structure information.
 *
 * If there is a problem with one of the sequences its name will be
 * assigned to offendingSeq
 */
void CFileReader::readSeqs(std::vector<CSequence>& vSequences, int vFirstSeq, std::string& vOffendingSeq)
{
    // Now we have the file name, we must open the file.
    m_pFileIn->open(m_sequenceFileName.c_str());
    if (m_pFileIn->fail())
        throw CAppException( DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "Can NOT open sequence file!\n" );

    // Check if the file exists!
	if (!m_pFileIn->is_open())
    {
		throw CAppException( DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "Can NOT open sequence file!\n" );
    }
   
    int nSeqs = 0;
	std::auto_ptr<CFileParser> fileParser; // Means we dont need to delete it!
    checkInfile(nSeqs, fileParser);

    if (nSeqs == 0)
    {
        throw CAppException( DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "NO sequences in file!\n" );
    }

	std::vector<CSequence> SeqRangeVector = fileParser->getSeqRange(1, nSeqs, vOffendingSeq);
	std::vector<int> outputIndex; // Local version of outputIndex.
	for (int i=0; i < (int)(SeqRangeVector.size()); i++) 
	{
        // Andreas Wilm (UCD): fixed wrong default output order
        // which is important when no alignment (convert!) is done
        // _outputIndex.push_back(i); // default output order
        outputIndex.push_back(i+1); // default output order

        CSequence tempSeq = SeqRangeVector[i];

             
        vSequences.push_back(tempSeq);

    } 
  
    if (m_pFileIn->is_open())
    {
        m_pFileIn->close();
    }
}

/*
 * The function profileInput is used to read profiles into the Alignment. If the first
 * profile is already there it will read in the second profile. It returns the number
 * of seqs. If it returns -1, couldnt open file.
 */
void CFileReader::profileInput(std::vector<CSequence>& vSequences)
{
}


/*
 * The function checkInfile is used to find out which format the file is in, and it
 * also calls the appropriate function to count the seqences.
 */
void CFileReader::checkInfile(int& vSeqNum, std::auto_ptr<CFileParser>& fileParser)
{
    char lineIn[MAXLINE + 1];
    int i;
    int lengthLine = 0;
    vSeqNum = 0;

    while (m_pFileIn->getline(lineIn, MAXLINE + 1))
    {
		if (!CGlobalSpace::m_sUtility.blankLine(lineIn))
        {
            break;
        }
    }
    lengthLine = strlen(lineIn) - 1;
    
    for (i = lengthLine; i >= 0; i--)
    {
        if (isgraph(lineIn[i]))
        {
            break;
        }
    }
    lineIn[i + 1] = '\0';

    // Put first 7 characters into upper case! 
    for (i = 0; i <= 6 && i <= lengthLine; i++)
    {
        lineIn[i] = toupper(lineIn[i]);
    }

    // Create the parser to read the file.
	// Create the parser to read the file.
	if (CGlobalSpace::m_sUtility.lineType(lineIn, "ID"))
	{
		// EMBL/Swiss-Prot format ?
	}
	else if (CGlobalSpace::m_sUtility.lineType(lineIn, "CLUSTAL"))
	{
	}
	else if (CGlobalSpace::m_sUtility.lineType(lineIn, "PILEUP")) // MSF
	{
		fileParser.reset(new CMSFFileParser(m_sequenceFileName));
		CGlobalSpace::m_sEventLog.writeEvent( "Sequence format is MSF" );
	}
	else if (CGlobalSpace::m_sUtility.lineType(lineIn, "!!AA_MULTIPLE_ALIGNMENT")) // MSF
	{
		fileParser.reset(new CMSFFileParser(m_sequenceFileName));
		CGlobalSpace::m_sAlignParams.setDNAFlag(false);
		CGlobalSpace::m_sEventLog.writeEvent( "Sequence format is MSF" );
	}
	else if (CGlobalSpace::m_sUtility.lineType(lineIn, "!!NA_MULTIPLE_ALIGNMENT")) // MSF
	{
		fileParser.reset(new CMSFFileParser(m_sequenceFileName));
		//userParameters->setDNAFlag(true);
		CGlobalSpace::m_sEventLog.writeEvent( "Sequence format is MSF" );
	}
	else if (strstr(lineIn, "MSF") && lineIn[strlen(lineIn) - 1] == '.' &&
		lineIn[strlen(lineIn) - 2] == '.') // MSF
	{
		fileParser.reset(new CMSFFileParser(m_sequenceFileName));
		CGlobalSpace::m_sEventLog.writeEvent( "Sequence format is MSF" );
	}
	else if (CGlobalSpace::m_sUtility.lineType(lineIn, "!!RICH_SEQUENCE")) // RSF
	{
	}
	else if (CGlobalSpace::m_sUtility.lineType(lineIn, "#NEXUS"))
	{
		//utilityObject->error("Cannot read NEXUS format\n");
		return;
	}
	else if (*lineIn == '>')
	{
		if (*lineIn == '>')
		{
			if ((lengthLine>=3) && lineIn[3] == ';') 
			{
			}
			else
			{
				// FASTA/PERSON
				fileParser.reset(new CFastaFileParser(m_sequenceFileName));
				CGlobalSpace::m_sEventLog.writeEvent("Sequence format is FASTA/Pearson"); 
			}
		}
	}
    else
    {
        throw CAppException(DEF_EXCEPTION_INVALID_FILE, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, "The sequence format is NOT right in the file!" );
    }
    
    try
    {
        // Get the number of sequences. This is passed back as a pointer!
        vSeqNum = fileParser->countSeqs();
        // no output in 1.83: cout << "number of seqs is: " << *nseqs << "\n";
    }
    catch( const std::exception& e)
    {
        vSeqNum = 0;
		CGlobalSpace::m_sEventLog.writeStdException(e);
    }
}

}

