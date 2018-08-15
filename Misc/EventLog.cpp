#pragma once

#include <sys/timeb.h>
#include <time.h>
#include <io.h>
#include<OpenThreads\mutex>
#include "EventLog.h"

namespace SeqAnsis
{

#define DEF_EVENT_BUF_DUMP_SIZE   10000  
#define DEF_NUM_EVENT_KEEP_IN_BUF 20

static FILE              *s_LogFile;
static FILE				*errorRW_LogFile;
static OpenThreads::Mutex s_OutputMutex;
static std::string        s_EventBuf;
static std::string        s_LockOwner;
static char               s_LogFileName[15] = "Events.log";
static char               errorRW_LogFileName[15] = "ErrorRW.log";

const long DEF_MAX_LOG_FILE_SIZE = 67108864; //64MB      2^32-1 是int在32位机器上面表示的最大值

CEventLog::CEventLog(void) : m_numEventInBuffer(0)
{
}

CEventLog::~CEventLog(void)
{
}

/****************************************************************************************************************************************************/
//FUNCTION: show the exception and write the exception information to event log file
void CEventLog::writeException(const CAppException& vException)
{
	char EventLog[500];
	sprintf_s(EventLog, 500, "The application throws an exception:\n    Description: %s\n    Source: %s\n    File: %s(%d)", 
		vException.m_Description.c_str(), vException.m_ExpFunc.c_str(), vException.m_ExpFile.c_str(), vException.m_ExpLine);

	s_OutputMutex.lock();
	s_LockOwner = "CEventLog::writeException";

	outputEvent(EventLog, DEF_OUTPUT_ERROR, false, true);
	//flushBufferEvents();

	s_LockOwner = "";
	s_OutputMutex.unlock();
}

/****************************************************************************************************************************************************/
//FUNCTION: display a warning message and also write the message to event log file
void CEventLog::writeWarning(int vLine, const std::string& vSource, const std::string& vFile, const std::string& vDescription)
{
	char EventLog[500];
	sprintf_s(EventLog, 500, "The application throws a warning:\n    Description: %s\n    Source: %s\n    File: %s(%d)", vDescription.c_str(), vSource.c_str(), vFile.c_str(), vLine);

	s_OutputMutex.lock();
	s_LockOwner = "CEventLog::writeWarning";

	outputEvent(EventLog, DEF_OUTPUT_WARNING, false, true);
	//flushBufferEvents();
	s_LockOwner = "";
	s_OutputMutex.unlock();
}

/****************************************************************************************************************************************************/
//FUNCTION: display a warning message and also write the message to event log file
//NOTE:     this function uses the callback function m_pMsgDisplayFunc to display message
void CEventLog::writeWarning(void *vOutputTarget, const std::string& vSource, const std::string& vFile, int vLine, const std::string& vDescription)
{
	if (vOutputTarget == NULL) 
	{
		writeWarning(vLine, vSource, vFile, vDescription);
	}
}

//*********************************************************************************************************************************************************************
//FUNCTION: 
void CEventLog::writeMessage(const std::string& vDescription, bool vWriteToEvent/* =true */)
{
	s_OutputMutex.lock();
	s_LockOwner = "CEventLog::displayMessage";

	if (vWriteToEvent)
	{
		outputEvent(vDescription.c_str(), DEF_OUTPUT_SUCCESS);
	}
	s_LockOwner = "";
	s_OutputMutex.unlock();
}

/****************************************************************************************************************************************************/
//FUNCTION: display a message and write it to event log file
//NOTE:     this function uses the callback function m_pMsgDisplayFunc to display message
void CEventLog::writeMessage(void *vOutputTarget, const std::string& vDescription, bool vWriteToEvent)
{
	if (vOutputTarget != NULL)
	{
		if (vWriteToEvent)
		{
			writeEvent(vDescription);
		}
	}
	else
	{
		writeMessage(vDescription, vWriteToEvent);
	}
}

//*********************************************************************************************************************************************************************
//FUNCTION: 
void CEventLog::writeStdException(const std::exception& vException)
{
	char EventLog[500];
	sprintf_s(EventLog, 500, "\nstd exception:   Description: %s\n    Type: %s\n", 
		      vException.what(), typeid(vException).name());

	s_OutputMutex.lock();
	s_LockOwner = "CEventLog::writeStdException";

	outputEvent(EventLog, DEF_OUTPUT_ERROR);
	//flushBufferEvents();

	s_LockOwner = "";
	s_OutputMutex.unlock();
}

/****************************************************************************************************************************************************/
//FUNCTION: create a string containing the current time
void CEventLog::assembleTimeString(char *voTimeString)
{
	struct _timeb timebuffer;
	char timeline[26];

	_ftime64_s( &timebuffer );

	ctime_s( timeline, 26, & ( timebuffer.time ) );
	for (int i=25; i>=0; i--)
	{
		if (timeline[i] == '\n')
		{
			timeline[i] = 0;
			break;
		}
	}
	sprintf_s(voTimeString, 128, "%.19s.%hu %s", timeline, timebuffer.millitm, &timeline[20] );
}

//*********************************************************************************************************************************************************************
//FUNCTION: 
void CEventLog::outputEvent(const char *vEvent, OutputType vEventType, bool vOutputToScreen/* =false */, bool vImmediateWrite/* =false */)
{
	assembleTimeString(m_TimeString);
	std::string Event;
	Event = std::string("EVENT at ") + m_TimeString + std::string("\n ") + std::string(vEvent);

	s_EventBuf += std::string("\n\n") + Event;
	m_numEventInBuffer++;
	if ((s_EventBuf.size() > DEF_EVENT_BUF_DUMP_SIZE) || vImmediateWrite || m_numEventInBuffer > DEF_NUM_EVENT_KEEP_IN_BUF)
	{
		bool bAppendOpen=true;
		if ( fopen_s(&s_LogFile, s_LogFileName, "r") == 0)
		{
			long Len = _filelength(_fileno(s_LogFile));
			bAppendOpen = (Len < DEF_MAX_LOG_FILE_SIZE);
			fclose(s_LogFile);
			if (bAppendOpen)
				fopen_s(&s_LogFile, s_LogFileName, "a+");
			else
				fopen_s(&s_LogFile, s_LogFileName, "w");
	    	fprintf(s_LogFile, "%s", s_EventBuf.c_str());
			fclose(s_LogFile);
			s_EventBuf.clear();
		}
		else
		{
			if (fopen_s(&errorRW_LogFile, errorRW_LogFileName, "r") != 0)
			{//for the first time to run application, the log file does not exist. thus we directly create the event log file
				fopen_s(&errorRW_LogFile, errorRW_LogFileName, "w");
				fprintf(errorRW_LogFile, "%s/n","Fail to open the log file for writing.");
				fprintf(errorRW_LogFile, "%s", s_EventBuf.c_str());
				fclose(errorRW_LogFile);
			}
			else
			{//check the size of the existed log file, clear the content of log file if its size is too big
				long Len = _filelength(_fileno(errorRW_LogFile));
				bAppendOpen = (Len < DEF_MAX_LOG_FILE_SIZE);
				fclose(errorRW_LogFile);

				if (bAppendOpen)
					fopen_s(&errorRW_LogFile, errorRW_LogFileName, "a+");
				else
					fopen_s(&errorRW_LogFile, errorRW_LogFileName, "w");

				fprintf(errorRW_LogFile, "%s/n","Fail to open the log file for writing.");
				fprintf(errorRW_LogFile, "%s", s_EventBuf.c_str());
				fclose(errorRW_LogFile);
				s_EventBuf.clear();
			}

			if(remove(s_LogFileName) == -1)
				std::cout << "Fail to Delete" << s_LogFileName <<std::endl;;
			std::cout << "Fail to open the log file for writing." << std::endl << std::endl;
		}
		m_numEventInBuffer = 0;
		s_EventBuf.clear();
	}

	if (vOutputToScreen)
	{
		std::cout << vEvent << std::endl;
	}
}

/****************************************************************************************************************************************************/
//FUNCTION: write the specified event to the event log file
void CEventLog::writeEvent(const char *vEvent, bool vImmediateWrite/* =true */, bool vOutputToScreen/* =false */)
{
	s_OutputMutex.lock();
	s_LockOwner = "CEventLog::writeEvent";
	outputEvent(vEvent, DEF_OUTPUT_SUCCESS, vOutputToScreen, vImmediateWrite);
	s_LockOwner = "";
	s_OutputMutex.unlock();
}

void CEventLog::writeEvent(const std::string& vEvent, bool vImmediateWrite/* =true */, bool vOutputToScreen/* =false */)
{
	s_OutputMutex.lock();
	s_LockOwner = "CEventLog::writeEvent";
	outputEvent(vEvent.c_str(), DEF_OUTPUT_SUCCESS, vOutputToScreen, vImmediateWrite);
	s_LockOwner = "";
	s_OutputMutex.unlock();
}

/****************************************************************************************************************************************************/
//FUNCTION: write the application start log
void CEventLog::writeAppStartEvent()
{
	bool bAppendOpen=true;

	if (fopen_s(&s_LogFile, s_LogFileName, "r") != 0)
	{//for the first time to run application, the log file does not exist. thus we directly create the event log file
		fopen_s(&s_LogFile, s_LogFileName, "w");
	}
	else
	{//check the size of the existed log file, clear the content of log file if its size is too big
		long Len = _filelength(_fileno(s_LogFile));
		bAppendOpen = (Len < DEF_MAX_LOG_FILE_SIZE);
		fclose(s_LogFile);

		if (bAppendOpen)
			fopen_s(&s_LogFile, s_LogFileName, "a+");
		else
			fopen_s(&s_LogFile, s_LogFileName, "w");

	}

	assembleTimeString(m_TimeString);
	fprintf(s_LogFile, "\n\n/*************************************************************************************************************/\n", m_TimeString);
	fprintf(s_LogFile, "/***********************************Application starts at %s***********************************/\n", m_TimeString);
	fclose(s_LogFile);
}

/****************************************************************************************************************************************************/
//FUNCTION: write the application end log
void CEventLog::writeAppEndEvent()
{
	if (fopen_s(&s_LogFile, s_LogFileName, "a+") == 0)
	{
		if (!s_EventBuf.empty())
		{
			fprintf(s_LogFile, "%s\n", s_EventBuf.c_str());
			s_EventBuf.clear();
		}
		assembleTimeString(m_TimeString);
		fprintf(s_LogFile, "\n\n/***********************************Application stops at %s***********************************/\n", m_TimeString);
		fprintf(s_LogFile, "/*******************************************************************************************************************/\n", m_TimeString);
		fclose(s_LogFile);
	}
}

/****************************************************************************************************************************************************/
//FUNCTION: force write the buffered content in s_EventBuf to log file 
void CEventLog::flushBufferEvents()
{
	if (s_EventBuf.empty()) return; 

	s_OutputMutex.lock();
	if (fopen_s(&s_LogFile, s_LogFileName, "a+") == 0)
	{
		fprintf(s_LogFile, "%s", s_EventBuf.c_str());
		fclose(s_LogFile);
		s_EventBuf.clear();
		m_numEventInBuffer = 0;
	}
	s_OutputMutex.unlock();
}

}