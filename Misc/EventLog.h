#pragma once

#include <string>
#include <iostream>
#include <stdio.h>
#include"AppException.h"

namespace SeqAnsis
{

enum OutputType
{
	DEF_OUTPUT_SUCCESS = 250,
	DEF_OUTPUT_WARNING = 251,
	DEF_OUTPUT_ERROR   = 252,
};

//FUNCTION: this class is used for handling exceptions and writing the application running log
//NOTE:     By default, the class outputs all the messages to std::cout. User can also redirect the output target by calling setMessageDisplayFunc()
class CEventLog
{
public:
	CEventLog(void);
	~CEventLog(void);

	void writeException(const CAppException& vException);
	void writeStdException(const std::exception& vException);
	//EVENT  buffered by program untile it is full and store to file
	void writeEvent(const std::string& vEvent, bool vImmediateWrite=true, bool vOutputToScreen=false);
	void writeEvent(const char *vEvent, bool vImmediateWrite=true, bool vOutputToScreen=false);
	//START and END EVENT immediately output to file
	void writeAppStartEvent();
	void writeAppEndEvent();
	//WARNING  immediately output to file
	void writeWarning(int vLine, const std::string& vSource, const std::string& vFile, const std::string& vDescription);
	void writeWarning(void *vOutputTarget, const std::string& vSource, const std::string& vFile, int vLine, const std::string& vDescription);
	//MESSAGE  buffered by program untile it is full and store to file
	void writeMessage(void *vOutputTarget, const std::string& vDescription, bool vWriteToEvent=true);
	void writeMessage(const std::string& vDescription, bool vWriteToEvent=true);

private:
	char              m_TimeString[128];
	int               m_numEventInBuffer;

	void assembleTimeString(char *voTimeString);
	void outputEvent(const char *vEvent, OutputType vEventType, bool vOutputToScreen=false, bool vImmediateWrite=false);
	void flushBufferEvents();
};

}