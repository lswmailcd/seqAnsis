#pragma once

#include<string>

#define __EXCEPTION_SITE__   __LINE__, __FUNCTION__, __FILE__ 

namespace SeqAnsis
{

//define the exception code
enum ExceptionCode
{
	DEF_EXCEPTION_UNEXPECTED         	= 100,
	DEF_EXCEPTION_INVALID_PARAMETER  	= 101,
	DEF_EXCEPTION_INVALID_FILE       	= 102,
	DEF_EXCEPTION_NULL_POINTER       	= 103,
	DEF_EXCEPTION_OPENGL_ERROR       	= 104,
	DEF_EXCEPTION_OBJECT_EXISTED     	= 105,
	DEF_EXCEPTION_INDEX_OUT_OF_RANGE 	= 106,
	DEF_EXCEPTION_UNINITIALIZED      	= 107,
	DEF_EXCEPTION_UNKNOWN_KEYWORD    	= 108,
	DEF_EXCEPTION_BAD_FORMAT         	= 109,
	DEF_EXCEPTION_UNSUPPORTED_FUNC   	= 110,
	DEF_EXCEPTION_MEM_ALLOC_FAILURE  	= 111,
	DEF_EXCEPTION_MEM_RELEASE_FAILURE	= 112,
	DEF_EXCEPTION_BAD_CONFIG            = 113,
	DEF_EXCEPTION_ASSERT_FAIL           = 114,
	DEF_EXCEPTION_STD                   = 115,
	DEF_EXCEPTION_OBJECT_NOT_EXIST      = 116,
	DEF_EXCEPTION_FILE_NOT_EXIST        = 117,
	DEF_EXCEPTION_WAIT_OUT_TIME         = 118,
	DEF_EXCEPTION_CONFIG_FILE_NOT_EXIST = 119,
	DEF_EXCEPTION_NETWORK               = 120,
	DEF_EXCEPTION_FAIL_CONNECT_DATABASE = 121,
	DEF_EXCEPTION_DATABASE_ERROR        = 122,
};

enum ExceptionLevel
{
	DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK = 200,
	DEF_EXCEPTION_LEVEL_EXIT_APP  = 201,
	DEF_EXCEPTION_LEVEL_WARNING   = 202,
};

//define an exception
class CAppException:public std::exception
{
public:
	CAppException(ExceptionCode vCode, ExceptionLevel vLevel, int vExpLine, const std::string& vExpFunc, const std::string& vExpFile, const std::string& vDescription) :
	m_Code(vCode), m_Level(vLevel), m_ExpLine(vExpLine), m_ExpFunc(vExpFunc), m_ExpFile(vExpFile), m_Description(vDescription) {}
	CAppException(const CAppException& e) {*this = e;}

	ExceptionCode  m_Code;
	ExceptionLevel m_Level;
	int            m_ExpLine;
	std::string    m_ExpFunc;
	std::string    m_Description;
	std::string    m_ExpFile;
};

}