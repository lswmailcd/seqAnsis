#include "Common.h"
#include "Utility.h"

namespace SeqAnsis
{

bool CUtility::blankLine( char* vLine )
{
    int i;

    for (i = 0; vLine[i] != '\n' && vLine[i] != '\0'; i++)
    {
        if (isdigit(vLine[i]) || isspace(vLine[i]) || (vLine[i] == '*') ||
            (vLine[i] == ':') || (vLine[i] == '.'))
            ;
        else
        {
            return false;
        }
    }
    return true;
}

bool CUtility::lineType(char *vpLine, const char *vpCode)
{
   int n;
   if (strlen(vpLine)<strlen(vpCode))
     n=strlen(vpLine);
   else
     n=strlen(vpCode);
   
   return (strncmp(vpLine, vpCode, strlen(vpCode)) == 0);
}

/**
 * Removes trailing blanks from a string
 */
void CUtility::rTrim(std::string *str)
{
   std::string::reverse_iterator rit = str->rbegin();
     
     while(rit != str->rend() && isspace(*rit))
     {
         str->erase((++rit).base());
     }
}

bool		CUtility::blankLineNumericLabel(char *line)
{
	int i;
	int dots = 0;
	bool isnumeric = false;

	for (i = 0; line[i] != '\n' && line[i] != EOS; i++)
	{
		if (isdigit(line[i]) || isspace(line[i]) || (line[i] == '*') ||
			(line[i] == ':') || (line[i] == '.'))
			;
		else
		{
			return false;
		}
		if(line[i] == '.')
			dots++;
		if(isdigit(line[i]))
			isnumeric = true;
	}
	if(isnumeric && dots > 10)
		return false;
	else
		return true;
}

/**
 * Removes trailing blanks from a string
 * @param str String to remove trailing blanks from.
 * @return Pointer to the processed string
 */
char * CUtility::rTrim(char *str)
{
    register int p;

    p = strlen(str) - 1;

    while (isspace(str[p]))
    {
        p--;
    }

    str[p + 1] = EOS;

    return str;
}

std::string CUtility::blankToUnderscore(std::string vStr)
{
    int i, p;

    p = vStr.size() - 1;

    for (i = 0; i <= p; i++)
        if ((vStr[i] == ' ') || (vStr[i] == ';') || (vStr[i] == ',') || (vStr[i] =='(') || (vStr[i] == ')') || (vStr[i] == ':'))
        {
            vStr[i] = '_';
        }

    return vStr;
}

long		CUtility::max3(long a, long b, long c)
{
	long result = a;	

	if(b > result) {	
		result = b;	
	}	
	if(c > result) {	
		result = c;	
	}	

	return result;	
}

}