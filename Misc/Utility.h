#pragma once

#include <string>
#include<time.h>

namespace SeqAnsis
{

class CUtility
{
public:
	void			initRandomNumber(){srand( (unsigned)time(0) );}
	double		getRandomNumber(){return	rand()/(double)RAND_MAX;}//get one random number at [0, 1)
	int			getRandomNumber( int mod ){	return	mod>0?rand()%mod:-1;	} //get one random number at [0, mod-1]

	bool blankLine( char* vLine );
	bool lineType(char *line, const char *code);
    char* rTrim(char *str);
     void rTrim(std::string *str);
	std::string blankToUnderscore(std::string vStr);
	long max3(long a, long b, long c);
    template <class T> T MIN(T x, T y){if (x <= y){return x;}return y;}
    template <class T> T MAX(T x, T y){if (x >= y){return x;}return y;}
	bool blankLineNumericLabel(char *line);
};

}

