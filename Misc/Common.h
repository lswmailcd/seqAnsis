#pragma once

#define EOS '\0'  
#define MAXLINE 5000            /* Max. line length */

#define SAFE_DELETE(x)       { delete   (x); (x)=NULL; }
#define SAFE_DELETE_ARRAY(x) { delete[] (x); (x)=NULL; }

const int MAXNAMES = 150;    /* Max chars read for seq. names */

const int  NONE = 0;

const float EPS=0.001f;

const int OK = -200;

