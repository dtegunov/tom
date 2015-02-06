#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* We use case-insensitive string comparison functions strcasecmp(), strncasecmp().
   These are a BSD addition and are also defined on Linux, but not on every OS,
   in particular Windows. The two "inline" declarations below fix this problem. If
   you get errors on other platforms, move the declarations outside the WIN32 block */

/*#define WIN32*/

#ifdef WIN32
#include <windows.h>
inline int strcasecmp(const char *s1,const char *s2) { return strcmp(s1,s2); }
inline int strncasecmp(const char *s1,const char *s2,size_t n) { return strncmp(s1,s2,n); }

#else
/* Wegen fehler bei mexErrMsgTxt unter linux */
//#define NO_MEXERRMSGTXT
#endif

#include <mysql.h>  //  Definitions for MySQL client API
#include <mex.h>    //  Definitions for Matlab API

const bool debug=false;  //  turn on information messages

/**********************************************************************
 **********************************************************************
 *
 * http://www.mmf.utoronto.ca/resrchres/mysql/
 *
 * Matlab interface to MySQL server
 * Documentation on how to use is in "mysql.m"
 *  Robert Almgren, 2000-2004; latest review and revision August 2004
 *
 * Here is a list of the routines in this file:
 *
 *  mexFunction()      Entry point from Matlab, presumably as "mysql"
 *   |
 *   +- getstring()    Extract string from Matlab array
 *   +- hostport()     Get host and port number from combined string
 *   |
 *   +- fix_types()    Some datatype fudges for MySQL C API
 *   +- can_convert()  Whether we can convert a MySQL data type to numeric
 *   +- field2num()    Convert a MySQL data string to a Matlab number
 *   |   |
 *   |   +- daynum()   Convert year,month,day to Matlab date number
 *   |   +- typestr()  Translate MySQL data type into a word, for errors
 *   |
 *   +- fancyprint()   Display query result in nice form
 *
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 *
 * hostport(s):  Given a host name s, possibly containing a port number
 *    separated by the port separation character (normally ':').
 * Modify the input string by putting a null at the end of
 * the host string, and return the integer of the port number.
 * Return zero in most special cases and error cases.
 * Modified string will not contain the port separation character.
 * Examples:  s="myhost:2515" modifies s to "myhost" and returns 2515.
 *            s="myhost"      leaves s unchanged and returns 0.
 *
 **********************************************************************/

const char portsep = ':';   //  separates host name from port number

static int hostport(char *s)
{
	//   Look for first portsep in s; return 0 if null or can't find
	if ( !s || !( s=strchr(s,portsep) ) ) return 0;

	//  If s points to portsep, then truncate and convert tail
	*s++ = 0;
	return atoi(s);   // Returns zero in most special cases
}

/**********************************************************************
 *
 * typestr(s):  Readable translation of MySQL field type specifier
 *              as listed in   mysql_com.h
 *
 **********************************************************************/

static const char *typestr( enum_field_types t )
{
	switch(t)
		{
		//  These are considered numeric by IS_NUM() macro
		case FIELD_TYPE_DECIMAL:      return "decimal";
		case FIELD_TYPE_TINY:         return "tiny";
		case FIELD_TYPE_SHORT:        return "short";
		case FIELD_TYPE_LONG:         return "long";
		case FIELD_TYPE_FLOAT:        return "float";
		case FIELD_TYPE_DOUBLE:       return "double";
		case FIELD_TYPE_NULL:         return "null";
		case FIELD_TYPE_LONGLONG:     return "longlong";
		case FIELD_TYPE_INT24:        return "int24";
		case FIELD_TYPE_YEAR:         return "year";
		case FIELD_TYPE_TIMESTAMP:    return "timestamp";

		//  These are not considered numeric by IS_NUM()
		case FIELD_TYPE_DATE:         return "date";
		case FIELD_TYPE_TIME:         return "time";
		case FIELD_TYPE_DATETIME:     return "datetime";
		case FIELD_TYPE_NEWDATE:      return "newdate";   // not in manual
		case FIELD_TYPE_ENUM:         return "enum";
		case FIELD_TYPE_SET:          return "set";
		case FIELD_TYPE_TINY_BLOB:    return "tiny_blob";
		case FIELD_TYPE_MEDIUM_BLOB:  return "medium_blob";
		case FIELD_TYPE_LONG_BLOB:    return "long_blob";
		case FIELD_TYPE_BLOB:         return "blob";
		case FIELD_TYPE_VAR_STRING:   return "var_string";
		case FIELD_TYPE_STRING:       return "string";
		case FIELD_TYPE_GEOMETRY:     return "geometry";   // not in manual 4.0

		default:                      return "unknown";
		}
}

/**********************************************************************
 *
 * fancyprint():  Print a nice display of a query result
 *     We assume the whole output set is already stored in memory,
 *     as from mysql_store_result(), just so that we can get the
 *     number of rows in case we need to clip the printing.
 *     In any case, we make only one pass through the data.
 *
 *     If the number of rows in the result is greater than NROWMAX,
 *     then we print only the first NHEAD and the last NTAIL.
 *     NROWMAX must be greater than NHEAD+NTAIL, normally at least
 *     2 greater to allow the the extra information
 *     lines printed when we clip (ellipses and total lines).
 *
 *     Display null elements as empty
 *
 **********************************************************************/

const char *contstr = "...";  //  representation of missing rows
const int contstrlen = 3;     //  length of above string
const int NROWMAX = 50;       //  max number of rows to print w/o clipping
const int NHEAD = 10;         //  number of rows to print at top if we clip
const int NTAIL = 10;         //  number of rows to print at end if we clip

static int fancyprint( MYSQL_RES *res ) {

    unsigned long nrow = mysql_num_rows(res);
	unsigned long nfield=mysql_num_fields(res);

	bool clip = ( nrow > NROWMAX );



	MYSQL_FIELD *f = mysql_fetch_fields(res);


	/************************************************************************/
	//  Determine column widths, and format strings for printing

	//  Find the longest entry in each column header,
	//    and over the rows, using MySQL's max_length



	int *len = (int *) mxMalloc( nfield * sizeof(int) );
	{
        for ( int j=0 ; j<nfield ; j++ ) {
            len[j] = strlen(f[j].name);
		    if ( f[j].max_length > len[j] )
                len[j] = f[j].max_length;
        }
    }




	//  Compare to the continuation string length if we will clip
	if (clip) {
	    for (int j=0; j<nfield; j++) {
		    if (contstrlen > len[j]) {
			    len[j] = contstrlen;
			}
		}
	}

	//  Construct the format specification for printing the strings
	char **fmt = (char**) mxMalloc(nfield * sizeof(char*));

	{
	    for (int j=0; j<nfield; j++ ) {
		    fmt[j] = (char *) mxCalloc( 10, sizeof(char) );
		    sprintf(fmt[j],"  %%-%ds ",len[j]);
		}
	}

	/************************************************************************/
	//  Now print the actual data

	mexPrintf("\n");

	//  Column headers
	{
	    for (int j=0; j<nfield; j++) {
		    mexPrintf(fmt[j], f[j].name);
		}
	}

	mexPrintf("\n");

	//  Fancy underlines
	{
	    for (int j=0; j<nfield; j++) {
		    mexPrintf(" +");
		    for (int k=0; k<len[j]; k++) {
			    mexPrintf("-");
			}
		    mexPrintf("+");
		}
	}
	mexPrintf("\n");

	//  Print the table entries
	if (nrow<=0) {
	    mexPrintf("(zero rows in result set)\n");
	} else {
	    if (!clip) {
		    mysql_data_seek(res,0);
		    for ( int i=0 ; i<nrow ; i++ ) {
    		    MYSQL_ROW row = mysql_fetch_row(res);
	    		if (!row) {
		    	    mexPrintf("Printing full table data from row %d\n",i+1);
#ifndef NO_MEXERRMSGTXT
			    	mexErrMsgTxt("Internal error:  Failed to get a row");
#else
				    mexPrintf("Internal error:  Failed to get a row");
				    return 1;
#endif
                }
                for ( int j=0 ; j<nfield ; j++ ) {
				    mexPrintf( fmt[j], ( row[j] ? row[j] : "" ) );
				}
	            mexPrintf("\n");
		    }
	    } else {         //  print half at beginning, half at end
			mysql_data_seek(res,0);
		    {
			    for ( int i=0 ; i<NHEAD ; i++ ) {
				    MYSQL_ROW row = mysql_fetch_row(res);
			        if (!row) {
					    mexPrintf("Printing head table data from row %d\n",i+1);
#ifndef NO_MEXERRMSGTXT
                        mexErrMsgTxt("Internal error:  Failed to get a row");
#else
                        mexPrintf("Internal error:  Failed to get a row");
                        return 1;
#endif
                    }
			        for ( int j=0 ; j<nfield ; j++ ) {
				        mexPrintf( fmt[j], ( row[j] ? row[j] : "" ) );
					}
			        mexPrintf("\n");
				}
			}
		    {
			    for ( int j=0 ; j<nfield ; j++ ) {
				    mexPrintf(fmt[j],contstr);
				}
			}
		    mexPrintf("\n");
		    mysql_data_seek( res, nrow - NTAIL );
		    {
    		    for ( int i=0 ; i<NTAIL ; i++ ) {
					MYSQL_ROW row = mysql_fetch_row(res);
					if (!row) {
						mexPrintf("Printing tail table data from row %d",nrow-NTAIL+i+1);
#ifndef NO_MEXERRMSGTXT
						mexErrMsgTxt("Internal error:  Failed to get a row");
#else
						mexPrintf("Internal error:  Failed to get a row");
                        return 1;
#endif
					}
					for ( int j=0 ; j<nfield ; j++ )
						mexPrintf( fmt[j], ( row[j] ? row[j] : "" ) );
					mexPrintf("\n");
                }
		    }
		    mexPrintf("(%d rows total)\n",nrow);
		}
	}
	mexPrintf("\n");

	// These should be automatically freed when we return to Matlab,
	//  but just in case ...
	mxFree(len);  mxFree(fmt);
    return 0;
}

/**********************************************************************
 *
 * field2num():  Convert field in string format to double number
 *
 **********************************************************************/

const double NaN = mxGetNaN();         //  Matlab NaN for null values
const double secinday = 24.*60.*60.;   //  seconds in one day

//==============================================
//  year,month,day --> serial date number, based on Matlab's own algorithm

const int cummonday[2][12] =
	{ {  0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 },
	  {  0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 } };

inline int cl( int y, int n )  { return ((y-1)/n)+1; }  // ceil(y/n)

static int daynum( int yr, int mon, int day )
{
	//  This is exactly Matlab's formula
	int leap = ( !(yr%4) && (yr%100) ) | !(yr%400);
	return 365*yr + cl(yr,4) - cl(yr,100) + cl(yr,400)
	        + cummonday[leap][mon-1] + day;
}

//================================================

static bool can_convert( enum_field_types t )
{
	return ( IS_NUM(t)
	   || ( t == FIELD_TYPE_DATE )
	   || ( t == FIELD_TYPE_TIME )
	   || ( t == FIELD_TYPE_DATETIME ) );
}

static double field2num( const char *s, enum_field_types t )
{
	if (!s) return NaN;  // MySQL null -- nothing there

	if ( IS_NUM(t) )
		{
		long double val=NaN;
		if ( sscanf(s,"%Lf",&val) != 1 )
			{ mexPrintf("Unreadable value \"%s\" of type %s\n",s,typestr(t));
			  return NaN; }
		return val;
		}
	else if ( t == FIELD_TYPE_DATE )
		{
		int yr, mon, day;
		if ( sscanf(s,"%4d-%2d-%2d",&yr,&mon,&day) != 3
		    || yr<1000 || yr>=3000 || mon<1 || mon>12 || day<1 || day>31 )
			{ mexPrintf("Unreadable value \"%s\" of type %s\n",s,typestr(t));
			  return NaN; }
		return (double) daynum(yr,mon,day);
		}
	else if ( t == FIELD_TYPE_TIME )
		{
		int hr, min, sec;
		if ( sscanf(s,"%2d:%2d:%2d",&hr,&min,&sec) != 3
		           || min>60 || sec>60 )
			{ mexPrintf("Unreadable value \"%s\" of type %s\n",s,typestr(t));
			  return NaN; }
		return (sec+60*(min+60*hr))/secinday;
		}
	else if ( t == FIELD_TYPE_DATETIME )
		{
		int yr, mon, day, hr, min, sec;
		if ( sscanf(s,"%4d-%2d-%2d %2d:%2d:%2d",&yr,&mon,&day,&hr,&min,&sec) != 6
		      || yr<1000 || yr>=3000 || mon<1 || mon>12 || day<1 || day>31
		      || min>60 || sec>60 )
			{ mexPrintf("Unreadable value \"%s\" of type %s\n",s,typestr(t));
			  return NaN; }
		return ((double) daynum(yr,mon,day))
		       + ((sec+60*(min+60*hr))/secinday );
		}
	else
		{
		mexPrintf("Tried to convert \"%s\" of type %s to numeric\n",s,typestr(t));
#ifndef NO_MEXERRMSGTXT
		mexErrMsgTxt("Internal inconsistency");
#else
		mexPrintf("Internal inconsistency");
#endif
		}
}

/**********************************************************************
 *
 * fix_types():   I have seen problems with the MySQL C API reporting
 *   types inconsistently. Returned data value from fields of type
 *   "date" or "time" has type set appropriately to DATE, TIME, etc.
 *   However, if any operation is performed on such data, even something
 *   simple like max(), then the result can sometimes be reported as the
 *   generic type STRING. The developers say this is forced by the SQL
 *   standard but to me it seems strange.
 *
 *   This function looks at the actual result data returned. For each
 *   column that is reported as type STRING, we see if we can make a
 *   more precise classification as DATE, TIME, or DATETIME.
 *
 *   For fields of type STRING, we look at the length, then content
 *    length = 8, 9, or 10:  can be time as HH:MM:SS, HHH:MM:SS, or HHHH:MM:SS
 *    length = 10:           can be date in form YYYY:MM:DD
 *    length = 19:           can be datetime in form YYYY:MM:DD HH:MM:SS
 *
 **********************************************************************/

static int fix_types( MYSQL_FIELD *f, MYSQL_RES *res ) {
	if (!res) {    //  This should never happen
#ifndef NO_MEXERRMSGTXT
		mexErrMsgTxt("Internal error:  fix_types called with res=NULL");
#else
		mexPrintf("Internal error:  fix_types called with res=NULL");
        return 1;
#endif
    }

	unsigned long nrow=mysql_num_rows(res), nfield=mysql_num_fields(res);

	if (nrow<1) {
	    return 0;  // nothing to look at, nothing to fix
	}

	bool *is_unknown = (bool *) mxMalloc( nfield * sizeof(bool) );
	int n_unknown=0;   //  count number of fields of unknown type
	for ( int j=0 ; j<nfield ; j++ ) {
	    is_unknown[j] = (f[j].type==FIELD_TYPE_STRING &&
		                 (f[j].length==8 || f[j].length==9 || f[j].length==10
		                                 || f[j].length==19));
		if (is_unknown[j]) {
		    n_unknown++;
		}
	}

	if (debug) {
	    mexPrintf("Starting types:");
		for ( int j=0 ; j<nfield ; j++ ) {
		    mexPrintf("  %s(%d)%s", typestr(f[j].type), f[j].length,
		              ( is_unknown[j] ? "?" : "" ) );
		}
	    mexPrintf("\n");
	}

	//  Look at successive rows as long as some columns are still unknown
	//  We go through columns only to find the first non-null data value.
	mysql_data_seek(res,0);
	for ( int i=0 ; n_unknown>0 && i<nrow ; i++ ) {
		MYSQL_ROW row = mysql_fetch_row(res);
		if (!row) {
            mexPrintf("Scanning row %d for type identification\n",i+1);
#ifndef NO_MEXERRMSGTXT
		    mexErrMsgTxt("Internal error in fix_types():  Failed to get a row");
#else
		    mexPrintf("Internal error in fix_types():  Failed to get a row");
            return 1;
#endif
        }

		if (debug) {
            mexPrintf("  row[%d]:",i+1);
			for ( int j=0 ; j<nfield ; j++ ) {
			   mexPrintf("  \"%s\"%s", ( row[j] ? row[j] : "NULL" ),
			             ( is_unknown[j] ? "?" : "" ) );
            }
		    mexPrintf("\n");
        }

		//  Look at each field to see if we can extract information
		for ( int j=0 ; j<nfield ; j++ ) {
			//  If this column is still a mystery, and if there is data here,
			//    then try extracting a date and/or time out of it
			if ( is_unknown[j] && row[j] ) {
				int yr, mon, day, hr, min, sec;

				if ( f[j].length==19
				      && sscanf(row[j],"%4d-%2d-%2d %2d:%2d:%2d",
				            &yr,&mon,&day,&hr,&min,&sec) == 6
				      && yr>=1000 && yr<3000 && mon>=1 && mon<=12 && day>=1 && day<=31
				      && min<=60 && sec<=60 ) {
				    f[j].type = FIELD_TYPE_DATETIME;
                } else if ( f[j].length==10
				      && sscanf(row[j],"%4d-%2d-%2d",&yr,&mon,&day) == 3
				      && yr>=1000 && yr<3000 && mon>=1 && mon<=12 && day>=1 && day<=31 ) {
				    f[j].type = FIELD_TYPE_DATE;
                } else if ( ( f[j].length==8 || f[j].length==9 || f[j].length==10 )
				      && sscanf(row[j],"%d:%2d:%2d",&hr,&min,&sec) == 3
				      && min<=60 && sec<=60 ) {
					f[j].type = FIELD_TYPE_TIME;
                }

				//  If the tests above failed, then the type is not date or time;
				//  it really is a string of unknown type.
				//  Whether the tests suceeded or failed, it is no longer a mystery.
				is_unknown[j]=false;  n_unknown--;
			}
		}
	}

	if (debug) {
        mexPrintf("  Ending types:");
		for ( int j=0 ; j<nfield ; j++ ) {
		   mexPrintf("  %s(%d)%s", typestr(f[j].type), f[j].length,
	                 ( is_unknown[j] ? "?" : "" ) );
        }
		mexPrintf("\n");
    }

	mxFree(is_unknown);  // should be automatically freed, but still...
    return 0;
}

/**********************************************************************
 *
 * getstring():   Extract string from a Matlab array
 *    (Space allocated by mxCalloc() should be freed by Matlab
 *     when control returns out of the MEX-function.)
 *   This is based on an original by Kimmo Uutela
 *
 **********************************************************************/

static char *getstring(const mxArray *a)
{
	int llen = mxGetM(a)*mxGetN(a)*sizeof(mxChar) + 1;
	char *c = (char *) mxCalloc(llen,sizeof(char));
	if (mxGetString(a,c,llen))
#ifndef NO_MEXERRMSGTXT
		mexErrMsgTxt("Can\'t copy string in getstring()");
#else
		mexPrintf("Can\'t copy string in getstring()");
#endif
	return c;
}

/**********************************************************************
 *
 * mysql():  Execute the actual action
 *
 *  Which action we perform is based on the first input argument,
 *  which must be present and must be a character string:
 *    'open', 'close', 'use', 'status', or a legitimate MySQL query.
 *
 *  This version does not permit binary characters in query string,
 *  since the query is converted to a C null-terminated string.
 *
 *  If no output argument is given, then information is displayed.
 *  If an output argument is given, then we operate silently, and
 *     return status information.
 *
 **********************************************************************/

/*********************************************************************/
//  Static variables that contain connection state
/*
 *  isopen gets set to true when we execute an "open"
 *  isopen gets set to false when either we execute a "close"
 *                        or when a ping or status fails
 *   We do not set it to false when a normal query fails;
 *   this might be due to the server having died, but is much
 *   more likely to be caused by an incorrect query.
 */

class conninfo {
public:
	MYSQL *conn;   //  MySQL connection information structure
	bool isopen;   //  whether we believe that connection is open
	conninfo()  { conn=NULL;  isopen=false; }
};

const int MAXCONN = 10;
static conninfo c[MAXCONN];   //  preserve state for MAXCONN different connections

extern "C" void mexFunction(int nlhs, mxArray *plhs[],
      int nrhs, const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[],
      int nrhs, const mxArray *prhs[]) {
	int cid=0;   // ID number of the current connection
	int jarg=0;  // Number of first string arg: becomes 1 if id is specified

    //mexPrintf("%s", mysql_get_client_info());
    //return;

	// Parse the first argument to see if it is a specific id number
	if (nrhs > 0 && mxIsNumeric(prhs[0])) {
        if (mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=1 ) {
            mexPrintf("Usage:  %s( [id], command, [options] )\n",
                mexFunctionName());
			mexPrintf("First argument is array %d x %d\n",
			    mxGetM(prhs[0]),mxGetN(prhs[0]) );
#ifndef NO_MEXERRMSGTXT
			mexErrMsgTxt("Invalid connection ID");
#else
			mexPrintf("Invalid connection ID");
            return;
#endif
        }
		double xid = *mxGetPr(prhs[0]);
		cid = int(xid);
		if ( double(cid)!=xid || cid<0 || cid>=MAXCONN ) {
            mexPrintf("Usage:  %s( [id], command, [options] )\n",
                mexFunctionName());
			mexPrintf("id = %g -- Must be integer between 0 and %d\n",
                xid, MAXCONN-1);
#ifndef NO_MEXERRMSGTXT
			mexErrMsgTxt("Invalid connection ID");
#else
			mexPrintf("Invalid connection ID");
            return;
#endif
        }
		jarg = 1;
    }

	if (debug) mexPrintf("cid = %d  jarg = %d\n", cid, jarg);

	// Shorthand notation now that we know which id we have
	typedef MYSQL *mp;
	mp &conn = c[cid].conn;
	bool &isopen = c[cid].isopen;

    char* command = NULL; // String with commandtype.

    //  Check the commandtype
    enum querytype { OPEN, CLOSE, USE, STATUS, CMD, ESCAPE } q;
    if (nrhs-jarg < 1) {
        q = STATUS;
        command = "status (implicit)";
    } else {
        if (!mxIsChar(prhs[jarg])) {
            mexPrintf("Usage:  %s( [id], command, [options] )\n",
                mexFunctionName());
#ifndef NO_MEXERRMSGTXT
            mexErrMsgTxt("Argument command must be String.");
#else
            mexPrintf("Argument command must be String.");
            return;
#endif
        } else {
            char* command = getstring(prhs[jarg]);

            // Check command
            if (!strcasecmp(command, "open")) {
                q = OPEN;
            } else if (!strcasecmp(command, "close")) {
                q = CLOSE;
            } else if (!strcasecmp(command, "use")) {
                q = USE;
            } else if (!strcasecmp(command, "status")) {
                q = STATUS;
            } else if (!strcasecmp(command, "escape")) {
                q = ESCAPE;
            } else if (!strcasecmp(command, "cmd")) {
                q = CMD;
            } else {
                mexPrintf("Usage:  %s( [id], command, [options] )\n",
                    mexFunctionName());
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("command must be one of [open, close, use, status, escape, cmd].");
#else
                mexPrintf("command must be one of [open, close, use, status, escape, cmd].");
                return;
#endif
            }
        }
    }

    // Check argumenttypes.
    switch (q) {
        case STATUS:
            if (nrhs-jarg >= 2) {
                mexPrintf("Usage:  %s( [id], command, [options] )\n",
                    mexFunctionName());
                mexPrintf("Syntax: %s([id], ['status'])\n",
                    mexFunctionName());
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("Too much arguments");
#else
                mexPrintf("Too much arguments");
                return;
#endif
            }
            break;
        case CLOSE:
            if (nrhs-jarg >= 2) {
                mexPrintf("Usage:  %s( [id], command, [options] )\n",
                    mexFunctionName());
                mexPrintf("Syntax: %s([id], 'close')\n",
                    mexFunctionName());
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("Too much arguments");
#else
                mexPrintf("Too much arguments");
                return;
#endif
            }
            break;
        case USE:
            if ((nrhs-jarg != 2) || !mxIsChar(prhs[jarg+1])) {
                mexPrintf("Usage:  %s( [id], command, [options] )\n",
                    mexFunctionName());
                mexPrintf("Syntax: %s([id], 'use', databasename)\n",
                    mexFunctionName());
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("Too much arguments or databasename is not a string");
#else
                mexPrintf("To much arguments or databasename is not a string");
                return;
#endif
            }
            break;
        case CMD:
            if ((nrhs-jarg < 2) || !mxIsChar(prhs[jarg+1])) {
                mexPrintf("Usage:  %s( [id], command, [options] )\n",
                    mexFunctionName());
                mexPrintf("Syntax: %s([id], 'cmd', sqlquery)\n",
                    mexFunctionName());
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("To much arguments or sqlquery is not a string");
#else
                mexPrintf("To much arguments or sqlquery is not a string");
                return;
#endif
            }
            break;
        case OPEN:
            if ((nrhs-jarg != 4) || !mxIsChar(prhs[jarg+1]) || !mxIsChar(prhs[jarg+2]) ||
                !mxIsChar(prhs[jarg+3])) {
                mexPrintf("Usage:  %s( [id], command, [options] )\n",
                    mexFunctionName());
                mexPrintf("Syntax: %s([id], 'open', hostname, username, password)\n",
                    mexFunctionName());
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("Wrong number of arguments or arguments are not strings");
#else
                mexPrintf("Wrong number of arguments or arguments are not strings");
                return;
#endif
            }
            break;
        case ESCAPE:
            if ((nrhs-jarg != 3) || !mxIsNumeric(prhs[jarg+1]) || !mxIsChar(prhs[jarg+2])) {
                mexPrintf("Usage:  %s( [id], command, [options] )\n",
                    mexFunctionName());
                mexPrintf("Syntax: %s([id], 'escape', length, string)\n",
                    mexFunctionName());
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("Wrong number of arguments or arguments are the wrong type");
#else
                mexPrintf("Wrong number of arguments or arguments are the wrong type");
                return;
#endif
            } else {
                // Is length an integer?
                if ( mxGetM(prhs[1+jarg])!=1 || mxGetN(prhs[1+jarg])!=1 ) {
                    mexPrintf("Usage:  %s( [id], command, [options] )\n",
                        mexFunctionName());
                    mexPrintf("Syntax: %s([id], 'escape', length, string)\n",
                        mexFunctionName());
                    mexPrintf("length argument is array %d x %d\n", mxGetM(prhs[1+jarg]),mxGetN(prhs[1+jarg]) );
#ifndef NO_MEXERRMSGTXT
                    mexErrMsgTxt("Wrong argumenttype for length");
#else
                    mexPrintf("Wrong argumenttype for length");
                    return;
#endif
                }
                double xlength = *mxGetPr(prhs[jarg+1]);
                int escape_length = int(xlength);
                if ( double(escape_length)!=xlength || escape_length<0) {
                    mexPrintf("Usage:  %s( [id], command, [options] )\n",
                        mexFunctionName());
                    mexPrintf("Syntax: %s([id], 'escape', length, string)\n",
                        mexFunctionName());
                    mexPrintf("length = %g -- Must be integer greater or equal 0\n", xlength);
#ifndef NO_MEXERRMSGTXT
                    mexErrMsgTxt("Wrong size for length or length is a double value");
#else
                    mexPrintf("Wrong size for length or length is a double value");
                    return;
#endif
                }
            }
            break;
        default:
            mexPrintf("Usage:  %s( [id], command, [options] )\n",
                mexFunctionName());
            mexPrintf("Unknown command %s\n", command);
#ifndef NO_MEXERRMSGTXT
            mexErrMsgTxt("Unknown command");
#else
            mexPrintf("Unknown command");
            return;
#endif
    }

	if (debug) {
        mexPrintf("command = %d (%s)\n", q, command);
    }


    switch (q) {
        case OPEN : {

            //  Close connection if it is open
            if (isopen) {
                mysql_close(conn);
                isopen = false;
                conn = NULL;
            }

            //  Extract information from input arguments
            char *host = getstring(prhs[jarg+1]);
            char *user = getstring(prhs[jarg+2]);
            char *pass = getstring(prhs[jarg+3]);
            int port = hostport(host);  // returns zero if there is no port

            if (nlhs<1) {
                mexPrintf("Connecting to  host=%s", host);
                mexPrintf("  user=%s", user);
                mexPrintf("  password=%s", pass);
                if (port) {
                    mexPrintf("  port=%d",port);
                }
                mexPrintf("\n");
            }

            //  Establish and test the connection
            //  If this fails, then conn is still set, but isopen stays false
            if (!(conn=mysql_init(conn))) {
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("Couldn\'t initialize MySQL connection object");
#else
                mexPrintf("Couldn\'t initialize MySQL connection object");
                return;
#endif
            }
            if (!mysql_real_connect( conn, host, user, pass, NULL, port, NULL,0 )) {
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt(mysql_error(conn));
#else
                mexPrintf(mysql_error(conn));
                return;
#endif
            }

            const char *c = mysql_stat(conn);

            if (c) {
                if (nlhs<1) {
                    mexPrintf("%s\n",c);
                } else {
#ifndef NO_MEXERRMSGTXT
                    mexErrMsgTxt(mysql_error(conn));
#else
                    mexPrintf(mysql_error(conn));
                    return;
#endif
                }
            }

            isopen = true;

            //  Now we are OK -- if he wants output, give him a 1
            if (nlhs>=1) {
                if (!(plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL))) {
#ifndef NO_MEXERRMSGTXT
                    mexErrMsgTxt("Unable to create matrix for output");
#else
                    mexPrintf("Unable to create matrix for output");
                    return;
#endif
                }
                *mxGetPr(plhs[0]) = 1.;
            }

        }
        break;
        case CLOSE : {
            if (isopen) {
                mysql_close(conn);
                isopen = false;
                conn = NULL;
            }
        }
        break;
        case ESCAPE : {
            if (nlhs > 2) {
                mexPrintf("You specified %i output arguments\n", nlhs);
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("Max. 2 output arguments");
#else
                mexPrintf("Max. 2 output arguments");
                return;
#endif
            } else {
                if (!isopen) {
#ifndef NO_MEXERRMSGTXT
                    mexErrMsgTxt("Not connected");
#else
                    mexPrintf("Not connected");
                    return;
#endif
                }
                if (mysql_ping(conn)) {
                    isopen=false;
#ifndef NO_MEXERRMSGTXT
                    mexErrMsgTxt(mysql_error(conn));
#else
                    mexPrintf(mysql_error(conn));
                    return;
#endif
                }

                char *esc_str = getstring(prhs[2+jarg]);

                int escape_length = int(*mxGetPr(prhs[jarg+1]));

                char* res = (char*)mxCalloc(escape_length*2+1, sizeof(char));
                unsigned long res_length = mysql_real_escape_string(conn, res, getstring(prhs[2+jarg]), escape_length);

                switch (nlhs) {
                    case 0:  mexPrintf("Neue Laenge ist %i\n\"%s\"\n", res_length, res);
                             break;
                    case 2:  plhs[1] = mxCreateDoubleScalar(res_length);
                    case 1:  plhs[0] = mxCreateString(res);
                    default:
                             ;
                }
            }
        }
        break;
        case USE : {
            if (!isopen) {
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("Not connected");
#else
                mexPrintf("Not connected");
                return;
#endif
            }

            if (mysql_ping(conn)) {
                isopen=false;
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt(mysql_error(conn));
#else
                mexPrintf(mysql_error(conn));
                return;
#endif
            }

            char *db = getstring(prhs[1]);

            if (mysql_select_db(conn,db)) {
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt(mysql_error(conn));
#else
                mexPrintf(mysql_error(conn));
                return;
#endif
            }

            if (nlhs<1) {
                mexPrintf("Current database is \"%s\"\n",db);
            } else {
    			if (!( plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL))) {
#ifndef NO_MEXERRMSGTXT
                    mexErrMsgTxt("Unable to create matrix for output");
#else
                    mexPrintf("Unable to create matrix for output");
                    return;
#endif
                }
                double *pr = mxGetPr(plhs[0]);
                *pr = 0.0;
            }
        }
        break;
        case STATUS : {
    		if (nlhs<1) {  //  He just wants a report
                if (jarg>0) {  //  mysql(cid,'status')  Give status of the specified id
    				if (!isopen) {
                        mexPrintf("%2d: Not connected\n", cid);
                        return;
                    }
                    if (mysql_ping(conn)) {
                        isopen=false;
#ifndef NO_MEXERRMSGTXT
                        mexErrMsgTxt(mysql_error(conn));
#else
                        mexPrintf(mysql_error(conn));
                        return;
#endif
                    }
                    mexPrintf("%2d:  %-30s   Server version %s\n",
                        cid, mysql_get_host_info(conn), mysql_get_server_info(conn) );
				}
                else {         //  mysql('status') with no specified connection
    				int nconn=0;
                    for (int j=0; j<MAXCONN; j++ ) {
                        if (c[j].isopen) {
                            nconn++;
                        }
                    }
                    if (debug) {
                        mexPrintf("%d connections open\n", nconn);
                    }
                    if (nconn==0) {
                        mexPrintf("No connections open\n");
                        return;
                    }
                    if (nconn==1 && c[0].isopen) { // Only connection number zero is open
                                                   //  Give simple report with no connection id #
                        if (mysql_ping(conn)) {
                            isopen=false;
#ifndef NO_MEXERRMSGTXT
                            mexErrMsgTxt(mysql_error(conn));
#else
                            mexPrintf(mysql_error(conn));
                            return;
#endif
                        }
                        mexPrintf("Connected to %s   Server version %s   Client %s\n",
                            mysql_get_host_info(conn), mysql_get_server_info(conn), mysql_get_client_info() );
					} else {       //  More than one connection is open
					               //     Give a detailed report of all open connections
                        for (int j=0; j<MAXCONN; j++ ) {
                            if (c[j].isopen) {
                                if (mysql_ping(c[j].conn)) {
                                    c[j].isopen = false;
                                    mexPrintf("%2d:  %s\n",mysql_error(c[j].conn));
                                } else {
                                    mexPrintf("%2d:  %-30s   Server version %s\n", j,
                                        mysql_get_host_info(c[j].conn), mysql_get_server_info(c[j].conn) );
                                }
                            }
                        }
					}
				}
			} else {        //  He wants a return value for this connection
    			if (!( plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL))) {
#ifndef NO_MEXERRMSGTXT
                    mexErrMsgTxt("Unable to create matrix for output");
#else
                    mexPrintf("Unable to create matrix for output");
                    return;
#endif
                }
                double *pr = mxGetPr(plhs[0]);
                if (!isopen) {
                    *pr=1.;
                    return;
                }
                if (mysql_ping(conn)) {
                    isopen = false;
                    *pr=2.0;
                    return;
                }
                if (!mysql_stat(conn)) {
                    isopen=false;
                    *pr = 3.0;
                    return;
                }
                *pr = 0.0;
			}
        }
        break;
        case CMD : {


            //  Check that we have a valid connection
            if (!isopen) {
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("No connection open");
#else
                mexPrintf("No connection open");
                return;
#endif
            }
            if (mysql_ping(conn)) {
                isopen = false;
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt(mysql_error(conn));
#else
                mexPrintf(mysql_error(conn));
                return;
#endif
            }

            //  Execute the query (data stays on server)
            if (mysql_query(conn, getstring(prhs[jarg + 1]))) {
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt(mysql_error(conn));
#else
                mexPrintf(mysql_error(conn));
                return;
#endif
            }


            //  Download the data from server into our memory
            //     We need to be careful to deallocate res before returning.
            //  Matlab's allocation routines return instantly if there is not
            //  enough free space, without giving us time to dealloc res.
            //  This is a potential memory leak but I don't see how to fix it.
            MYSQL_RES *res = mysql_store_result(conn);

            //  As recommended in Paul DuBois' MySQL book (New Riders, 1999):
            //  A NULL result set after the query can indicate either
            //    (1) the query was an INSERT, DELETE, REPLACE, or UPDATE, that
            //        affect rows in the table but do not return a result set; or
            //    (2) an error, if the query was a SELECT, SHOW, or EXPLAIN
            //        that should return a result set but didn't.
            //  Distinguish between the two by checking mysql_field_count()
            //  We return in either case, either correctly or with an error
            if (!res) {
                if (!mysql_field_count(conn)) {
                    unsigned long nrows = mysql_affected_rows(conn);
                    if (nlhs<1) {
                        mexPrintf("%u rows affected\n", nrows);
                        return;
                    } else {
                        if (!(plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL))) {
#ifndef NO_MEXERRMSGTXT
                            mexErrMsgTxt("Unable to create numeric matrix for output");
#else
                            mexPrintf("Unable to create numeric matrix for output");
                            return;
#endif
                        }
                        *(mxGetPr(plhs[0])) = (double) nrows;
                        return;
                    }
                } else {
#ifndef NO_MEXERRMSGTXT
                    mexErrMsgTxt(mysql_error(conn));
#else
                    mexPrintf(mysql_error(conn));
                    return;
#endif
                }
            }

            unsigned long nrow = mysql_num_rows(res);
            unsigned long nfield = mysql_num_fields(res);

            //  If he didn't ask for any output (nlhs=0),
            //       then display the output and return
            if (nlhs<1) {
                if (fancyprint(res) != 0) {
                    return;
                }
                mysql_free_result(res);
                return;
            }


            //  If we are here, he wants output
            //  He must give exactly the right number of output arguments
            if (nlhs != nfield) {
                mysql_free_result(res);
                mexPrintf("You specified %d output arguments, and got %d columns of data\n",
                    nlhs, nfield);
#ifndef NO_MEXERRMSGTXT
                mexErrMsgTxt("Must give one output argument for each column");
#else
                mexPrintf("Must give one output argument for each column");
                return;
#endif
            }

            //  Fix the column types to fix MySQL C API sloppiness
            MYSQL_FIELD *f = mysql_fetch_fields(res);
            if (fix_types(f, res) != 0) {
                return;
            }

            //  Create the Matlab arrays for output
            double** pr = (double**) mxMalloc(nfield * sizeof(double*));
            for (int j=0; j<nfield; j++ ) {
                if (can_convert(f[j].type)) {
                    if (!(plhs[j] = mxCreateDoubleMatrix(nrow, 1, mxREAL))) {
                        mysql_free_result(res);
#ifndef NO_MEXERRMSGTXT
                        mexErrMsgTxt("Unable to create numeric matrix for output");
#else
                        mexPrintf("Unable to create numeric matrix for output");
                        return;
#endif
                    }
                    pr[j] = mxGetPr(plhs[j]);
                } else {
                    if (!(plhs[j] = mxCreateCellMatrix(nrow, 1))) {
                        mysql_free_result(res);
#ifndef NO_MEXERRMSGTXT
                        mexErrMsgTxt("Unable to create cell matrix for output");
#else
                        mexPrintf("Unable to create cell matrix for output");
                        return;
#endif
                    }
                    pr[j] = NULL;
                }
            }


            //  Load the data into the cells
            mysql_data_seek(res, 0);
            for (int i=0; i<nrow; i++) {
                MYSQL_ROW row = mysql_fetch_row(res);
                if (!row) {
                    mexPrintf("Scanning row %d for data extraction\n",i+1);
#ifndef NO_MEXERRMSGTXT
                    mexErrMsgTxt("Internal error:  Failed to get a row");
#else
                    mexPrintf("Internal error:  Failed to get a row");
                    return;
#endif
                }
                for (int j=0; j<nfield; j++) {
                    if (can_convert(f[j].type)) {
                        pr[j][i] = field2num(row[j],f[j].type);
                    } else {
                        mxArray *c = mxCreateString(row[j]);
                        mxSetCell(plhs[j],i,c);
                    }
                }
            }
            mysql_free_result(res);
        }
        break;
        default : {
            mexPrintf("Unknown query type q = %d (%s)\n", q, command);
#ifndef NO_MEXERRMSGTXT
            mexErrMsgTxt("Internal code error");
#else
            mexPrintf("Internal code error");
            return;
#endif
        }
        break;
    }
}
