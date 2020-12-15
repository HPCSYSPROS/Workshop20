/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   strlib contains a bunch of functions that are useful for basic 
   file input and string manipulation.  See strlib.h for a list and 
   description of the functions that are available.
*/

#include "strlib.h"

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_read_line				*/
/*									*/
/*   INPUTS:								*/
/*	fd - FILE to read line from					*/
/*	buf - buffer to put line into					*/
/*									*/
/*   OUTPUTS:								*/
/*	this function returns a 0 if the line is read successfully or   */
/*   a -1 if an EOF is encountered					*/
/*									*/
/*	NAMD_read_line reads in a line from an open file.  It places    */
/*   the line that is read in into the buffer passed in via buf.  If    */
/*   an EOF is encountered, a -1 is returned.  If an EOF is encountered */
/*   in the middle of a line, the program will terminate abnormally.    */
/*   Also, any comments of the form {blah blah blah} will automatically */
/*   be skipped by NAMD_read_line, even if usch comments span lines.    */
/*   Also, the string will be left justified.  That is, any spaces at   */
/*   the begining of the line will be removed.				*/
/*									*/
/************************************************************************/

int NAMD_read_line(FILE *fd, char *buf, int bufsize)

{
	int i=0;	//  Current position in buf
	int c;		//  Character read in from file

	/*  Loop and read characters until we get either an EOF or a    */
	/*  newline							*/
	while ( ((c=fgetc(fd)) != EOF) && (c != '\n') )
	{
		/*  If we encounter a bracketed comment, skip it.  This */
		/*  basically means read EVERYTHING until the next } and*/
		/*  throw it into the big bit bucket			*/
		if (c == '{')
		{
			while ( ((c=fgetc(fd)) != EOF) && (c!='}') )
			{
			}

			if (c==EOF)
			{
				buf[i]=STRINGNULL;
				char *err_msg = new char[128 + strlen(buf)];
				sprintf(err_msg, "ABNORMAL EOF FOUND - buffer=*%s*\n", buf);
				NAMD_die(err_msg);
			}

			continue;
		}

		/*  Also, use this little piece of logic to remove      */
		/*  any leading spaces from the line			*/
		if ((i>0) || !isspace(c))
		{
			buf[i] = c;
			i++;
			if ( i == bufsize ) {
				i--;
	                        buf[i]=STRINGNULL;
                                char *err_msg = new char[128 + strlen(buf)];
				sprintf(err_msg, "BUFFER OVERRUN - buffer=*%s*\n", 
				   buf);
				NAMD_die(err_msg);
			}
		}
	}

	/*  NULL terminate the string					*/
	buf[i]=STRINGNULL;

	/*  Check for an EOF in the middle of a line			*/
	if ((c==EOF) && (i!=0))
	{
		buf[i]=STRINGNULL;
		char *err_msg = new char[128 + strlen(buf)];
		sprintf(err_msg, "ABNORMAL EOF FOUND - buffer=*%s*\n", buf);
		NAMD_die(err_msg);
	}

	/*  Return the appropriate value				*/
	if (c==EOF)
		return(-1);
	else
		return(0);
}
/*			END OF FUNCTION NAMD_read_line			*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_remove_comment			*/
/*									*/
/*   INPUTS:								*/
/*	str - String to remove comment from				*/
/*									*/
/*	This function removes comments from the end of a line that	*/
/*   are of the form:							*/
/*									*/
/*	sample line		! This is a comment			*/
/*									*/
/************************************************************************/

void NAMD_remove_comment(char *str)

{
	int i=0;

	while ( (str[i] != STRINGNULL) && (str[i] != '!') )
	{
		i++;
	}

	str[i] = STRINGNULL;
}
/*			END OF FUNCTION NAMD_remove_comment		*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_truncate				*/
/*									*/
/*   INPUTS:								*/
/*	str - string to truncate					*/
/*									*/
/*      NAMD_truncate will remove any trailing spaces from a string.    */
/*   i.e.  "AB DF FG     "  would be truncated to "AB DF FG".		*/
/*									*/
/************************************************************************/

void NAMD_truncate(char *str)

{
	int i;		//  Current position in str

	i=strlen(str);

	/*  Loop from the back of the string until we find a non-space  */
	for (i--; i>=0 && isspace(str[i]); i--)
	{
	}
	
	str[i+1]=STRINGNULL;
}
/*			END OF FUNCTION NAMD_truncate			*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_find_word				*/
/*									*/
/*   INPUTS:								*/
/*	source - the string to be searched in				*/
/*	search - the string to be searched for				*/
/*									*/
/*   OUTPUTS:								*/
/*	a 1 is returned if the string search is found within the string */
/*   source, otherwise a 0 is returned.					*/
/*									*/
/*	NAMD_find_word searches for one string within another.  It is   */
/*   usually used to determine if a word appears within a given line.   */
/*   If the word is found, a 1 is returned.  Otherwise, 0 is returned.  */
/*   Case is ignored while doing this search.				*/
/*									*/
/************************************************************************/

int NAMD_find_word(const char *source, const char *search)

{
	int i=0;		//  Position inside source
	int search_length;	//  Length of string search
	int source_length;	//  Length of string source
	int found=0;		//  Flag 1-> found the value

	search_length=strlen(search);
	source_length=strlen(source);

	/*  While we haven't found it and we haven't readched the      */
	/*  point where our current position plus the length of the    */
	/*  string we are looking for is longer than the string itself,*/
	/*  check the next search_length characters for a match	       */
	while (!found && ((i+search_length)<=(source_length)))
	{
		found = (strncasecmp(source+i, search, search_length)==0);

		i++;
	}

	return(found);
}
/*			END OF FUNCTION NAMD_find_word			*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_blank_str				*/
/*									*/
/*   INPUTS:								*/
/*	str - string to test						*/
/*									*/
/*   OUTPUTS:								*/
/*	a 1 is returned if the string is blank, otherwise a 0 is        */
/*   returned								*/
/*									*/
/*	NAMD_blank_str tests to see if a string is blank.  That is,     */
/*   contains only characters where isspace() is true			*/
/*									*/
/************************************************************************/

int NAMD_blank_string(char *str)
{
	int i;		// Current position in str
	int blank=1;	// Flag 1-> string is blank
	int len;	// length of the string str

	len=strlen(str);

	for (i=0; i<len && blank; i++)
	{
		blank = isspace(str[i]);
	}

	return(blank);
}
/*			END OF FUNCTION NAMD_blank_string		*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_find_first_word			*/
/*									*/
/*   INPUTS:								*/
/*	source - string to obtain word from				*/
/*	word - buffer to place word into				*/
/*									*/
/*   OUTPUTS:								*/
/*	word is returned containing the first word of source		*/
/*									*/
/*	This function finds the first word in a string.  The first word */
/*   is defined to be the first set of continuous non-space charaters   */
/*   in a string.  So in the string "   AB14^  FDGFD GFDG"  the first   */
/*   word would be "AB14^".  The word is returned in the string pointed */
/*   to by word.							*/
/*									*/
/************************************************************************/

void NAMD_find_first_word(char *source, char *word)

{
	int i=0;	//  Position within source
	int j=0;	//  Position within word

	/*  Skip leading white space					*/
	while ( (source[i] != STRINGNULL) && isspace(source[i]))
		i++;

	/*  Copy the word						*/
	while ( (source[i] != STRINGNULL) && !isspace(source[i]))
	{
		word[j]=source[i];
		i++;
		j++;
	}

	word[j]=STRINGNULL;

	return;
}
/*			END OF FUNCTION NAMD_find_first_word		*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_read_int 				*/
/*									*/
/*   INPUTS:								*/
/*	fd - file descriptor to read from				*/
/*	msg - string indicating what we are trying to read		*/
/*									*/
/*   OUTPUTS:								*/
/*	the value of the next integer in the file is returned		*/
/*									*/
/*	NAMD_read_int is used to read in integer lists from the .psf    */
/*   file.  It will read the next integer it finds in the file passed   */
/*   to it.  If an alpha character is encountered, the program          */
/*   terminates.  If an EOF is encountered, the program terminates.     */
/*   The string msg is used to indicate what we were trying to read     */
/*   in any error messages.						*/
/*									*/
/************************************************************************/

int NAMD_read_int(FILE *fd, const char *msg)

{
	int i;			//  Loop counter
	int c;			//  Character read in from file
	char tmp_string[11];	//  Temporary string for integer
	int isNeg;
    
	/*  Skip white space				*/
	while ( ((c=fgetc(fd)) == '\n') || isspace(c) )
	{
	}

	/*  Check to make sure we didn't hit EOF	*/
	if (c==EOF)
	{
		char err_msg[128];

		sprintf(err_msg, "EOF ENCOUNTERED READING %s FROM PSF FILE", msg);
		NAMD_die(err_msg);
	}

	/*  Now read in the integer itself		*/
	i=0;

	/* Modified to read an integer with '-' or '+' sign --Chao Mei */
	isNeg = 0;
	if(c=='-'){
	    c = fgetc(fd);
	    isNeg = 1;
	}
	if(c=='+')
	    c = fgetc(fd);
		

	while (!isspace(c))
	{
		/*  Check to make sure we only get #'s  */
		if (!isdigit(c))
		{
			char err_msg[128];

			sprintf(err_msg, "ALPHA CHARCTER ENCOUNTERED WHILE READING %s FROM PSF FILE", msg);
			NAMD_die(err_msg);
		}

		tmp_string[i] = c;
		i++;

		c=fgetc(fd);

		/*  Check to make sure we didn't hit EOF*/
		if (c==EOF)
		{
			char err_msg[128];

			sprintf(err_msg, "EOF ENCOUNTERED WHILE READING %s FROM PSF FILE", msg);
			NAMD_die(err_msg);
		}
	}

	tmp_string[i]=STRINGNULL;

	/*  Convert the string to an integer and return its value	*/
	if(isNeg)
	    return(-atoi(tmp_string));
	else
	   return(atoi(tmp_string));
}
/*			END OF FUNCTION NAMD_read_int			*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_pad				*/
/*									*/
/*	This function pads a string with leading spaces to a specified  */
/*   length.								*/
/*									*/
/************************************************************************/

void NAMD_pad(char *str, size_t length)

{
	char tmp_string[128];
	size_t i;

	if (strlen(str) >= length)
		return;

	for (i=0; i<(length-strlen(str)); i++)
	{
		tmp_string[i] = ' ';
	}

	tmp_string[i] = STRINGNULL;

	strcat(str, tmp_string);
}
/*			END OF FUNCTION NAMD_pad			*/

