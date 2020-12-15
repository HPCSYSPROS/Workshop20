 /*@@
   @file      ParseFile.c
   @date      Tue Jan 12 15:58:31 1999
   @author    Tom Goodale
   @desc
              Routines to read in a parameter file and pass the resulting data
              to a user-supplied subroutine.
              Currently taken from the old cactus ones and slightly modifed.
   @enddesc
 @@*/

/*#define DEBUG*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>

#include <cctk.h>

#include "cctk_Flesh.h"
#include "cctk_CommandLine.h"
#include "util_String.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(util_ParseFile_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static void CheckBuf(int, int);
static void removeSpaces(char *stripMe);
static char *ReadFile(FILE *file, unsigned long *filesize);
static char *ParseDefines(char *buffer, unsigned long *buffersize);
static void convert_crlf_to_lf(char *buffer);
int ParseBuffer(char *buffer,
                int (*set_function)(const char *, const char *, int),
                tFleshConfig *ConfigData);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

int ParseFile(FILE *ifp,
              int (*set_function)(const char *, const char *, int),
              tFleshConfig *ConfigData);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/* parse buffer size */
#define BUF_SZ   (8 * 1024)

/* line number */
static int lineno = 1;

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine ParseFile
   @author Paul Walker, Frank Loeffler
   @desc
   This routine actually parses the parameter file. The
   syntax we allow is
   <ul>
     <li>a = b
         <li>a,b,c = d,e,f
         <li># rest of the line is ignored
         <li>x = "string value"
   </ul>
   So it is easy to parse
   <p>
   We go through the file looking for stuff and then set
   it in the global database using calls to the passed in set_function.
   @enddesc
   @history
   @hdate Tue Jan 12 16:41:36 1999 @hauthor Tom Goodale
   @hdesc Moved to CCTK.
          Changed to pass data to arbitrary function.
          Changed to take a file descriptor rather than a filename.
          Use a buffer to parse which can be preprocessed before main parsing
   @endhistory
   @var     ifp
   @vdesc   The filestream to parse
   @vtype   FILE *
   @vio     in
   @vcomment

   @endvar
   @var     set_function
   @vdesc   The function to call to set the value of a parameter
   @vtype   int (*)(const char *, const char *)
   @vio     in
   @vcomment

   @endvar
   @var     ConfigData
   @vdesc   Flesh configuration data
   @vtype   tFleshConfig *
   @vio     in
   @vcomment

   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/

int cctk_PirahaParser(const char *buffer,unsigned long buffersize,int (*set_function)(const char *, const char *, int));

int ParseFile(FILE *ifp,
              int (*set_function)(const char *, const char *, int),
              tFleshConfig *ConfigData)
{
  int retval=1;
  unsigned long buffersize;
  char *buffer = ReadFile(ifp, &buffersize);
  if (!buffer)
    return 1;

  const int piraha_active = 1;
  if(piraha_active)
  {
    /* the new way */
    buffersize = strlen(buffer);

    retval = cctk_PirahaParser(buffer, buffersize, set_function);
  }
  else
  {
    /* The old way */
    /* Ensure Unix line endings */
    convert_crlf_to_lf(buffer);
    buffersize = strlen(buffer);

    buffer = ParseDefines(buffer, &buffersize);
    /* ParseBuffer can get confused with detecting the end of the buffer
       (in comment or in a string), and may overrun.  Therefore
       we te a buffer that is a bit longer.  */
    {
      buffer = realloc (buffer, strlen(buffer) + 10);
      memset (buffer+strlen(buffer), '\0', 10);
    }
    retval = ParseBuffer(buffer, set_function, ConfigData);
    free(buffer);
  }
  return retval;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine CheckBuf
   @author Paul Walker
   @desc
   A simple description and warning message in case of
   a fixed parse buffer overflow.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
   @var     p
   @vdesc   buffer location
   @vtype   int
   @vio     in
   @vcomment

   @endvar
   @var     l
   @vdesc   Line number
   @vtype   int
   @vio     in
   @vcomment

   @endvar

@@*/

static void CheckBuf(int p, int l)
{
  if (p >= BUF_SZ)
  {
    fprintf(stderr,"WARNING: Parser buffer overflow on line %d\n",
            l);
    fprintf(stderr,"This indicates either an incorrect parm file or\n");
    fprintf(stderr,"the need to recompile " __FILE__ " with a bigger\n");
    fprintf(stderr,"BUF_SZ parm.\n");

    assert(0);
    exit(1);
  }
}


 /*@@
   @routine removeSpaces
   @author Paul Walker
   @desc
   removes the spaces from a char * <b>in place</b>. Beware
   that this function will change the input value and if you
   want to keep a copy you have to do so yourself!
   @enddesc
   @calls
   @calledby
   @history

   @endvar
   @var     stripMe
   @vdesc   String to strip
   @vtype   char *
   @vio     inout
   @vcomment

   @endvar

@@*/
static void removeSpaces(char *stripMe)
{
  char *to = stripMe;
  for (char *from = stripMe; *from; ++from) {
    if (!isspace(*from)) {
      *to++ = *from;
    }
  }
  *to = '\0';
}

/* Convert string CRLF line endings to LF */
static void convert_crlf_to_lf(char *buffer)
{
  char *to = buffer;
  for (char *from = buffer; *from; ++from) {
    if (*from == '\r' && *(from+1) == '\n') {
      // do nothing -- skip the \r
    } else {
      *to++ = *from;
    }
  }
  *to = '\0';
}


/* #define TEST_ParseFile */

#ifdef TEST_ParseFile

int parameter_printer(const char *param, const char *val)
{
  printf("Parameter %s has value %s\n", param, val);

  return 0;
}

int main(int argc, char *argv[])
{
  int retval;
  FILE *parameter_file;

  if(argc > 1)
  {
    parameter_file = fopen(argv[1], "r");
    if(parameter_file)
    {
      ParseFile(parameter_file, parameter_printer);
      fclose(parameter_file);
      retval = 0;
    }
    else
    {
      retval=2;
    }
  }
  else
  {
    printf("Usage: %s <filename>\n", argv[0]);
    retval = 1;
  };

  return 0;
}

#endif


 /*@@
   @routine ReadFile
   @author Frank Loeffler
   @desc
   This routine reads a file into a buffer
   @enddesc
   @history
   @var     file
   @vdesc   The filestream to read
   @vtype   FILE *
   @vio     in
   @vcomment

   @endvar
   @var     filesize
   @vdesc   The size of the file
   @vtype   *unsigned long
   @vio     out
   @vcomment

   @endvar

   @returntype char *
   @returndesc
   NULL - failure
   !NULL allocated buffer
   @endreturndesc
@@*/
static char *ReadFile(FILE *file, unsigned long *filesize)
{
  char *buffer;

  if (!file)
  {
    fprintf(stderr, "Could not use file for reading.\n");
    return NULL;
  }
  /* Get the file size */
  fseek(file, 0, SEEK_END);
  *filesize = ftell(file);
  fseek(file, 0, SEEK_SET);
  /* Allocate buffer */
  buffer = (char *)malloc(*filesize+1);
  if (!buffer)
  {
    fprintf(stderr, "Could not allocate memory.\n");
    return NULL;
  }
  /* Read file into buffer and return */
  fread(buffer, *filesize, 1, file);
  /* Protect buffer for string operations */
  buffer[*filesize] = '\0';
  return buffer;
}

 /*@@
   @routine ParseDefines
   @author Frank Loeffler
   @desc
   This routine parses a buffer for defines of the form \$[^ \n] and
   replaces them with some predefined defines, returning the new buffer,
   which might be the same as the old buffer. Otherwise the old buffer
   gets freed by this function. In case of parsing errors, the whole
   parsing is stopped and the old buffer is returned
   (FIXME: Should we return NULL?)
   @enddesc
   @history
   @var     buffer
   @vdesc   The buffer to parse
   @vtype   char *
   @vio     in
   @vcomment

   @endvar
   @var     buffersize
   @vdesc   The size of the buffer
   @vtype   *unsigned long
   @vio     out
   @vcomment

   @endvar

   @returntype char *
   @returndesc
   !NULL - new buffer, might be == buffer
   @endreturndesc
@@*/
static char *ParseDefines(char *buffer, unsigned long *buffersize)
{
  /* define name */
  char define[1024];
  /* Position in define name */
  size_t defpos = 0;
  /* Current position in buffer */
  size_t pos = 0;
  /* Character at current position */
  char c;
  /* Position of start of found define */
  size_t def_start = 0;
  /* Flag to indicate if we are inside a definition name */
  int indef = 0;
  if (!buffer)
    return buffer;
  /* Main loop over the buffer */
  while (((c=buffer[pos]) != '\0') && (pos < *buffersize) )
  {
    /* Check if we found a define */
    if (c == '$')
    {
      /* Mark the state and location and go on */
      indef = 1;
      def_start = pos;
    }
    /* If we are inside a define */
    else if (indef)
    {
      /* Look for the end of the define name */
      if (c == '}' || c == '\0' || c == '\n')
      {
        /* make the define name a proper string */
        define[defpos] = '\0';
        /* this holds the value of the define to be filled in */
        char *value = NULL;
        /* a new buffer for the whole thing */
        char * new_buffer;
        unsigned long new_size;
        /* Check for different kinds of defines */
        /* Environment variables */
        if ((defpos > 6) &&
            (strncmp(define, "ENV{'", 5) == 0) &&
            (define[defpos-1] == '\''))
        {
          /* overwrite the trailing ' */
          define[defpos-1] = '\0';
          value = getenv(define+5);
          if (!value)
          {
            CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "No environment variable %s found\n", define+5);
            /* TODO: Should we abort here? (return NULL) */
          }
          /* increase pos to jump over the trailing } */
          if (pos+1 < *buffersize) { pos++; }
        }
        /* Parameter file name ? */
        else if ((strcmp(define, "parfile") == 0) ||
                 (strcmp(define, "{parfile") == 0))
        {
          char path[500];
          CCTK_ParameterFilename(500, path);
          value = strrchr (path, '/');
          if (value == NULL)
          {
            value = path;
          }
          else
          {
            value++;
          }
          /* skip the parameter file extension */
          if (strcmp (value + strlen (value) - 4, ".par") == 0)
          {
            value[strlen (value) - 4] = '\0';
          }
          if (define[0] == '{')
          {
            /* increase pos to jump over the trailing } */
            if (pos+1 < *buffersize) { pos++; }
          }
        }
        /* Else: unknown define - or no define at all: ignore */
        if (value)
        {
          /* Replace the old buffer with the new, copying in all the data */
          new_size = *buffersize - strlen(define) + strlen(value) + 1;
          if (new_size < def_start)
          {
            CCTK_WARN(0, "Something is wrong with me, HELP!");
            return buffer;
          }
          new_buffer = (char *)malloc(new_size);
          if (!new_buffer)
          {
            CCTK_WARN(0, "I am out of memory and give up parsing for defines.");
            return buffer;
          }
          new_buffer[0]='\0';
          strncat(new_buffer, buffer, def_start);
          strncat(new_buffer, value, new_size-1);
          strncat(new_buffer, buffer+pos, new_size-1);
          new_buffer[new_size-1] = '\0';
          /* free the old buffer */
          free(buffer);
          /* Start further processing of new buffer */
          pos = def_start + strlen(define);
          buffer = new_buffer;
          *buffersize = new_size;
        }
        /* General define cleanup */
        indef = 0;
        defpos = 0;
        define[0] = '\0';
        def_start = 0;
      }
      else if (defpos > 1023)
      {
        /* We only print this warning and ignore the possible define problem.
         * It might not be a define after all but a valid, long parameter
         * value containing a '$'.*/
        define[1023] = '\0';
        CCTK_VWarn(CCTK_WARN_PICKY, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Possible define too long: %s", define);
        indef = 0;
        defpos = 0;
        define[0] = '\0';
        def_start = 0;
      }
      else
        define[defpos++] = c;
    }
    pos++;
  }
  return buffer;
}

 /*@@
   @routine ParseBuffer
   @author Paul Walker
   @desc
   This routine actually parses the parameter file. The
   syntax we allow is
   <ul>
     <li>a = b
         <li>a,b,c = d,e,f
         <li># rest of the line is ignored
         <li>x = "string value"
   </ul>
   So it is easy to parse
   <p>
   We go through the file looking for stuff and then set
   it in the global database using calls to the passed in set_function.
   @enddesc
   @history
   @hdate Tue Jan 12 16:41:36 1999 @hauthor Tom Goodale
   @hdesc Moved to CCTK.
          Changed to pass data to arbitrary function.
          Changed to take a file descriptor rather than a filename.
          Use a buffer to parse which can be preprocessed before main parsing
   @endhistory
   @var     buffer
   @vdesc   The buffer to parse
   @vtype   char *
   @vio     in
   @vcomment

   @endvar
   @var     set_function
   @vdesc   The function to call to set the value of a parameter
   @vtype   int (*)(const char *, const char *)
   @vio     in
   @vcomment

   @endvar
   @var     ConfigData
   @vdesc   Flesh configuration data
   @vtype   tFleshConfig *
   @vio     in
   @vcomment

   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
int ParseBuffer(char *buffer,
              int (*set_function)(const char *, const char *, int),
              tFleshConfig *ConfigData)
{
  /* position in buffer */
  unsigned int pos = 0;
  /* Buffers for parsing from the file */
  char *tokens, *value;
  char *subtoken, *subvalue;
  /* Positions in the buffers */
  int ntokens;
  /* Status flags */
  int intoken, inval;
  /* Current char */
  char c;
  int num_errors; /* number of errors in file parsing */

  num_errors = 0;

  /* avoid compiler warning about unused parameter */
  (void) (ConfigData + 0);

  /* allocate parse buffers */
  tokens = (char *) malloc (4 * BUF_SZ);
  value    = tokens + 1*BUF_SZ;
  subtoken = tokens + 2*BUF_SZ;
  subvalue = tokens + 3*BUF_SZ;

  intoken = 0; inval = 0;

  while ((c=buffer[pos++]) != '\0')
  {
#ifdef DEBUG
    printf("%c",c);
#endif
    /* Main Loop */
    while (c == '#' || c == '!' )
    {
      /* Comment line.  So forget rest of line */
      while ((c=buffer[pos++]) != '\n' && c != '\r' && c != '\0')
      {
#ifdef DEBUG
        printf("%c",c);
#endif
      }
      if (c == '\n')
      {
        lineno++;
      }
      c = buffer[pos++];
#ifdef DEBUG
      printf("%c",c);
#endif
    }

    /* End of line */
    if (c == '\n')
    {
      if(intoken)
      {
        fprintf(stderr, "Parse error at line %d.  No value supplied.\n", lineno);
        num_errors++;
        intoken = 0;
      }

#ifdef DEBUG
      printf ("LINE %d\n",lineno);
#endif
      lineno ++;

    }

    /* Token character */
    if (intoken && c != '=')
    {
      tokens[intoken++] = (char)c;
      CheckBuf(intoken,lineno);
    }


    /* Start of a new token */
    if (c != ' ' && c != '\t' && c != '\n' && c != '\r' && !inval && !intoken)
    {
      intoken = 0;
      tokens[intoken++] = (char)c;
    }

    /* End of a token signified by an = */
    if (c == '=')
    {
      if (intoken)
      {
        unsigned int ll;
        tokens[intoken] = '\0';  /* Very important! */
        intoken = 0;
        inval = 0;
        removeSpaces(tokens);
        ntokens = 1;
        for (ll=0;ll < strlen(tokens); ll++)
          if (tokens[ll] == ',') ntokens++;
#ifdef DEBUG
        printf ("\nNew token! >>%s<<\n",tokens);
        printf ("%d token elements\n",ntokens);
#endif

        /* Scan ahead to the beginning of the value
         * and check if the value is a string or not.
         * This parser DOES strip quotes off of the strings.
         */
        while ((c = buffer[pos++]) == ' ' || c == '\n' || c == '\r' || c == '\t')
        {
#ifdef DEBUG
          printf("%c",c);
#endif
          if (c == '\n')
          {
#ifdef DEBUG
            printf ("LINE %d\n",lineno);
#endif
            lineno++;
          }
        }

        if (c == '"')
        {
          int startline = lineno;

          /* Just handle the thing. */
          int p = 0;
          if (ntokens > 1)
          {
            fprintf (stderr, "%s%s%s\n",
                     "WARNING: Multiple string ",
                     "tokens not supported for ",
                     tokens);
            fprintf(stderr, "This is a fatal error");
            /* deallocate parse buffers */
            free (tokens);
            return 1;
          }
          while ((c = buffer[pos++]) != '"')
          {
#ifdef DEBUG
            printf("%c",c);
#endif
            if (c == '#' && lineno != startline)
            {
              /* found a comment inside a multi-line string */
              while ((c = buffer[pos++]) != '\n')
              {
                if (c == '\0')
                {
                  break;
                }
              }
              if (c == '\0')
              {
                value[p++] = (char)c;
                CheckBuf(p,lineno);
                break;
              }
#ifdef DEBUG
              printf ("LINE %d\n",lineno);
#endif
              lineno++;
            }
            else
            {
              value[p++] = (char)c;
              CheckBuf(p,lineno);
              if (c == '\n')
              {
#ifdef DEBUG
                printf ("LINE %d\n",lineno);
#endif
                lineno++;
              }
              else if (c == '\0')
              {
                break;
              }
            }
          }
          /* ignore all extra spaces or a trailing comment in this line */
          if (c == '"')
          {
            int is_comment = 0;
            while ((c = buffer[pos++]) != '\n')
            {
              if (c == '\0')
              {
                break;
              }
              if (c == '#')
              {
                is_comment = 1;
              }
              else if (! (is_comment || c == ' ' || c == '\t' || c == '\r'))
              {
                fprintf (stderr, "ERROR: extra characters found after closing "
                                 "quote for value of parameter '%s' starting "
                                 "in line %d\n",
                         tokens, startline);
                fprintf(stderr, "This is a fatal error.\n");
                free (tokens);
                return 1;
              }
            }
            if (c == '\n' || c == '\0')
            {
              if (c == '\n') lineno++;
              /* mark it as properly quoted string value */
              c = '"';
            }
          }
          /* fail if the quoted string did not terminate with a quote */
          if (c != '"')
          {
            fprintf (stderr, "ERROR: no closing quote found for quoted value "
                             "of parameter '%s' starting in line %d\n",
                     tokens, startline);
            fprintf(stderr, "This is a fatal error.\n");
            free (tokens);
            return 1;
          }
          value[p] = '\0';
#ifdef DEBUG
          printf ("\nString %s -> %s\n",
                  tokens,value);
#endif
          set_function(tokens,value, lineno);
        }
        else
        {

          int p = 0;
          value[p++] = (char)c;
          if (ntokens == 1)
          {
            /* Simple case. We have an int
               or a double which contain no
               spaces! */
            c = buffer[pos++];
#ifdef DEBUG
            printf("%c",c);
#endif
            while (!(c==' ' || c=='\t' || c == '\n' || c == '\r' || c == '\0'))
            {
              value[p++] = (char)c;
              CheckBuf(p,lineno);
              c = buffer[pos++];
#ifdef DEBUG
              printf("%c",c);
#endif
            }
            value[p] = '\0';
#ifdef DEBUG
            printf ("Parsed %d characters\n", p);
            printf("\nFloat/Int: %s -> %s\n", tokens,value);
#endif
            set_function(tokens,value,lineno);
            if (c=='\n')
            {
#ifdef DEBUG
              printf ("LINE %d\n",lineno);
#endif
              lineno++;
            }
          }
          else
          {
            /* Harder case of multiple tokens */
            int ncommas = 0;
            int pp=0, i;
            int pt, pv;

            value[pp++] = (char)c;
            /* OK, since we only have numbers in the
               old input stream, we can go along getting
               ntokens-1 commas, stripping spaces, and
               make a nice little string.
               */
            c = buffer[pos++];
#ifdef DEBUG
            printf("%c",c);
#endif
            while (ncommas < ntokens-1 && c != '\0')
            {
              if (!(c == ' ' || c == '\t' || c == '\n' || c == '\r'))
              {
                value[pp++] = (char)c;
                CheckBuf(pp,lineno);
              }
              if (c == ',') ncommas ++;
              c = buffer[pos++];
#ifdef DEBUG
              printf("%c",c);
#endif
            }
            if (c == ' ' || c == '\t')
            {
              /* Great now strip out the spaces */
              while((c = buffer[pos++]) == ' ' || c=='\t' || c == '\n' || c == '\r')
              {
#ifdef DEBUG
                printf("%c",c);
#endif
                if (c =='\n')
                {
#ifdef DEBUG
                  printf ("LINE %d\n",lineno);
#endif
                  lineno++;
                }
              }
            }

            /* And tack the rest on */
            value[pp++] = (char)c;
            CheckBuf(p,lineno);

            c = buffer[pos++];
#ifdef DEBUG
            printf("%c",c);
#endif
            while (c != ' ' && c != '\t' && c != '\n' && c != '\r' && c != '\0')
            {
              value[pp++] = (char)c;
              CheckBuf(pp,lineno);
              c = buffer[pos++];
#ifdef DEBUG
              printf("%c",c);
#endif
            }
            value[pp] = '\0';
#ifdef DEBUG
            printf("Comma list: %s -> %s\n",
                   tokens,value);
#endif
            /* So parse out the tokens */
            pt = 0;
            pv = 0;
            for (i=0;i<ncommas;i++)
            {
              pp = 0;
              while (tokens[pt] != ',')
              {
                subtoken[pp++] = tokens[pt++];
                CheckBuf(p,lineno);
              }
              subtoken[pp] = '\0';
              pp = 0;
              while (value[pv] != ',')
              {
                subvalue[pp++] = value[pv++];
                CheckBuf(pp,lineno);
              }
              subvalue[pp] = '\0';

              set_function(subtoken,subvalue,lineno);
#ifdef DEBUG
              printf("Setting sub-token %s -> %s\n",
                     subtoken, subvalue);
#endif
              /* Now remember we are sitting on a comma
               * in both our input strings, so bump by one
               */
              pv ++; pt ++;
            }
            /* And OK, so now we have one parameter left
             * so lets handle that
             */
            pp = 0;
            while (tokens[pt] != '\0')
            {
              subtoken[pp++] = tokens[pt++];
              CheckBuf(pp,lineno);
            }
            subtoken[pp] = '\0';
            pp = 0;
            while (value[pv] != '\0')
            {
              subvalue[pp++] = value[pv++];
              CheckBuf(pp,lineno);
            }
            subvalue[pp] = '\0';

            set_function(subtoken,subvalue,lineno);
          }
        }
      }
      else
      {
        fprintf (stderr, "Parser failed at = on line %d\n",
                 lineno);
      }
    }
  }

  /* deallocate parse buffers */
  free (tokens);

  return num_errors;
}

