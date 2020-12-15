 /*@@
   @file      ActiveThorns.c
   @date      Sun Jul  4 16:15:36 1999
   @author    Tom Goodale
   @desc
              Stuff to deal with activethorns.
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SKBinTree.h"

#include "cctk_Flesh.h"
#include "cctk_FortranString.h"
#include "cctk_WarnLevel.h"

#include "util_String.h"
#include "util_StringList.h"

#include "cctk_ActiveThorns.h"
#include "cctki_ActiveThorns.h"
#include "cctki_Parameter.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_ActiveThorns_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

struct THORN
{
  int active;
  char *implementation;
  const char **requires_thorns;
  const char **activates_thorns;
};

struct IMPLEMENTATION
{
  int active;
  t_sktree *thornlist;
  char *activating_thorn;

  int n_ancestors;
  char **ancestors;

  int n_friends;
  char **friends;
};

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int RegisterImp(const char *name,
                       const char *thorn,
                       const char **ancestors,
                       const char **friends);

static int ActivateThorn(const char *name);
static int ActivateImp(const char *implementation, const char *thorn);

static int CompareStrings(const void *string1, const void *string2);
static int JustPrintThornName(const char *key,void *input, void *dummy);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

int CCTK_FCALL CCTK_FNAME (CCTK_IsThornActive)
                          (ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(CCTK_IsThornCompiled)
     (int *retval, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_IsImplementationCompiled)
                           (int *retval, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_IsImplementationActive)
                           (int *retval, ONE_FORTSTRING_ARG);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static t_sktree *thornlist = NULL;
static t_sktree *implist   = NULL;

static int n_thorns = 0;
static int n_imps   = 0;

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTKi_RegisterThorn
   @date       Sun Jul  4 17:44:14 1999
   @author     Tom Goodale
   @desc
   Registers a thorn with the flesh.
   @enddesc

   @var     attributes
   @vdesc   Thorn attributes
   @vtype   const struct iAttributeList
   @vio     in
   @endvar

   @returntype int
   @returndesc
    0 - success
   -1 - Duplicate thorn
   -2 - memory failure storing thorn
   -3 - memory failure storing implementation name
   -4 - failed to store thorn in tree
   @endreturndesc
@@*/
int CCTKi_RegisterThorn(const struct iAttributeList *attributes)
{
  int i, retval;
  t_sktree *node, *temp;
  struct THORN *thorn;
  const char *name, *imp;
  const char **ancestors, **friends, **requires_thorns, **activates_thorns;


  name = imp = NULL;
  ancestors = friends = requires_thorns = activates_thorns = NULL;

  for(i=0; attributes[i].attribute; i++)
  {
    if(!strcmp(attributes[i].attribute, "name"))
    {
      if(attributes[i].AttributeData.StringList)
      {
        name = attributes[i].AttributeData.StringList[0];
      }
    }
    else if(!strcmp(attributes[i].attribute, "implementation"))
    {
      if(attributes[i].AttributeData.StringList)
      {
        imp = attributes[i].AttributeData.StringList[0];
      }
    }
    else if(!strcmp(attributes[i].attribute, "ancestors"))
    {
      ancestors = attributes[i].AttributeData.StringList;
    }
    else if(!strcmp(attributes[i].attribute, "friends"))
    {
      friends = attributes[i].AttributeData.StringList;
    }
    else if(!strcmp(attributes[i].attribute, "requires thorns"))
    {
      requires_thorns = attributes[i].AttributeData.StringList;
    }
    else if(!strcmp(attributes[i].attribute, "activates thorns"))
    {
      activates_thorns = attributes[i].AttributeData.StringList;
    }
    else
    {
      fprintf(stderr, "Unknown/unimplemented thorn attribute %s\n", attributes[i].attribute);
    }
  }

  /*  printf("Registering thorn %s, which provides %s\n", name, imp);*/

  /* Does the thorn already exist ? */
  node = SKTreeFindNode(thornlist, name);
  if(!node)
  {
    n_thorns++;

    /* Create the structure to hold thorn info. */
    thorn = malloc (sizeof(struct THORN));
    if(thorn)
    {
      thorn->requires_thorns = NULL;
      if (requires_thorns)
      {
        /* count the number of thorns */
        for (i = 0; requires_thorns[i]; i++);

        thorn->requires_thorns = malloc ((i+1) * sizeof(char *));
        thorn->requires_thorns[i] = NULL;
        while (--i >= 0)
        {
          thorn->requires_thorns[i] = Util_Strdup (requires_thorns[i]);
        }
      }

      thorn->activates_thorns = NULL;
      if (activates_thorns)
      {
        /* count the number of thorns */
        for (i = 0; activates_thorns[i]; i++);

        thorn->activates_thorns = malloc ((i+1) * sizeof(char *));
        thorn->activates_thorns[i] = NULL;
        while (--i >= 0)
        {
          thorn->activates_thorns[i] = Util_Strdup (activates_thorns[i]);
        }
      }

      thorn->implementation = Util_Strdup(imp);
      if(thorn->implementation)
      {
        /* Fill out data for the thorn. */
        thorn->active = 0;

        /* Store the data in the tree */
        temp = SKTreeStoreData(thornlist, thornlist, name, thorn);

        if(!thornlist)
        {
          thornlist = temp;
        }

        if(temp)
        {
          /* Register the implementation */
          RegisterImp(imp, name, ancestors, friends);

          retval = 0;
        }
        else
        {
          retval = -4;
        }
      }
      else
      {
        retval = -3;
      }
    }
    else
    {
      retval = -2;
    }
  }
  else
  {
    retval = -1;
  }

  return retval;
}

 /*@@
   @routine    CCTKi_ActivateThorn
   @date       Sun Jul  4 17:46:15 1999
   @author     Tom Goodale
   @desc
   Activates a thorn and the associated implementation assuming
   the implementation isn't already active.
   @enddesc

   @var     name
   @vdesc   Name of thorn to activate
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
    0 - success
   -1 - non-existent thorn
   -2 - internal error
   -3 - thorn already active
   -4 - implementation already active
   @endreturndesc
@@*/
int CCTKi_ActivateThorn(const char *name)
{
  int retval;
  t_sktree *thornnode;
  t_sktree *impnode;

  struct THORN *thorn;
  struct IMPLEMENTATION *imp;

  printf("Activating thorn %s...", name);

  /* Find the thorn */
  thornnode = SKTreeFindNode(thornlist, name);
  if(thornnode)
  {
    thorn = thornnode->data;

    /* Find the implementation */
    impnode = SKTreeFindNode(implist, thorn->implementation);
    if(impnode)
    {
      imp = impnode->data;

      if(!thorn->active)
      {
        if(!imp->active)
        {
          /* Activate the thorn. */
          printf("Success -> active implementation %s\n", thorn->implementation);
          thorn->active = imp->active = 1;

          /* Remember which thorn activated this imp. */
          imp->activating_thorn = Util_Strdup(name);
          retval = 0;
        }
        else
        {
          printf("Failure -> Implementation %s already activated by %s\n", thorn->implementation, imp->activating_thorn);
          retval = -4;
        }
      }
      else
      {
        printf("Failure -> Thorn %s already active\n", name);
        retval = -3;
      }
    }
    else
    {
      printf("Internal error - can't find imp %s from thorn %s\n", thorn->implementation, name);
      retval = -2;
    }
  }
  else
  {
    printf("Failure -> non-existent thorn.\n");
    retval = -1;
  }

  return retval;
}


/*@@
   @routine    CCTK_IsThornActive
   @date       Sun Jul  4 17:46:56 1999
   @author     Tom Goodale
   @desc
   Checks if a thorn is active.
   @enddesc

   @var     name
   @vdesc   Name of thorn
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - inactive
   1 - active
   @endreturndesc
@@*/
int CCTK_IsThornActive(const char *name)
{
  const t_sktree *node;


  node = SKTreeFindNode(thornlist, name);

  return (node ? ((struct THORN *) node->data)->active : 0);
}

int CCTK_FCALL CCTK_FNAME (CCTK_IsThornActive)
                          (ONE_FORTSTRING_ARG)
{
  int retval;
  ONE_FORTSTRING_CREATE(name)
  retval = CCTK_IsThornActive(name);
  free(name);
  return retval;
}

 /*@@
   @routine    CCTK_ThornImplementation
   @date       Sun Oct 17 21:10:19 1999
   @author     Tom Goodale
   @desc
   Returns the implementation provided by the thorn.
   @enddesc

   @var     name
   @vdesc   Name of the thorn
   @vtype   const char *
   @vio     in
   @endvar

   @returntype const char *
   @returndesc
   Name of the implementation or NULL
   @endreturndesc
@@*/
const char *CCTK_ThornImplementation(const char *name)
{
  const t_sktree *node;


  node = SKTreeFindNode(thornlist, name);

  return (node ? ((struct THORN *) node->data)->implementation : NULL);
}

 /*@@
   @routine    CCTK_ImplementationThorn
   @date       Sun Oct 17 22:04:13 1999
   @author     Tom Goodale
   @desc
   Returns the name of one thorn providing an implementation.
   @enddesc

   @var     name
   @vdesc   Name of the implementation
   @vtype   const char *
   @vio     in
   @endvar

   @returntype const char *
   @returndesc
   Name of the thorn or NULL
   @endreturndesc
@@*/
const char *CCTK_ImplementationThorn(const char *name)
{
  const t_sktree *node;


  node = SKTreeFindNode(implist, name);

  return (node ? ((struct IMPLEMENTATION *) node->data)->thornlist->key : NULL);
}


/*@@
   @routine    CCTK_IsThornCompiled
   @date       Sun Jul  4 17:46:56 1999
   @author     Tom Goodale
   @desc
   Checks if a thorn is compiled in.
   @enddesc

   @var     name
   @vdesc   Name of thorn
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - not compiled
   1 - compiled
   @endreturndesc
@@*/
int CCTK_IsThornCompiled(const char *name)
{
  const t_sktree *node;


  node = SKTreeFindNode(thornlist, name);

  return (node != NULL);
}

void CCTK_FCALL CCTK_FNAME(CCTK_IsThornCompiled)
     (int *retval, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(name)
  *retval = CCTK_IsThornCompiled(name);
  free(name);
}


/*@@
   @routine    CCTK_IsImplementationCompiled
   @date       Sun June 3 2001
   @author     Gabrielle Allen
   @desc
   Checks if a implementation is compiled in.
   @enddesc

   @var     name
   @vdesc   Name of implementation
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - not compiled
   1 - compiled
   @endreturndesc
@@*/
int CCTK_IsImplementationCompiled(const char *name)
{
  const t_sktree *node;


  node = SKTreeFindNode(implist, name);

  return (node != NULL);
}

void CCTK_FCALL CCTK_FNAME (CCTK_IsImplementationCompiled)
                           (int *retval, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(name)
  *retval = CCTK_IsImplementationCompiled(name);
  free(name);
}


/*@@
   @routine    CCTK_IsImplementationActive
   @date       Sun Jul  4 17:46:56 1999
   @author     Tom Goodale
   @desc
   Checks if an implementation is active.
   @enddesc

   @var     name
   @vdesc   Name of implementation
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - inactive
   1 - active
   @endreturndesc
@@*/
int CCTK_IsImplementationActive(const char *name)
{
  const t_sktree *node;


  node = SKTreeFindNode(implist, name);

  return (node ? ((struct IMPLEMENTATION *) node->data)->active : 0);
}

void CCTK_FCALL CCTK_FNAME (CCTK_IsImplementationActive)
                           (int *retval, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(name)
  *retval = CCTK_IsImplementationActive(name);
  free(name);
}

 /*@@
   @routine    CCTKi_PrintThorns
   @date       Mon Jul  5 10:02:15 1999
   @author     Tom Goodale
   @desc
   Prints a list of thorns.
   Only lists active ones if the 'active' parameter is true.
   @enddesc

   @var     file
   @vdesc   File stream to print to
   @vtype   FILE *
   @vio     in
   @endvar
   @var     format
   @vdesc   format string for file
   @vtype   const char *
   @vio     in
   @endvar
   @var     active
   @vdesc   Just print active thorns ?
   @vtype   int
   @vio     in
   @vcomment
   Set to 0 to print all thorns.
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
int CCTKi_PrintThorns(FILE *file, const char *format, int active)
{
  int retval;
  const t_sktree *node;


  for(node = SKTreeFindFirst(thornlist), retval = 0;
      node;
      node = node->next, retval++)
  {
    if(((struct THORN *) node->data)->active || !active)
    {
      fprintf(file, format, node->key);
    }
  }

  return retval;
}


 /*@@
   @routine    CCTKi_PrintImps
   @date       Mon Jul  5 10:08:19 1999
   @author     Tom Goodale
   @desc
   Prints a list of implementations.
   Only lists active ones if the 'active' parameter is true.
   @enddesc

   @var     file
   @vdesc   File stream to print to
   @vtype   FILE *
   @vio     in
   @endvar
   @var     format
   @vdesc   format string for file
   @vtype   const char *
   @vio     in
   @endvar
   @var     active
   @vdesc   Just print active implementations ?
   @vtype   int
   @vio     in
   @vcomment
   Set to 0 to print all implementations.
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
int CCTKi_PrintImps(FILE *file, const char *format, int active)
{
  int retval;
  const t_sktree *node;


  for(node = SKTreeFindFirst(implist), retval = 0;
      node;
      node = node->next, retval++)
  {
    if(((struct IMPLEMENTATION *) node->data)->active || !active)
    {
      fprintf(file, format, node->key);
    }
  }

  return retval;
}


 /*@@
   @routine    CCTK_ActivatingThorn
   @date       Thu Oct 14 16:08:42 1999
   @author     Tom Goodale
   @desc
   Finds the thorn which activated a particular implementation
   @enddesc

   @var     name
   @vdesc   implementation name
   @vtype   const char *
   @vio     in
   @endvar

   @returntype const char *
   @returndesc
   Name of activating thorn or NULL if inactive
   @endreturndesc
@@*/
const char *CCTK_ActivatingThorn(const char *name)
{
  const char *retval;
  const t_sktree *node;


  node = SKTreeFindNode(implist, name);
  retval = node && ((struct IMPLEMENTATION *) node->data)->active ?
           ((struct IMPLEMENTATION *) node->data)->activating_thorn : NULL;

  return (retval);
}


 /*@@
   @routine    CCTK_ImpThornList
   @date       Tue Jul 27 09:15:58 1999
   @author     Tom Goodale
   @desc
   Return the thorns for an implementation.
   For now return an sktree - FIXME
   @enddesc

   @var     name
   @vdesc   Name of implementation
   @vtype   const char *
   @vio     in
   @endvar

   @returntype t_sktree *
   @returndesc
   Thornlist
   @endreturndesc
@@*/
t_sktree *CCTK_ImpThornList(const char *name)
{
  const t_sktree *node;


  node = SKTreeFindNode(implist, name);

  return (node ? ((struct IMPLEMENTATION *) node->data)->thornlist : NULL);
}


 /*@@
   @routine    CCTK_NumCompiledThorns
   @date       Tue Feb 02 2000
   @author     Thomas Radke
   @desc
   Return the number of thorns compiled in.
   @enddesc

   @returntype int
   @returndesc
   Number of thorns compiled in
   @endreturndesc
@@*/
int CCTK_NumCompiledThorns(void)
{
  return n_thorns;
}


 /*@@
   @routine    CCTK_CompiledThorn
   @date       Tue Feb 02 2000
   @author     Thomas Radke
   @desc
   Return the name of the compiled thorn with given index.
   @enddesc

   @var     tindex
   @vdesc   thorn index
   @vtype   int
   @vio     in
   @endvar

   @returntype const char *
   @returndesc
   Name of thorn
   @endreturndesc
@@*/
const char *CCTK_CompiledThorn(int tindex)
{
  int i;
  t_sktree *node;


  for(node = SKTreeFindFirst(thornlist), i = 0;
      node;
      node = node->next, i++)
  {
    if (i == tindex)
    {
      return (node->key);
    }
  }

  return (NULL);
}


 /*@@
   @routine    CCTK_NumCompiledImplementations
   @date       Tue Feb 02 2000
   @author     Thomas Radke
   @desc
   Return the number of implementations compiled in.
   @enddesc

   @returntype int
   @returndesc
   Number of implementations compiled in
   @endreturndesc
@@*/
int CCTK_NumCompiledImplementations(void)
{
  return n_imps;
}


 /*@@
   @routine    CCTK_CompiledImplementation
   @date       Tue Feb 02 2000
   @author     Thomas Radke
   @desc
   Return the name of the compiled implementation with given index.
   @enddesc

   @var     tindex
   @vdesc   implementation index
   @vtype   int
   @vio     in
   @endvar

   @returntype const char *
   @returndesc
   Name of implementation
   @endreturndesc
@@*/
const char *CCTK_CompiledImplementation(int tindex)
{
  int i;
  t_sktree *node;


  for(node = SKTreeFindFirst(implist), i = 0;
      node;
      node = node->next, i++)
  {
    if (i == tindex)
    {
      return (node->key);
    }
  }

  return (NULL);
}


 /*@@
   @routine    CCTK_ImplementationRequires
   @date       Sat Oct 20 2001
   @author     Gabrielle Allen
   @desc
   Return the ancestors for an implementation
   @enddesc

   @returntype int
   @returndesc
   @endreturndesc
@@*/

uStringList *CCTK_ImplementationRequires(const char *imp)
{
  int i;
  struct IMPLEMENTATION *impdata;
  uStringList *ancestors;


  impdata = SKTreeFindNode(implist, imp)->data;

  ancestors = Util_StringListCreate(n_thorns);

  /* Get ancestors */
  for(i=0; impdata->ancestors[i]; i++)
  {
    Util_StringListAdd(ancestors,impdata->ancestors[i]);
  }

  /* Get friends */
  for(i=0; impdata->friends[i]; i++)
  {
    Util_StringListAdd(ancestors,impdata->ancestors[i]);
  }

  return ancestors;
}

 /*@@
   @routine    CCTKi_ActivateThorns
   @date       Mon May 21 22:06:37 2001
   @author     Tom Goodale
   @desc
   Activates a list of thorns if they are self consistent.
   @enddesc

   @var     activethornlist
   @vdesc   The list of thorns to activate.
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   -ve Number of errors encountered.
   @endreturndesc
@@*/
int CCTKi_ActivateThorns(const char *activethornlist)
{
  int retval;
  char *local_list;
  uStringList *activated_thorns;
  uStringList *new_thorns;
  uStringList *required_thorns;
  uStringList *requested_imps;
  uStringList *required_imps;
  char *token;
  const char *this_imp;
  int n_warnings;
  int n_errors;
  t_sktree *thornnode;
  t_sktree *impnode;
  t_sktree *impthornlist;

  struct IMPLEMENTATION *imp;
  struct THORN *this_thorn;
  int did_add_thorns;
  int i, j;

  const char *imp1, *imp2;
  const char *this, *thorn, *new_thorn;
  struct iInternalStringList *current;

  local_list = Util_Strdup(activethornlist);

  activated_thorns = Util_StringListCreate(n_thorns);
  required_thorns  = Util_StringListCreate(n_thorns);
  required_imps    = Util_StringListCreate(n_imps);
  requested_imps   = Util_StringListCreate(n_imps);

  printf("Activation requested for \n--->%s<---\n", activethornlist);

  n_errors = 0;
  n_warnings = 0;

  /* Parse list of activated thorns */
  for(token = strtok(local_list, " \r\t\n"); token; token = strtok(NULL," \r\t\n"))
  {
    switch (Util_StringListAdd(activated_thorns, token))
    {
    case 1:
      /* Thorn was added; do nothing */
      break;
    case 0:
      printf("Warning: thorn %s already scheduled for activation\n", token);
      n_warnings++;
      break;
    default:
      CCTK_Warn(0, __LINE__, __FILE__, "Cactus", "Internal error");
    }
  }

  /* Auto-activate some thorns */
  did_add_thorns = 1;
  while(did_add_thorns) {
    new_thorns = Util_StringListCreate(n_thorns);
    /* Copy existing thorns */
    for(thorn = Util_StringListNext(activated_thorns,1); thorn;
        thorn = Util_StringListNext(activated_thorns,0))
    {
      Util_StringListAdd(new_thorns, thorn);
    }
    /* Add all thorns activated by these thorns */
    did_add_thorns = 0;
    for(thorn = Util_StringListNext(activated_thorns,1); thorn;
        thorn = Util_StringListNext(activated_thorns,0))
    {
      thornnode = SKTreeFindNode(thornlist, thorn);
      if(thornnode) {
        this_thorn = thornnode->data;
        if(this_thorn->activates_thorns)
        {
          for(i = 0; this_thorn->activates_thorns[i]; i++)
          {
            new_thorn = this_thorn->activates_thorns[i];
            if (! CCTK_IsThornActive(new_thorn))
            {
              const int ierr = Util_StringListAdd(new_thorns, new_thorn);
              switch (ierr)
              {
              case 0:
                /* Thorn already scheduled for activation */
                break;
              case 1:
                printf("Thorn %s requests automatic activation of %s\n",
                       thorn, new_thorn);
                did_add_thorns = 1;
                break;
              default:
                CCTK_VError(__LINE__, __FILE__, "Cactus", "Internal error: %d",
                            ierr);
              }
            }
          }
        }
      }
    }
    Util_StringListDestroy(activated_thorns);
    activated_thorns = new_thorns;
  }

  for(thorn = Util_StringListNext(activated_thorns,1); thorn;
      thorn = Util_StringListNext(activated_thorns,0))
  {
    if(CCTK_IsThornActive(thorn))
    {
      printf("Warning: thorn %s already active\n", thorn);
      n_warnings++;
    }
    else if (! (this_imp = CCTK_ThornImplementation(thorn)))
    {
      printf("Error: Thorn %s not found\n", thorn);
      n_errors++;
      /*  Give some more help */
      if (CCTK_IsImplementationCompiled(thorn))
      {
        impthornlist = CCTK_ImpThornList(thorn);

        printf("       However, implementation %s was found and is\n",thorn);
        printf("       provided by thorn(s):");
        SKTreeTraverseInorder(impthornlist, JustPrintThornName, NULL);
        printf("\n");
      }
    }
    else if (CCTK_IsImplementationActive(this_imp))
    {
      printf("Error: thorn %s provides implementation %s - already active\n",
             thorn, this_imp);
      n_errors++;
    }
    else if (! Util_StringListAdd(required_thorns,thorn))
    {
      printf("Warning: thorn %s already scheduled for activation\n", thorn);
      n_warnings++;
    }
    else if (! Util_StringListAdd(requested_imps,this_imp))
    {
      printf("Error: thorn %s provides implementation %s which is already "
             "scheduled for activation\n", thorn, this_imp);
      n_errors++;
    }
    else if ((impnode = SKTreeFindNode(implist, this_imp)))
    {
      /* Ok, this thorn exists, and isn't active, a duplicate, or provide
       * the same imp as another thorn which is active or has just been
       * schedule for activation, so can get on with cataloging dependencies.
       */

      Util_StringListAdd(required_imps,this_imp);

      imp = (struct IMPLEMENTATION *)(impnode->data);

      /* Look at ancestors */
      for(i=0; imp->ancestors[i]; i++)
      {
        if(!CCTK_IsImplementationActive(imp->ancestors[i]))
        {
          /* We need this imp */
          Util_StringListAdd(required_imps, imp->ancestors[i]);
        }
      }

      /* Look at friends */
      for(i=0; imp->friends[i]; i++)
      {
        if(!CCTK_IsImplementationActive(imp->friends[i]))
        {
          /* We need this imp */
          Util_StringListAdd(required_imps, imp->friends[i]);
        }
      }
    }
    else
    {
      CCTK_Warn(0, __LINE__, __FILE__, "Cactus", "Internal error :- please "
                "report this to cactusmaint@cactuscode.org");
    }
  }

  /* No longer need the local list */
  free(local_list);

  if(! n_errors)
  {
    /* So, let's see if we are requesting all the imps we need */

    for(imp1=Util_StringListNext(requested_imps,1),
        imp2=Util_StringListNext(required_imps,1);
        imp1&&imp2;
        imp1=Util_StringListNext(requested_imps,0),
        imp2=Util_StringListNext(required_imps,0))
    {
      do
      {
        if(Util_StrCmpi(imp1,imp2))
        {
          printf("Error: Implementation '%s' not activated.\n", imp2);
          printf("       This implementation is required by activated "
                 "thorn(s):\n");

          for (thorn = Util_StringListNext (required_thorns, 1);
               thorn;
               thorn = Util_StringListNext (required_thorns, 0))
          {
            this_imp = CCTK_ThornImplementation (thorn);
            impnode = SKTreeFindNode (implist, this_imp);
            imp = (struct IMPLEMENTATION *) impnode->data;

            /* check ancestors and friends */
            j = 1;
            for (i = 0; imp->ancestors[i] && j; i++)
            {
              j = strcmp (imp2, imp->ancestors[i]);
            }
            for (i = 0; imp->friends[i] && j; i++)
            {
              j = strcmp (imp2, imp->friends[i]);
            }

            if (j == 0)
            {
              printf ("           %s (implementing %s)\n",
                      CCTK_ImplementationThorn (this_imp), this_imp);
            }
          }

          /*  Give some more help */
          if (CCTK_IsImplementationCompiled(imp2))
          {
            impthornlist = CCTK_ImpThornList(imp2);

            printf("       This implementation is provided by compiled "
                   "thorn(s):\n");
            printf("          ");
            SKTreeTraverseInorder(impthornlist, JustPrintThornName, NULL);
            printf("\n");
          }
          else
          {
            printf("       This implementation is not provided by any "
                   "compiled thorn\n");
          }

          printf("       Add a thorn providing this implementation to the "
                 "ActiveThorns parameter.\n");
          n_errors++;
        }
        else
        {
          break;
        }
      } while((imp2=Util_StringListNext(required_imps,0)));
    }
    /* Since the requested imps is a subset of the required imps,
     * we may still have some required imps to go through.
     */
    while(imp2)
    {
      printf("Error: Required implementation %s was not requested\n", imp2);
      printf("       Add a thorn providing this implementation to the ActiveThorns "
             "parameter.\n");
      n_errors++;
      /*  Give some more help */
      if (CCTK_IsImplementationCompiled(imp2))
      {
        impthornlist = CCTK_ImpThornList(imp2);

        printf("       For example, this implementation is provided by compiled thorns:\n");
        printf("          ");
        SKTreeTraverseInorder(impthornlist, JustPrintThornName, NULL);
        printf("\n");
      }
      else
      {
        printf("       This implementation is not provided by any "
               "compiled thorn\n");
      }
      imp2=Util_StringListNext(required_imps,0);
    }
  }

  /* check that all thorns that a thorn requires will also be activated */
  if(! n_errors)
  {
    for (this = Util_StringListNext (required_thorns, 1);
         this;
         this = Util_StringListNext (required_thorns, 0))
    {
      /* cannot recursively browse through string list
         so we have to remember the current entry for the outer loop */
      current = required_thorns->current;
      this_thorn = ((t_sktree *) SKTreeFindNode (thornlist, this))->data;
      for (i = 0; this_thorn->requires_thorns &&
                  this_thorn->requires_thorns[i]; i++)
      {
        for (thorn = Util_StringListNext (required_thorns, 1);
             thorn;
             thorn = Util_StringListNext (required_thorns, 0))
        {
          if (! Util_StrCmpi (thorn, this_thorn->requires_thorns[i]) ||
              CCTK_IsThornActive (this_thorn->requires_thorns[i]))
          {
            break;
          }
        }
        if (! thorn)
        {
          printf ("Error: Thorn %s requires thorn %s to be active.\n"
                  "       Please add this thorn to the ActiveThorns "
                  "parameter.\n", this, this_thorn->requires_thorns[i]);
          n_errors++;
        }
      }
      required_thorns->current = current;
    }
  }

  if(! n_errors)
  {
    /* Ok, so we have all required imps, so can activate the thorns, finally */

    for(thorn = Util_StringListNext(required_thorns, 1);
        thorn;
        thorn = Util_StringListNext(required_thorns, 0))
    {
      ActivateThorn(thorn);
    }

    /* Now do any necessary parameter activation */
    for(thorn = Util_StringListNext(required_thorns, 1);
        thorn;
        thorn = Util_StringListNext(required_thorns, 0))
    {
      CCTKi_ParameterActivateThornParameters(thorn);
    }

    retval = 0;
  }
  else
  {
    printf("Activation failed - %d errors in activation sequence\n", n_errors);
    retval = -n_errors;
  }

  Util_StringListDestroy(required_thorns);
  Util_StringListDestroy(required_imps);
  Util_StringListDestroy(requested_imps);

  return retval;
}



/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    RegisterImp
   @date       Sun Jul  4 17:44:42 1999
   @author     Tom Goodale
   @desc
   Registers an implementation.
   @enddesc

   @var     name
   @vdesc   name of the implementation
   @vtype   const char *
   @vio     in
   @endvar
   @var     thorn
   @vdesc   name of the thorn
   @vtype   const char *
   @vio     in
   @endvar
   @var     ancestors
   @vdesc   ancestors of the implementation
   @vtype   const char **
   @vio     in
   @endvar
   @var     friends
   @vdesc   friends of the implementation
   @vtype   const char **
   @vio     in
   @endvar

   @returntype int
   @returndesc
    0 - success
   -1 - failed to store thorn in implementation
   -2 - memory failure creating implementation
   -3 - failed to store implementtion in tree
   @endreturndesc
@@*/
static int RegisterImp(const char *name,
                       const char *thorn,
                       const char **ancestors,
                       const char **friends)
{
  int retval;
  int count;
  t_sktree *node;
  t_sktree *temp;

  struct IMPLEMENTATION *imp;

  /* Does the implementation already exist ? */
  node = SKTreeFindNode(implist, name);
  if(!node)
  {
    n_imps++;

    /* Create the structure to hold info about it. */
    imp = malloc(sizeof(struct IMPLEMENTATION));
    if(imp)
    {
      imp->active = 0;

      /* Store the name of this thorn in a tree */
      imp->thornlist = SKTreeStoreData(NULL,NULL, thorn, NULL);

      /* Store the info in the tree. */
      temp = SKTreeStoreData(implist, implist, name, imp);

      if(!implist) implist = temp;

      retval = temp ? 0 : -3;
      if(!retval)
      {
        /* Count the ancestors */
        for(count=0; ancestors && ancestors[count];count++);

        imp->n_ancestors = count;
        imp->ancestors = malloc((count+1)*sizeof(char *));

        if(imp->ancestors)
        {
          for(count=0; count < imp->n_ancestors; count++)
          {
            imp->ancestors[count] = Util_Strdup(ancestors[count]);
          }
          imp->ancestors[count] = NULL;

          qsort(imp->ancestors, count, sizeof(char *), CompareStrings);

        }

        /* Count the friends */
        for(count=0; friends && friends[count];count++);

        imp->n_friends = count;
        imp->friends = malloc((count+1)*sizeof(char *));

        if(imp->friends)
        {
          for(count=0; count < imp->n_friends; count++)
          {
            imp->friends[count] = Util_Strdup(friends[count]);
          }
          imp->friends[count] = NULL;

          qsort(imp->friends, count, sizeof(char *), CompareStrings);
        }
      }
    }
    else
    {
      retval = -2;
    }
  }
  else
  {
    imp = node->data;
    SKTreeStoreData(imp->thornlist,imp->thornlist, thorn, NULL);

    retval = -1;
  }

  return retval;
}


 /*@@
   @routine    ActivateThorn
   @date       Mon May 21 22:09:47 2001
   @author     Tom Goodale
   @desc
   Activate one thorn - assumes all error checking done by calling routine.
   @enddesc

   @var     name
   @vdesc   Name of thorn to activate
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0  - success
   -1 - can't find thorn
   @endreturndesc
@@*/
static int ActivateThorn(const char *name)
{
  int retval;
  t_sktree *thornnode;

  struct THORN *thorn;

  printf("Activating thorn %s...", name);

  /* Find the thorn */
  thornnode = SKTreeFindNode(thornlist, name);
  if(thornnode)
  {
    thorn = thornnode->data;

    thorn->active = 1;

    printf("Success -> active implementation %s\n", thorn->implementation);

    retval = ActivateImp(thorn->implementation, name);
  }
  else
  {
    retval = -1;
  }

  return retval;
}


 /*@@
   @routine    ActivateImp
   @date       Mon May 21 22:09:47 2001
   @author     Tom Goodale
   @desc
   Activate one implementation - assumes all error checking done by calling routine.
   @enddesc

   @var     implementation
   @vdesc   Name of implementation to activate
   @vtype   const char *
   @vio     in
   @endvar
   @var     thorn
   @vdesc   Name of thorn activating this imp
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0  - success
   -1 - can't find implementation
   @endreturndesc
@@*/
static int ActivateImp(const char *implementation, const char *thorn)
{
  t_sktree *impnode;
  struct IMPLEMENTATION *imp;

  /* Find the implementation */
  impnode = SKTreeFindNode(implist, implementation);
  if(impnode)
  {
    imp = impnode->data;

    imp->active = 1;
    /* Remember which thorn activated this imp. */
    imp->activating_thorn = Util_Strdup(thorn);
  }

  return (impnode ? 0 : -1);
}

/*@@
   @routine    CompareStrings
   @date       Thu Sep 14 18:57:52 2000
   @author     Tom Goodale
   @desc
   Case independent string comparison to pass to qsort.
   @enddesc
@@*/
static int CompareStrings(const void *string1, const void *string2)
{
  return Util_StrCmpi(*(const char *const *)string1, *(const char *const *)string2);
}

 /*@@
   @routine    JustPrintThornName
   @date       Mon Jun  4 19:05:45 2001
   @author     Tom Goodale
   @desc
   Print the name of a thorn if it is passed from an sktree.
   @enddesc
@@*/
static int JustPrintThornName(const char *key, void *input, void *dummy)
{
  input = input;
  dummy = dummy;

  printf(" %s", key);

  return 0;
}
