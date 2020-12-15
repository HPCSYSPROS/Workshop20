 /*@@
   @header    util_StringList.h
   @date      Mon May 21 16:58:15 2001
   @author    Tom Goodale
   @desc 
   Prototypes and data structures for fixed-size lists of unique strings.
   @enddesc
   @version $Header$
 @@*/

#ifndef _UTIL_STRINGLIST_H_
#define _UTIL_STRINGLIST_H_ 1

#ifdef __cplusplus
extern "C" 
{
#endif
struct iInternalStringList
{
  struct iInternalStringList *next;
  char *string;
};

typedef struct
{
  int max_size;
  int fill;
  struct iInternalStringList *list;
  struct iInternalStringList *head;
  struct iInternalStringList *current;
} uStringList;

uStringList *Util_StringListCreate(int size);
int Util_StringListAdd(uStringList *list, const char *item);
void Util_StringListDestroy(uStringList *list);
const char *Util_StringListNext(uStringList *list, int flag);


#ifdef __cplusplus
}
#endif

#endif /* _UTIL_STRINGLIST_H_ */
