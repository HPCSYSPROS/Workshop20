 /*@@
   @header    util_Expression.h
   @date      Tue Sep 19 22:02:45 2000
   @author    Tom Goodale
   @desc 
   Header for expression stuff.
   @enddesc
   @version $Header$
 @@*/

#ifndef __UTIL_EXPRESSION_H__
#define __UTIL_EXPRESSION_H__ 1

#ifdef __cplusplus
extern "C" 
{
#endif

  /* Structure to hold values. */
typedef struct 
{
  enum {rval,ival,sval} type;
  
  union 
  {
    CCTK_REAL rval;
    CCTK_INT  ival;
    const char *sval;
  } value;
} uExpressionValue;

#ifndef __UTILI_EXPRESSION_H__
  /* Externally visible representation of the expression. */
typedef void *uExpression;
#endif /*__UTIL_EXPRESSION_H__ */

uExpression Util_ExpressionParse(const char *expression);

int Util_ExpressionEvaluate(const uExpression buffer,
                            uExpressionValue *retval,
                            int (*eval)(int, const char * const *, uExpressionValue *, const void *),
                            const void *data);

void Util_ExpressionFree(uExpression buffer);
typedef int (*uExpressionEvaluator)(int, const char * const *, uExpressionValue *, const void *);

#ifdef __cplusplus
}
#endif

#endif /* __UTIL_EXPRESSION_H__ */
