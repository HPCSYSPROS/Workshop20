 /*@@
   @file      ParameterBindings.h
   @date      Mon Feb 22 15:01:14 1999
   @author    Tom Goodale
   @desc 
   Defines for parameter stuff
   @enddesc 
   @version $Header$
 @@*/

#ifndef _PARAMETERBINDINGS_H_
#define _PARAMETERBINDINGS_H_


#ifdef __cplusplus
extern "C" {
#endif

int CCTKi_ParameterCreate(const char *name,        /* The parameter name */
                          const char *thorn,       /* The thorn          */ 
                          const char *type,        /* The parameter type */
                          const char *scope,       /* The scoping block  */
                          int steerable,           /* Is it steerable ?  */
                          const char *description, /* The description    */ 
                          const char *defval,      /* The default value  */ 
                          void *datapointer,       /* The actual data    */
                          int array,               /* Is it an array, if so what size */ 
                          const char *accumulator, /* Is it an accumultor parameter */
                          int n_ranges,            /* How many allowed ranges it has */
                          ...);

int CCTKi_ParameterAddRange(const char *implementation, 
                      const char *name,
                      const char *range_origin,
                      const char *range,
                      const char *range_description);

void CCTKi_ParameterAccumulatorBase(const char *thorn,
                                    const char *parameter,
                                    const char *importhorn,
                                    const char *acc_base);

#ifdef __cplusplus
}
#endif

/* These definitions must be identical to the ones in cctk_Parameter.h */
#define PARAMETER_FIRST    701
#define PARAMETER_KEYWORD  701
#define PARAMETER_STRING   702
#define PARAMETER_SENTENCE 703
#define PARAMETER_INT      704
#define PARAMETER_INTEGER  704
#define PARAMETER_REAL     705
#define PARAMETER_BOOLEAN  706

#endif
