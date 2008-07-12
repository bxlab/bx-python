#include <stdlib.h>
/*
  Author: Ian N Schenck
  Version: 7/21/2006

  Most of this was ripped out of James Taylor's never-released code,
  and plugged in here for use in Python.  Slight modifications were
  made where I saw fit.

  It looks as if CpG's are typically not next to gaps.
*/

static inline int is_cpg( char * sp1, char * sp2, int pos)
{
  if ( pos < 1 ) return 0;
  if ( sp1[pos + 1] == '\0' ) return 0;
  if ( sp1[pos - 1] != 'C' && sp2[pos - 1] != 'C' &&
       sp1[pos + 1] == 'G' && sp2[pos + 1] == 'G' &&
       (sp1[pos] == 'C' || sp2[pos] == 'C') ) return 1;
  if ( sp1[pos + 1] != 'G' && sp2[pos + 1] != 'G' &&
       sp1[pos - 1] == 'C' && sp2[pos - 1] == 'C' &&
       (sp1[pos] == 'G' || sp2[pos] == 'G') ) return 1;
  return 0;
}

static inline int is_non_cpg( char * sp1, char * sp2, int pos)
{
  // first one can't assuredly be cpg
  if ( pos < 1 ) return 1;
  if ( sp1[pos + 1] == '\0' ) return 0;
  return
    ( sp1[pos - 1] != 'C' && sp2[pos - 1] != 'C' &&
      sp1[pos + 1] != 'G' && sp2[pos + 1] != 'G' );
}

static inline int is_cpg_restricted( char * sp1, char * sp2, int pos )
{
  return !is_non_cpg( sp1, sp2, pos );
}

int next( char * sp1, char * sp2, int start, int (*func)(char*,char*,int))
{
  while( sp1[start+1] != '\0')
    {
      if( func(sp1, sp2, start) )
	return start;
      start++;
    }
  // nothing found
  return -1;
}

int next_cpg( char * sp1, char * sp2, int start)
{
  return next( sp1, sp2, start, &is_cpg);
}

int next_cpg_restricted( char * sp1, char *sp2, int start)
{
  return next( sp1, sp2, start, &is_cpg_restricted );
}

int next_non_cpg( char * sp1, char * sp2, int start)
{
  return next( sp1, sp2, start, &is_non_cpg);
}
