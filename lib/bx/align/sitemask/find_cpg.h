#ifndef __find_cpg__
#define __find_cpg__
inline int is_cpg( char * sp1, char * sp2, int pos);
inline int is_non_cpg( char * sp1, char * sp2, int pos);
inline int is_cpg_restricted( char * sp1, char * sp2, int pos );
int next( char * sp1, char * sp2, int start, int (*func)(char*,char*,int));
int next_cpg( char * sp1, char * sp2, int start);
int next_cpg_restricted( char * sp1, char *sp2, int start);
int next_non_cpg( char * sp1, char * sp2, int start);
#endif
