#ifndef CLUSTER_H
#define CLUSTER_H

struct ClusterNode
{
  int start;
  int end;
  int priority;
  int regions;
  struct ClusterNode *left;
  struct ClusterNode *right;
  struct linelist *linenums;
};

struct listitem
{
  int value;
  struct listitem *next;
};

struct linelist
{
  struct listitem *head;
  struct listitem *tail;
};

struct treeitr
{
  struct treeitr * next;
  struct ClusterNode * value;
};

struct ClusterNode* clusterNodeAlloc( int start, int end );
struct ClusterNode* clusterNodeInsert( struct ClusterNode** cn, int start, int end, int linenum, int mincols);
void clusterPushUp( struct ClusterNode **ln, struct ClusterNode **cn, int mincols );
void clusterRotateRight( struct ClusterNode **cn );
void clusterRotateLeft( struct ClusterNode **cn );
struct ClusterNode* clusterNodeFind( struct ClusterNode *cn, int position );

inline int max(int a, int b);
inline int min(int a, int b);

void get_itr_pre(struct ClusterNode*, struct treeitr**);
void get_itr_in(struct ClusterNode*, struct treeitr**);
void get_itr_post(struct ClusterNode*, struct treeitr**);
struct ClusterNode* next(struct treeitr**);
int has_next(struct treeitr**);

void append(struct linelist** ptr, int value);
void merge(struct linelist* left, struct linelist* right);
void freelist(struct linelist** ptr);
void freetree(struct ClusterNode** cn);

void dumpTree(struct ClusterNode *cn);

#endif
