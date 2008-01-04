#include "common.h"
#include "cluster.h"
#include "stdlib.h"
#include "math.h"
#include "stdio.h"

/* Allocates a new clusternode.  Lines is an arbitrary pointer, and
   will hold a python list object most likely */
struct ClusterNode* clusterNodeAlloc( int start, int end)
{
  struct ClusterNode* cn;
  AllocVar( cn );
  cn->start = start;
  cn->end = end;
  cn->left = NULL;
  cn->right = NULL;
  cn->linenums = NULL;
  cn->regions = 0;

  // calculate coin-flip priority
  double uniform = ((double)rand()) / (RAND_MAX);
  if (uniform == 1.0)
    uniform = 0;
  cn->priority = (int) ceil( (-1.0 / log(.5)) * log( -1.0 / (uniform - 1)));
  return cn;
}

/* Insert a new range into node cn, return reference to node inserted
   at to allow cn->lines to be changed. */
struct ClusterNode* clusterNodeInsert( struct ClusterNode** cn, int start, int end, int linenum, int mincols)
{
  if( (*cn) == NULL ) {
    (*cn) = clusterNodeAlloc( start, end );
  }
  struct ClusterNode* root = *cn;
  if( start - mincols > root->end ) {
    clusterNodeInsert( &(root->right), start, end, linenum, mincols );
    if( root->priority < root->right->priority ) {
      clusterRotateLeft( cn );
    }
    return root->right;
  } else if( end + mincols < root->start) {
    clusterNodeInsert( &(root->left), start, end, linenum, mincols );
    if( root->priority < root->left->priority ) {
      clusterRotateRight( cn );
    }
    return root->left;
  } else {
    root->regions++;
    int oldstart = root->start;
    int oldend = root->end;
    root->start = min(start, root->start);
    root->end = max(end, root->end);
    append(&(root->linenums), linenum);
    if( oldstart != root->start && root->left != NULL ) {
      clusterPushUp( &(root->left), cn, mincols );
    }
    if( oldend != root->end && root->right != NULL ) {
      clusterPushUp( &(root->right), cn, mincols );
    }
    return root;
  }
}

struct ClusterNode* clusterNodeFind(struct ClusterNode *cn, int position)
{
  if(cn == NULL) return NULL;
  if(position < cn->start) {
    if(cn->left)
      return clusterNodeFind(cn->left, position);
    else
      return NULL;
  } else if(position > cn->end) {
    if(cn->right)
      return clusterNodeFind(cn->right, position);
    else
      return NULL;
  } else
    return cn;
}

void clusterRotateRight( struct ClusterNode **cn )
{
  struct ClusterNode* root = (*cn)->left;
  (*cn)->left = (*cn)->left->right;
  root->right = (*cn);
  *cn = root;
}

void clusterRotateLeft( struct ClusterNode **cn )
{
  struct ClusterNode* root = (*cn)->right;
  (*cn)->right = (*cn)->right->left;
  root->left = (*cn);
  *cn = root;
}

void clusterPushUp( struct ClusterNode **ln, struct ClusterNode **cn, int mincols )
{
  struct ClusterNode* local = *ln;
  struct ClusterNode* root = *cn;
  int maxstart = max(root->start, local->start);
  int maxend = max(local->end, root->end);
  int minstart = min(root->start, local->start);
  int minend = min(root->end, local->end);
  int distance = maxstart - minend;

  if( distance <= mincols ) {
    // merge nodes, free merged children
    root->start = minstart;
    root->end = maxend;
    merge(root->linenums, local->linenums);
    root->regions += local->regions;
    if( local->right) clusterPushUp(&(local->right), cn, mincols);
    if( local->left) clusterPushUp(&(local->left), cn, mincols);
    if((local->right == NULL) && (local->left == NULL)) {
      free(local->linenums);
      free(local);
      *ln = NULL;
    } else if(local->right) {
      *ln = local->right;
      free(local->linenums);
      free(local);
    } else if (local->left) {
      *ln = local->left;
      free(local->linenums);
      free(local);
    } else {
      // Impossible
    }
    return;
  }

  if( local->end < root->start && local->left ) {
    clusterPushUp(&(local->left), cn, mincols );
  }
  if( local->start > root->end && local->right ) {
    clusterPushUp(&(local->right), cn, mincols );
  }
}

inline int min(int a, int b)
{
  if( a < b )
    return a;
  else
    return b;
}

inline int max(int a, int b)
{
  if( a > b )
    return a;
  else
    return b;
}

/* Tree iterator.  Dump references into a single-linked list on
   get_itr(ClusterNode), walk the list with next(treeitr) and
   has_next(treeitr). */

/* Infix, in reverse because we stack. */
void get_itr_in(struct ClusterNode* root, struct treeitr** itr)
{
  if( root == NULL ) return;
  struct treeitr* temp;
  get_itr_in(root->right, itr);
  AllocVar( temp );
  temp->next = (*itr);
  temp->value = root;
  (*itr) = temp;
  get_itr_in(root->left, itr);
}

/* Postfix, looks like mirrored prefix. */
void get_itr_post(struct ClusterNode* root, struct treeitr** itr)
{
  if( root == NULL ) return;
  struct treeitr* temp;
  AllocVar( temp );
  temp->next = (*itr);
  temp->value = root;
  (*itr) = temp;
  get_itr_post(root->right, itr);
  get_itr_post(root->left, itr);
}

/* Prefix, looks like mirrored postfix. */
void get_itr_pre(struct ClusterNode* root, struct treeitr** itr)
{
  if( root == NULL) return;
  struct treeitr* temp;
  get_itr_pre(root->right, itr);
  get_itr_pre(root->left, itr);
  AllocVar( temp );
  temp->next = (*itr);
  temp->value = root;
  (*itr) = temp;
}

/* next( itr )*/
struct ClusterNode* next(struct treeitr** itr)
{
  if( (*itr) == NULL) return NULL;
  // get pointer to clusternode
  struct ClusterNode* temp = (*itr)->value;
  // get pointer to next
  struct treeitr* nextitr = (*itr)->next;
  // free old
  free(*itr);
  // repoint
  (*itr) = nextitr;
  
  return temp;
}

/* has_next(itr).  Returns a boolean value indicating if there is
   another item to be iterated.  Mind you it fakes it by evaluating
   whether or not the pointer has been deallocated and zeroed.*/
int has_next(struct treeitr **itr) 
{
  return ((*itr) != NULL);
}

/* Linked list functions. Our linked list has the advantage of only
   needing to accept writes, reads, and clear.  No removal is
   needed.*/
void append(struct linelist** ptr, int value) {
  struct linelist *l = *ptr;
  if( l == NULL ) {
    // Allocate a new list
    AllocVar( l );
    *ptr = l;
    l->head = NULL;
    l->tail = NULL;
  }
  struct listitem* newitem;
  AllocVar( newitem );
  newitem->value = value;
  newitem->next = NULL;
  if(l->head==NULL && l->tail==NULL) {
    l->head = newitem;
    l->tail = newitem;
  } else {
    l->tail->next = newitem;
    l->tail = newitem;
  }
}

void merge(struct linelist* left, struct linelist* right)
{
  if(left->tail && right->head) {
    left->tail->next = right->head;
    left->tail = right->tail;
  } else if(left->tail == NULL) {
    left->head = right->head;
    left->tail = right->tail;
  }
}

void freelist(struct linelist** ptr)
{
  if( (*ptr) == NULL) return;
  struct linelist* l = *ptr;
  struct listitem* itemptr = l->head;
  struct listitem* nextptr;
  while(itemptr != NULL) {
    nextptr = itemptr->next;
    free(itemptr);
    itemptr = nextptr;
  }
  l->head = NULL;
  l->tail = NULL;
  l = NULL;
}

void freetree(struct ClusterNode** cn) {
  if( (*cn) == NULL) return;
  struct ClusterNode* root = *cn;
  if(root) {
    if(root->left) {
      freetree(&(root->left));
    }
    if(root->right) {
      freetree(&(root->right));
    }
    if(root->linenums)
      freelist(&(root->linenums));
    free(root);
  }
}

void dumpTree(struct ClusterNode *cn)
{
  if(cn->left != NULL) {
    dumpTree(cn->left);
  }
  printf("Start: %d End: %d \n",cn->start, cn->end);
  if(cn->linenums && cn->linenums->head)
    printf("Head linenum: %d\n",cn->linenums->head->value);
  if(cn->right != NULL) {
    dumpTree(cn->right);
  }
}

int main();
int main() {
  struct ClusterNode* root = NULL;
  int x;
  for(x=1; x<50; x++) {
    int start = ((double) rand()) / RAND_MAX * 512;
    int end = start + 10;
    clusterNodeInsert(&root, start, end, x, 1);
  }
  dumpTree(root);
  struct listitem* listptr = root->linenums->head;
  while(listptr) {
    printf("%d ",listptr->value);
    listptr=listptr->next;
  }
  struct treeitr* myitr = NULL;
  get_itr_in(root, &myitr);
  while( has_next(&myitr) )
    printf("%d\n", (next(&myitr))->start);
  freetree(&root);
  return 0;
}
