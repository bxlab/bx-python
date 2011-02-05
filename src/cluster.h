typedef struct struct_interval {
    int start;
    int end;
    int id;
    
    struct struct_interval *next;
} interval;

typedef struct struct_clusternode {
    int start;
    int end;
    int priority;
    
    struct struct_interval *interval_head;
    struct struct_interval *interval_tail;
    int num_ivals;
    
    struct struct_clusternode *left;
    struct struct_clusternode *right;
} clusternode;

typedef struct {
    int max_dist;
    int min_intervals;
    
    clusternode *root;
} clustertree;

typedef struct struct_treeitr {
    struct struct_treeitr *next;
    struct struct_clusternode *node;
} treeitr;


clusternode* clusternode_insert(clustertree *tree, clusternode *node, int start, int end, int id);
clustertree* create_clustertree(int max_dist, int min_intervals);
treeitr* clusteritr(clustertree *tree);
void freeclusteritr(treeitr *itr);
void free_tree(clustertree *tree);
