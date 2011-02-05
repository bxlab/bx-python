/*
    Kanwei Li, 2009
    Inspired by previous ClusterTree
    
    This clustering algorithm uses a binary tree structure. Nodes correspond to 
    non-overlapping intervals, where overlapping means that the distance between
    two intervals is less or equal to max_dist, which is the max separation.
    
    The tree self-balances using rotations based on the binomial sequence. Merges
    among nodes are performed whenever a node is changed/added that will cause other
    nodes to form a new cluster.
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cluster.h"

#define ALLOC(pt) (malloc(sizeof(pt)))

static int min(int a, int b) {
    if( a < b )
        return a;
    else
        return b;
}

static int max(int a, int b) {
    if( a > b )
        return a;
    else
        return b;
}

/* Create new tree with given max_dist (max distance between intervals to be
    considered a cluster), and min_intervals, the minimum number of intervals
    needed for a cluster to be considered significant */
clustertree* create_clustertree(int max_dist, int min_intervals) {
    clustertree *tree = ALLOC(clustertree);
    tree->max_dist = max_dist;
    tree->min_intervals = min_intervals;
    tree->root = NULL;
    return tree;
}

static interval* create_interval(int start, int end, int id) {
    interval *ival = ALLOC(interval);
    
    ival->start = start;
    ival->end = end;
    ival->id = id;
    ival->next = NULL;
    return ival;
}

static clusternode* create_node(int start, int end, int id) {
    clusternode *new_node = ALLOC(clusternode);
    
    new_node->start     = start;
    new_node->end       = end;
    new_node->interval_head = create_interval(start, end, id);
    new_node->interval_tail = new_node->interval_head;
    new_node->num_ivals = 1;
    new_node->left      = NULL;
    new_node->right     = NULL;
    
    double uniform = ((double)rand()) / (RAND_MAX);
    if (uniform == 1.0)
        uniform = 0;
    new_node->priority = (int)ceil( (-1.0 / log(.5)) * log( -1.0 / (uniform - 1)));
    
    return new_node;
}

static void recursively_free_intervals(interval *ival) {
    interval *next;
    if(ival) {
        next = ival->next;
        free(ival);
        recursively_free_intervals(next);
    }
}

static void recursively_free_nodes(clusternode *node) {
    if(node) {
        recursively_free_nodes(node->left);
        recursively_free_nodes(node->right);
        recursively_free_intervals(node->interval_head);
        free(node);
    }
}

void free_tree(clustertree *tree) {
    recursively_free_nodes(tree->root);
    free(tree);
}

void cluster_rotateright(clusternode **node) {
    clusternode* root = (*node)->left;
    (*node)->left = (*node)->left->right;
    root->right = (*node);
    *node = root;
}

void cluster_rotateleft(clusternode **node) {
    clusternode* root = (*node)->right;
    (*node)->right = (*node)->right->left;
    root->left = (*node);
    *node = root;
}

/* Go down the tree and merge nodes if necessary */
void cluster_fixup(clustertree *tree, clusternode **ln, clusternode **rn) {
    clusternode* local = *ln;
    clusternode* root = *rn;
    int maxstart = max(root->start, local->start);
    int maxend = max(local->end, root->end);
    int minstart = min(root->start, local->start);
    int minend = min(root->end, local->end);

    if( maxstart - minend <= tree->max_dist ) {
        /* Have to merge this node and children */
        root->start = minstart;
        root->end = maxend;
        root->interval_tail->next = local->interval_head;
        root->interval_tail = local->interval_tail;
        root->num_ivals += local->num_ivals;
        if( local->right) cluster_fixup(tree, &(local->right), rn);
        if( local->left) cluster_fixup(tree, &(local->left), rn);
        if((local->right == NULL) && (local->left == NULL)) {
            free(local);
            *ln = NULL;
        } else if(local->right) {
            *ln = local->right;
            free(local);
        } else if (local->left) {
            *ln = local->left;
            free(local);
        }
        return;
    }
    // Even if we miss, we still have to check children
    if(local->left) {
        cluster_fixup(tree, &(local->left), rn);
    }
    if(local->right) {
        cluster_fixup(tree, &(local->right), rn);
    }
}

/* Pyrex "getregions" implements this. Only used for C debugging */
void clustereach(clustertree *tree, clusternode *node) {
    interval* ival;
    if (node == NULL) {
        exit(1); /* Shouldn't happen */
    }
    if (node->left != NULL) {
        clustereach(tree, node->left);
    }
    printf("Node: %d\t%d\n", node->start, node->end);
    ival = node->interval_head;
    while(ival) {
        printf("\tInterval %d: %d\t%d\n", ival->id, ival->start, ival->end);
        ival = ival->next;
    }
    
    if (node->right != NULL) {
        clustereach(tree, node->right);
    }
}

void clusteritr_recursive(clustertree *tree, clusternode *node, treeitr* *itr) {
    treeitr *newitr;

    if (node == NULL) {
        return;
    }
    if (node->right != NULL) {
        clusteritr_recursive(tree, node->right, itr);
    }
    if (node->num_ivals >= tree->min_intervals) {
        newitr = ALLOC(treeitr);
        newitr->next = *itr;
        newitr->node = node;
        *itr = newitr;
    }
    if (node->left != NULL) {
        clusteritr_recursive(tree, node->left, itr);
    }
}

/* Create an infix iterator */
treeitr* clusteritr(clustertree *tree) {
    treeitr *itr = NULL;
    
    clusteritr_recursive(tree, tree->root, &itr);
    if (itr != NULL) {
        return itr;
    }
    return NULL;
}

/* Free iterator (tail recursive) */
void freeclusteritr(treeitr *itr) {
    treeitr *next;
    if (itr == NULL) {
        return;
    }
    
    next = itr->next;
    free(itr);
    freeclusteritr(next);
}

/* Insert based on the start position of intervals */
clusternode* clusternode_insert(clustertree *tree, clusternode *node, int start, int end, int id) {
    int oldstart;
    int oldend;
    interval* ival;
    
    // printf("Inserting %d %d %d\n", start, end, id);
    if (node == NULL) {
        node = create_node(start, end, id);
        
    } else if ( (start - tree->max_dist) > node->end ) { /* We're to the right of this cluster */
        node->right = clusternode_insert(tree, node->right, start, end, id);
        if (node->priority < node->right->priority) cluster_rotateleft(&node);
    
    } else if ( (end + tree->max_dist) < node->start) { /* We're to the left of this cluster */
        node->left = clusternode_insert(tree, node->left, start, end, id);
        if (node->priority < node->left->priority) cluster_rotateright(&node);
                        
    } else { /* We're in the range of this cluster */
        /* Update the start and end to match to new values */
        oldstart    = node->start;
        oldend      = node->end;
        node->start = min(start, node->start);
        node->end   = max(end, node->end);
        ival = create_interval(start, end, id);
        ival->next = node->interval_head; /* Add this interval as the head of the interval list */
        node->interval_head = ival;
        node->num_ivals += 1;
                
        if ( oldstart > node->start && node->left != NULL ) { /* New interval added to the start, and there's a left child */
            cluster_fixup(tree, &(node->left), &node);
        }
        if ( oldend < node->end && node->right != NULL ) { /* New interval added to the end, and there's a right child */
            cluster_fixup(tree, &(node->right), &node);
        }
    }
    return node;
}

int main() {
    
    // Simple test
    clustertree* tree = create_clustertree(0, 1);
    
    tree->root = clusternode_insert(tree, tree->root, 3, 4, 0);
    tree->root = clusternode_insert(tree, tree->root, 6, 7, 1);
    tree->root = clusternode_insert(tree, tree->root, 9, 10, 2);
    tree->root = clusternode_insert(tree, tree->root, 1, 2, 3);
    tree->root = clusternode_insert(tree, tree->root, 3, 8, 4);
    
    clustereach(tree, tree->root);
    return 0;
    
}
