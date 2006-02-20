#include "common.h"

void *needMem(size_t size)
        /* Need mem calls abort if the memory allocation fails. The memory
         *  * is initialized to zero. */
{
        void *pt;
        if ((pt = malloc(size)) == NULL)
        {        
                    fprintf( stderr, "Out of memory needMem - request size %llu bytes\n",
                                                 (unsigned long long)size);
                    exit(1);
        }            
        memset(pt, 0, size);
        return pt;
}

void freeMem(void *pt)
        /* Free memory will check for null before freeing. */
{
        if (pt != NULL)
                    free(pt);
}


void *needLargeZeroedMem(size_t size)
        /* Request a large block of memory and zero it. */
{
        void *v;
        /*v = needLargeMem(size);*/
        v = malloc(size);
        memset(v, 0, size);
        return v;
}

void freez(void *vpt)
        /* Pass address of pointer.  Will free pointer and set it 
         *  * to NULL. */
{
        void **ppt = (void **)vpt;
        void *pt = *ppt;
        *ppt = NULL;
        freeMem(pt);
}

/* fill a specified area of memory with zeroes */
void zeroBytes(void *vpt, int count)
{
        char *pt = (char*)vpt;
        while (--count>=0)
                    *pt++=0;
}
