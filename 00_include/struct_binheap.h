#ifndef _BinHeap_H
#define _BinHeap_H

#ifdef	__cplusplus
extern "C" {
#endif

typedef struct {
    double rmsd;
    int triple_x[3];
    int triple_y[3];
    double quat[4];
} Triple;

typedef Triple ElementType;


struct HeapStruct;
typedef struct HeapStruct *PriorityQueue;

PriorityQueue Initialize(int MaxElements);
void Destroy(PriorityQueue H);
void MakeEmpty(PriorityQueue H);
void Insert(ElementType X, PriorityQueue H);
ElementType DeleteMax(PriorityQueue H);
ElementType FindMax(PriorityQueue H);
int IsEmpty(PriorityQueue H);
int IsFull(PriorityQueue H);

#ifdef	__cplusplus
}
#endif

#endif

