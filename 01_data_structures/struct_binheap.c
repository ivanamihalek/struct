#include "struct.h"

#define MinPQSize (10)
#define MinData (-32767) // napraviti za nasu potrebu

struct HeapStruct {
    int Capacity;
    int Size;
    ElementType *Elements;
};

PriorityQueue Initialize(int MaxElements) {
    PriorityQueue H;

    /* 1*/ if (MaxElements < MinPQSize) {
        /* 2*/ printf("Priority queue size is too small");
               exit(1);
           }
    /* 3*/ H = malloc(sizeof ( struct HeapStruct));
    /* 4*/ if (H == NULL){
        /* 5*/ printf("Out of space!!!");
               exit(1);
           }
    /* Allocate the array plus one extra for sentinel */
    /* 6*/ H->Elements = malloc((MaxElements + 1)
            * sizeof ( ElementType));
    /* 7*/ if (H->Elements == NULL){
        /* 8*/ printf("Out of space!!!");
               exit(1);
           }
    /* 9*/ H->Capacity = MaxElements;
    /*10*/ H->Size = 0;
    /*11*/ H->Elements[ 0 ].rmsd = 10000;
           H->Elements[0].triple_x[0] = -1;
           

    /*12*/ return H;
}

/* END */

void MakeEmpty(PriorityQueue H) {
    H->Size = 0;
}

/* H->Element[ 0 ] is a sentinel */

void Insert(ElementType X, PriorityQueue H) {
    int i;

    
    
    if (X.rmsd > FindMax(H).rmsd) return;
    
    if (IsFull(H)) DeleteMax(H);
    

    for (i = ++H->Size; H->Elements[ i / 2 ].rmsd < X.rmsd; i /= 2)
        H->Elements[ i ] = H->Elements[ i / 2 ];
    H->Elements[ i ] = X;
}

/* END */


ElementType DeleteMax(PriorityQueue H) {
    int i, Child;
    ElementType MaxElement, LastElement;

    /* 1*/ if (IsEmpty(H)) {
        /* 2*/ printf("Priority queue is empty");
        /* 3*/ return H->Elements[ 0 ];
    }
    /* 4*/ MaxElement = H->Elements[ 1 ];
    /* 5*/ LastElement = H->Elements[ H->Size-- ];

    /* 6*/ for (i = 1; i * 2 <= H->Size; i = Child) {
        /* Find bigger child */
        /* 7*/ Child = i * 2;
        /* 8*/ if (Child != H->Size && H->Elements[ Child + 1 ].rmsd
                /* 9*/ > H->Elements[ Child ].rmsd)
            /*10*/ Child++;

        /* Percolate one level */
        /*11*/ if (LastElement.rmsd < H->Elements[ Child ].rmsd)
            /*12*/ H->Elements[ i ] = H->Elements[ Child ];
        else
            /*13*/ break;
    }
    /*14*/ H->Elements[ i ] = LastElement;
    /*15*/ return MaxElement;
}

ElementType FindMax(PriorityQueue H) {
    if (!IsEmpty(H))
        return H->Elements[ 1 ];
    printf("Priority Queue is Empty");
    return H->Elements[ 0 ];
}

int IsEmpty(PriorityQueue H) {
    return H->Size == 0;
}

int IsFull(PriorityQueue H) {
    return H->Size == H->Capacity;
}

void Destroy(PriorityQueue H) {
    free(H->Elements);
    free(H);
}

#if 0

for (i = N / 2; i > 0; i--)
    PercolateDown(i);

#endif
