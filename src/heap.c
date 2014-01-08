/*
 *
 * heap_t        : Type to declare new heaps.
 * heap_init     : Initializes the heap.
 * heap_destroy  : Destroys the heap.
 * heap_isempty  : Function that returns true if the heap is empty.
 * HEAP_ISEMPTY  : Macro that returns true if the heap is empty.
 * heap_insert   : To insert a new element into the heap.
 * heap_top      : Function that returns the data at the top of the heap.
 * HEAP_TOP      : Macro that returns the data at the top of the heap.
 * heap_removetop: Same than heap_top but it also removes the top element.
 * heap_rm_element: Returns the element, at position, value and removes it from the heap.
 *
 * This code can be compiled without modification for tree different
 * environments: (1) Linux kernel, (2) RTLinux and (3) a normal Linux
 * process.
 * 
 *
 * Released under the GPL V2 (GNU Public License, Version 2)
 * Copyright (C) 2000 Ismael Ripoll  (Valencia, Spain)
 * 
 *
 *
 *
 */

/* Print debug info */
#define DEBUG_HEAP

/* Include code to perform sany checks... but it is not for free */
/* Comment this define to improve the performance */
//#define DEBUG_SANITY_CHECKS


#ifdef __KERNEL__

#include <linux/kernel.h>
#include <linux/version.h>
#include <linux/malloc.h> 

#ifdef __RTL__
#include <rtl.h>
#define printf rtl_printf
#else
#define printf printk
#endif

#else

#include <stdio.h>
//#include <malloc.h>

#endif


#ifdef CRITICAL_SEC 
#include <pthread.h> 
#endif

#include "heap.h"



/*
 * heap_init - initializes the heap
 * 
 * @h: heap to initialize.
 * @size: max number of elements that may hold the heap
 * 
 * This function uses dynamic memory to allocate the heap. When
 * compiled in the Linux Kernel, it uses kmalloc to get the
 * mem... which means that the maximum size is 128K/sizeof(heapdata_t)
 * which is around 10000.
 *
 * Returns true if the heap has been created correctly.
 *
 */
int heap_init(heap_t *h, unsigned int size){
    int i;

#ifdef DEBUG_SANITY_CHECKS
    if (h->heap != NULL){
        printf("Heap: trying to initialize a heap already in use\n");	
        return -1;
    }
#endif

#ifdef FIFO_PRIO
    h->actual =0;
#endif  
    h->heap_len=0;
    h->max_len=size;

#ifdef __KERNEL__
    h->heap=(heapdata_t *) kmalloc((size)* sizeof(heapdata_t), GFP_KERNEL);
    h->data_p=(item_data_t *) kmalloc((size)* sizeof(item_data_t), GFP_KERNEL);
#else
    h->heap=(heapdata_t *) malloc((size)* sizeof(heapdata_t));
    h->data_p=(item_data_t *) malloc((size)* sizeof(item_data_t));
#endif

    h->first_free=h->data_p;
    for(i=0;i<size-1;i++){
        h->data_p[i].p_ctrl.next_free = &h->data_p[i+1];
    };

    /* The last one never could be used 
    h->data_p[size].p_ctrl.next_free = NULL;
    */

    /* J.C. Perez advice. Tricky code to use the first element of the
       heap array. heap[0] is never used. */
    if (h->heap){
        h->heap--;
        return 1;
    } else {
        return 0;
    }
}



/*
 * heap_destroy - delete the heap
 * 
 * @h: heap to destroy.
 *
 * The heap structure is free. The real data is not destroyed since the
 * heap do not have the real data, just pointers to it.
 */
int heap_destroy(heap_t *h){

#ifdef FIFO_PRIO
  h->actual =0;
#endif
  h->heap_len=0;
  h->max_len=0;
  h->first_free=0;

#ifdef DEBUG_SANITY_CHECKS
  if (h->heap == NULL){
    printf("Heap: trying to free an already destroyed heap\n");
    return -1;
  }
#endif

  h->heap++;
#ifdef __KERNEL__
  kfree((void *) h->heap);
  kfree((void *) h->data_p);
#else
  free((void *) h->heap);
  free((void *) h->data_p);
#endif

  h->heap   = NULL;
  h->data_p = NULL;
  return 1;
}



/*
 * lessthan - true if the first parameter is smaller than the second
 *
 * You may want to modify this function to define the order relation
 * used in the heap structure.
 */
  /* return (n>m) if you want a greater values first */
#ifdef FIFO_PRIO
inline int lessthan_2Obj(heapdata_t *n,heapdata_t *m) {
  return((n->val->key < m->val->key)||((n->val->key == m->val->key) && (n->val->order < m->val->order)));
}  /* return (n>m) if you want a greater values first */
inline int lessthan_ValObj(HEAP_KEY_SORT_TYPE n,ORDER_TYPE n_order, heapdata_t *m) {
  return ((n < m->val->key)||((n == m->val->key) && (n_order < m->val->order)));
}  /* return (n>m) if you want a greater values first */
inline int lessthan_2Val(HEAP_KEY_SORT_TYPE n,ORDER_TYPE n_order,HEAP_KEY_SORT_TYPE m,ORDER_TYPE m_order) {
  return ((n < m)||((n == m) && (n_order < m_order)));
}  /* return (n>m) if you want a greater values first */
#else
inline int lessthan(HEAP_KEY_SORT_TYPE n, HEAP_KEY_SORT_TYPE m){
  return (n < m); /* This way the top value of the heap will be the smallest */
  /* return (n>m) if you want a greater values first */
}
#endif


/*
 * heap_isempty - returns true if the heap is empty
 * @h: heap
 */
inline int heap_isempty(heap_t *h){
    return (h->heap_len==0);
}


/*
 * heap_top - returns the data associated to the top of the heap
 *
 * @h: heap to get the top
 *
 * This function do not modify the heap structure. See also heap_removetop.
 */
inline HEAP_ITEM_TYPE *heap_top(heap_t *h){

#ifdef DEBUG_SANITY_CHECKS
    if (h->heap_len==0){
        printf("Heap: call to 'head_top' function in an empty heap\n");
        return (void *) -1;
    }
#endif 
    return (&h->heap[SMALLEST].val->user_data);
}


/*
 * heap_topkey - returns the key associated to the top of the heap
 *
 * @h: heap to get the top key 
 *
 * This function do not modify the heap structure.
 */
inline HEAP_KEY_SORT_TYPE heap_topkey(heap_t *h){

#ifdef DEBUG_SANITY_CHECKS
    if (h->heap_len==0){
        printf("Heap: call to 'head_topkey' function in an empty heap\n");
        return -1.0;
    }
#endif 
    return (h->heap[SMALLEST].val->key);
}




/*
 * heapify - push a data item down the heap until in the proper place
 *
 * @h: the heap structure to work with
 * @top: heap element to start pushing
 *
 * This is an "internal" function used by heap_remove.
 * 
 */
void heapify(heap_t *h, int top) {

    heapdata_t tempvalue = h->heap[top];  /* Data to pushdown */
    int child;                /* left son of top */

    for( child= top<<1; child <= h->heap_len;top = child, child <<= 1 ) {
#ifdef FIFO_PRIO
        if (child < h->heap_len && lessthan_2Obj(&h->heap[child+1], &h->heap[child]))
#else
        if (child < h->heap_len && lessthan(H_KEY(h,child+1),H_KEY(h,child)))
#endif
            child++;
#ifdef FIFO_PRIO
        if (lessthan_2Obj(&tempvalue, &h->heap[child])) 
#else
        if (lessthan(tempvalue.val->key,H_KEY(h,child))) 
#endif
            break;
	
        h->heap[top].val->p_ctrl.heapdata = child ;
        h->heap[top] = h->heap[child];
        h->heap[top].val->p_ctrl.heapdata = top ;
	h->heap[child] = tempvalue;
    }
/*    if (h->heap_len != 0) {
        h->heap[top].val->p_ctrl.heapdata = top;
	h->heap[top] = tempvalue;
    }
*/
}


/* 
 * heap_removetop - returns the heap top value and removes it from the heap.
 * @h: heap
 *
 * Remove the smallest element from the heap and recreate the heap with
 * one less element.
 */
inline void  *heap_removetop(heap_t *h){

    HEAP_ITEM_TYPE *val;

#ifdef DEBUG_SANITY_CHECKS
    if (h->heap_len==0){
        printf("Heap: call to 'head_top' function in an empty heap\n");
        return (void *) -1;
    }
#endif 

    // Storing the user data to return it
    val = &h->heap[SMALLEST].val->user_data;

    // Deleting the root of the heap
    h->heap[SMALLEST].val->p_ctrl.next_free = h->first_free;
    h->first_free = h->heap[SMALLEST].val;

    if (h->heap_len > 1) {
        // Copying the las data at the root 
        h->heap[h->heap_len].val->p_ctrl.heapdata = SMALLEST;
	h->heap[SMALLEST] = h->heap[h->heap_len--];

	// Reordering the heap
	heapify(h, SMALLEST);
    } else {
        h->heap_len = 0;
#ifdef FIFO_PRIO
	h->actual = 0;
#endif
    }
    return (val);
}


/* 
 * heap_insert - insert a new iten into the heap
 * @h: heap
 * @val: data to be inserted into the heap
 * @key: priority of the data. It must be of type HEAP_KEY_SORT_TYPE
 *
 * Inserts a new item in the heap using "key" to order it.
 * "val" is a pointer to the real data.
 */
//int heap_insert(heap_t *h, void *val, HEAP_KEY_SORT_TYPE key) {
//int heap_insert(heap_t *h, HEAP_ITEM_TYPE val, HEAP_KEY_SORT_TYPE key) {
//item_data_t *heap_insert(heap_t *h, HEAP_ITEM_TYPE val, HEAP_KEY_SORT_TYPE key) {
heap_it_handler heap_insert(heap_t *h, HEAP_ITEM_TYPE val, HEAP_KEY_SORT_TYPE key) {
    int position;

    item_data_t  *aux_data;

    if (h->heap_len >= h->max_len){
        /* Heap maximum capacity reached */
#ifdef DEBUG_HEAP
        printf("Heap: Warning the heap is full\n");
#endif
        return NULL;
    }

    position = ++h->heap_len;

#ifdef FIFO_PRIO
    h->actual++ ;
    while ( (position > SMALLEST) && lessthan_ValObj(key,h->actual ,&h->heap[position>>1]))
#else
      while ( (position > SMALLEST) && lessthan(key,H_KEY(h,position>>1)))
#endif
    {
        /* If the new son is less than the father then move down the father */
        h->heap[position] = h->heap[position>>1];
	h->heap[position>>1].val->p_ctrl.heapdata = position;
	//h->heap[position].val->p_ctrl.heapdata = position>>1;
	position = position >> 1;
    }

    aux_data = h->first_free;
    h->first_free = h->first_free->p_ctrl.next_free;
    h->heap[position].val = aux_data;
    h->heap[position].val->user_data = val;
    h->heap[position].val->p_ctrl.heapdata = position;

    H_KEY(h,position) = key;
#ifdef FIFO_PRIO
    h->heap[position].val->order=h->actual;
#endif
    return (h->heap[position].val);
}

/////////////////////////////////////////////////////////////////////////
// NEW IN VERSION 1.1
/////////////////////////////////////////////////////////////////////////


/*
 * float - push a data item up the heap until in the proper place
 *
 * @h: the heap structure to work with
 * @element: heap element to start pushing
 *
 * This is an "internal" function used by float_heapify.
 * 
 */
void heap_float(heap_t *h, INDEX_TYPE position) {

     heapdata_t tempvalue;
     INDEX_TYPE father;
     father = position>>1;
#ifdef FIFO_PRIO
    while ( (position > SMALLEST) && lessthan_2Obj(&h->heap[position],&h->heap[father]))
#else
    while ( (position > SMALLEST) && lessthan(H_KEY(h,position),H_KEY(h,father)))
#endif
    {
        tempvalue = h->heap[position];
        h->heap[position]=h->heap[father];
        h->heap[father] = tempvalue;
	h->heap[position].val->p_ctrl.heapdata = position;
	h->heap[father].val->p_ctrl.heapdata   = father;

	position = father;
	father = position>>1;
    }
}    

/*
 * heapfloat_heapify - push a data item up or down the heap until in the proper place
 *
 * @h: the heap structure to work with
 * @element: heap element to start pushing
 * 
 */
void heapfloat_heapify(heap_t *h, int position) {

    INDEX_TYPE father;

    father = position>>1;

    if (position > SMALLEST) {
#ifdef FIFO_PRIO
        if (lessthan_2Obj(&h->heap[position],&h->heap[father])) 
#else
	if (lessthan(H_KEY(h,position),H_KEY(h,father)))
#endif
	    {
	        heap_float(h,position);
		return;
	    } else {
	        heapify(h,position);
		return;
	    }
    } else {
        if (position == SMALLEST) {
	    heapify(h,position);
	    return;
	} else {
	    return;
	}
    }

}

HEAP_ITEM_TYPE heap_element_from_handler(heap_it_handler element){

  return (element->user_data);
  //return (element->p_ctrl.heapdata);
}

/* 
 * heap_rm_element -returns the element value and removes it from the heap.
 * @h: heap
 * @element: element's index to remove, (othes implementations used pointers)
 *
 * Remove the element from the heap and recreate the heap with
 * one less element.
 */
HEAP_ITEM_TYPE heap_rm_element(heap_t *h, heap_it_handler element){

    HEAP_ITEM_TYPE *val;
    INDEX_TYPE ele;
    ele = element->p_ctrl.heapdata;
#ifdef DEBUG_SANITY_CHECKS
    if ((h->heap_len< ele)||(ele < 1)){
        printf("Heap: the element it is trying to remove doesn't exist.\n");
        return NULL;
    }
#endif 

/*
    val = h->data_p[h->heap[element].val];
    h->data_p[h->heap[element].val]= h->first_free;
    h->first_free = h->heap[element].value;
*/

    if ((h->heap_len< ele)||(ele < 1)){
        printf("Heap: the element it is trying to remove doesn't exist.PAU\n");
        return NULL;
    }


    val = &h->heap[ele].val->user_data;
    h->heap[ele].val->p_ctrl.next_free= h->first_free;
    h->first_free = h->heap[ele].val;
    
    if (ele < h->heap_len){
        h->heap[ele] = h->heap[h->heap_len--];
	h->heap[ele].val->p_ctrl.heapdata = ele;
	heapfloat_heapify(h, ele);
    } else {
        if (ele == h->heap_len){
	    h->heap_len--;
	} else {
	  printf ("Error, when try to delete\n");
	}

    }
#ifdef FIFO_PRIO
    if (h->heap_len==0){
      h->actual = 0;
    }
#endif
    return (*val);
}


#ifdef DEBUG_AND_TEST
void print_heap(heap_t *h) {
    int i;


    for (i=1;i<=h->heap_len; i++){
        printf ("%3d -%3d, %3d, index:%3d\n",i,h->heap[i].val->key,h->heap[i].val->user_data,h->heap[i].val->p_ctrl.heapdata);
    }
    printf ("\n");
}

int test_heap(heap_t *h) {
    int son;
    INDEX_TYPE father;

    for (son=2;son<=(h->heap_len); son++){
        father = son>>1;
        if (h->heap[father].val->key > h->heap[son].val->key)
	    {
  	        printf ("Error\n");
  	        printf ("\n");
	        print_heap(h);
		return -1;
	    }
    }
    for (son=1;son<=(h->heap_len); son++){
        if (h->heap[son].val->p_ctrl.heapdata != son){
	    printf ("Error\n");
	    printf ("\n");
	}
  
    }
    return 0;
}
#endif
