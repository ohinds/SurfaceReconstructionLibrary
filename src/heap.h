/*
 * Binary HEAP.
 *
 * This code implements:
 *
 * heap_t        : Type to delare new heaps.
 * heap_init     : Initilizes the heap.
 * heap_destroy  : Destroys the heap.
 * heap_isempty  : Function that returns true if the heap is empty.
 * heap_insert   : To insert a new element into the heap.
 * heap_top      : Returns the data at the top of the heap.
 * heap_removetop: Same than heap_top but it also removes the top element.
 * heap_rm_element: Returns the element, at position, value and removes it from the heap.
 *
 *
 * Released under the GPL V2 (GNU Public License, Version 2)
 * Copyright (C) 2000 Ismael Ripoll  (Valencia, Spain)
 * Copyright (C) 2006 Oliver Hinds <oph@bu.edu>
 *
 *
 *
 *
 */

#include"libsrTypes.h"

#ifndef __BINARY_HEAP__
#define __BINARY_HEAP__

/* Print the heap and test if it is correct*/
//#define DEBUG_AND_TEST

//#define INDEX_TYPE     unsigned char
//#define INDEX_TYPE     unsigned short
#define INDEX_TYPE     unsigned long
//#define INDEX_TYPE     unsigned long long

//#include "mqueue.h"
//#define HEAP_ITEM_TYPE      long
//#define HEAP_ITEM_TYPE        mp_item_handler
#define HEAP_ITEM_TYPE  edge*

/* key's TYPE to sort the heap */
//#define HEAP_KEY_SORT_TYPE  unsigned long long
#define HEAP_KEY_SORT_TYPE  double

/* Among same priority, keep FIFO order. */
#define FIFO_PRIO

/*There is a posible error of overflow in the order values
 *  So to minimize this posibility redefine the ORDER_TYPE
 *    as you think.
 * (unsigned long long) is practicaly impossible to overflow,
 *    but is takes a lot of space.
 */
#ifdef FIFO_PRIO
//#define ORDER_TYPE     unsigned short
#define ORDER_TYPE     unsigned long
//#define ORDER_TYPE     unsigned long long
#endif


#ifndef NULL
#define NULL ((void *)0)
#endif

#define SMALLEST 1

struct heapdt_t_t;
struct item_dt_t_t;

typedef struct item_dt_t_t{
  HEAP_ITEM_TYPE   user_data;
  union {
    struct item_dt_t_t *next_free;
    //struct heapdt_t_t  *heapdata_p;
    INDEX_TYPE heapdata;
  }p_ctrl;
  HEAP_KEY_SORT_TYPE key;
#ifdef FIFO_PRIO
  ORDER_TYPE order;
#endif
} item_data_t;

typedef struct  heapdt_t_t{
  item_data_t *val;
} heapdata_t;

//typedef INDEX_TYPE heap_it_handler;
typedef item_data_t * heap_it_handler;

#define H_KEY(__HEAP_ ,__IND__)  (__HEAP_)->heap[(__IND__)].val->key


/* This is a heap */
typedef struct {
  heapdata_t   *heap; /* The first element of the array is unused: heap[0] */
  unsigned int heap_len;   /* number of elements in the heap */
  unsigned int max_len;
  item_data_t  *data_p;
  item_data_t  *first_free;
#ifdef FIFO_PRIO
  ORDER_TYPE actual;
#define _ACTUAL_ ,((ORDER_TYPE)0)
#else
#define _ACTUAL_
#endif
#define _POINT_ ,NULL,((INDEX_TYPE)0)
} heap_t;


#define HEAP_INITIALIZER {NULL, 0, 0  _ACTUAL_ _POINT_}

#define HEAP_ISEMPTY(_h_) ((_h_)->heap_len==0)
#define HEAP_TOP(_h_)     ((_h_)->heap[SMALLEST].val->user_data)
#define HEAP_TOP_KEY(_h_) ((_h_)->heap[SMALLEST].val->key)

/* heap_init - initializes the heap */
extern int heap_init(heap_t *, unsigned int );


/* heap_destroy - delete the heap */
extern int heap_destroy(heap_t *);


/* heap_isempty - returns true if the heap is empty */
extern int heap_isempty(heap_t *h);

/* heap_top - returns the data associated to the top of the heap */
extern  HEAP_ITEM_TYPE  *heap_top(heap_t *h);


/* heap_topkey - returns the key associated to the top of the heap */
extern  HEAP_KEY_SORT_TYPE  heap_topkey(heap_t *h);


/* heap_removetop - returns the heap top value and removes it from the heap */
extern inline void *heap_removetop(heap_t *);


/* heap_insert - insert a new iten into the heap */
//extern int heap_insert(heap_t *, void *, HEAP_KEY_SORT_TYPE );
//extern item_data_t *heap_insert(heap_t *h, HEAP_ITEM_TYPE val, HEAP_KEY_SORT_TYPE key);
extern heap_it_handler heap_insert(heap_t *h, HEAP_ITEM_TYPE val, HEAP_KEY_SORT_TYPE key);

/* push a data item up or down the heap until in the proper place */
extern void heapfloat_heapify(heap_t *h, int position);

/* heap_rm_element - returns the element value and removes it from the heap */
extern HEAP_ITEM_TYPE heap_rm_element(heap_t *h, heap_it_handler element);


extern HEAP_ITEM_TYPE heap_element_from_handler(heap_it_handler element);

#ifdef DEBUG_AND_TEST
/* print_heap - prints the array, that stores the heap */
extern void print_heap(heap_t *h);

/* test_heap - returns 0 if all the elements in the heap are less than its parent, else return -1  */
extern int test_heap(heap_t *h);
#endif

#endif
