#ifndef H_MULTIDX
#define H_MULTIDX

#include "utils.h"
#include "htslib/kstring.h"

struct multi_index {
    int *mi;
    union {
        void *data;
        uint32_t count;
    };
};

struct index_list {
    struct multi_index *mi;
    int l_idx;
    int l, m;
};

struct multi_index *index_push(struct index_list *list, int l, int *idx);
struct multi_index *index_push1(struct index_list *list, int l, int *idx, void *data);

struct multi_index *index_query(struct index_list *list, int l, int *idx);
struct index_list *index_list_init();
void index_list_destroy(struct index_list *list);

#endif
