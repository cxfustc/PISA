#include "multi_index_list.h"

int lstcmp(int *idx1, int *idx2, int l)
{
    while (l-- > 0) {
        if (*idx1++ != *idx2++) 
            return idx1[-1] - idx2[-1];
    }
    return 0;
}
void enlarge_list(struct index_list *list)
{
    if (list->l >= list->m) {
        list->m = list->l == 0 ? 4 : list->l*2;
        list->mi = realloc(list->mi, sizeof(struct multi_index)*list->m);
    }
}
int index_add(struct index_list *list, int i, int *idx)
{
    list->l++;
    enlarge_list(list);
    memmove(list->mi+i+1, list->mi+i, (list->l-i)*sizeof(struct multi_index));
    memset(&list->mi[i], 0, sizeof(struct multi_index));
    list->mi[i].mi = idx;
    return i;
}

static int index_push_core(struct index_list *list, int l, int *idx)
{
    enlarge_list(list);
    if (list->l_idx == 0) list->l_idx = l;

    if (list->l == 0) {
        return index_add(list, 0, idx);
    }

    if (l != list->l_idx) error("Try to push an unequal length idx? %d vs %d", l, list->l_idx);
    int ret;
    if (list->l == 1) {
        ret = lstcmp(list->mi[0].mi, idx, l);
        if (ret == 0) return 0;
        
        if (ret > 0) return index_add(list, 0, idx);
        else return index_add(list, 1, idx);
    }
    
    int i= 0, j = list->l-1;
    
    for (;;) {
        ret = lstcmp(list->mi[i].mi, idx, l);
        if (ret == 0) return i;
        if (ret > 0) return index_add(list, i, idx);

        ret = lstcmp(list->mi[j].mi, idx, l);
        if (ret == 0) return j;
        if (ret < 0) return index_add(list, j+1, idx);

        int m;
        m = (i+j)/2;
        if (m == i||m==j) return index_add(list, j, idx);
        ret = lstcmp(list->mi[m].mi, idx, l);
        if (ret == 0) return m;
        if (ret > 0) j = m;
        else i = m;
    }
}
struct multi_index *index_push(struct index_list *list, int l, int *idx)
{
    int i = index_push_core(list, l, idx);
    list->mi[i].count++;
    return &list->mi[i];
}
struct multi_index *index_push1(struct index_list *list, int l, int *idx, void *data)
{
    int i = index_push_core(list, l, idx);
    list->mi[i].data = data;
    return &list->mi[i];
}
struct multi_index *index_query(struct index_list *list, int l, int *idx)
{
    if (list->l == 0) return NULL;
    if (list->l == 1) {
        if (lstcmp(list->mi[0].mi, idx, l) == 0) return &list->mi[0];
        return NULL;
    }

    if (l != list->l_idx) error("Try to query an unequal length idx? %d vs %d", l, list->l_idx);
    
    int i= 0, j = list->l-1;
    int ret;
    for (;;) {
        ret = lstcmp(list->mi[i].mi, idx, l);
        if (ret == 0) return &list->mi[i];
        if (ret > 0) return NULL;

        ret = lstcmp(list->mi[j].mi, idx, l);
        if (ret == 0) return &list->mi[j];
        if (ret < 0) return NULL;

        int m;
        m = (i+j)/2;
        if (m == i || m == j) return NULL;

        ret = lstcmp(list->mi[m].mi, idx, l);
        if (ret == 0) return &list->mi[m];
        if (ret > 0) j = m;
        else i = m;
    }
}

struct index_list *index_list_init()
{
    struct index_list *list = malloc(sizeof(*list));
    memset(list, 0, sizeof(*list));
    return list;
}

void index_list_destroy(struct index_list *list)
{
    int i;
    for (i = 0; i < list->l; ++i)
        free(list->mi[i].mi);
    free(list->mi);
    free(list);
}
