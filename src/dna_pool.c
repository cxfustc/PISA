#include "dna_pool.h"
#include "htslib/kstring.h"

static const unsigned char seq_nt16_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

static const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

#define MAX_DISTANCE_UMI 1

int dna_dist(const uint8_t *a, const uint8_t *b, size_t l)
{
    int dist = 0;
    while (l-- > 0) {
        if (*a++ != *b++) {
            if ((a[-1] & 0xf) != (b[-1]&0xf)) dist++; 
            if (a[-1]>>4 != b[-1]>>4) dist++;
        }
    }
    return dist;
}
// correct similar pool::dna by pool::count
void dna_pool_corr(struct PISA_dna_pool *p, int e)
{
    int i, j;
    int *idx = malloc(p->l*sizeof(int));
    for (i = 0; i < p->l; ++i) idx[i] = i;

    if (e <0) e = 1;
    
    if (e == 0) { // disable correction
    }
    else {    

        for (i = 0; i < p->l; ++i) {
            struct PISA_dna *a = &p->data[i];
            if (idx[i] != i) continue; // already be updated
            int best = i;
            int best_cnt = a->count;
            for (j = i+1; j < p->l; ++j) {
                struct PISA_dna *b = &p->data[j];
                int dist = dna_dist(a->dna, b->dna, p->len/2);
                if (dist <= e) {                
                    if (best_cnt >= b->count) idx[j] = best;
                    else { // refresh index
                        int k;
                        for (k = 0; k < p->l; ++k)
                            if (idx[k] == best) idx[k] = j;
                        best = j;
                        best_cnt = b->count;
                    }
                }
            }
        }
    }

    // update alias
    for (i = 0; i < p->l; ++i) {
        struct PISA_dna *a = &p->data[i];
        a->alias = idx[i];
    }
    
    free(idx);
}
char *dna_decode_str(struct PISA_dna *b, int l)
{
    kstring_t str = {0,0,0};
    int i;
    for (i = 0; i < l; ++i) {
        kputc(seq_nt16_str[b->dna[i]>>4], &str);
        kputc(seq_nt16_str[b->dna[i]&0xf], &str);
    }
    return str.s;
    
}
// use dna_pool_corr() before dna_corr(); 
char *dna_corr(struct PISA_dna_pool *p, const char *seq)
{
    struct PISA_dna *a = PISA_dna_query(p, seq);
    if (a == NULL) {
        LOG_print("%s", seq);
        return NULL;
    }
    //assert(a);
    if (a->alias < 0) return NULL;
    struct PISA_dna *b = &p->data[a->alias];
    return dna_decode_str(b, p->len/2);
}
//static const int seq_nt16_int[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

void PISA_dna_destroy(struct PISA_dna_pool *p)
{
    int i;
    for (i = 0; i < p->l; ++i) 
        free(p->data[i].dna);
    free(p->data);
    free(p);
}
void PISA_idx_destroy(struct PISA_dna_pool *p)
{
    free(p->data); free(p);
}
                           
struct PISA_dna_pool *PISA_dna_pool_init()
{
    struct PISA_dna_pool *p = malloc(sizeof(*p));
    memset(p, 0, sizeof(*p));
    return p;
}

static uint8_t *PISA_dna_pack(const char *seq, int l)
{
    int l2 = l&~1, i;
    uint8_t *p = malloc(l2);
    for (i = 0; i < l2; i+=2) 
        p[i>>1] = seq_nt16_table[(unsigned char)seq[i]] << 4 | seq_nt16_table[(unsigned char)seq[i+1]];
    return p;
}
// return index of query sequence, because all the sequence save in the pool by order, the index is not stable
int PISA_dna_query0(const struct PISA_dna_pool *p, const char *seq)
{
    if (p->len == 0 || p->l== 0) return -1;
    int l = strlen(seq);    
    if (p->len != l) error("Try to insert an unequal length sequence. %s, %d vs %d.", seq, l, p->len);
    uint8_t *a = PISA_dna_pack(seq, l);
    int i = 0, j = p->l > 0 ? p->l-1 : 0;
    int l2 = l&~1;
    for (;;) {
        if (j < i) break;
        int ret;
        ret = memcmp(p->data[i].dna, a, l2);
        if (ret == 0) return i;
        if (ret > 0) return -1;
        
        ret = memcmp(p->data[j].dna, a, l2);
        if (ret == 0) return j;
        if (ret < 0) return -1;
        
        int m = (i+j)/2;
        
        ret = memcmp(p->data[m].dna, a, l2);
        if (ret == 0) return m;

        if (i == m||j==m) return -1;
        
        if (ret > 0) j = m;
        else i = m;
    }

    return -1;
}
int PISA_idx_query0(struct PISA_dna_pool *p, const int idx)
{
    if (p->l == 0) return -1;
    int i = 0, j = p->l > 0 ? p->l-1 : 0;
   
    for (;;) {
        assert(i>=0 && j>=0);
        if (j < i) break;
        if (p->data[i].idx == idx) return i;
        if (p->data[i].idx > idx) return -1;
        if (p->data[j].idx == idx) return j;
        if (p->data[j].idx < idx) return -1;

        int m = (i+j)/2;
        if (p->data[m].idx == idx) return m;

        if (i == m||j==m) return -1;
        
        if (p->data[m].idx > idx) j = m;
        else i = m;
    }
    return -1;
}
static void* enlarge_pool(struct PISA_dna_pool *p)
{
    if (p->l >= p->m) {
        p->m = p->l==0 ? 2 : p->l*2;
        p->data = realloc(p->data, p->m*sizeof(struct PISA_dna));
    }
    return p->data;
}
static int PISA_idx_pool_add(struct PISA_dna_pool *p, int i, const int idx)
{
    p->l++;
    enlarge_pool(p);

    memmove(p->data+i+1, p->data+i, (p->l-i)*sizeof(struct PISA_dna));
    memset(&p->data[i], 0, sizeof(struct PISA_dna));
    p->data[i].idx = idx;
    return i;
}
static int PISA_idx_push_core(struct PISA_dna_pool *p, const int idx)
{
    enlarge_pool(p);
    if (p->l == 0) {
        p->data[0].idx = idx;
        return p->l++;
    }
    int i = 0, j = p->l > 0 ? p->l-1 : 0;
    for (;;) {
        if (p->data[i].idx == idx) return i;
        if (p->data[i].idx > idx) return PISA_idx_pool_add(p, i, idx);
        if (p->data[j].idx == idx) return j;
        if (p->data[j].idx < idx) return PISA_idx_pool_add(p, j+1, idx);

        int m = (i+j)/2;
        //if (m==i) return PISA_idx_pool_add(p, i, idx);
        if (m==i||m==j) return PISA_idx_pool_add(p, j, idx); // add to j
        
        if (p->data[m].idx == idx) return m;
        if (p->data[m].idx > idx) j = m;
        else i = m;
    }
    return -1;
}
static int PISA_dna_pool_add(struct PISA_dna_pool *p, int i, uint8_t *a, int l)
{
    p->l++;
    enlarge_pool(p);
    memmove(p->data+i+1, p->data+i, (p->l-i)*sizeof(struct PISA_dna));
    memset(&p->data[i], 0, sizeof(struct PISA_dna));
    p->data[i].dna = a;
    return i;
}
void PISA_dna_pool_print(struct PISA_dna_pool *p) {
    int i, j;
    kstring_t str = {0,0,0};
    debug_print("l : %d", p->l);
    for (i = 0; i < p->l; ++i) {
        uint8_t *a = p->data[i].dna;
        str.l = 0;
        for (j = 0; j < p->len/2; ++j) {
            kputc(seq_nt16_str[a[j]>>4], &str);
            kputc(seq_nt16_str[a[j]&0xf], &str);
        }
        debug_print("%s\t%d", str.s, p->data[i].count);
    }
    free(str.s);
}

static int PISA_dna_push_core(struct PISA_dna_pool *p, const char *seq)
{
    //PISA_dna_pool_print(p);
    int l = strlen(seq);
    if (p->len == 0) p->len = l;
    if (p->len != l) error("Try to insert an unequal length sequence. %s, %d vs %d.", seq, l, p->len);
    uint8_t *a = PISA_dna_pack(seq, l);
    enlarge_pool(p);

    if (p->l ==0) {
        p->data[0].dna = a;
        p->data[0].count = 0;
        return p->l++;
    }

    int ret;
    int i= 0, j=p->l>0?p->l-1:0;
    int l2 = l/2;
    for (;;) {
        
        ret = memcmp(p->data[i].dna, a, l2);
        if (ret == 0) {
            free(a);
            return i;
        }
        if (ret > 0) // add to front of data[i]
            return PISA_dna_pool_add(p, i, a, l2);
            
        ret = memcmp(p->data[j].dna, a, l2);
        if (ret == 0) {
            free(a);
            return j;
        }
        if (ret < 0) // add to end of data[j]
            return PISA_dna_pool_add(p, j+1, a, l2);
        
        int m = (i+j)/2;

        //if (m == i) 
            //return PISA_dna_pool_add(p, i, a, l2);
        
        if (m==i||m == j)
            return PISA_dna_pool_add(p, j, a, l2);
        
        ret = memcmp(p->data[m].dna, a, l2);
        if (ret == 0) {
            free(a);
            return m;
        }

        if (ret < 0) i = m;
        else j = m;
    }
    return -1;
}
struct PISA_dna *PISA_dna_push(struct PISA_dna_pool *p, const char *seq)
{
    int idx = PISA_dna_push_core(p, seq);
    p->data[idx].count++;
    return &p->data[idx];
}
struct PISA_dna *PISA_dna_push1(struct PISA_dna_pool *p, const char *seq, void *data)
{
    int idx = PISA_dna_push_core(p, seq);
    p->data[idx].data = data;
    return &p->data[idx];
}
struct PISA_dna *PISA_dna_query(struct PISA_dna_pool *p, const char *seq)
{
    int idx = PISA_dna_query0(p, seq);
    if (idx == -1) return NULL;
    return &p->data[idx];
}

struct PISA_dna *PISA_idx_push(struct PISA_dna_pool *p, const int idx)
{
    int ret = PISA_idx_push_core(p, idx);
    p->data[ret].count++;
    return &p->data[ret];
}
struct PISA_dna *PISA_idx_push1(struct PISA_dna_pool *p, const int idx, void *data)
{
    int ret = PISA_idx_push_core(p, idx);
    p->data[ret].data = data;
    return &p->data[ret];
}
struct PISA_dna *PISA_idx_query(struct PISA_dna_pool *p, const int idx)
{
    int ret = PISA_idx_query0(p, idx);
    if (ret == -1) return NULL;
    return &p->data[ret];
}



#ifdef _DNA_POOL_MAIN
int main()
{
    char *a = "AAGATGGCCTTACGAAGTGC";
    char *b = "AAGATGGCCTTACGAAGTGG";
    char *c = "AAGATGGCCTTACGAAGTGA";
    char *d = "GAGATGGCCTTACGAAGTGA";

    struct PISA_dna_pool *p = PISA_dna_pool_init();
    PISA_dna_push(p, a);
    PISA_dna_pool_print(p);
    PISA_dna_push(p, a);
    PISA_dna_pool_print(p);
    PISA_dna_push(p, a);
    PISA_dna_pool_print(p);
    PISA_dna_push(p, b);
    PISA_dna_pool_print(p);
    PISA_dna_push(p, c);
    PISA_dna_push(p, c);        
    PISA_dna_pool_print(p);

    PISA_dna_push(p, a);
    PISA_dna_pool_print(p);

    PISA_dna_push(p, d);
    PISA_dna_pool_print(p);
        
    PISA_dna_destroy(p);
    return 0;
}
#endif
