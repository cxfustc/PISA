/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2021-03-11 14:40:16
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kson.h"
#include "util.h"
#include "config.h"
#include "number.h"
#include "json_config.h"
#include "sim_search.h"

void bcode_reg_clean(struct bcode_reg *br)
{
    if (br->wl) ss_destroy(br->wl);
}

void bcode_reg_destroy(struct bcode_reg *br)
{
    bcode_reg_clean(br);
    free(br);
}

struct BRstat {
    int q30_bases;
    int bases;
    int exact_match;
    int filter;
};

struct BRstat *BRstat_alloc()
{
    struct BRstat *b = malloc(sizeof(*b));
    memset(b, 0, sizeof(*b));
    return b;
}

void
config_destroy2 (config_t * config)
{
    if (config->cell_barcode_tag) free(config->cell_barcode_tag);
    if (config->sample_barcode_tag) free(config->sample_barcode_tag);
    if (config->raw_cell_barcode_tag) free(config->raw_cell_barcode_tag);
    if (config->raw_cell_barcode_qual_tag) free(config->raw_cell_barcode_qual_tag);
    int i;
    for (i = 0; i < config->n_cell_barcode; ++i) {
        struct bcode_reg *br = &config->cell_barcodes[i];
        bcode_reg_clean(br);
    }
    free(config->cell_barcodes);
    if (config->sample_barcodes) bcode_reg_destroy(config->sample_barcodes);
    if (config->UMI) bcode_reg_destroy(config->UMI);
    if (config->read_1) bcode_reg_destroy(config->read_1);
    if (config->read_2) bcode_reg_destroy(config->read_2);
    if (config->raw_sample_barcode_tag) free(config->raw_sample_barcode_tag);
    if (config->raw_sample_barcode_qual_tag) free(config->raw_sample_barcode_qual_tag);
    if (config->umi_tag) free(config->umi_tag);
    if (config->umi_qual_tag) free(config->umi_qual_tag);

    free (config);
}

config_t *
config_init2(const char *fn, int nth)
{
    config_t * config = (config_t *) ckalloc (1, sizeof(config_t));

    char *config_str = json_config_open(fn);
    if (config_str == NULL) err_mesg("Empty configure file.");
    kson_t *json = kson_parse(config_str);
    free(config_str);

    const kson_node_t *root = json->root;
    if (root == NULL) err_mesg("Format error. Root node is empty.");
    int i;
    for (i = 0; i < root->n; ++i) {
        const kson_node_t *node = kson_by_index(root,i);
        if (node == NULL) continue;
        if (node->key == NULL) {
            warn_mesg("Format error. Node key is empty. skip..");
            continue;
        }
        if (strcmp(node->key, "platform") == 0) {
            continue;
        }
        else if (strcmp(node->key, "version") == 0) {
            continue;
        }
        else if (strcmp(node->key, "cell barcode tag") == 0) {
            if (node->v.str) config->cell_barcode_tag = strdup(node->v.str);
        }
        else if (strcmp(node->key, "cell barcode raw tag") == 0) {
            if (node->v.str) config->raw_cell_barcode_tag = strdup(node->v.str);
        }
        else if (strcmp(node->key, "cell barcode raw qual tag") == 0) { 
            if (node->v.str) config->raw_cell_barcode_qual_tag = strdup(node->v.str);
        }        
        else if (strcmp(node->key, "cell barcode") == 0) {
            if (node->type != KSON_TYPE_BRACKET) err_mesg("Format error. \"cell barcode\":[{},{}]");
            config->n_cell_barcode = node->n;
            config->cell_barcodes = calloc(node->n,sizeof(struct bcode_reg));

            int j;
            for (j = 0; j < node->n; ++j) {
                const kson_node_t *n1 = kson_by_index(node, j);                
                if (n1 == NULL) err_mesg("cell barcode is empty.");
                if (n1->type != KSON_TYPE_BRACE) err_mesg("Format error. \"cell barcode\":[{},{}]");
                if (n1->n == 0) continue; // empty record
                struct bcode_reg *br = &config->cell_barcodes[j];
                memset(br, 0, sizeof(struct bcode_reg));
                int k;
                for (k = 0; k < n1->n; ++k) {
                    const kson_node_t *n2 = kson_by_index(n1, k);
                    if (n2 == NULL) err_mesg("cell barcode record is empty.");
                    if (strcmp(n2->key, "location") == 0) {
                        // assume format R1:1-2
                        char *p = n2->v.str;
                        if (strlen(p) < 6) err_mesg("Unknown location format, should be like \"R1:1-2\"");
                        if (p[0] == 'R' && p[2] == ':') {
                            if (p[1] == '1') br->rd = 1;
                            else if (p[1] == '2') br->rd = 2;                            
                            else  err_mesg("Unknown location format, should be like \"R[12]:1-2\"");
                        }
                        else {
                            err_mesg("Unknown location format, should be like \"R1:1-2\"");
                        }
                        p = p + 3;
                        char *s = p;
                        int c = 0;
                        for (; check_char_num(*s); s++,c++);
                        br->start = str2int_l(p, c);
                        if (*s != '-') err_mesg("Unknown location format, should be like \"R1:1-2\"");
                        p = ++s;
                        c = 0;
                        for (; check_char_num(*s); s++,c++);
                        br->end = str2int_l(p, c);
                        br->len = br->end - br->start +1;
                    }
                    else if (strcmp(n2->key, "distance") == 0) {
                        br->dist = str2int(n2->v.str);
                    }
                    else if (strcmp(n2->key, "white list") == 0) {
                        if (n2->type != KSON_TYPE_BRACKET) err_mesg("Format error. \"white list\":[]");
                        br->n_wl = n2->n;
                        br->white_list = malloc(n2->n*sizeof(char*));
                        int l;
                        for (l = 0; l < n2->n; ++l) {
                            const kson_node_t *n3 = kson_by_index(n2, l);
                            br->white_list[l] = strdup(n3->v.str);                            
                        }                                        
                    }
                    else err_mesg("Unknown key : \"%s\"", n2->key);
                }
            }
        }
        else if (strcmp(node->key, "sample barcode tag") == 0) {
            if (node->v.str) config->sample_barcode_tag = strdup(node->v.str);
        }
        else if (strcmp(node->key, "sample barcode") == 0) {
            
        }
        else if (strcmp(node->key, "read 1") == 0) {
            if (node->type != KSON_TYPE_BRACE) err_mesg("Format error. \"read 1\":{}");
            config->read_1 = malloc(sizeof(struct bcode_reg));
            struct bcode_reg *br = config->read_1;
            memset(br, 0, sizeof(struct bcode_reg));
            const kson_node_t *n1 = kson_by_index(node, 0);
            if (n1 == NULL) err_mesg("read 1 location at config is empty.");
            if (strcmp(n1->key, "location") == 0) {
                // assume format R1:1-2
                char *p = n1->v.str;
                if (strlen(p) < 6) err_mesg("Unknown location format, should be like \"R1:1-2\"");
                if (p[0] == 'R' && p[2] == ':') {
                    if (p[1] == '1') br->rd = 1;
                    else if (p[1] == '2') br->rd = 2;                            
                    else  err_mesg("Unknown location format, should be like \"R[12]:1-2\"");
                }
                else {
                    err_mesg("Unknown location format, should be like \"R1:1-2\"");
                }
                p = p + 3;
                char *s = p;
                int c = 0;
                for (; check_char_num(*s); s++,c++);
                br->start = str2int_l(p, c);
                if (*s != '-') err_mesg("Unknown location format, should be like \"R1:1-2\"");
                p = ++s;
                c = 0;
                for (; check_char_num(*s); s++,c++);
                br->end = str2int_l(p, c);
                br->len = br->end - br->start +1;
            }
        }
        else if (strcmp(node->key, "read 2") == 0) {
            if (node->type != KSON_TYPE_BRACE) err_mesg("Format error. \"read 1\":{}");
            config->read_2 = malloc(sizeof(struct bcode_reg));
            struct bcode_reg *br = config->read_2;
            memset(br, 0, sizeof(struct bcode_reg));
            const kson_node_t *n1 = kson_by_index(node, 0);
            if (n1 == NULL) err_mesg("read 2 location at config is empty.");
            if (strcmp(n1->key, "location") == 0) {
                // assume format R1:1-2
                char *p = n1->v.str;
                if (strlen(p) < 6) err_mesg("Unknown location format, should be like \"R1:1-2\"");
                if (p[0] == 'R' && p[2] == ':') {
                    if (p[1] == '1') br->rd = 1;
                    else if (p[1] == '2') br->rd = 2;                            
                    else  err_mesg("Unknown location format, should be like \"R[12]:1-2\"");
                }
                else {
                    err_mesg("Unknown location format, should be like \"R1:1-2\"");
                }
                p = p + 3;
                char *s = p;
                int c = 0;
                for (; check_char_num(*s); s++,c++);
                br->start = str2int_l(p, c);
                if (*s != '-') err_mesg("Unknown location format, should be like \"R1:1-2\"");
                p = ++s;
                c = 0;
                for (; check_char_num(*s); s++,c++);
                br->end = str2int_l(p, c);
                br->len = br->end - br->start +1;
            }
        }
        else if (strcmp(node->key, "UMI tag") == 0) {
            if (node->v.str) config->umi_tag = strdup(node->v.str);
        }
        else if (strcmp(node->key, "UMI qual tag") == 0) {
            if (node->v.str) config->umi_qual_tag = strdup(node->v.str);
        }
        else if (strcmp(node->key, "UMI") == 0) {
            if (node->type != KSON_TYPE_BRACE) err_mesg("Format error. \"UMI\":{}");
            config->UMI = malloc(sizeof(struct bcode_reg));
            struct bcode_reg *br = config->UMI;
            memset(br, 0, sizeof(*br));
            const kson_node_t *n1 = kson_by_index(node, 0);
            if (n1 == NULL) err_mesg("UMI location at config is empty.");
            if (strcmp(n1->key, "location") == 0) {
                // assume format R1:1-2
                char *p = n1->v.str;
                if (strlen(p) < 6) err_mesg("Unknown location format, should be like \"R1:1-2\"");
                if (p[0] == 'R' && p[2] == ':') {
                    if (p[1] == '1') br->rd = 1;
                    else if (p[1] == '2') br->rd = 2;                            
                    else  err_mesg("Unknown location format, should be like \"R[12]:1-2\"");
                }
                else {
                    err_mesg("Unknown location format, should be like \"R1:1-2\"");
                }
                p = p + 3;
                char *s = p;
                int c = 0;
                for (; check_char_num(*s); s++,c++);
                br->start = str2int_l(p, c);
                if (*s != '-') err_mesg("Unknown location format, should be like \"R1:1-2\"");
                p = ++s;
                c = 0;
                for (; check_char_num(*s); s++,c++);
                br->end = str2int_l(p, c);
                br->len = br->end - br->start +1;
            }
        }
        else {
            err_mesg("Unknown key \"%s.\"", node->key);
        }            
    }

    // init white list hash
    for (i = 0; i < config->n_cell_barcode; ++i) {
        struct bcode_reg *br = &config->cell_barcodes[i];
        if (br->n_wl == 0) continue;
        br->wl = ss_init(nth);
        int j;
        for (j = 0; j < br->n_wl; j++) {
            int len = strlen(br->white_list[j]);
            if (len != br->len) err_mesg("Inconsistance white list length. %d vs %d, %s", br->len, len, br->white_list[j]);
            if (br->dist>3) err_mesg("Set too much distance for cell barcode, allow 3 distance at max.");
            if (br->dist > len/2) err_mesg("Allow distance greater than half of barcode! Try to reduce distance.");
            ss_push(br->wl, br->white_list[j], len);
            free(br->white_list[j]);
        }
        free(br->white_list);
    }
    kson_destroy(json);

    return config;
}
