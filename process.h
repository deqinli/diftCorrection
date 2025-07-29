#ifndef PROCESS_H
#define PROCESS_H

#include "accord.h"

#define FN_LEN 800
#define PATH_LEN 500

typedef struct file_list{
    char filename[FN_LEN];
    struct file_list *next;
    t_im_struct image;
    char invalid;
} t_file_list;

typedef struct conf_struct{
    char xcorrs; // no xcorrs images is default
    char traditional; // no traditional averaging is default
    char all_images; // only final image is default
    char output_path[PATH_LEN];
    char output_prefix[FN_LEN]; //output-file base name
    char full_weight_normalize; //max and min values are counted from only those pixels that have full information weight (all frames contributed to them)
    ACC_DIM_TYPE cr_left, cr_right, cr_top, cr_bottom; // crop values - default is 0 for all
} t_conf_struct;

extern char (*iteration_callback)(void *data);

void parse_config_file(t_conf_struct *conf, t_file_list **file_list_o, char *filename);
void write_config_file(t_conf_struct *conf, t_file_list *file_list, char *filename);
int process_images(t_file_list *file_list, t_conf_struct *conf);
void destroy_file_list(t_file_list *file_list);

#endif // PROCESS_H
