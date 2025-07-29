#include "accord.h"
#include "process.h"


char (*iteration_callback)(void *data);

void trim_whites(char *str){/*{{{*/
    size_t pos;
    while(str[pos=strlen(str)-1] == ' ') str[pos]=0;
    while(str[pos=strlen(str)-1] == '\t') str[pos]=0;
    while(str[0] == ' ') memcpy(str, str+1, strlen(str));
    while(str[0] == '\t') memcpy(str, str+1, strlen(str));
}/*}}}*/

void parse_config_file(t_conf_struct *conf, t_file_list **file_list_o, char *filename){/*{{{*/


    t_file_list *file_list = NULL;
    // set the config defaults
    conf->xcorrs = 0;
    conf->traditional = 0;
    conf->all_images = 0;
    conf->full_weight_normalize = 1;
    conf->cr_left = 0;
    conf->cr_right = 0;
    conf->cr_top = 0;
    conf->cr_bottom = 0;

    strncpy_s(conf->output_prefix, "out", FN_LEN);
    strncpy_s(conf->output_path, ".", FN_LEN);



    char line[FN_LEN];
    FILE *infile = fopen(filename, "rt");
    if (!infile) {
        fprintf(stderr, "Fatal error: Could not open file %s\n", filename);
        exit(-1);
    }

    char section = '0'; // indicator of the section being parsed
    while (fgets(line, FN_LEN, infile)) {

        char *helper;
        // strip newline
        helper = strrchr(line,'\n');
        if (helper) *helper=0;
        helper = strrchr(line,'\r');
        if (helper) *helper=0;

        // whatever follows # is a comment
        helper = strchr(line,'#');
        if (helper) *helper= 0;

        trim_whites(line);

        if (!strcmp(line,"")) continue; // ignore empty lines

        // reset section indicator, if line contains [ and ]
        if ((strchr(line,'[')) && (strchr(line,']'))) section = '0'; // end of the previous section
        if (strstr(line,"[input files]")) {
            section = 'i';
            continue;
        }
        if (strstr(line,"[config]")) {
            section = 'c';
            continue;
        }

        // process input files
        if (section == 'i'){
            t_file_list *newitem = (t_file_list *) malloc(sizeof(t_file_list));
            newitem->next=NULL;
            strncpy_s(newitem->filename, line, FN_LEN);
            newitem->invalid = 1; // image is not open - it is not valid.
            if (!file_list) { // we are adding the first item;
                file_list = newitem;
                continue;
            }
            t_file_list *h = file_list;
            while (h->next) h = h->next; // looking for the last item
            h->next = newitem; // adding the new item to the end of the list.
        }

        if (section == 'c'){
            char name[100], value[100];
            sscanf_s(line, "%99[^:]:%99s", name, value);
            trim_whites(name);
            trim_whites(value);
            if ((!strcmp(name,"xcorrs")) && (!strcmp(value,"yes"))) conf->xcorrs = 1;
            if ((!strcmp(name,"traditional")) && (!strcmp(value,"yes"))) conf->traditional = 1;
            if ((!strcmp(name,"all images")) && (!strcmp(value,"yes"))) conf->all_images = 1;
            if ((!strcmp(name,"full weight normalize")) && (strcmp(value,"yes"))) conf->full_weight_normalize = 0;
            if (!strcmp(name,"output path")) strncpy_s(conf->output_path, value, PATH_LEN);
            if (!strcmp(name,"output prefix")) strncpy_s(conf->output_prefix, value, FN_LEN);
            if (!strcmp(name,"crop left")) sscanf_s(value,"%ld",&conf->cr_left);
            if (!strcmp(name,"crop right")) sscanf_s(value,"%ld",&conf->cr_right);
            if (!strcmp(name,"crop top")) sscanf_s(value,"%ld",&conf->cr_top);
            if (!strcmp(name,"crop bottom")) sscanf_s(value,"%ld",&conf->cr_bottom);
        }
    }
    fclose(infile);
    *file_list_o = file_list;

}/*}}}*/

void write_config_file(t_conf_struct *conf, t_file_list *file_list, char *filename){/*{{{*/
    FILE *cfout; //config-file output
    cfout = fopen(filename, "wt");
    if (!cfout) {
        fprintf(stderr, "Fatal error: Could not open file %s\n", filename);
        exit(-1);
    }


    fprintf(cfout, "#!/usr/bin/env caccord -c\n\n#This is an ACCORD configuration file.\n# see http://accord.sf.net.\n\n"); // print header

    // Write configuration

    fprintf(cfout, "\n\n[config]\n");
    fprintf(cfout, "output path: %s\n", conf->output_path);
    fprintf(cfout, "output prefix: %s\n", conf->output_prefix);

    fprintf(cfout, "traditional: %s      #yes for plain averaging. Default: no\n", (conf->traditional)? "yes": "no");
    fprintf(cfout, "all images: %s      #yes if all sum images should be saved. Default: no\n", (conf->all_images)? "yes": "no");
    fprintf(cfout, "full weight normalize: %s       #yes if only full-weight pixels determine the norm. Default: yes\n", (conf->full_weight_normalize)? "yes": "no");

    fprintf(cfout,"\n");
    fprintf(cfout, "crop top: %ld\n", conf->cr_top);
    fprintf(cfout, "crop bottom: %ld\n", conf->cr_bottom);
    fprintf(cfout, "crop left: %ld\n", conf->cr_left);
    fprintf(cfout, "crop right: %ld\n", conf->cr_right);

    // Write file list

    fprintf(cfout, "\n\n[input files]\n");
    while(file_list) {
        fprintf(cfout, "%s\n", file_list->filename);
        file_list = file_list->next;
    }
    fclose(cfout);

}/*}}}*/

int process_images(t_file_list *file_list, t_conf_struct *conf){/*{{{*/
    unsigned int i;
    t_file_list *h = file_list;


    char out_fn_base[FN_LEN];
#ifdef WIN32
    snprintf(out_fn_base, FN_LEN, "%s\\%s", conf->output_path, conf->output_prefix);
#else
    snprintf(out_fn_base, FN_LEN, "%s/%s", conf->output_path, conf->output_prefix);
#endif

    FILE *displ_file;
    char outname[FN_LEN];
    snprintf(outname, FN_LEN, "%s_disp.txt", out_fn_base);
    displ_file = fopen(outname, "wt");
    if (!displ_file) {
        fprintf(stderr, "Fatal error: Could not open file %s\n", outname);
        exit(-1);
    }

    t_im_struct final;
    h = file_list;

    if (!h) return -1;
    while (load_image(&h->image, h->filename , IM_NOT_INITIALIZED)) h = h->next;
    if (conf->cr_left || conf->cr_right || conf->cr_top || conf->cr_bottom){
        crop_image(&h->image, conf->cr_left, conf->cr_right, conf->cr_top, conf->cr_bottom);
    }
    copy_image(&final, &h->image, IM_NOT_INITIALIZED);
    initialize_weight_counting(&final);
    h = h->next;
    i=0;
    while(h){
        char err = ACC_ERR_NO_ERROR;
        char msg[500];
        snprintf(msg, 500, "Status: Reading filename: %s",h->filename);
        status_callback(0,msg);

        char res = load_image(&h->image, h->filename , IM_NOT_INITIALIZED);
        if (res) {
            fprintf(stderr, "Cannot open file: \"%s\". Ignoring it.", h->filename);
            h->invalid = 1;
            h = h->next;
            continue;
        }
        h->invalid = 0;

        if (conf->cr_left || conf->cr_right || conf->cr_top || conf->cr_bottom){
            crop_image(&h->image, conf->cr_left, conf->cr_right, conf->cr_top, conf->cr_bottom);
        }

#ifdef WIN32
        char *fnh = strrchr(h->filename, '\\');
#else
        char *fnh = strrchr(h->filename, '/');
#endif
        if (!fnh) fnh = h->filename;
        snprintf(msg, 500, "Status: Processing %s", fnh);
        status_callback(i,msg);
        double shx = 0;
        double shy = 0;

        if (!conf->traditional) err = add_with_correction(&final, &h->image, IM_DESTRUCTIVE, &shx, &shy);
        else err = add_images(&final, &h->image);
        if (err) return err;

        snprintf(msg, 500, "Displacement: [%0.3f, %0.3f]", shx, shy);
        status_callback(i,msg);
        fprintf(displ_file, "%d; %0.5f; %0.5f\n", i, shx, shy);

        if (conf->all_images) {
            t_im_struct final_copy;

            copy_image_with_weight(&final_copy, &final, IM_NOT_INITIALIZED);
            normalize_image_w(&final_copy, conf->full_weight_normalize);
            char outname[FN_LEN];
            snprintf(outname, FN_LEN, "%s%04un.tiff", out_fn_base, i);
            save_image(&final_copy, outname);
            destroy_image(&final_copy);
        }

        destroy_image(&h->image);
        h->invalid = 1;
        h = h->next;
        i++;
        if (iteration_callback((void *) &final)) break;

    }

    normalize_image_w(&final, conf->full_weight_normalize);
    snprintf(outname, FN_LEN, "%s.tiff", out_fn_base);

    save_image(&final, outname);
    iteration_callback(NULL); // this tells the main program that the processing ended.

    /*  CLEANUP */

    destroy_image(&final);
    fclose(displ_file);
    return ACC_ERR_NO_ERROR;
}/*}}}*/

void destroy_file_list(t_file_list *file_list){/*{{{*/
    t_file_list *h = file_list;
    while(h){
        if (!h->invalid) destroy_image(&h->image);
        file_list = file_list->next;
        free(h);
        h = file_list;
    }
}/*}}}*/
