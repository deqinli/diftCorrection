#include "imageutils.h"

// CALC THREAD **********

// void qCalcThread::run(){/*{{{*/
//     process_images(file_list, conf);
// }/*}}}*/

// void qCalcThread::set_task(t_file_list *file_list, t_conf_struct *conf){/*{{{*/
//     this->file_list = file_list;
//     this->conf = conf;
// }/*}}}*/

// void qCalcThread::iteration_callback(void *data){/*{{{*/
//     if (!data) { // this is a message from the process.c that the processing has finished. We can clean up.
//         emit processing_finished();
//         return;
//     }
//     mutex.lock();
//     t_im_struct *image = (t_im_struct *)data;
//     t_im_struct *image_copy = new t_im_struct;
//     copy_image_with_weight(image_copy, image, IM_NOT_INITIALIZED);
//     emit image_ready(image_copy);
//     mutex.unlock();
// }/*}}}*/

// void qCalcThread::status_callback(int i, char *message){/*{{{*/
//     char *msg = new char[500];
//     strncpy(msg, message, 500);
//     emit status_output(i, msg);
// }/*}}}*/
