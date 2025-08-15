#include "myprocess.h"

#include <QImage>
#include <QDebug>
#include <QString>

void qCalcThread::run(){
    //process_images(file_list);
    QImage img;

    while(file_list->filename)
    {
        std::string strFilename = file_list->filename;
        qDebug() <<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<"file_list->filename = "<<QString::fromStdString(strFilename);
    }
    process.process_image(file_list, img);
    emit imgToGUI(img);
    this->quit();
}/*}}}*/

void qCalcThread::set_task(t_file_list *file_list){
    this->file_list = file_list;
}

void qCalcThread::iteration_callback(void *data){
    if (!data) { // this is a message from the process.c that the processing has finished. We can clean up.
        emit processing_finished();
        return;
    }
    mutex.lock();
    t_im_struct *image = (t_im_struct *)data;
    t_im_struct *image_copy = new t_im_struct;
    copy_image_with_weight(image_copy, image, IM_NOT_INITIALIZED);
    emit image_ready(image_copy);
    mutex.unlock();
}


myProcess::myProcess() {}

myProcess::~myProcess() {}

char myProcess::process_image(t_file_list* img_list, QImage& qImg)
{
    t_file_list* h = img_list;
    t_im_struct final;

    // while (load_image(&h->image, h->filename , IM_NOT_INITIALIZED))
    // {
    //     h = h->next;
    // }

    while(h)
    {
        char err = ACC_ERR_NO_ERROR; // 程序执行状态，错误类型
        char res = load_image(&h->image, h->filename , IM_NOT_INITIALIZED);
        if (res)
        {
            fprintf(stderr, "Cannot open file: \"%s\". Ignoring it.", h->filename);
            h->invalid = 1;
            h = h->next;
            continue;
        }
        h->invalid = 0;
        double shx = 0;
        double shy = 0;
        err = add_with_correction(&final, &h->image, IM_DESTRUCTIVE, &shx, &shy);
        if (err)
        {
            return err;
        }


        destroy_image(&h->image);
        h->invalid = 1;
        h = h->next;
    }

    //t_im_struct *image = &final;

    if (image_buffer) delete image_buffer;
    image_buffer = image_to_32bit(&final);
    QImage im(image_buffer, final.sizex, final.sizey, QImage::Format_RGB32);
    qImg = im.copy();
    //destroy_image(&final);
    return ACC_ERR_NO_ERROR;
}


unsigned char* myProcess::image_to_32bit(t_im_struct* image)
{
    unsigned char *buffer = new unsigned char[4*image->sizex * image->sizey];
    for (unsigned long int i=0; i< image->sizex * image->sizey; i++) {
        buffer[4*i+3] = 0xff;
        for (int k=0; k<3; k++)
            buffer[4*i+k] = (unsigned char)(image->buffer[i] * 255);
    }
    return buffer;
}

void myProcess::normalize_image_w(t_im_struct *image, char full_weight_only){/*{{{*/
    unsigned long int i;
    unsigned int max_weight;

    ACC_IM_TYPE max_value = image->buffer[0]; // initialization
    ACC_IM_TYPE min_value = image->buffer[0]; // initialization

    if (image->weight_buffer) {
        max_weight = image->weight_buffer[0]; // initialization as well
        for (i=0; i<IM_SIZE; i++) if (image->weight_buffer[i] > max_weight) max_weight = image->weight_buffer[i];
        for (i=0; i<IM_SIZE; i++) image->buffer[i] *= (ACC_IM_TYPE) max_weight / image->weight_buffer[i];
    }


    for (i=0; i<IM_SIZE; i++)
        if ((!image->weight_buffer) || (!full_weight_only) || (image->weight_buffer[i] == max_weight)) {
            if (image->buffer[i] > max_value) max_value = image->buffer[i];
            if (image->buffer[i] < min_value) min_value = image->buffer[i];
        }

    if (image->weight_buffer) {
        free(image->weight_buffer);
        image->weight_buffer = NULL;
    }

    if (max_value - min_value)
        for (i=0; i<IM_SIZE; i++) {
            image->buffer[i] = (image->buffer[i] - min_value) / (max_value - min_value);
            if (image->buffer[i] > 1) image->buffer[i] = 1; // overflow protection
            if (image->buffer[i] < 0) image->buffer[i] = 0; // underflow protection
        }

    image->fourier_ok = 0;
}/*}}}*/
