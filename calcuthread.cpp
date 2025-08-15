#include "calcuthread.h"
#include <QDebug>


void CalcuThread::setTask(QList<QString> fileList, QString qstrFilesDir, int nNumOfCompositionImg)
{
    this->file_list_widget = fileList;
    this->nNumOfCompositionImg = nNumOfCompositionImg;
    this->qstrFilesDir = qstrFilesDir;
}

void CalcuThread::normalize_image_w(t_im_struct *image, char full_weight_only){
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
}

unsigned char* CalcuThread::image_to_32bit(t_im_struct* image)
{
    unsigned char *buffer = new unsigned char[4*image->sizex * image->sizey];
    for (unsigned long int i=0; i< image->sizex * image->sizey; i++) {
        buffer[4*i+3] = 0xff;
        for (int k=0; k<3; k++)
            buffer[4*i+k] = (unsigned char)(image->buffer[i] * 255);
    }
    return buffer;
}

void CalcuThread::run()
{
    int nCount = file_list_widget.count();
    fileList_queue.clear();
    int nRecord=0;
    while((!isStopped())&&(nRecord<nCount))
    {
        QString FileName = file_list_widget.at(nRecord);
        QString absFileName = qstrFilesDir + FileName;
        if(fileList_queue.size() >= nNumOfCompositionImg)
        {
            fileList_queue.dequeue();
        }

        fileList_queue.enqueue(absFileName);

        if(fileList_queue.size()<nNumOfCompositionImg)
        {
            QImage img(fileList_queue.back());
            emit sg_resultImgToUIDisplay(img,img);
            ++nRecord;
            continue;
        }

        for(int i=0; i<fileList_queue.size();i++)
        {
            QString qstrFilename = fileList_queue.at(i);
        }

        // 3-将容器中所有图像进行漂移校正和融合
        t_file_list* im_list_queue = NULL;
        t_file_list* h = NULL;
        for(int i=0;i<fileList_queue.size();i++)
        {
            t_file_list* newitem;
            newitem = (t_file_list*)malloc(sizeof(t_file_list));
            newitem->next = nullptr;
            std::string srnm = fileList_queue.at(i).toStdString();
            strncpy_s(newitem->filename,srnm.c_str(),srnm.length());
            newitem->invalid = 1;
            if(!h)
            {
                im_list_queue = newitem;
                h = newitem;
            }
            else
            {
                h->next = newitem;
                h = h->next;
            }
        }
        //file_list = im_list_queue;
        //t_file_list* h3 = im_list_queue;
        h = im_list_queue;

        if(!h)
        {
            return;
        }
        auto start_loadImg = std::chrono::high_resolution_clock::now();

        while(load_image(&h->image, h->filename , IM_NOT_INITIALIZED))
        {
            h = h->next;
        }
        auto end_loadImg = std::chrono::high_resolution_clock::now();
        auto duration_loadImg = std::chrono::duration_cast<std::chrono::milliseconds>(end_loadImg - start_loadImg);
        qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<"load images cost time: "<<duration_loadImg.count() << "ms";


        //test wether the image data is empty or not

        t_im_struct final;

        copy_image(&final, &h->image, IM_NOT_INITIALIZED);
        initialize_weight_counting(&final);

        auto start_deal = std::chrono::high_resolution_clock::now();

        while(h&&(!isStopped()))
        {
            char err = ACC_ERR_NO_ERROR; // 程序执行状态，错误类型
            char res = load_image(&h->image, h->filename , IM_NOT_INITIALIZED);
            qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<"filename: "<<QString::fromStdString(h->filename);
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
                return;
            }

            destroy_image(&h->image);
            h->invalid = 1;
            h = h->next;
        }

        normalize_image_w(&final, 1);
        if (image_buffer) delete image_buffer;
        image_buffer = image_to_32bit(&final);
        QImage imgDest2(image_buffer, final.sizex, final.sizey, QImage::Format_RGB32);
        QImage imgDest = imgDest2.copy();

        auto end_deal = std::chrono::high_resolution_clock::now();

        // 计算时间差（毫秒）
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_deal - start_deal);
        qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<"deal images cost time: "<<duration.count()<<"ms";
        qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<"---------------------------------------------------------------";

        QImage img(fileList_queue.back());

        emit sg_resultImgToUIDisplay(img,imgDest);


        delete [] image_buffer;
        image_buffer = nullptr;

        destroy_image(&final);
        destroy_file_list(im_list_queue);

        ++nRecord;

    }

}

void CalcuThread::stop()
{
    QMutexLocker locker(&m_mutex);
    m_stopped = true;
}

bool CalcuThread::isStopped()
{
    QMutexLocker locker(&m_mutex);
    return m_stopped;
}

void CalcuThread::restart()
{
    QMutexLocker locker(&m_mutex);
    m_stopped = false;
}
