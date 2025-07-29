#ifndef MYPROCESS_H
#define MYPROCESS_H
#include <opencv2/opencv.hpp>
#include <QThread>
#include <QMutex>
#include <QObject>
#include <QImage>
#include "accord.h"

#define MAX_FILENAME_LEN 800

class myProcess
{
public:
    myProcess();
    ~myProcess();

    char process_image(t_file_list* img_list, QImage& qImg);

    unsigned char * image_to_32bit(t_im_struct* image);
    void normalize_image_w(t_im_struct *image, char full_weight_only);
private:

    unsigned char* image_buffer;

};

class qCalcThread : public QThread
{
    Q_OBJECT

public:
    void run();
    void set_task(t_file_list* file_list);
    void iteration_callback(void* data);

signals:
    void image_ready(t_im_struct* image);
    void status_output(int i,char* message);

    void processing_finished();
    void imgToGUI(QImage img);

private:
    myProcess process;
    void setConnection();

protected:
    QMutex mutex;
    t_file_list* file_list;
    t_conf_struct* conf;
};



#endif // MYPROCESS_H
