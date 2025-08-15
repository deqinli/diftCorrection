#ifndef CALCUTHREAD_H
#define CALCUTHREAD_H

#include <QObject>
#include <QThread>
#include <QMutex>
#include <QList>
#include <QQueue>
#include "process.h"
#include <QImage>
#include <chrono>

class CalcuThread : public QThread
{
    Q_OBJECT
public:
    void setTask(QList<QString> fileList, QString qstrFilesDir, int nNumOfCompositionImg=8);
    void stop();
    bool isStopped();
    void restart();

signals:
    void sg_resultImgToUIDisplay(QImage imgOrg, QImage imgResult);

private:
    void normalize_image_w(t_im_struct *image, char full_weight_only);
    unsigned char* image_to_32bit(t_im_struct* image);



    t_file_list* file_list;
    QQueue<QString> fileList_queue;
    QList<QString> file_list_widget;
    int nNumOfCompositionImg;
    QString qstrFilesDir;
    unsigned char *image_buffer = nullptr;

    QMutex m_mutex;
    bool m_stopped = false;
protected:
    void run() override;
};

#endif // CALCUTHREAD_H
