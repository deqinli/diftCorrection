#ifndef IMAGEUTILS_H
#define IMAGEUTILS_H

#include <QObject>
#include <QWidget>
#include <QThread>
#include "process.h"
#include "accord.h"
#include <QMutex>


// class qCalcThread : public QThread{
//     Q_OBJECT

// public:
//     void run();
//     void set_task(t_file_list* file_list, t_conf_struct* conf);
//     void iteration_callback(void* data);
//     void status_callback(int i,char* message);

// signals:
//     void image_ready(t_im_struct* image);
//     void status_output(int i,char* message);

//     void processing_finished();

// protected:
//     QMutex mutex;
//     t_file_list* file_list;
//     t_conf_struct* conf;
// };


#endif // IMAGEUTILS_H
