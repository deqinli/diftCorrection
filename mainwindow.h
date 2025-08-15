#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QImage>
#include <QPixmap>
#include <QList>
#include "process.h"
#include <QQueue>
#include "myprocess.h"
#include "calcuthread.h"
#include <QCloseEvent>


QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void setConnection();

    void mySleep(int ms);

    void screen2conf();

signals:
    void sg_dispImg(QImage imgOrg, QImage imgDest);

private slots:
    void update_image(QImage image);
    void sl_dispImg(QImage imgOrg, QImage imgDest);

private:
    Ui::MainWindow *ui;

    QString qstrCurrentFileName;
    QString qstrFilesDir;

    QImage qImg;
    QList<QString>  fileNameList_abs;
    t_file_list* file_list;
    int nNumOfCompositionImg;

    QQueue<QString> fileList_queue;
    qCalcThread calc_thread;
    t_conf_struct conf;
    char process_running;
    char stop_process;
    unsigned char *image_buffer = nullptr;

    myProcess processImg;

    CalcuThread calThread;


private slots:
    void on_btnOpenFiles_clicked();
    void currentSelectFileNameChanged();
    void on_btnStartProcess_clicked();
    void on_btnStartProcessMultiThread_clicked();
    void on_btnStop_clicked();
    void setNumberOfCompositionImg(int n);

protected:
    void closeEvent(QCloseEvent* event)override;

};
#endif // MAINWINDOW_H
