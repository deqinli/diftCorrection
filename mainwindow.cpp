#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QDebug>
#include <QListWidget>
#include <Windows.h>
#include <opencv2/opencv.hpp>
#include <QTime>
#include <QCoreApplication>
#include <QImage>
#include <chrono>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->lineEdit_filePath->setReadOnly(true);
    ui->label_orgDisp->setScaledContents(true);
    ui->label_destDisp->setScaledContents(true);
    qstrFilesDir = "./";
    nNumOfCompositionImg = ui->spinBox_frameCount->value();
    setConnection();
    getImagePathDir();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::closeEvent(QCloseEvent* event)
{
    calThread.stop();
    calThread.wait();
}

void MainWindow::setConnection()
{
    connect(ui->list_files,SIGNAL(itemSelectionChanged()),this,SLOT(currentSelectFileNameChanged()));
    connect(ui->spinBox_frameCount,SIGNAL(valueChanged(int)),this,SLOT(setNumberOfCompositionImg(int)));
    connect(&calc_thread,SIGNAL(imgToGUI(QImage)), this, SLOT(update_image(QImage)));
    connect(this,&MainWindow::sg_dispImg,this,&MainWindow::sl_dispImg);
    connect(&calThread,SIGNAL(sg_resultImgToUIDisplay(QImage,QImage)),this,SLOT(sl_dispImg(QImage,QImage)));
}

void MainWindow::update_image(QImage image)
{
    //qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<"imageSize=("<<image.width()<<","<<image.height()<<")";
    //qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<"fileList_queue.last()=("<<fileList_queue.last();
    QImage img(fileList_queue.last());
    ui->label_orgDisp->setPixmap(QPixmap::fromImage(img));
    ui->label_destDisp->setPixmap(QPixmap::fromImage(image));
}

void MainWindow::sl_dispImg(QImage imgOrg, QImage imgDest)
{
    //qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<"::"<<__LINE__<<"imgOrg.width()="<<imgOrg.width()<<",imgDest.width()="<<imgDest.width();
    ui->label_orgDisp->setPixmap(QPixmap::fromImage(imgOrg));
    ui->label_destDisp->setPixmap(QPixmap::fromImage(imgDest));
}

void MainWindow::on_btnOpenFiles_clicked()
{

    QStringList qstrFileList = QFileDialog::getOpenFileNames(0,"打开图像",qstrFilesDir,"*.tif *tiff *.bmp *.jpg *.png");
    if(qstrFileList.size()>0)
    {
        ui->list_files->clear();
        qstrFilesDir = qstrFileList[0].mid(0,qstrFileList[0].lastIndexOf("/")+1);
        ui->lineEdit_filePath->setText(qstrFilesDir);
        qImg.load(qstrFileList[0]);
        ui->label_orgDisp->setPixmap(QPixmap::fromImage(qImg));
        ui->label_destDisp->setPixmap(QPixmap::fromImage(qImg));
        setImagePathDir();
    }
    for(int i=0;i<qstrFileList.size();i++)
    {
        QString qstrFileName = qstrFileList[i].mid(qstrFileList[0].lastIndexOf("/")+1,-1);

        ui->list_files->addItem(qstrFileName);
    }
}

void MainWindow::screen2conf()
{
    t_file_list *im_list = NULL;
    t_file_list *h = NULL;
    for( int i=0; i<ui->list_files->count(); i++){
        t_file_list *newitem;
        newitem = (t_file_list *) malloc(sizeof(t_file_list));
        newitem->next = NULL;
        strncpy(newitem->filename, ui->list_files->item(i)->text().toLocal8Bit().data(), FN_LEN);
        newitem->invalid = 1;
        if (!h) {
            im_list = newitem;
            h = newitem;
        } else {
            h->next = newitem;
            h = h->next;
        }
    }
    file_list = im_list;
}

void MainWindow::currentSelectFileNameChanged()
{
    QList fileList = ui->list_files->selectedItems();
    if(fileList.size()==0)
    {
        return;
    }
    qstrCurrentFileName = fileList[0]->text();
    qDebug() << qstrCurrentFileName;
    QString absFileName = qstrFilesDir + qstrCurrentFileName;
    qImg.load(absFileName);
    ui->label_orgDisp->setPixmap(QPixmap::fromImage(qImg));
    ui->label_destDisp->setPixmap(QPixmap::fromImage(qImg));
}

void MainWindow::on_btnStartProcess_clicked()
{

    int nCount = ui->list_files->count();

    fileList_queue.clear();
    for(int i=0;i<nCount;i++)
    {
        // 1-获取所有图像文件名
        QString FileName = ui->list_files->item(i)->text();
        QString absFileName = qstrFilesDir + FileName;
        if(fileList_queue.size() >= nNumOfCompositionImg)
        {
            fileList_queue.dequeue();

        }
        // 2-按照设置图像数量将图像名装入固定大小的队列容器
        fileList_queue.enqueue(absFileName);

        if(fileList_queue.size()<=0)
        {
            return;
        }

        // 如果队列中只有一幅图，则直接显示图像。
        if(fileList_queue.size()==1)
        {
            QImage img(fileList_queue.head());
            ui->label_orgDisp->setPixmap(QPixmap::fromImage(img));
            ui->label_destDisp->setPixmap(QPixmap::fromImage(img));
            continue;
        }

        for(int i=0; i<fileList_queue.size();i++)
        {
            QString qstrFilename = fileList_queue.at(i);
            qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<qstrFilename;
        }
        qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<"-----------------------------------------------";

#if 1
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
        file_list = im_list_queue;
        t_file_list* h3 = im_list_queue;

        if(!h3)
        {
            return;
        }
        auto start_loadImg = std::chrono::high_resolution_clock::now();

        while(load_image(&h3->image, h3->filename , IM_NOT_INITIALIZED))
        {
            h3 = h3->next;
        }

        auto end_loadImg = std::chrono::high_resolution_clock::now();
        auto duration_loadImg = std::chrono::duration_cast<std::chrono::milliseconds>(end_loadImg - start_loadImg);
        qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<"load images cost time: "<<duration_loadImg.count() << "ms";


        //test wether the image data is empty or not

        t_im_struct final;

        copy_image(&final, &h3->image, IM_NOT_INITIALIZED);
        initialize_weight_counting(&final);

        auto start_deal = std::chrono::high_resolution_clock::now();

        while(h3)
        {
            char err = ACC_ERR_NO_ERROR; // 程序执行状态，错误类型
            char res = load_image(&h3->image, h3->filename , IM_NOT_INITIALIZED);
            if (res)
            {
                fprintf(stderr, "Cannot open file: \"%s\". Ignoring it.", h3->filename);
                h3->invalid = 1;
                h3 = h3->next;
                continue;
            }
            h3->invalid = 0;
            double shx = 0;
            double shy = 0;
            err = add_with_correction(&final, &h3->image, IM_DESTRUCTIVE, &shx, &shy);
            if (err)
            {
                return;
            }


            destroy_image(&h3->image);
            h3->invalid = 1;
            h3 = h3->next;
        }

        processImg.normalize_image_w(&final, 1);
        if (image_buffer) delete image_buffer;
        image_buffer = processImg.image_to_32bit(&final);
        QImage imgDest(image_buffer, final.sizex, final.sizey, QImage::Format_RGB32);

        auto end_deal = std::chrono::high_resolution_clock::now();

        // 计算时间差（毫秒）
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_deal - start_deal);
        qDebug()<<"ldq "<<__FILE__<<":"<<__FUNCTION__<<" "<<__LINE__<<"deal images cost time: "<<duration.count()<<"ms";
#endif

        QImage img(fileList_queue.back());

        emit sg_dispImg(img,imgDest);
        int nMs = ui->spinBox_delayTime->value();
        mySleep(nMs);

        destroy_file_list(file_list);
        file_list = nullptr;

    }

}

void MainWindow::on_btnStartProcessMultiThread_clicked()
{

    if(calThread.isStopped())
    {
        calThread.restart();
    }

    int nCount = ui->list_files->count();

    fileNameList_abs.clear();
    for(int i=0;i<nCount;i++)
    {
        // 1-获取所有图像文件名
        fileNameList_abs.append(ui->list_files->item(i)->text());
    }

    calThread.setTask(fileNameList_abs,qstrFilesDir,nNumOfCompositionImg);
    calThread.start();


}

void MainWindow::on_btnStop_clicked()
{
    calThread.stop();
}

void MainWindow::mySleep(int ms)
{
    QTime dtime = QTime::currentTime().addMSecs(ms);
    while(QTime::currentTime() < dtime)
    {
        QCoreApplication::processEvents(QEventLoop::AllEvents,100);
    }
}

void MainWindow::setNumberOfCompositionImg(int n)
{
    nNumOfCompositionImg = n;
    qDebug() << __FUNCTION__<<"::"<<__LINE__<<"::nNumOfCompositionImg="<<nNumOfCompositionImg;
}

void MainWindow::setImagePathDir()
{
    QSettings config(CONFIG_INI,QSettings::IniFormat);
    config.setValue("PARAM/ImageDir",qstrFilesDir);
}

void MainWindow::getImagePathDir()
{
    QSettings config(CONFIG_INI,QSettings::IniFormat);
    qstrFilesDir = config.value("PARAM/ImageDir","").toString();
    ui->lineEdit_filePath->setText(qstrFilesDir);
}
