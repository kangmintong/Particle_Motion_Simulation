#include "mainwindow.h"
#include "worker.h"
#include <QtWidgets>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    workerThread(new QThread),
    isWorking(false),
    resultText(new QPlainTextEdit)
{
    createFormGroupBox();
    createButtons();

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(formGroupBox);
    mainLayout->addWidget(buttons);
    mainLayout->addWidget(resultText);
    mainLayout->setAlignment(Qt::AlignCenter);

    QWidget *widget = new QWidget;
    widget->setLayout(mainLayout);

    setWindowTitle(tr("宇宙射线μ子蒙特卡洛模拟 (浙江大学)"));
    setCentralWidget(widget);
    setMinimumWidth(1000);
    setMinimumHeight(700);
    setStyleSheet("background-color:rgb(247,238,214);");

    Worker *worker = new Worker;
    worker->moveToThread(workerThread);
    connect(this, &MainWindow::operate, worker, &Worker::doWork);
    connect(worker, &Worker::resultReady, this, &MainWindow::handleResult);
    connect(worker, &Worker::printText, this, &MainWindow::updateResultText);
    workerThread->start();

    qRegisterMetaType<Result>();
}

void MainWindow::createFormGroupBox()
{
    chanel2 = new QLabel(tr("智能算法参数设定："));
    temperatureLineEdit = new QLineEdit;
    colding_factorLineEdit = new QLineEdit;

    chanel1 = new QLabel(tr("实验参数设定："));
    chanel3 = new QLabel(tr("实验精度设定："));
    chanel4 = new QLabel(tr("说明：当前程序只适用于长方体容器的探测模拟。"));

    densityLineEdit = new QLineEdit;
    molecularMassLineEdit = new QLineEdit;
    containerVolumeLineEdit = new QLineEdit;
    precisionLineEdit = new QLineEdit;
    averageTimesSpinBox = new QLineEdit;
    durationSpinBox = new QLineEdit;
//    image = new QImage("./zheda.jpg");
//    image->load("./zheda.jpg");
//    chanel5->setPixmap(QPixmap::fromImage(*image));



    averageTimesSpinBox->setText("5");
    durationSpinBox->setText("1.0");
    densityLineEdit->setText("1.06");
    molecularMassLineEdit->setText("104");
    containerVolumeLineEdit->setText("4000");
    precisionLineEdit->setText("0.01");
    temperatureLineEdit->setText("100");
    colding_factorLineEdit->setText("1.5");

    QFormLayout *layout = new QFormLayout;
//    chanel1->setGeometry(rect().x()+50,rect().y()+30,20,20);
//    chanel1->setParent(this);

    //layout->addWidget(chanel5);
    layout->addRow(new QLabel(tr("")),chanel1);
    layout->addRow(new QLabel(tr("模拟退火初始温度")),temperatureLineEdit);
    layout->addRow(new QLabel(tr("冷却系数")),colding_factorLineEdit);
    layout->addRow(new QLabel(tr("")),chanel2);
    layout->addRow(new QLabel(tr("靶物质密度（g/cm³）")), densityLineEdit);
    layout->addRow(new QLabel(tr("相对分子质量")), molecularMassLineEdit);
    layout->addRow(new QLabel(tr("容器体积（cm³）")), containerVolumeLineEdit);
    layout->addRow(new QLabel(tr("持续小时(h)")), durationSpinBox);
    layout->addRow(new QLabel(tr("")),chanel3);
    layout->addRow(new QLabel(tr("平均次数")), averageTimesSpinBox);
    layout->addRow(new QLabel(tr("精度(辛普森积分的精度)")), precisionLineEdit);
    layout->addRow(new QLabel(tr("")),chanel4);


    formGroupBox = new QGroupBox;
    formGroupBox->setLayout(layout);
}

void MainWindow::resetFormGroupBox()
{
    averageTimesSpinBox->setText("0");
    durationSpinBox->setText("0");
    densityLineEdit->setText("");
    molecularMassLineEdit->setText("");
    containerVolumeLineEdit->setText("");
    precisionLineEdit->setText("");
}

void MainWindow::createButtons()
{
    startButton = new QPushButton(tr("开始模拟"));
    connect(startButton, &QPushButton::clicked, this, &MainWindow::startSimulation);

    resetButton = new QPushButton(tr("重置参数"));
    connect(resetButton, &QPushButton::clicked, this, &MainWindow::resetFormGroupBox);

    quitButton = new QPushButton(tr("退出程序"));
    connect(quitButton, &QPushButton::clicked, qApp, &QApplication::quit);

    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(startButton);
    layout->addWidget(resetButton);
    layout->addWidget(quitButton);

    buttons = new QWidget;
    buttons->setLayout(layout);
}

void MainWindow::startSimulation()
{
    int averageTimes = averageTimesSpinBox->text().toDouble();
    double density = densityLineEdit->text().toDouble();
    int molecularMass = molecularMassLineEdit->text().toInt();
    int containerVolume = containerVolumeLineEdit->text().toInt();
    double duration = durationSpinBox->text().toDouble();
    double precision = precisionLineEdit->text().toDouble();
    double temperature = temperatureLineEdit->text().toDouble();
    double colding_factor = colding_factorLineEdit->text().toDouble();


    if (averageTimes == 0 || density == 0 || molecularMass == 0 || containerVolume == 0 || duration == 0 || precision == 0) {
        QMessageBox::warning(this, tr("参数错误"), tr("请输入正确的参数"), QMessageBox::Yes, QMessageBox::Yes);
    } else if (isWorking) {
        QMessageBox::warning(this, tr("正在模拟"), tr("本次模拟正在运行中……"), QMessageBox::Yes, QMessageBox::Yes);
    } else {
        isWorking = true;
        resultText->setPlainText("开始模拟……");
        operate(averageTimes, density, molecularMass, containerVolume, duration, precision,temperature,colding_factor);
    }
}

void MainWindow::handleResult(Result result)
{
    QMessageBox::information(this, tr("模拟结束"), tr("最佳长宽高为 %1, %2, %3\n此时沉积在容器中的粒子数为 %4").arg(result.length).arg(result.width).arg(result.height).arg(result.particleNumber), QMessageBox::Yes, QMessageBox::Yes);
    isWorking = false;
}

void MainWindow::updateResultText(QString string)
{
    resultText->appendPlainText(string);
}
