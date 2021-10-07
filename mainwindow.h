#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QThread>
#include <QMetaType>
#include <QLabel>
#include <qlabel.h>
#include <qmovie.h>

QT_BEGIN_NAMESPACE
class QWidget;
class QSpinBox;
class QDoubleSpinBox;
class QLineEdit;
class QGroupBox;
class QPushButton;
class QPlainTextEdit;
QT_END_NAMESPACE

struct Result {
    int length, width, height, particleNumber;

    Result() :
        length(0), width(0), height(0), particleNumber(0) {}

    Result(int l, int w, int h, int n) :
        length(l), width(w), height(h), particleNumber(n) {}
};
Q_DECLARE_METATYPE(Result)

class MainWindow : public QMainWindow
{
    Q_OBJECT

    QThread *workerThread;
    bool isWorking;

    QLineEdit *averageTimesSpinBox;
    QLineEdit *durationSpinBox;
    QLineEdit *densityLineEdit;
    QLineEdit *molecularMassLineEdit;
    QLineEdit *containerVolumeLineEdit;
    QLineEdit *precisionLineEdit;
    QLineEdit *temperatureLineEdit;
    QLineEdit *colding_factorLineEdit;
    QLabel *chanel1;
    QLabel *chanel2;
    QLabel *chanel3;
    QLabel *chanel4;
    QLabel *chanel5;
    QImage *image;

    QGroupBox *formGroupBox;

    QPushButton *startButton;
    QPushButton *resetButton;
    QPushButton *quitButton;
    QWidget *buttons;

    QPlainTextEdit *resultText;

    void createFormGroupBox();
    void createButtons();

public:
    MainWindow(QWidget *parent = 0);

signals:
    void operate(int averageTimes, double density, int molecularMass, int containerVolume, double duration, double precision,double Temperature,double colding_factor);

private slots:
    void startSimulation();
    void resetFormGroupBox();

public slots:
    void handleResult(Result result);
    void updateResultText(QString string);
};

#endif // MAINWINDOW_H
