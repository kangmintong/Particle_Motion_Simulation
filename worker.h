#ifndef WORKER_H
#define WORKER_H

#include "mainwindow.h"

class Worker : public QObject
{
    Q_OBJECT

signals:
    void resultReady(Result result);
    void printText(QString string);

public slots:
    void doWork(int averageTimes, double density, int molecularMass, int containerVolume, double duration, double precision,double temperature,double colding_factor);
    void process(int averageTimes, double density, int molecularMass, int containerVolume, double duration, double precision,int AA,int BB,double CC);

};

#endif // WORKER_H
