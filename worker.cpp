#include "worker.h"
#include "mainwindow.h"
#include <QPlainTextEdit>
#include <cmath>
#include <ctime>
#include <random>

#define PI 3.1415926536

int sum1, sum2, sum3;
int best1, best2, best3;
int besta, bestb, bestc;

#define PI 3.1415926536

void Worker::doWork(int averageTimes, double density, int molecularMass, int containerVolume, double duration, double precision,double temperature,double colding_factor)
{
    srand((int)time(0)); //时间作为随机数种子
       double T = temperature;
       double v = containerVolume;

       int a, b;
       double c;
       a = (1.0 * rand() / (RAND_MAX)) * pow(v,1.0/3) ;
       while (a < 6)  a = (1.0 * rand() / (RAND_MAX)) * pow(v, 1.0 / 3);
       b = (1.0 * rand() / (RAND_MAX)) * pow(v,1.0/3);
       while (b < 6)  b = (1.0 * rand() / (RAND_MAX)) * pow(v, 1.0 / 3);
       c = v / a / b;
       double bound = pow(v,1.0/5);
       int step = pow(v,1.0/3) / 10;
       if(step <= 0) step = 1;

       while (T > 1)
       {
           printText(tr("当前长宽高为 %1 %2 %3").arg(a).arg(b).arg(c));
           Worker::process(averageTimes, density, molecularMass, v, duration , precision,a,b,c);
           if (sum2 > best2)
           {
               best1 = sum1;
               best2 = sum2;
               best3 = sum3;
               besta = a;
               bestb = b;
               bestc = c;
               int t = (1.0 * rand() / (RAND_MAX)) * step;
               if (rand() % 2) a += 1;
               else a -= 1;
               if (a < bound) a = bound;

               t = (1.0 * rand() / (RAND_MAX)) * b / step;
               if (rand() % 2) b += 1;
               else b -= 1;
               if (b < bound) b = bound;
               c = 1.0 * v / a / b;
               if (a > b)
               {
                   int temp = a;
                   a = b;
                   b = temp;
               }
           }else if ((1.0 * rand() / (RAND_MAX)) <= pow(2.718, -(best2 - sum2) / T))
           {
               best1 = sum1;
               best2 = sum2;
               best3 = sum3;
               besta = a;
               bestb = b;
               bestc = c;

               int t = (1.0 * rand() / (RAND_MAX)) * step;
               if (rand() % 2) a += 1;
               else a -= 1;
               if (a < bound) a = bound;

               t = (1.0 * rand() / (RAND_MAX)) * b / step;
               if (rand() % 2) b += 1;
               else b -= 1;
               if (b < bound) b = bound;
               c = 1.0 * v / a / b;
               if (a > b)
               {
                   int temp = a;
                   a = b;
                   b = temp;
               }
           }
           else
           {
               int t = (1.0 * rand() / (RAND_MAX)) * step;
               if (rand() % 2) a += 1;
               else a -= 1;
               if (a < bound) a = bound;

               t = (1.0 * rand() / (RAND_MAX)) * b / step;
               if (rand() % 2) b += 1;
               else b -= 1;
               if (b < bound) b = bound;
               c = 1.0 * v / a / b;
               if (a > b)
               {
                   int temp = a;
                   a = b;
                   b = temp;
               }
           }
           T = T / colding_factor;
       }

       printText(tr("\n当前模拟状态下最佳长宽高为 %1 %2 %3\n此时沉积在容器中的粒子数为 %4").arg(besta).arg(bestb).arg(bestc).arg(best2));
       Result result(besta, bestb, bestc, best2);
       emit resultReady(result);


    }

void Worker::process(int calculate_times, double density, int molecular_mass, int volume, double duration, double precision, int AA, int BB, double CC)
{

    int n1=0,n2=0,n3=0;
    for (int times = 0; times < calculate_times; times++)
    {
        //计算半径
        double R = sqrt((AA * AA + BB * BB + CC * CC)) / 2;

        //计算生成的粒子数
        int n = (int)(18.3275 * duration * R * R + 1);


        for (int t = 0; t < n; t++)
        {
            //随机粒子的初始位置
            double r1 = (1.0 * rand() / (RAND_MAX));
            double theta1 = acos(r1);/*粒子位于球面上任意一点(等概率),对于天顶角，能量的分布,下面再表示，0<=cos(theta1)<=1*/
            double r2 = (1.0 * rand() / (RAND_MAX));
            double theta2 = 2 * PI * r2;/*方位角*/
            double x1 = R * sin(theta1) * cos(theta2);
            double x2 = R * sin(theta1) * sin(theta2);
            double x3 = R * cos(theta1);

            /*蒙特卡罗方法选出一个符合分布的粒子*/
            double p1 = 0.102573, p2 = -0.068287, p3 = 0.958633, p4 = 0.0407253, p5 = 0.817285;
            double e, theta_temp;
            for (int j = 1; j <= 1000; j++) {
                double r3 = (1.0 * rand() / (RAND_MAX));
                double r4 = (1.0 * rand() / (RAND_MAX));
                double r5 = (1.0 * rand() / (RAND_MAX));
                double a = r3;/*角度*/
                double d = (1 - 0.105658) * r4 + 0.105658;/*能量*/  /*最小的静止能量为0.105658*/
                double b = pow((a * a + p1 * p1 + p2 * pow(a, p3) + p4 * pow(a, p5)) / (1 + p1 * p1 + p2 + p4), 0.5);/*角度换算*/
                double x = (3 * d + 7 / b) / 10;
                double y = 2.06 / 1000 * (950.0 / b - 90) + x;
                double AA = 1.1 * pow(90.0 / 1030 * pow(0.001 + a, 0.5), 4.5 / (y * b));
                double u = pow(10.0, 3) * AA * 0.14 * pow(x, -2.70) * ((1 / (1 + 1.1 * y * b / 115) + 0.054 / (1 + 1.1 * y * b / 850)) + 0.0001);
                if (r5 * 4.61681 <= u)
                {
                    /*4.61681为最大值,选出在曲线之下的点*/
                    e = d;/*能量:以Gev为单位*/
                    theta_temp = acos(a);/*选出天顶角*/
                    break;
                }
            }
            double M = 105.658;
            double energy = (e - 0.105658) * 1000.0 / M;/*1000表示化为以Mev为单位*/
            double h = pow(energy * (energy + 2), 0.5); /*h与能量之间的转换关系*/
            double theta_sky = theta_temp;
            double r6 = (1.0 * rand() / (RAND_MAX));
            double theta_direction = 2 * PI * r6;


            //改进:当粒子没有进入长方体时,不必逐步逼近,可以直接计算.

            //由对称性可知 转化坐标到第一象限等价原问题

//            x1 = fabs(x1);
//            x2 = fabs(x2);
//            x3 = fabs(x3);
//            theta_sky = theta_sky < (PI - theta_sky) ? theta_sky : (PI - theta_sky);


//            double RR1 = (x3 - CC/2) / cos(theta_sky);

//                x1 = x1 - RR1 * sin(theta_sky) * cos(theta_direction);
//                x2 = x2 - RR1 * sin(theta_sky) * sin(theta_direction);
//                x3 = x3 - RR1 * cos(theta_sky);

//             double RR2 = (x1 - AA/2) / sin(theta_sky) / cos(theta_direction);
//                double x33 = x3 - RR2 * cos(theta_sky);
//                double  x11 = x1 - RR2 * sin(theta_sky) * cos(theta_direction);
//                double x22 = x2 - RR2 * sin(theta_sky) * sin(theta_direction);

//                double RR3 = (x2 - BB/2) / sin(theta_sky) / sin(theta_direction);
//                   double x333 = x3 - RR2 * cos(theta_sky);
//                   double  x111 = x1 - RR2 * sin(theta_sky) * cos(theta_direction);
//                   double x222 = x2 - RR2 * sin(theta_sky) * sin(theta_direction);

            for (int qq = 0; qq < (int)R * 2; qq++)
                       {
                           x1 = x1 - 1 * sin(theta_sky) * cos(theta_direction);/*未进入探测器前，不妨将步长取长一些*/
                           x2 = x2 - 1 * sin(theta_sky) * sin(theta_direction);
                           x3 = x3 - 1 * cos(theta_sky);
                           if ((fabs(x1) <= AA / 2) && (fabs(x2) <= BB / 2) && (fabs(x3) <= CC / 2))
                           {
                               n1++;
                               x1 = x1 + 1 * sin(theta_sky) * cos(theta_direction);
                               x2 = x2 + 1 * sin(theta_sky) * sin(theta_direction);
                               x3 = x3 + 1 * cos(theta_sky);

                               for (int l = 0; l <= 1000; l++)
                               {
                                   x1 = x1 - 0.001 * sin(theta_sky) * cos(theta_direction);
                                   x2 = x2 - 0.001 * sin(theta_sky) * sin(theta_direction);
                                   x3 = x3 - 0.001 * cos(theta_sky);
                                   if ((fabs(x1) <= AA / 2) && (fabs(x2) <= BB / 2) && (fabs(x3) <= CC / 2))
                                       break;
                               }

                    //变步长辛普森改为自适应算法
                    auto calculate = [](double f, int z, double h) {
                        return (-0.30705 * f * z / 1 + 0.153525 * f * z / 1 * (1 + pow(h, -2)) *       /*能量损失*/  /*除以1，1代表阿伏伽德罗常数*/
                            log(1.0404 * pow(10.0, 10) * pow(h, 4) / pow((double)z, 2) / (1 + 0.009680542 * pow(1 + pow(h, 2), 0.5) + 2.34282 * pow(10.0, -5)
                                ))) * pow(1 + pow(h, 2), 0.5) / h;
                    };
                    int z1 = 6, z2 = 1, m1 = 8, m2 = 8;
                    double f1 = m1 * density * 1 / molecular_mass;/*每立方厘米含有的原子个数 以后可以修改 阿伏伽德罗常数用1代替*/
                    double f2 = m2 * density * 1 / molecular_mass;
                    double s1;
                    double a = 0, b = h, A = (b - a) / 2, B = (b + a) / 2;
                    double S = 0, S1 = 0;/*S1记录前一次的积分值，S记录下一次的积分值*/
                    double RP;
                    double s = calculate(f1, z1, -A + B) + calculate(f2, z2, -A + B);
                    s1 = calculate(f1, z1, A + B) + calculate(f2, z2, A + B);
                    s = M / s;
                    s1 = M / s1;
                    RP = s + s1;
                    double RC = 0, H = 2 * A;
                    int M1, M2;
                    for (M1 = 1; M1 < 40; M1++)
                    {
                        RC = 0;
                        S1 = S;
                        for (M2 = 1; M2 <= pow(2.0, M1 - 1); M2++)
                        {
                            s = calculate(f1, z1, -A - H / 2 + M2 * H + B) + calculate(f2, z2, -A - H / 2 + M2 * H + B);
                            s = M / s;
                            RC += s;
                        }
                        S = H / 6 * (RP + 4 * RC);
                        H = H / 2;
                        RP += 2 * RC;

                        if (fabs(S1 - S) <= e)
                        {

                            break;
                        }
                    }
                    //下面三行计算粒子动能为0时的位置坐标
                    x1 = x1 - S * sin(theta_sky) * cos(theta_direction);
                    x2 = x2 - S * sin(theta_sky) * sin(theta_direction);
                    x3 = x3 - S * cos(theta_sky);
                    if ((fabs(x1) <= AA / 2) && (fabs(x2) <= BB / 2) && (fabs(x3) <= CC / 2))
                        n2++;
                    else
                        n3++;
                    break;
                }
            }

        }
    }
        sum1 = n1 / calculate_times;
        sum2 = n2 / calculate_times;
        sum3 = n3 / calculate_times;

         printText(tr("进入探测器的粒子数为 %1\n沉积在探测器中的粒子数为 %2\n进入探测器后又逸出的粒子数为 %3").arg(sum1).arg(sum2).arg(sum3));


}
