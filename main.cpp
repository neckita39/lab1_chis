#include <iostream>
#include "math.h"
const double eps= 0.0000001;
const double npoints=100.0;
double getFunc(double x){
    return 2.0*log(5.0*x)-1.0/(3.0*x);
}
double getTable(int a, int b){
    for (int i=(int)a; i<(int)b; i+=(int)((b-a)/npoints)){
        std::cout<<"point "<<i<<"; value "<<getFunc(i)<<std::endl;
    }
}
double getFuncFirstDer(double x){
    return 2.0/x+1/(3.0*x*x);
}

double getFuncSecondDer(double x){
    return -2.0/(x*x)-2.0/(3.0*x*x*x);
}
double f(double x){
    return exp(1.0/(6.0*x))/5.0;
}

double getBorder(double a, double b){
    if (getFunc(a)* getFuncSecondDer(a)>0)
    {
        return a;
    }
    return b;
}

void solveNewton(double a, double b)
{
    double fault=0;
    int cnt=0;
    double x[2];
    x[0]=getBorder(a, b);
    x[1]=x[0]- getFunc(x[0])/ getFuncFirstDer(x[0]);
    while (abs(x[0]-x[1])>eps)
    {
        cnt++;
        x[0]=x[1];
        x[1]=x[0]- getFunc(x[0])/ getFuncFirstDer(x[0]);
        fault= getFunc(x[0])/ getFuncFirstDer(b);
    }
    std::cout<<"solveNewton"<<std::endl;
    std::cout<<"iter="<<cnt<<std::endl;
    std::cout<<"x="<<x[1]<<std::endl;
    std::cout<<"fault="<<abs(fault)<<std::endl;
}
void solveChrod(double a, double b)
{
    double fault=0;
    int cnt=0;
    double x[2];
    x[0]=b;
    x[1]=x[0]- getFunc(x[0])*(a-x[0])/(getFunc(a)- getFunc(x[0]));
    while (abs(x[0]-x[1])>eps)
    {
        cnt++;
        x[0]=x[1];
        x[1]=x[0]- getFunc(x[0])*(a-x[0])/(getFunc(a)- getFunc(x[0]));
        fault= getFunc(x[0])/ getFuncFirstDer(b);
    }
    std::cout<<"solveChrod"<<std::endl;
    std::cout<<"iter="<<cnt<<std::endl;
    std::cout<<"x="<<x[1]<<std::endl;
    std::cout<<"fault="<<abs(fault)<<std::endl;
}
void solveSecant(double a, double b){
    double fault=0;
    int cnt=0;
    double x[3];
    x[0]=a;
    x[1]=b;
    x[2]=x[1]- getFunc(x[1])*(x[1]-x[0])/(getFunc(x[1])- getFunc(x[0]));
    while (abs(x[0]-x[1])>eps)
    {
        cnt++;
        x[0]=x[1];
        x[1]=x[2];
        x[2]=x[1]- getFunc(x[1])*(x[1]-x[0])/(getFunc(x[1])- getFunc(x[0]));
        fault= getFunc(x[0])/ getFuncFirstDer(b);
    }
    std::cout<<"solveSecant"<<std::endl;
    std::cout<<"iter="<<cnt<<std::endl;
    std::cout<<"x="<<x[1]<<std::endl;
    std::cout<<"fault="<<abs(fault)<<std::endl;
}
void solveDifNewton(double a, double b){
    double fault=0;
    double h=0.001;
    int cnt=0;
    double x[2];
    x[0]=a;
    x[1]=x[0]-h* getFunc(x[0])/(getFunc(x[0]+h)- getFunc(x[0]));
    while (abs(x[0]-x[1])>eps)
    {
        cnt++;
        x[0]=x[1];
        x[1]=x[0]-h* getFunc(x[0])/(getFunc(x[0]+h)- getFunc(x[0]));
        fault= getFunc(x[0])/ getFuncFirstDer(b);
    }
    std::cout<<"solveDifNewton"<<std::endl;
    std::cout<<"iter="<<cnt<<std::endl;
    std::cout<<"x="<<x[1]<<std::endl;
    std::cout<<"fault="<<abs(fault)<<std::endl;
}

void solveSteffanson(double b, double a=0.3001){
    double fault=0;
    int cnt=0;
    double x[2];
    x[0]=a;
    x[1]=x[0]- getFunc(x[0])* getFunc(x[0])/ (getFunc(x[0]+ getFunc(x[0]))- getFunc(x[0]));
    while (abs(x[0]-x[1])>eps)
    {
        cnt++;
        x[0]=x[1];
        x[1]=x[0]- getFunc(x[0])* getFunc(x[0])/ (getFunc(x[0]+ getFunc(x[0]))- getFunc(x[0]));
        fault= getFunc(x[0])/ getFuncFirstDer(b);
    }
    std::cout<<"solveSteffanson"<<std::endl;
    std::cout<<"iter="<<cnt<<std::endl;
    std::cout<<"x="<<x[1]<<std::endl;
    std::cout<<"fault="<<abs(fault)<<std::endl;
}
void solveIterations(double a, double b){
    double fault=0;
    int cnt=0;
    double x[2];
    x[0]=a;
    x[1]=b;
    while (abs(x[0]-x[1])>eps)
    {
        cnt++;
        x[1]=x[0];
        x[0]=f(x[0]);
        fault= getFunc(x[0])/ getFuncFirstDer(b);
    }
    std::cout<<"solveIterations"<<std::endl;
    std::cout<<"iter="<<cnt<<std::endl;
    std::cout<<"x="<<x[1]<<std::endl;
    std::cout<<"fault="<<abs(fault)<<std::endl;
}
int main() {
    std::cout.precision(8);
    double a=0.1;
    double b=1.0;
    solveNewton(a, b);
    solveChrod(a, b);
    solveSecant(a, b);
    solveDifNewton(a, b);
    solveSteffanson(b);
    solveIterations(a, b);

}
