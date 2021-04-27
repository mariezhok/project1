#include <iostream>
#include <math.h>
#include <cmath>
#include <locale.h>
using namespace std;



double q = 2, w = 0.25;
double z = 6; //координационное число решетки
double l = 0.08e-3;//средняя длина капилляров
double k = 5e-13;//проницаемость
double rm=pow(8*l*l*k/3.14,0.25);
double a = 0.0, b =10.0, eps = 0.001;
double m = 0.3;

double f(double r) //функция распределения
{
	return (1/ (sqrt(2 * 3.14) * r * w)*exp(-(pow((log(r) - q),2)) / (2 * pow(w,2))) );
	
}


double Fr(double r)
{
	return (pow(rm, 4) - pow(r, 4) / (pow(r, 4) + (z / 2 - 1) * pow(rm, 4)) * f(r));
}

double dFdq(double r)
{
    return (pow(rm, 4) - pow(r, 4) / (pow(r, 4) + (z / 2 - 1) * pow(rm, 4)) * (log(r) - q) / pow(w, 2) * f(r));
}

double dFdw(double r)
{
    return (pow(rm, 4) - pow(r, 4) / (pow(r, 4) + (z / 2 - 1) * pow(rm, 4)) * f(r)*((-(pow(log(r) - q, 2) / (4 * w * w))) - 1) / w);
}


double Gr(double r)
{
	return (pow(r,2) * f(r));
}

double dGdq(double r)
{
    return (pow(r, 2) * f(r)* (log(r) - q) / pow(w, 2));
}

double dGdw(double r)
{
    return (pow(r, 2) * f(r) * ((-(pow(log(r)-q,2)/(4*w*w)))-1) / w);
}

double integral(double a, double b, int n, double (*foo)(double)) //функция,вычисляющая интеграл - прав.прямоуг-ки
{
    double r, h;
    double sum = 0.0;
       h = (b - a) / n;
    for (int i = 1; i <= n; i++)
    {   
        r = a + i * h;
        sum += foo(r);
    }
    return (sum * h);
}

double A11() /**********************производная F по мю ***************/
{
    double a11, s22;
   int n = 1;
    a11 = integral(a, b, n, dFdq);

    do
    {
        s22 = a11;
        n = 2 * n;
        a11 = integral(a, b, n, dFdq);

    } while (fabs(a11 - s22) > eps);
    return a11;
}

double A12() {
    /********************************** производна F по сигма ********************************************************************/

    double a12, y;
    int n = 1;
    a12 = integral(a, b, n, dFdw);
    do
    {
        y = a12;
        n = 2 * n;
        a12 = integral(a, b, n, dFdw);

    } while (fabs(a12 - y) > eps);
    return a12;
}

double A21()
{
    //*********************************************** dG/ dq ********************************************************************/

    double a21, z;
    int n = 1;
    a21 = integral(a, b, n, dGdq);
    do
    {
        z = a21;
        n = 2 * n;
        a21 = integral(a, b, n, dGdq);

    } while (fabs(a21 - z) > eps);
    return a21;
}

double A22() { //*********************************************** dG/ dw ********************************************************************/
    double a22, p;
   int n = 1;
    a22 = integral(a, b, n, dGdw);
    do
    {
        p = a22;
        n = 2 * n;
        a22 = integral(a, b, n, dGdw);

    } while (fabs(a22 - p) > eps);
    return a22;
}

double B1() 
{
    /*************************************integral F0 *****************************************************/

    double b1, r;
    int n = 1;
    b1 = integral(a, b, n, Fr);
    do
    {
        r = b1;
        n = 2 * n;
        b1 = integral(a, b, n, Fr);

    } while (fabs(b1 - r) > eps);
    return b1;
}

double B2()
{
    /*************************************integral G0 *****************************************************/

    double b2, c;
    int n = 1;
    b2 = integral(a, b, n, Gr);
    do
    {
        c = b2;
        n = 2 * n;
        b2 = integral(a, b, n, Gr);

    } while (fabs(b2 - c) > eps);
    return b2;
}


/**********************************************        MAIN     *******************************************************/

int main()
{
          
    
    double q1, w1;
    unsigned short int i;
    double E = 0.0001;
    double delta, deltaq, deltaw;
    
    double m1 = m * 2 * l * l / (z * 3.14);

    q1 = q;
    w1 = w;
    cout << "qi wi dqi dwi" << endl;
    do
    {     
        q = q1; 
        w = w1;
        delta = A11() * A22() - A21() * A12();
        deltaq = -B1() * A22() - (-B2()+m1) * A12();
        deltaw = A11() * (-B2()+m1) + A21() * B1();
        q1 = deltaq / delta;
        w1 = deltaw / delta;
        printf("%10.6f %10.6f %10.6f %10.6f\n", q1, w1, fabs(q1 - q), fabs(w1 - w));
    } while (fabs(q1 - q)> E && (fabs(w1 - w)>E));
    cout << "q =" << q1 <<endl<< "w = " << w1<<endl;



    return 0;
}
