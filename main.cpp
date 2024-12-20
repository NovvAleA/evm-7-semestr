#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
const double eps       = 0.00001;  // epsilon
const double alpha     = 3.0;
const double a         = 0.9;     // левая граница отрезка
const double b         = 1.45;     // правая граница отрезка (тоже подобрать как и а)
const double q         = 0.02;     // константа в сжимающем отображении
const double residual_ = 0.00005;
const double x0        = 1.5;
double M = 6.8;                   //ограничение модуля производной сверху на отрезке
double m = 2.4;                   //ограничение модуля производной снизу на отрезке


double phi(double x){
    return cbrt(alpha*sin(sqrt(1 + 2*sin(x))));
}

double f(double x){
    return alpha*sin(sqrt(1 + 2*sin(x))) - pow(x,3);
}

double first_derivative_f(double x){
    return alpha*cos(sqrt(1 + 2*sin(x)))*cos(x)/sqrt(1 + 2*sin(x)) - 3*pow(x,2);
}

double simple_iteration_method(ostream& res_out){
    //ofstream res_out;
    //res_out.open("result.txt", ios::app);
    int N = 0;
    double x = x0, x1 = x0;
    double delta_x = b - a; //delta_x - разность между n и n+1 результатами итераций
    cout << endl << "simple iteration method" << endl;
    while (delta_x >= eps * (1 - q)/q) {
        x1 = phi(x);
        delta_x = fabs(x1 - x);
        x = x1;
        N++;
        cout << "x = " << x1 << endl;
    }
    cout << "N = " << N << endl;
    res_out << "|" << "Простых итераций" << "|" << setw(16) << x1 << "|" << setw(12) << f(x1) << "|" << " [" << setw(4) << a << ", " << setw(4) << b << "] " << "|";
    res_out << setw(5) << x0 << "|" << setw(7) << N << "|" << setw(5) << M << "|" << setw(5) << m<< "|" << setw(5) << q << "|" << endl;
    res_out << "+----------------+----------------+------------+--------------+-----+-------+-----+-----+-----+" << endl;
    return x1;
}

double chord_method(ostream& res_out){
    //ofstream res_out;
    //res_out.open("result.txt", ios::app);
    int N = 0;
    double x = x0, x1 = x0;
    double delta_x = b - a; //delta_x - разность между n и n+1 результатами итераций
    cout << endl << "Chord method" << endl;
    while (delta_x >= eps * m / (M - m)) {
        x1 = x - f(x)*(x - a)/(f(x) - f(a));
        delta_x = fabs(x1 - x);
        x = x1;
        N++;
        cout << "x = " << x1 << endl;
    }
    cout << "N = " << N << endl;
    res_out << "|" << "Метод хорд      " << "|" << setw(16) << x1 << "|" << setw(12) << f(x1) << "|" << " [" << setw(4) << a << ", " << setw(4) << b << "] " << "|";
    res_out << setw(5) << x0 << "|" << setw(7) << N << "|" << setw(5) << M << "|" << setw(5) << m<< "|" << setw(5) << q << "|" << endl;
    res_out << "+----------------+----------------+------------+--------------+-----+-------+-----+-----+-----+" << endl;
    return x1;

}

double bisect_method(ostream& res_out){
    double temp_a = a, temp_b = b; //концы отрезков в методе бисекции
    int N = log2((b - a)/eps) + 1; //число итераций для достижения необходимой точности
    double ksi = (temp_a + temp_b) / 2; //середина отрезка
    cout << endl << "Bisect method" << endl;

    for (int i = 1; i <= N && fabs(f(ksi))>=residual_; i++){
        if (f(ksi) * f(temp_a) < 0){
            temp_b = ksi;
        }
        else temp_a = ksi;
        cout << "i = " << i << " [a,b]=[" << temp_a << ", " << temp_b << "]" << endl;
        ksi = (temp_a + temp_b) / 2;
    }

    cout << "x = " << ksi << endl;
    cout << "N = " << N << endl;
    res_out << "|" << "Метод бисекции  " << "|" << setw(16) << ksi << "|" << setw(12) << f(ksi) << "|" << " [" << setw(4) << a << ", " << setw(4) << b << "] " << "|";
    res_out << setw(5) << x0 << "|" << setw(7) << N << "|" << setw(5) << M << "|" << setw(5) << m<< "|" << setw(5) << q << "|" << endl;
    res_out << "+----------------+----------------+------------+--------------+-----+-------+-----+-----+-----+" << endl;
    return ksi;

}

double newton_method(ostream& res_out){
    int N = 0;
    double x = x0, x1 = x0;
    double delta_x = b - a; //delta_x - разность между n и n+1 результатами итераций
    cout << endl << "Newton method" << endl;
    while (delta_x >= eps * m / (M - m)) {
        x1 = x - f(x)/first_derivative_f(x);
        delta_x = fabs(x1 - x);
        x = x1;
        N++;
        cout << "x = " << x1 << endl;
    }
    cout << "N = " << N << endl;
    res_out << "|" << "Метод Ньютона   " << "|" << setw(16) << x1 << "|" << setw(12) << f(x1) << "|" << " [" << setw(4) << a << ", " << setw(4) << b << "] " << "|";
    res_out << setw(5) << x0 << "|" << setw(7) << N << "|" << setw(5) << M << "|" << setw(5) << m<< "|" << setw(5) << q << "|" << endl;
    res_out << "+----------------+----------------+------------+--------------+-----+-------+-----+-----+-----+" << endl;
    return x1;

}

double etken_method(ostream& res_out){
    int N = 0;
    double x = x0, x1 = phi(x0), x1_5 = phi(x1), x2 = 0;//x = x'n-1, x1 = x'n, x1_5 = x'n+1/2, x2 = x'n+1
    double delta_x = b - a; //delta_x - разность между n и n+1 результатами итераций
    cout << endl << "Etken method" << endl;
    while (delta_x >= eps) {
        if((x - 2*x1 + x1_5) == 0){
            cout << "x = " << x1 << endl;
            cout << "zero in the denominator, exiting while " << x1 << endl;
            break;
        }
        x2 = (x*x1_5 - x1*x1)/(x - 2*x1 + x1_5);
        delta_x = fabs(x2 - x1);
        x = x1;
        x1_5 = phi(x1);
        x1 = x2;
        N++;
        cout << "x = " << x2 << endl;
    }
    cout << "N = " << N << endl;
    res_out << "|" << "Метод Еткена    " << "|" << setw(16) << x2 << "|" << setw(12) << f(x2) << "|" << " [" << setw(4) << a << ", " << setw(4) << b << "] " << "|";
    res_out << setw(5) << x0 << "|" << setw(7) << N << "|" << setw(5) << M << "|" << setw(5) << m<< "|" << setw(5) << q << "|" << endl;
    res_out << "+----------------+----------------+------------+--------------+-----+-------+-----+-----+-----+" << endl;
    return x1;

}


int main() {
    ofstream res_out;
    res_out.open("result.txt");
    res_out << "eps   = " << eps << endl;
    res_out << "alpha = " << alpha << endl << endl;
    res_out << "+================+================+============+==============+=====+=======+=====+=====+=====+" << endl;
    res_out << "|  Метод         | Корень         |Невязка     | Отрезок      | x0  | N + 1 | M   | m   | q   |" << endl;
    res_out << "+================+================+============+==============+=====+=======+=====+=====+=====+" << endl;
    simple_iteration_method(res_out);
    newton_method(res_out);
    chord_method(res_out);
    bisect_method(res_out);
    etken_method(res_out);
    res_out.close();
    return 0;
}

