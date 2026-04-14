#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <iomanip>

using namespace std;

typedef function<double(double, double)> Function2D;

/**
 * Численное интегрирование по треугольной области методом Римана
 * 
 * @param a, b, c - вершины треугольника
 * @param n, m - степени x и y в подынтегральной функции
 * @param n_alpha - количество шагов по alpha (внутренний цикл)
 * @param m_beta - количество шагов по beta (внешний цикл)
 * @param q - дополнительная функция q(x, y)
 * @return значение интеграла
 */
double intRim(vector<double>& a, vector<double>& b, vector<double>& c,
              int n, int m, int n_alpha, int m_beta,
              Function2D q = [](double x, double y) { return 1.0; }) {
    
    double ax = a[0], ay = a[1];
    double bx = b[0], by = b[1];
    double cx = c[0], cy = c[1];
    
    double Ax = ax - cx;
    double Bx = bx - cx;
    double Ay = ay - cy;
    double By = by - cy;
    
    // Якобиан преобразования координат
    double Jacob = Ax * By - Ay * Bx;
    
    double dbeta = 1.0 / m_beta;
    double dalpha = 1.0 / n_alpha;  // ИСПРАВЛЕНО: должна быть константой
    
    double summ = 0.0;
    
    // Внешний цикл по beta
    for (int i = 0; i < m_beta; i++) {
        double beta = dbeta * (i + 0.5);
        double suma = 0.0;
        
        // Внутренний цикл по alpha
        for (int j = 0; j < n_alpha; j++) {
            double alpha = dalpha * (j + 0.5);
            
            // Вычисление координат точки в треугольнике
            double x = Ax * alpha + Bx * beta + cx;
            double y = Ay * alpha + By * beta + cy;
            
            // Подынтегральная функция: x^n * y^m * q(x, y)
            suma += pow(x, n) * pow(y, m) * q(x, y);
        }
        
        summ += suma * dalpha;
    }
    
    // Результат с учётом якобиана и шага по beta
    double result = summ * dbeta * Jacob;
    
    cout << "Якобиан: " << fixed << setprecision(6) << Jacob << "\n";
    
    return result;
}

//int main() {
//    // Параметры интегрирования
//    int n = 5;
//    int m = 3;
//    int an = 500 * 2;   // количество шагов по alpha
//    int bm = 500 * 2;   // количество шагов по beta
//    
//    // Вершины треугольника
//    vector<double> A = {-3.0, 3.0};
//    vector<double> B = {1.0, 1.0};
//    vector<double> C = {-2.0, 3.0};
//    
//    // Подынтегральная функция: x^n * y^m * x*sin(y)
//    Function2D q = [](double x, double y) {
//        return x * sin(y);
//    };
//    
//    // Вычисление интеграла
//    double result = intRim(A, B, C, n, m, an, bm, q);
//    
//    cout << "Результат интеграла: " << fixed << setprecision(10) << result << "\n";
//    
//    return 0;
//}
