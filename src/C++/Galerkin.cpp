#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

struct Node { double x, y; };
struct Triangle { int v[3]; };

// Правая часть
double q(double x, double y) {
    double x2 = x * x, y2 = y * y;
    double x1 = x - 1.0, y1 = y - 1.0;
    return 24.0 * x2 * x1 * x1 + 24.0 * y2 * y1 * y1 + 2.0 * (12 * x2 - 12 * x + 2) * (12 * y2 - 12 * y + 2);
}

// Локальные матрицы: stiffness ∫∇φ·∇φ и mass ∫φ·φ
void localStiffnessMass(const Node& A, const Node& B, const Node& C,
    double K[3][3], double M[3][3]) {
    double x1 = A.x, y1 = A.y, x2 = B.x, y2 = B.y, x3 = C.x, y3 = C.y;
    double area = 0.5 * fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

    double b[3] = { y2 - y3, y3 - y1, y1 - y2 };
    double c[3] = { x3 - x2, x1 - x3, x2 - x1 };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            K[i][j] = (b[i] * b[j] + c[i] * c[j]) / (4.0 * area);
            M[i][j] = (i == j) ? area / 6.0 : area / 12.0;
        }
    }
}

// Сборка глобальной матрицы
void assemble(const vector<Node>& nodes, const vector<Triangle>& tris,
    vector<vector<double>>& K, vector<vector<double>>& M) {
    int N = nodes.size();
    K.assign(N, vector<double>(N, 0.0));
    M.assign(N, vector<double>(N, 0.0));

    for (auto& t : tris) {
        double Ke[3][3], Me[3][3];
        localStiffnessMass(nodes[t.v[0]], nodes[t.v[1]], nodes[t.v[2]], Ke, Me);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                int I = t.v[i], J = t.v[j];
                K[I][J] += Ke[i][j];
                M[I][J] += Me[i][j];
            }
        }
    }
}

// Формируем локальный вектор нагрузки
void assembleLoad(const vector<Node>& nodes, const vector<Triangle>& tris,
    vector<double>& F) {
    int N = nodes.size();
    F.assign(N, 0.0);
    for (auto& t : tris) {
        double area = 0.5 * fabs(
            (nodes[t.v[1]].x - nodes[t.v[0]].x) * (nodes[t.v[2]].y - nodes[t.v[0]].y) -
            (nodes[t.v[2]].x - nodes[t.v[0]].x) * (nodes[t.v[1]].y - nodes[t.v[0]].y)
        );
        // точка в центре треугольника
        double xc = (nodes[t.v[0]].x + nodes[t.v[1]].x + nodes[t.v[2]].x) / 3.0;
        double yc = (nodes[t.v[0]].y + nodes[t.v[1]].y + nodes[t.v[2]].y) / 3.0;
        double val = q(xc, yc) * area / 3.0;
        for (int i = 0; i < 3; i++) F[t.v[i]] += val;
    }
}

// Простое решение СЛАУ методом Гаусса (для учебного примера)
vector<double> solveLinearSystem(vector<vector<double>> A, vector<double> b) {
    int n = b.size();
    for (int i = 0; i < n; i++) {
        // поиск главного элемента
        int pivot = i;
        for (int j = i + 1; j < n; j++)
            if (fabs(A[j][i]) > fabs(A[pivot][i])) pivot = j;
        swap(A[i], A[pivot]);
        swap(b[i], b[pivot]);
        // исключение
        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) A[j][k] -= factor * A[i][k];
            b[j] -= factor * b[i];
        }
    }
    // обратный ход
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }
    return x;
}

// Пример сетки: единичный квадрат 2x2 линейная сетка
void createMesh(vector<Node>& nodes, vector<Triangle>& tris) {
    int nx = 5, ny = 5;
    double dx = 1.0 / (nx - 1), dy = 1.0 / (ny - 1);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            nodes.push_back({ i * dx,j * dy });
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            int n0 = j * nx + i, n1 = n0 + 1, n2 = n0 + nx, n3 = n2 + 1;
            tris.push_back({ n0,n1,n3 });
            tris.push_back({ n0,n3,n2 });
        }
    }
}

// Граничные условия u=0 на границе
void applyDirichlet(vector<vector<double>>& K, vector<double>& b, const vector<Node>& nodes) {
    int N = nodes.size();
    for (int i = 0; i < N; i++) {
        if (nodes[i].x == 0 || nodes[i].x == 1 || nodes[i].y == 0 || nodes[i].y == 1) {
            for (int j = 0; j < N; j++) K[i][j] = 0.0;
            K[i][i] = 1.0;
            b[i] = 0.0;
        }
    }
}

int main() {
    vector<Node> nodes;
    vector<Triangle> tris;
    createMesh(nodes, tris);

    int N = nodes.size();

    vector<vector<double>> K, M;
    vector<double> F;

    // 1. Сборка матриц
    assemble(nodes, tris, K, M);
    assembleLoad(nodes, tris, F);

    // 2. Решаем w: K w = F
    vector<double> w = solveLinearSystem(K, F);

    // 3. Решаем u: K u = M*w
    vector<double> Mw(N, 0.0);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            Mw[i] += M[i][j] * w[j];

    vector<double> u = solveLinearSystem(K, Mw);

    // 4. Вывод решения
    cout << "u at nodes:\n";
    for (int i = 0; i < N; i++)
        cout << fixed << setprecision(6) << nodes[i].x << " " << nodes[i].y << " " << u[i] << "\n";

    return 0;
}