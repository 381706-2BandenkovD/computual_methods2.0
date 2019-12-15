#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <ctime>

using namespace std;
double nuton(double x, double eps)
{
  double f, df; int iter = 0;
  do {
    f = sin(M_PI * x / 180) - 1 / x;
    df = M_PI / 180 * cos(M_PI * x / 180) + 1 / (x * x);
    x = x - f / df;
    iter++;
  } while (fabs(f) > eps&& iter < 20000);
  cout << "Решение по методу Ньютона:" << std::endl;
  cout << "1" << endl; cout << "2" << endl;
  cout << "Итераций по методу Ньютона: " << iter << endl;
  return x;
}

double hordi(double x0, double x1, double eps)
{
  double rez = x1, f0, f;
  int iter = 0;
  do {
    f = sin(M_PI * rez / 180) - 1 / rez;
    f0 = sin(M_PI * x0 / 180) - 1 / x0;
    rez = rez - f / (f - f0) * (rez - x0);
    iter++;
  } while (fabs(f) > eps&& iter < 20000);
  cout << "Решение по методу хорд:" << std::endl;
  cout << "1" << endl; cout << "2" << endl;
  cout << "Итераций по методу хорд: " << iter << endl;
  return rez;
}

typedef std::vector<std::vector<double>> Matrix;
template <class T>
bool equal(const T& a, const T& b) {
  if (abs(a - b) < 0.0000001) return true;
  else return false;
}
//------------------GAUS---------------------
void Gauss(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double> x) {
  int i, j, k;
  int size = b.size();
  double alfa;
  for (j = 0; j < size; j++)
  {
    for (i = j + 1; i < size; i++)
    {
      alfa = A[i][j] / A[j][j];
      for (k = j; k < size; k++)
      {
        A[i][k] -= alfa * A[j][k];
      }
      b[i] -= alfa * b[j];
    }

  }
  x[size - 1] = b[size - 1] / A[size - 1][size - 1];
  for (i = size - 2; i >= 0; i--)
  {
    double sum = 0;
    for (j = i + 1; j < size; j++)
    {
      sum += A[i][j] * x[j];
    }
    x[i] = (b[i] - sum) / A[i][i];
  }
  std::cout << "Решение по методу Гаусса:" << std::endl;
  for (i = 0; i < size; i++)
  {
    std::cout << x[i] << std::endl;
  }
}

//------------------------------------------------
std::vector<double> subtr(const std::vector<double>& v1, const std::vector<double>& v2) {
  std::vector<double> result(v1.size());
  if (v1.size() != v2.size()) {
    std::cout << "subtraction error!" << std::endl;
  }
  else {
    for (int i = 0; i < v1.size(); i++) {
      result[i] = v1[i] - v2[i];
    }
  }
  return result;
}

void printMatrix(Matrix const& mat) {
  for (int row = 0; row < mat.size(); ++row) {
    for (int col = 0; col < mat[0].size(); ++col) {
      cout.setf(ios::left);
      cout.width(8);
      std::cout << mat[row][col] << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

void printVector(const std::vector<double>& v) {
  for (double d : v) {
    std::cout << d << "\n";
  }
}

double det(Matrix a) {
  double determinant = 0;
  if (a.size() == 2) {
    determinant = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    return  determinant;
  }
  else {
    for (int i = 0; i < a.size(); i++) {
      //prepare tmp matrix
      Matrix tmp = a;
      tmp.erase(tmp.begin());
      for (int j = 0; j < tmp.size(); j++) {
        tmp[j].erase(tmp[j].begin() + i);
      }

      //calc det
      determinant += pow(-1, i + 2) * a[0][i] * det(tmp);
    }
  }
}

std::vector<double> mul_Z(Matrix A, std::vector<double> v) {
  std::vector<double> res(v.size());
  for (int i = 0; i < v.size(); i++) { //matrix row
    double sum = 0;
    for (int j = 0; j < v.size(); j++) { //matrix column and vector column
      sum += A[i][j] * v[j];
    }
    res[i] = sum;
  }
  return res;
}
std::vector<double> mul_R(Matrix A, std::vector<double> v, const double& coeff) {
  std::vector<double> res(v.size());
  for (int i = 0; i < v.size(); i++) { //matrix row
    double sum = 0;
    for (int j = 0; j < v.size(); j++) { //matrix column and vector column
      sum += A[i][j] * v[j];
    }
    res[i] = coeff * sum;
  }
  return res;
}

std::vector<double> sum(std::vector<double> v1, std::vector<double> v2) {
  int size = v1.size();
  std::vector<double> res(size);
  if (v1.size() != v2.size()) {
    std::cout << "sum error! " << "v1 is " << v1.size() << " v2 is " << v2.size() << std::endl;
    return res;
  }
  for (int i = 0; i < size; i++) {
    res[i] = v1[i] + v2[i];
  }
  return res;
}

double secondVectorNorm(const std::vector<double>& v) {
  double tmpSum = 0;
  for (double val : v) {
    tmpSum += pow(val, 2);
  }
  return sqrt(tmpSum);
}

//------------------EasyIterations--------------
void easyIterations(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double> x, const double eps) {
  int size = b.size();
  std::vector<double> x0(size);
  int iterCounter = 0;
  double currentEps = secondVectorNorm(subtr(b, mul_Z(A, x)));

  while (currentEps > eps) {

    for (int i = 0; i < size; i++) {
      x0 = x;
      double tmpX = 0;
      for (int j = 0; j < size; j++) {
        if (j == i) continue;
        tmpX -= A[i][j] * x0[j];
      }
      x[i] = (tmpX + b[i]) / A[i][i];
    }
    currentEps = secondVectorNorm(subtr(b, mul_Z(A, x)));
    iterCounter++;
  }
  std::cout << "Решение по методу простых итераций:" << std::endl;
  for (int i = 0; i < size; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << "Итераций по методу простых итераций: " << iterCounter << std::endl;
}
bool chekingDiag(const std::vector<std::vector<double>>& A) {
  bool returnval = true;
  bool flag = false;
  int size = A.size();
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) {
      if (A[i][i] < A[j][i]) {
        returnval = false;
        flag = true;
        break;
      }
      if (flag) break;

    }
  return returnval;
}


//----------------------------Zeidel-------------------------
bool isConvergenceZeidel(const Matrix& A) {
  double sum = 0;
  for (std::vector<double> v : A) {
    for (double val : v) {
      sum += pow(val, 2);
    }
  }
  double eNorm = sqrt(sum);
  return (eNorm < 1);
}

void Zeidel(Matrix A, std::vector<double> b, std::vector<double> x, double eps) {
  Matrix alpha;
  std::vector<double> beta;
  int size = A.size();

  for (int i = 0; i < b.size(); i++) {
    beta.push_back(b[i] / A[i][i]);
  }

  for (const std::vector<double>& v : A) {
    alpha.push_back(v);
  }
  for (int i = 0; i < size; i++) { //select row
    for (int j = 0; j < size; j++) { //select column
      if (i != j) alpha[i][j] = -(A[i][j] / A[i][i]);
      else alpha[i][j] = 0;
    }
  }

  if (isConvergenceZeidel(alpha)) {
    int iterCounter = 0;
    double currentEps = secondVectorNorm(subtr(b, mul_Z(A, x)));
    while (currentEps > eps) {
      x = sum(beta, mul_Z(alpha, x));
      iterCounter++;
      currentEps = secondVectorNorm(subtr(b, mul_Z(A, x)));
    }
    //std::cout << x.size();
    std::cout << "Решение по методу Зейделя:" << std::endl;
    for (int i = 0; i < size; i++)
    {
      std::cout << x[i] << std::endl;
    }
    std::cout << "Итераций по методу Зейделя: " << iterCounter << std::endl;
  }
  else std::cout << "Матрица на сходится для метода Зейделя\n";

}
//-------------------Relaxing--------------------------------
const std::vector<double> operator*(const double& left, const std::vector<double>& right) {
  std::vector<double> tmp(right.size());
  tmp = right;
  for (int i = 0; i < tmp.size(); i++)
    tmp[i] = tmp[i] * left;
  return tmp;

}
void relaxing(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double> x, const double eps) {
  double coeff;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(1, 2);
  coeff = dist(gen);

  Matrix alpha;
  std::vector<double> beta;
  int size = A.size();

  for (int i = 0; i < b.size(); i++) {
    beta.push_back(coeff * (b[i] / A[i][i]));
  }

  for (const std::vector<double>& v : A) {
    alpha.push_back(v);
  }
  for (int i = 0; i < size; i++) { //select row
    for (int j = 0; j < size; j++) { //select column
      if (i != j) alpha[i][j] = -(A[i][j] / A[i][i]);
      else alpha[i][j] = 0;
    }
  }

  if (isConvergenceZeidel(alpha)) {
    int iterCounter = 0;
    double currentEps = secondVectorNorm(subtr(b, mul_Z(A, x)));
    while (currentEps > eps) {
      x = sum((1 - coeff) * x, sum(beta, mul_R(alpha, x, coeff)));
      iterCounter++;
      currentEps = secondVectorNorm(subtr(b, mul_Z(A, x)));
    }
    std::cout << "Решение по методу верхней релаксации:" << std::endl;
    for (int i = 0; i < size; i++)
    {
      std::cout << x[i] << std::endl;
    }
    std::cout << "Итераций по методу верхней релаксации: " << iterCounter << std::endl;
    std::cout << "coeff=" << coeff << std::endl;
  }
  else std::cout << "Матрица не сходится для метода верхней релаксации\n";

}

//----------------------------Kramer----------------------------
void Kramer(Matrix A, std::vector<double> b, std::vector<double> x) {
  int size = A.size();
  x.clear();
  std::vector<double> dets;
  double Adet = det(A);
  if (Adet != 0) {
    for (int i = 0; i < size; i++) {
      Matrix tmp;
      for (const std::vector<double>& v : A) {
        tmp.push_back(v);
      }
      for (int j = 0; j < tmp.size(); j++) {
        tmp[j][i] = b[j];
      }
      dets.push_back(det(tmp));
    }
    for (double d : dets) {
      x.push_back(d / Adet);
    }
  }
  else {
    cout << "Определитель 0";
    
  }

  std::cout << "Решение по методу Крамера:" << std::endl;
  for (int i = 0; i < size; i++)
  {
    std::cout << x[i] << std::endl;
  }
}

void printMenu() {
  cout << "_________РЕШЕНИЕ СИСТЕМ ЛИНЕЙНЫХ УРАВНЕНИЙ________" << std::endl;
  std::cout << "1. Ввод системы уравнений\n"
    << "2. Генерация системы уравнений\n"
    << "3. Выход\n";
  cout << "___________________________________________________" << std::endl;
}

Matrix inputMatrix(int elementsCount) {
  std::vector<std::vector<double>> matrix(elementsCount, std::vector<double>(elementsCount));
  for (std::vector<double>& row : matrix) {
    for (double& column : row) {
      std::cin >> column;
    }
  }
  return matrix;
}

std::vector<double> inputVector(int elementsCount) {
  std::vector<double> v(elementsCount);
  for (int i = 0; i < elementsCount; i++) {
    std::cin >> v[i];
  }
  return v;
}

Matrix generateMatrix(int elementsCount) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(1, 2);

  std::vector<std::vector<double>> matrix(elementsCount, std::vector<double>(elementsCount));
  for (std::vector<double>& row : matrix) {
    for (double& column : row) {
      column = 100 - dist(gen) * 41;
    }
  }
  return matrix;
}

std::vector<double> generateVector(int elementsCount) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(1, 2);

  std::vector<double> v(elementsCount);
  for (double& val : v) {
    val = 100 - dist(gen) * 50;
  }
  return v;
}

void useMenu() {
  int number = 0;
  bool hasMatrix = false;
  int elementsCount = 0;

  Matrix A;
  std::vector<double> b;
  std::vector<double> x;

  printMenu();
  while (number != 3) {
    std::cin >> number;
    switch (number) {
    case 1: {
      std::cout << "Ведите порядок системы\n";
      std::cin >> elementsCount;
      if (elementsCount < 2) std::cout << "ERROR";
      else {
        std::cout << "Ведите матрицу\n";
        A = inputMatrix(elementsCount);
        std::cout << "Ведите вектор свободных членов\n";
        b = inputVector(elementsCount);
        std::cout << "Ведите начальное приближение\n";
        x = inputVector(elementsCount);
      }
      hasMatrix = true;
      int flad = 0;
      while (flad != 6) {
        cout << "_________Выберите способ решения________" << std::endl;
        cout << "1.Решение методом Гауса\n";
        cout << "2.Решение методом Крамера\n";
        cout << "3.Решение методом простых итераций\n";
        cout << "4.Решение методом Зейделя\n";
        cout << "5.Решение методом релаксации\n";
        cout << "6.Выход\n";
        cout << "__________________________________________" << std::endl;
        cin >> flad;
        if (flad == 1) {
          unsigned int start_time = clock();
          Gauss(A, b, x);
          unsigned int end_time = clock();
          unsigned int search_time = end_time - start_time;
          cout << "время выполнения " << search_time << "мс" << endl;
        }
        if (flad == 2) {
          unsigned int start_time = clock();
          Kramer(A, b, x);
          unsigned int end_time = clock();
          unsigned int search_time = end_time - start_time;
          cout << "время выполнения " << search_time << "мс" << endl;
        }
        if (flad == 3) {
          unsigned int start_time = clock();
          easyIterations(A, b, x, 0.00001);
          unsigned int end_time = clock();
          unsigned int search_time = end_time - start_time;
          cout << "время выполнения " << search_time << "мс" << endl;
        }
        if (flad == 4) {
          unsigned int start_time = clock();
          Zeidel(A, b, x, 0.00001);
          unsigned int end_time = clock();
          unsigned int search_time = end_time - start_time;
          cout << "время выполнения " << search_time << "мс" << endl;
        }
        if (flad == 5) {
          unsigned int start_time = clock();
          relaxing(A, b, x, 0.00001);
          unsigned int end_time = clock();
          unsigned int search_time = end_time - start_time;
          cout << "время выполнения " << search_time << "мс" << endl;
        }
        if (flad == 6)
          break;
      }
      printMenu();
      break;
    }
    case 2: {
      std::cout << "Ведите порядок системы\n";
      std::cin >> elementsCount;
      if (elementsCount < 2) std::cout << "ERROR";
      else {
        A = generateMatrix(elementsCount);
        printMatrix(A);
        b = generateVector(elementsCount);
        x = generateVector(elementsCount);
        int flad = 0;
        while (flad != 6) {
          cout << "_________Выберите способ решения________" << std::endl;
          cout << "1.Решение методом Гауса\n";
          cout << "2.Решение методом Крамера\n";
          cout << "3.Решение методом простых итераций\n";
          cout << "4.Решение методом Зейделя\n";
          cout << "5.Решение методом релаксации\n";
          cout << "6.Выход\n";
          cout << "__________________________________________" << std::endl;
          cin >> flad;
          if (flad == 1) {
            unsigned int start_time = clock();
            Gauss(A, b, x);
            unsigned int end_time = clock();
            unsigned int search_time = end_time - start_time;
            cout << "время выполнения " << search_time << "мс" << endl;
          }
          if (flad == 2) {
            unsigned int start_time = clock();
            Kramer(A, b, x);
            unsigned int end_time = clock();
            unsigned int search_time = end_time - start_time;
            cout << "время выполнения " << search_time << "мс" << endl;
          }
          if (flad == 3) {
            unsigned int start_time = clock();
            easyIterations(A, b, x, 0.00001);
            unsigned int end_time = clock();
            unsigned int search_time = end_time - start_time;
            cout << "время выполнения " << search_time << "мс" << endl;
          }
          if (flad == 4) {
            unsigned int start_time = clock();
            Zeidel(A, b, x, 0.00001);
            unsigned int end_time = clock();
            unsigned int search_time = end_time - start_time;
            cout << "время выполнения " << search_time << "мс" << endl;
          }
          if (flad == 5) {
            unsigned int start_time = clock();
            relaxing(A, b, x, 0.00001);
            unsigned int end_time = clock();
            unsigned int search_time = end_time - start_time;
            cout << "время выполнения " << search_time << "мс" << endl;
          }
          if (flad == 6)
            break;
        }
        printMenu();
      }
      hasMatrix = true;
      break;
    }
    }
  }
}


int main() {
  setlocale(LC_ALL, "Russian");
 // nuton(1, 0.00001);
  //hordi(1.0, 10.0, 0.000001);
  useMenu();

}