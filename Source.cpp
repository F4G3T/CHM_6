#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <iomanip>

using namespace std;
using namespace std::chrono;

using Complex = complex<double>;
const double PI = 3.14159265358979323846;

//  DFT/IDFT 
vector<Complex> DFT(const vector<Complex>& x) {
    int N = x.size();
    vector<Complex> X(N, 0);
    for (int k = 0; k < N; ++k) {
        for (int n = 0; n < N; ++n) {
            double angle = -2 * PI * k * n / N;
            X[k] += x[n] * Complex(cos(angle), sin(angle));
        }
    }
    return X;
}

vector<Complex> IDFT(const vector<Complex>& X) {
    int N = X.size();
    vector<Complex> x(N, 0);
    for (int n = 0; n < N; ++n) {
        for (int k = 0; k < N; ++k) {
            double angle = 2 * PI * k * n / N;
            x[n] += X[k] * Complex(cos(angle), sin(angle));
        }
        x[n] /= N;
    }
    return x;
}

//   FFT/IFFT 
void bitReverseReorder(vector<Complex>& a) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            swap(a[i], a[j]);
    }
}

void FFT_iterative(vector<Complex>& a, bool invert = false) {
    int n = a.size();
    bitReverseReorder(a);

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        Complex wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            Complex w(1);
            for (int j = 0; j < len / 2; ++j) {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (Complex& x : a)
            x /= n;
    }
}

void IFFT_iterative(vector<Complex>& a) {
    FFT_iterative(a, true);
}

int main() {
    setlocale(LC_ALL, "");
    
    // Параметры 
    int n = 10;          
    int N = 1 << n;      
    double A = 2.56;
    double B = 0.13;
    double omega1 = 3;
    double omega2 = 195;
    double phi = PI / 3;  

        // Генерация сигнала
    vector<Complex> z(N);
    for (int j = 0; j < N; ++j) {
        double val = A * cos(2 * PI * omega1 * j / N + phi) +
            B * cos(2 * PI * omega2 * j / N);
        z[j] = Complex(val, 0);
    }

    // DFT и замер времени
    auto start = high_resolution_clock::now();
    vector<Complex> Z_dft = DFT(z);
    auto stop = high_resolution_clock::now();
    auto dft_time = duration_cast<microseconds>(stop - start);

    // FFT и замер времени
    vector<Complex> z_fft_copy = z;
    start = high_resolution_clock::now();
    FFT_iterative(z_fft_copy);
    stop = high_resolution_clock::now();
    auto fft_time = duration_cast<microseconds>(stop - start);
    vector<Complex> Z_fft = z_fft_copy;

    cout << fixed << setprecision(6);
    cout << "Время выполнения:\n";
    cout << "DFT: " << dft_time.count() / 1000000.0 << " сек ("
        << dft_time.count() << " мкс)\n";
    cout << "FFT: " << fft_time.count() / 1000000.0 << " сек ("
        << fft_time.count() << " мкс)\n";

    double acceleration = (double)dft_time.count() / fft_time.count();
    cout << "Ускорение: " << fixed << setprecision(2)
        << acceleration << " раз\n\n";

    // Таблица DFT только частоты с ненулевыми амплитудами
    cout << "Таблица DFT (только частоты с ненулевыми амплитудами):\n";
    cout << "   m      Re z        Re z^       Im z^     Амплитуда      Фаза    \n"; 

    vector<int> significant_freqs;
    for (int m = 0; m < N; ++m) {
        double amp = abs(Z_fft[m]);
        if (amp > 10.0) {  
            significant_freqs.push_back(m);
        }
    }

    sort(significant_freqs.begin(), significant_freqs.end());

    cout << fixed << setprecision(6);
    for (int m : significant_freqs) {
        double amp = abs(Z_fft[m]);
        double phase = arg(Z_fft[m]);
        cout << " " << setw(4) << m << "  "
            << setw(10) << z[m].real() << "  "
            << setw(10) << Z_fft[m].real() << "  "
            << setw(10) << Z_fft[m].imag() << "  "
            << setw(10) << amp << "  "
            << setw(11) << phase << " \n";
    }
 
    // Фильтрация гармонического сигнала
    cout << "Фильтрация: обнуляем компоненты при m = " << omega2
        << " (высокочастотный шум) и m = " << (N - omega2) << "\n\n";

    vector<Complex> Z_filtered = Z_fft;

    // Обнуляем высокочастотный шум и его эрмитово-сопряженную частоту
    Z_filtered[omega2] = 0;
    Z_filtered[N - omega2] = 0;

    // Также обнуляем соседние частоты для полного удаления шума
    for (int offset = 1; offset <= 2; offset++) {
        int m1 = omega2 + offset;
        int m2 = omega2 - offset;
        int sym_m1 = N - m1;
        int sym_m2 = N - m2;

        if (m1 < N) Z_filtered[m1] = 0;
        if (m2 >= 0) Z_filtered[m2] = 0;
        if (sym_m1 < N && sym_m1 >= 0) Z_filtered[sym_m1] = 0;
        if (sym_m2 < N && sym_m2 >= 0) Z_filtered[sym_m2] = 0;
    }

    // Таблица после фильтрации
    cout << "Таблица после фильтрации (включая обнуленные компоненты):\n";
   
    cout << "   m      Re z^       Im z^   \n";
    

    for (int m : significant_freqs) {
        cout << " " << setw(4) << m << "  "
            << setw(10) << Z_filtered[m].real() << "  "
            << setw(10) << Z_filtered[m].imag() << " \n";
    }
    

    // Восстановление сигнала после фильтрации с использованием IFFT
    vector<Complex> z_filtered = Z_filtered;
    IFFT_iterative(z_filtered);

    // Проверка, что сигнал остался вещественным
    for (int i = 0; i < N; ++i) {
        z_filtered[i] = Complex(z_filtered[i].real(), 0);
    }

    // Сохранение данных для гармонического сигнала
    ofstream out1("harmonic_signal.txt");
    out1 << "index original filtered\n";
    for (int i = 0; i < N; ++i) {
        out1 << i << " " << z[i].real() << " " << z_filtered[i].real() << "\n";
    }
    out1.close();

     // Генерация к-п сигнала
    vector<Complex> z_pw(N);
    for (int j = 0; j < N; ++j) {
        double val = 0;
        if (j >= N / 4 && j <= N / 2) {
            val = A + B * cos(2 * PI * omega2 * j / N);
        }
        else if (j > 3 * N / 4 && j < N) {
            val = A + B * cos(2 * PI * omega2 * j / N);
        }
        z_pw[j] = Complex(val, 0);
    }

    // FFT кусочно-постоянного сигнала
    vector<Complex> Z_pw = z_pw;
    FFT_iterative(Z_pw);

    // Фильтрация к-п сигнала
    vector<Complex> Z_pw_filtered = Z_pw;
    // Низкочастотная фильтрация
    int cutoff = 20;
    for (int m = 0; m < N; ++m) {
        if (m > cutoff && m < N - cutoff) {
            Z_pw_filtered[m] = 0;
        }
    }

    // Восстанавливаем эрмитову симметрию
    for (int k = 1; k <= cutoff; ++k) {
        Z_pw_filtered[N - k] = conj(Z_pw_filtered[k]);
    }

    // Обратное преобразование
    vector<Complex> z_pw_filtered = Z_pw_filtered;
    IFFT_iterative(z_pw_filtered);

    // Убираем мнимую часть
    for (int i = 0; i < N; ++i) {
        z_pw_filtered[i] = Complex(z_pw_filtered[i].real(), 0);
    }

    // Сохранение данных кусочно-постоянного сигнала
    ofstream out2("piecewise_signal.txt");
    out2 << "index original filtered\n";
    for (int i = 0; i < N; ++i) {
        out2 << i << " " << z_pw[i].real() << " " << z_pw_filtered[i].real() << "\n";
    }
    out2.close();
       return 0;
}