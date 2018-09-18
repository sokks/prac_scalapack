# Тривиальная эволюция с логгированием состояний

**Вход:**  
Размерность системы N  
R - изначальная матрица плотности, в виде либо матрицы N*N (читать из файла), либо вектора размера N (читать из файла), либо числа (номера базисного состояния системы).  
dT - шаг по времени  
H - гамильтониан системы (читать из бинарного файла)  
количество шагов — n  
  
**Выход:** 
последовательность диагоналей матрицы плотности, с шагом по времени dT в виде:  
Abs(Diag(R(0))  
Abs(Diag(R(dT)))  
Abs(Diag(R(2*dT)))  
…  
Abs(Diag(R(n*dT)))  
  
Код оформить в виде библиотеки на C/C++, матричные операции — через Scalapack.  
В примере можно либо загружать матрицы из файла, либо генерировать их динамически.  
Тип чисел: complex<double>.  
  
Программа должна по максимуму использовать выделенные при запуске распределённые процессы.  

## Алгоритм
1. Ввод матриц и тд, параметры как аргументы командной строки
2. Вычистить <a href="https://www.codecogs.com/eqnedit.php?latex=U_{\Delta&space;t}&space;=&space;e^{-1\frac{1}{i}H\Delta&space;t}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_{\Delta&space;t}&space;=&space;e^{-1\frac{1}{i}H\Delta&space;t}" title="U_{\Delta t} = e^{-1\frac{1}{i}H\Delta t}" /></a> 
   Для взятия экспоненты от матрицы нужно разложить ее по собственным векторам <a href="https://www.codecogs.com/eqnedit.php?latex=H&space;=&space;V^*DV" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H&space;=&space;V^*DV" title="H = V^*DV" /></a>
   Можно без скалапака (даже удобнее)
3. 
for int i = 1; i < n; i++ {  
    <a href="https://www.codecogs.com/eqnedit.php?latex=\rho_i&space;=&space;U^*_{\Delta&space;t}&space;\rho_{i-1}&space;U_{\Delta&space;t}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho_i&space;=&space;U^*_{\Delta&space;t}&space;\rho_{i-1}&space;U_{\Delta&space;t}" title="\rho_i = U^*_{\Delta t} \rho_{i-1} U_{\Delta t}" /></a> // со скалапаком  
}

## Замечания
Начать с того чтобы разложить матрицу по собственным векторам и восстановить обратно.  
В скалапаке перевернуты массивы, порядок векторов важен. Самая большая проблема в том, чтобы преобразовать матрицу во что-то, с чем будет работать фортран.  
Header-файлы нужно написать по примеру либо найти в интернете.