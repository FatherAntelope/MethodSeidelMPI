#include "mpi.h"
#include <stdio.h>
#include "FileTransfer.h"
#include <math.h>

/// Debug x86       | Debug x64     | Release x64
/// ----------------|---------------|---------------
/// 1: 0.12752s     | 0.126089s     | 1: 0.0379844s
/// 2: 0.130218s    | 0.133961s     | 2: 0.0587723s
/// 4: 0.142875s    | 0.143954s     | 4: 0.0876631s
/// 6: 27.8739s     | 25.5916s      | 6: 24.956s
/// 8: 67.33s       | 60.299s       | 8: 56.0289s

using namespace std;

int nSize, nRank;
double** A;
double* B;
double* X;
double* tempX;
double epsilon, norm;
int N, countIterations = 0, /*номер файла:*/ numberFile = 8;
FileTransfer* fileInit;


int main(int argc, char* argv[]) {
    stringstream strNumberFile;
    strNumberFile << numberFile;
    double time;
    MPI_Init(&argc, &argv);

    //Общее количество процессов
    MPI_Comm_size(MPI_COMM_WORLD, &nSize);
    //Номер текущего процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &nRank);


    //Инициализация данных с файла
    if (nRank == 0) {


        //Директория текстовых документов
        fileInit = new FileTransfer("D:\\cpp\\lr3\\lr3\\lr3\\");
        fileInit->initData(N, epsilon, A, B, "input2_" + strNumberFile.str() + ".txt");

        X = new double[N];
        tempX = new double[N];

        for (int i = 0; i < N; i++) {
            X[i] = tempX[i] = 0;
        }

        cout << "START!" << endl;
        time = MPI_Wtime();
    }


    /// Передаем по всем потокам epSilon и N
    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //Выделяем память для остальных потоков, кроме 0
    if (nRank > 0)
    {
        B = new double[N];
        X = new double[N];
        tempX = new double[N];

        A = new double* [N];
        for (int i = 0; i < N; i++) {
            A[i] = new double[N];
        }
    }
        
    //Передаем массивы по всем потокам
    MPI_Bcast(X, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tempX, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < N; i++) {
        MPI_Bcast(A[i], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

        
    ///Алгоритм
    int i, j;
    double pieceSum = 0, resSum = 0;

    int i1 = (N / nSize) * nRank;
    int i2 = (N / nSize) * (nRank + 1);

    if (nRank == nSize - 1) {
        i2 = N;
    }


    while (true) {
        
        for (i = 0; i < N; i++)
            tempX[i] = X[i];

        //Объединяем результаты двух сумм из формулы в одну
        for (i = 0; i < N; i++) {
            pieceSum = 0;

            for (j = 0; j < i; j++) {
                if ((i1 <= j) && (j < i2) && i != j) {
                    pieceSum += (-A[i][j] / A[i][i]) * X[j];
                }
            }

            for (j = i; j < N; j++) {
                if ((i1 <= j) && (j < i2) && i != j) {
                    pieceSum += (-A[i][j] / A[i][i]) * tempX[j];
                }
            }

            //MPI_Barrier(MPI_COMM_WORLD);
            resSum = 0;
            //Объединяем pieceSum из всех потоков в resSum
            MPI_Reduce(&pieceSum, &resSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Bcast(&resSum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


            //Формула
            X[i] = resSum + B[i] / A[i][i];
        }
        
        //Считаем норму и количество итераций
        if (nRank == 0) {
            norm = 0;
            for (i = 0; i < N; i++) {
                norm += (X[i] - tempX[i]) * (X[i] - tempX[i]);
            }
            norm = sqrt(norm);
            countIterations++;
        }

        MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (norm < epsilon)
            break;
    }
    
    //cout << "Sum: " << resSum << endl;



    if (nRank == 0) {
        double timeRes = MPI_Wtime() - time;
        cout << "Time: " << (timeRes) << "s." << endl;
        cout << "Count iterations: " << countIterations << endl;
        
        cout << "X: ";
        for (int i = 0; i < N; i++) {
            cout << X[i] << " ";
        }
        
        fileInit->outData(X, N, countIterations, timeRes, "output2_" + strNumberFile.str() + ".txt");
    }

    MPI_Finalize();
}

