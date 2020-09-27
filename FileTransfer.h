//
// Created by gorbu on 12.09.2020.
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#ifndef LR2_FILETRANSFER_H
#define LR2_FILETRANSFER_H

using namespace std;

class FileTransfer {
    private: string filePath;

    public:FileTransfer(string filePath) {
        this->filePath = filePath;
    }

    public: void initData(int & N, double & epsilon, double ** &A, double * &B, string fileName) {
        string textLine;
        ifstream init(filePath + fileName);

        if (init.is_open()) {
            int lineNumber = 0;
            while (getline(init, textLine)) {
                lineNumber++;
                stringstream textPiece(textLine);
                if (lineNumber == 1) {
                    textPiece >> N;
                    textPiece >> epsilon;
                    A = new double *[N];
                    B = new double[N];
                } else if (lineNumber == 2) {
                    int i = 0;
                    while (!textPiece.eof()) {
                        textPiece >> B[i++];
                    }
                } else {
                    int i = 0;
                    A[lineNumber - 3] = new double[N];
                    double number = 0;
                    while (!textPiece.eof()) {
                        textPiece >> number;
                        A[lineNumber - 3][i] = number;
                        i++;
                    }
                }
            }
        }
    }

    public: void outData(double * numberIterations, int N,  int countApproximate, double times, string fileName) {

        ofstream out;
        out.open(filePath + fileName);
        if(out.is_open()) {
            out << countApproximate << endl;
            for(int i = 0; i < N; i++) {
                out << fixed << setprecision(5) << numberIterations[i] << " ";
            }
            out << endl;
            //for (int i = 0; i < 4; ++i) {
                out << setprecision(5) << times << " ";
            //}
        }
    }
};


#endif //LR2_FILETRANSFER_H
