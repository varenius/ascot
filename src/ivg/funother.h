#ifndef funother_H
#define funother_H

#include <iterator>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

// Hilfsfunktion2 fuer die verbesserte und matlabaehnlichere Bedienung der Matrix-Klasse
// 18.06.2012 von Malwin Eichborn



class Double
{
        double value;
    public:
        Double(double value_) : value(value_) {}
        operator double() const
        {
            return value;
        }
};
double operator^(const Double  &o1,
                 const double &o2); // globale Funktion mit class Double
double operator^(const Double  &o1,
                 const int &o2); // globale Funktion mit class Double

void show(const vector<int> &vec);

template <typename T>
void showT(const vector<T> &vec)
{
    cout <<"[ " ;
    //copy ( vec. begin () , vec.end () , ostream_iterator <short >( cout , ", " ) );
    for (unsigned int i=0; i< vec.size() ; i++)
        cout << vec.at(i) << ", " ;
    cout <<"]^T" << endl;
    return;
};

// modulo nach MATLAB =! fmod! (fmod erhaelt das Vorzeichen!)
double modulo(double x, double y);

double maxVec(const vector<double> &data);
double minVec(const vector<double> &data);
int maxVec(const vector<int> &data);
int minVec(const vector<int> &data);
double roundD(double d);


#endif /* funother */
