#include <funother.h>

// z.B. Indexvektoren anschauen


double operator^(const Double  &o1,
                 const double &o2) // globale Funktion // laueft nicht!
{
    return pow(double(o1),o2);
}

double operator^(const Double  &o1,
                 const int &o2) // globale Funktion // laueft nicht!
{
    return pow(double(o1),o2);
}

void show(const vector<int> &vec)
{
    cout <<"[ " ;
    copy ( vec. begin () , vec.end () , ostream_iterator <short >( cout ,
            ", " ) );
    cout <<"]^T" << endl;
    return;
}


// modulo nach MATLAB =! fmod! (fmod erhaelt das Vorzeichen!)
double modulo(double x, double y)
{
    return x - y * floor(x / y);
}


double maxVec(const vector<double> &data)
{
    return *max_element(data.begin(), data.end());
}

double minVec(const vector<double> &data)
{
    return *min_element(data.begin(), data.end());
}

int maxVec(const vector<int> &data)
{
    return *max_element(data.begin(), data.end());
}

int minVec(const vector<int> &data)
{
    return *min_element(data.begin(), data.end());
}

double roundD(double d)
{
    return std::floor(d + 0.5);
}
