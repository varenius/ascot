#include "tsa.h"

namespace ivgat
{

Tsa::Tsa() {
}

ivg::Matrix Tsa::autocorr(ivg::Matrix m1, bool all_lags) {
    ivg::Matrix out = autocov(m1, all_lags);
    return out / out(0);
}

//ivg::Matrix Tsa::crosscorr(ivg::Matrix m1, ivg::Matrix m2, bool all_lags) {
//    ivg::Matrix out = crosscov(m1, m2, all_lags);
//    return out ;/// out(0);
//}

// hier nicht erwartungstreuer Schätzer: xcorr mit 'biased'
ivg::Matrix Tsa::autocov(ivg::Matrix m1, bool all_lags) {
    
    int N = m1.size(1);
    int maxlag;
    if(all_lags){
        maxlag = N-1;
    }else{
        maxlag = floor(((double) N) / 10);
    }

    ivg::Matrix out(maxlag+1, 1, 0);
    for (int k = 0; k <= maxlag; ++k) {
        double v = 0;
        for (int i = 0; i < N - k; ++i) {
            v = v + m1(i) * m1(i + k);
        }
        out(k,0) = v;
    }

    return out/N;
}

// entspricht octave xcorr mit erwartungstreuem Schätzer 'unbiased'
ivg::Matrix Tsa::crosscorr(ivg::Matrix m1, ivg::Matrix m2 , bool all_lags , ivg::Matrix &lags) {
    // use only if m1 and m2 have identical time vectors
    
    int N = max(m1.size(1),m2.size(1));
    int maxlag;
    if(all_lags){
        maxlag = N-1;
    }else{
        maxlag = floor(((double) N) / 10);
    }
    
    lags = ivg::Matrix(2*maxlag+1,1,0.0);

    ivg::Matrix out(maxlag+1,1, 0);
    for (int k = 0; k <= maxlag; ++k) {
        double v = 0;
        for (int i = 0; i < N - k; ++i) {
            v = v + m2(i) * m1(i + k);
        }
        out(k,0) = v/(N-k);
        lags(k+maxlag,0) = k;
    }

    ivg::Matrix out2(maxlag,1, 0);
    for (int k = 1; k <= maxlag; ++k) {
        double v = 0;
        for (int i = N-1; i >= k; --i) {
            v = v + m2(i) * m1(i - k);
        }
        out2(maxlag-k,0) = v/(N-k);
        lags(maxlag-k,0) = -k;
    }
    
    out2.append_rows(out);
    
    return out2;
}



void Tsa::autocorrNEQD(ivg::Matrix t1 , ivg::Matrix m1,int numClasses,ivg::Matrix &tout , ivg::Matrix &mout) {
    crosscovNEQD( t1 , m1, t1, m1, numClasses,tout , mout);
    mout = mout/mout(0);
}
void Tsa::autocovNEQD(ivg::Matrix t1 , ivg::Matrix m1,int numClasses,ivg::Matrix &tout , ivg::Matrix &mout) {
    
    crosscovNEQD( t1 , m1, t1, m1, numClasses,tout , mout);

}

void Tsa::crosscorrNEQD(ivg::Matrix t1 , ivg::Matrix m1,ivg::Matrix t2, ivg::Matrix m2,int numClasses,ivg::Matrix &tout , ivg::Matrix &mout) {
    crosscovNEQD(t1 , m1, t2, m2,numClasses, tout , mout);
    mout = mout/mout.max();
}

void Tsa::crosscorrNEQD2(const ivg::Matrix &t1 , const ivg::Matrix &m1,const ivg::Matrix &t2, const ivg::Matrix &m2,const int numClasses,ivg::Matrix &tout , ivg::Matrix &mout) {
    cerr<<"aa"<<endl;
    crosscovNEQD2(t1 , m1, t2, m2,numClasses, tout , mout);
    cerr<<"bb"<<endl;
    mout = mout/mout.max();
}
void Tsa::crosscovNEQD2(const ivg::Matrix &t1 , const ivg::Matrix &m1,const ivg::Matrix &t2, const ivg::Matrix &m2,const int numClasses,ivg::Matrix &tout , ivg::Matrix &mout) {
    
    crosscovNEQD_Calc(t1 , m1, t2, m2,numClasses, tout , mout, 0);

}


//void Tsa::crosscovNEQD2(const ivg::Matrix &t1 , const ivg::Matrix &m1,const ivg::Matrix &t2, const ivg::Matrix &m2,const int numClasses,ivg::Matrix &tout , ivg::Matrix &mout) {
//    cerr<<"a1"<<endl;
//    double maxlag = std::max(abs(t1(t1.size(1)-1,0)-t2(0,0)),abs(t2(t2.size(1)-1,0)-t1(0,0)))/10;
//    cerr<<"a2"<<maxlag<<endl;
//    double schrittweite = maxlag/numClasses;
//
//cerr<<"a"<<endl;
////    ivg::Matrix mult = m1 * m2.transpose();
//    cerr<<"b"<<endl;
//    ivg::Matrix lag = (t1 * ivg::Matrix(1, t2.size(1), 1) - (t2 * ivg::Matrix(1, t1.size(1), 1)).transpose()).abs();
////    ivg::Matrix lag = t1 * ivg::Matrix(1, t2.size(1), 1);
////    ivg::Matrix lag2 = (t2 * ivg::Matrix(1, t1.size(1), 1)).transpose();
//    
////    plus_product_of(const Matrix & A,const Matrix & B,
////                             CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB)
//    
//    double minimum = lag.min();
//    cerr<<"c "<<minimum<<" "<<schrittweite<<" "<<minimum+schrittweite* (numClasses-1)<<endl;
//    tout = ivg::Matrix(minimum, schrittweite, minimum+schrittweite * (numClasses-1), 1);
//    mout = ivg::Matrix(tout.size(1), 1, 0);
//cerr<<"d"<<endl;
//    double sh = schrittweite/2;
//    
//    for (int k = 0; k < numClasses; ++k) {
//        vector<int> idx;
//
//        if (k == numClasses-1) {
//            lag.find_elem(lt, tout(k)+schrittweite-sh,ge, tout(k)-sh, idx);
//        } else {
//            lag.find_elem(lt, tout(k + 1)-sh, ge, tout(k)-sh, idx);
//        }
//        double summe = 0;
//        for (int kk=0;kk<idx.size();kk++){
////            i+j*_rows)
//            int iii = idx.at(kk)%lag.rows();
//            int jjj = (idx.at(kk)-iii)/lag.rows();
//            summe += m1(iii)*m2(jjj);
//        }
//        mout(k) = (summe) / idx.size();
////        ivg::Matrix m = mult(idx); 
////        mout(k) = (m.sum_col())(0) / m.size(1);
//    }
//   cerr<<"e"<<endl; 
//}


void Tsa::crosscovNEQD_Calc(const ivg::Matrix &t1 , const ivg::Matrix &m1,const ivg::Matrix &t2, const ivg::Matrix &m2,const int numClasses,ivg::Matrix &tout , ivg::Matrix &mout,int flag) {

    double maxlag;
    if(flag==1){
        maxlag = std::max(abs(t1(t1.size(1)-1,0)-t2(0,0)),abs(t2(t2.size(1)-1,0)-t1(0,0)))/10;
    }else{
        maxlag = std::max((t1(t1.size(1)-1,0)-t2(0,0)),(t2(t2.size(1)-1,0)-t1(0,0)))/10;
    }

    double schrittweite = maxlag/numClasses;


//    ivg::Matrix mult = m1 * m2.transpose();

    ivg::Matrix lag;
    if(flag==1){
        lag = (t1 * ivg::Matrix(1, t2.size(1), 1) - (t2 * ivg::Matrix(1, t1.size(1), 1)).transpose()).abs();
    }else{
        lag = (t1 * ivg::Matrix(1, t2.size(1), 1) - (t2 * ivg::Matrix(1, t1.size(1), 1)).transpose());
    }

    double maximum = lag.max();
    double minimum = lag.min();

    maxlag =  (maximum-minimum)/10;
    schrittweite = maxlag/(numClasses);
    minimum = -schrittweite*floor((numClasses-1)/2);
    maximum = schrittweite*ceil((numClasses-1)/2);

    tout = ivg::Matrix(minimum, schrittweite, minimum+schrittweite * (numClasses-1), 1);
    mout = ivg::Matrix(tout.size(1), 1, 0);

    double sh = schrittweite/2;
    
    for (int k = 0; k < numClasses; ++k) {
        vector<int> idx;

        if (k == numClasses-1) {
            lag.find_elem(lt, tout(k)+schrittweite-sh,ge, tout(k)-sh, idx);
        } else {
            lag.find_elem(lt, tout(k + 1)-sh, ge, tout(k)-sh, idx);
        }
        double summe = 0;
        for (int kk=0;kk<idx.size();kk++){
            int iii = idx.at(kk)%lag.rows();
            int jjj = (idx.at(kk)-iii)/lag.rows();
            summe += m1(iii)*m2(jjj);
        }
        mout(k) = (summe) / idx.size();
    }
}

void Tsa::crosscovNEQD(ivg::Matrix t1 , ivg::Matrix m1,ivg::Matrix t2, ivg::Matrix m2,int numClasses,ivg::Matrix &tout , ivg::Matrix &mout) {

    double maxlag = std::max(t1(t1.size(1)-1)-t1(0),t2(t2.size(1)-1)-t2(0))/10;
    double schrittweite = maxlag/numClasses;


    ivg::Matrix mult = m1 * m2.transpose();
    ivg::Matrix lag = (t1 * ivg::Matrix(1, t2.size(1), 1) - (t2 * ivg::Matrix(1, t1.size(1), 1)).transpose()).abs();
    
    double minimum = lag.min();
    
    tout = ivg::Matrix(minimum, schrittweite, schrittweite * (numClasses-1), 1);
    mout = ivg::Matrix(tout.size(1), 1, 0);

    double sh = schrittweite/2;
    
    for (int k = 0; k < numClasses; ++k) {
        vector<int> idx;

        if (k == numClasses-1) {
            lag.find_elem(lt, tout(k)+schrittweite-sh,ge, tout(k)-sh, idx);
        } else {
            lag.find_elem(lt, tout(k + 1)-sh, ge, tout(k)-sh, idx);
        }

        ivg::Matrix m = mult(idx);
        mout(k) = (m.sum_col())(0) / m.size(1);
    }
}

ivg::Matrix Tsa::movingMedian(ivg::Matrix m1, int window) { // für  Anzahl Werte             //funktioniert noch nicht

    int min_Iter;
    if ((double) (window / 2) != ((double) window) / 2) {
        min_Iter = (window - 1) / 2;
    }

    ivg::Matrix out(m1.size(1) - 2 * min_Iter, 1, 0);

    for (int k = min_Iter; k < m1.size(1) - min_Iter; ++k) {
        out(k - min_Iter, 0) = m1.get_sub(k - min_Iter, 0, k + min_Iter, 0).median();
    }

    return out;
}

void Tsa::movingMedian(ivg::Matrix t1, ivg::Matrix m1, double window, double stepsize, ivg::Matrix &t_neu, ivg::Matrix &out) { // für Zeit statt Anzahl Werte

    double diff = t1(t1.rows()-1) - t1(0);
    if(diff < window)
    {
        cerr << "void Tsa::movingMedian(ivg::Matrix t1, ivg::Matrix m1, double window, double stepsize, ivg::Matrix &t_neu, ivg::Matrix &fout): Selected window bigger than timescale." << endl;
        t_neu = t1;
        out = m1;
    }
    else
    {
        int min_Index = -1;
        double valA;
        int max_Index = -1;
        double valB;
        double add = 0;
        for (int k = 0; k < t1.size(1); ++k) {
            add += t1(k + 1, 0) - t1(k, 0);
            if (add > window / 2.0) {
                min_Index = k;
                valA = t1(k, 0);
                break;
            }
        }

        add = 0;
        for (int k = 0; k < t1.size(1); ++k) {
//            cerr << m1.size(1) - k - 2 << " " << m1.size(1) - k - 1 << endl;
            add += t1(m1.size(1) - k - 1) - t1(m1.size(1) - k - 2);
            if (add > window / 2.0) {
                max_Index = m1.size(1) - k - 1;
                valB = t1(m1.size(1) - k - 1);
                break;
            }
        }

//        ivg::Matrix out;
        vector<int> int_idx;
        if (stepsize == 0) {
            out = ivg::Matrix(max_Index - min_Index + 1, 1, 0);
            for (int k = min_Index; k <= max_Index; ++k) {
                int_idx = t1.find_idx(ge, t1(k) - window / 2, le, t1(k) + window / 2);
                out(k - min_Index, 0) = m1.get_sub(int_idx,{0}) .median();
            }
            t_neu = t1.get_sub(min_Index, 0, max_Index, 0);
        } else {
            
            int anzahlWerte = (valB-valA)/stepsize + 1;
            out = ivg::Matrix(anzahlWerte, 1, 0);
            t_neu = ivg::Matrix(anzahlWerte, 1, 0);
//            int k=0;
            double val = valA;
//            for (double val = valA; val <= valB; val+=stepsize) {
            for (int k = 0; k < anzahlWerte; k++) {
                int_idx = t1.find_idx(ge, val - window / 2, le, val + window / 2);
                out(k, 0) = m1.get_sub(int_idx,{0}) .median();
                t_neu(k, 0) = val;
//                k++;
                val+=stepsize;
            }
        }
        
    }
}


void Tsa::intervalMeanBiases(ivg::Matrix t1, ivg::Matrix m1, ivg::Matrix t2, ivg::Matrix m2, double window, ivg::Matrix &t_neu, ivg::Matrix &out, stepMode flag) { 
    
    ivg::Matrix t1_neu, out1,t2_neu, out2;
    this->intervalMeanValues(t1, m1,  window, t1_neu, out1, flag);
    this->intervalMeanValues(t2, m2,  window, t2_neu, out2, flag);
    
    
    sort(t1_neu.begin(),t1_neu.end());
    sort(t2_neu.begin(),t2_neu.end());
    
    vector<double> vector3;
    
    set_intersection(t1_neu.begin(),t1_neu.end(),t2_neu.begin(),t2_neu.end(),back_inserter(vector3));
    
    
    out = ivg::Matrix(vector3.size(),1);
    
    int k = 0;
    for (double d: vector3){

        int p1 = std::find (t1_neu.begin(), t1_neu.end(), d) - t1_neu.begin();

        int p2 = std::find (t2_neu.begin(), t2_neu.end(), d) - t2_neu.begin();


        out(k,0) =  out2(p2) - out1(p1);
        k++;

    }
    
    t_neu = ivg::Matrix(vector3);
    
}


void Tsa::intervalMeanValues(ivg::Matrix t1, ivg::Matrix m1, double window, ivg::Matrix &t_neu, ivg::Matrix &out, stepMode flag) { 

    double stepsize = window;
    
    double diff = t1(t1.rows()-1) - t1(0);
    
    if(diff < window)
    {
        cerr << "void Tsa::movingMedian(ivg::Matrix t1, ivg::Matrix m1, double window, double stepsize, ivg::Matrix &t_neu, ivg::Matrix &fout): Selected window bigger than timescale." << endl;
        t_neu = t1;
        out = m1;
    }
    else
    {
        int min_Index = -1;
        double valA;
        int max_Index = -1;
        double valB;
        double add = 0;
        for (int k = 0; k < t1.size(1); ++k) {
            add += t1(k + 1, 0) - t1(k, 0);
            if (add > window / 2.0) {
                min_Index = k;
                valA = t1(k, 0);
                break;
            }
        }

        add = 0;
        for (int k = 0; k < t1.size(1); ++k) {
//            cerr << m1.size(1) - k - 2 << " " << m1.size(1) - k - 1 << endl;
            add += t1(m1.size(1) - k - 1) - t1(m1.size(1) - k - 2);
            if (add > window / 2.0) {
                max_Index = m1.size(1) - k - 1;
                valB = t1(m1.size(1) - k - 1);
                break;
            }
        }
        if(flag==stepMode::daily){
            valA = floor(valA);
        }

            vector<int> int_idx;
            
            int anzahlWerte = (valB-valA)/stepsize + 1;
            out = ivg::Matrix(anzahlWerte, 1, 0);
            t_neu = ivg::Matrix(anzahlWerte, 1, 0);
//            int k=0;
            double val = valA;
//            for (double val = valA; val <= valB; val+=stepsize) {
            for (int k = 0; k < anzahlWerte; k++) {

                
//                int_idx = t1.find_idx(ge, val - window / 2, le, val + window / 2);
                int_idx = t1.find_idx(ge, val , le, val + window );
//                std::copy (int_idx.begin(), int_idx.end(), std::ostream_iterator<int>(std::cerr,", "));
                if(int_idx.size()>0){
//                out(k, 0) = m1.get_sub(int_idx,{0}) .median();
//                t_neu(k, 0) = val/*+ window / 2.0*/;
                t_neu(k, 0) = val+ window / 2.0;
//                m1.get_sub(int_idx,{0}). 
                        ivg::Matrix A(int_idx.size(),1,1);
                        A.append_cols(t1.get_sub(int_idx,{0}));

                      ivg::Matrix param = m1.get_sub(int_idx,{0});  
                      param.ols(A);
                        
                        out(k, 0) = param(1)*(val+ window / 2.0)+param(0);
//                        out(k, 0) = param(1)*(val/*+ window / 2.0*/)+param(0);
                }
                val+=stepsize;
            }
        
        
    }
}



double Tsa::WRMS(ivg::Matrix v1, ivg::Matrix std) {
    ivg::Matrix p = std^(-2);
    ivg::Matrix summe = (((v1^2).transpose()*p)*v1.size(1))/
                        (p.sum_col()(0)*(v1.size(1) - 1));
    return sqrt(summe(0, 0));
};

double Tsa::RMS(ivg::Matrix v1) {
    ivg::Matrix summe = (((v1).transpose()) * (v1)) / (v1.size(1) - 1);
    return sqrt(summe(0, 0));
};

void Tsa::fourierAmplitudeSpectrumEQD(ivg::Matrix t, ivg::Matrix f, ivg::Matrix &nu, ivg::Matrix &A) {

    int n_points = f.size(1);




    complex<double>* hfafreq = new std::complex<double> [n_points];
    fftw_plan p1 = fftw_plan_dft_r2c_1d(f.length(), f.data_ptr(),
            reinterpret_cast<fftw_complex*> (hfafreq),
            FFTW_ESTIMATE);
    fftw_execute(p1);



    int n_spectrum_points = (n_points + 1) / 2;
    A = ivg::Matrix(1, n_spectrum_points, 0);
    //   ivg::Matrix A(1, n_spectrum_points, 0);
    A(0) = hfafreq[0].real() * hfafreq[0].real();
    for (int i = 1; i < n_spectrum_points; ++i) {
        A(i) = hfafreq[i].real() * hfafreq[i].real() + hfafreq[n_points - i].real() * hfafreq[n_points - i].real();
    }


    double schrittweite = t(2) - t(1);

    double nu_Nyquist = 1 / (2 * schrittweite)*2 * M_PI;
    double delta_nu = nu_Nyquist / (n_spectrum_points - 1);
    nu = ivg::Matrix(0, delta_nu, nu_Nyquist, 1);

    return;
}

void Tsa::fourierSpectrumEQD(ivg::Matrix t, ivg::Matrix f, ivg::Matrix &nu, complex<double>* hfafreq) {

    int n_points = f.size(1);


    hfafreq = new std::complex<double> [n_points];
    fftw_plan p1 = fftw_plan_dft_r2c_1d(f.length(), f.data_ptr(),
            reinterpret_cast<fftw_complex*> (hfafreq),
            FFTW_ESTIMATE);
    fftw_execute(p1);



    int n_spectrum_points = (n_points + 1) / 2;


    double schrittweite = t(2) - t(1);

    double nu_Nyquist = 1 / (2 * schrittweite)*2 * M_PI;
    double delta_nu = nu_Nyquist / (n_spectrum_points - 1);
    nu = ivg::Matrix(0, delta_nu, nu_Nyquist, 1);

    return;
}

//void Tsa::fourierSpectrumEQD2(ivg::Matrix t, ivg::Matrix f, ivg::Matrix &nu, fftw_complex *out) {
//
//    int N = f.size(1);
//    double *input;
//    input = new double [N];
//
//    for (int i = 0; i < N; i++) {
//
//        input[i] = f(i);
//    }
//
//
//    out = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * (N / 2) + 1);
//    fftw_plan p;
//    p = fftw_plan_dft_r2c_1d(N, input, out, FFTW_ESTIMATE);
//    fftw_execute(p);
//
//
//    ivg::Matrix ampl(1, N / 2 + 1, 0);
//    ivg::Matrix ampl_x(1, N / 2 + 1, 0);
//
//
//    for (int i = 0; i < (N / 2) + 1; i++) {
//
//        ampl(i) = (out[i][0]);
//
//        ampl_x(i) = (input[i]);
//
//    }
//}
//
//void Tsa::inverseFourierTransform(ivg::Matrix nu, complex<double> *hfafreq, ivg::Matrix &t, ivg::Matrix &x) {
//
//    int n_pts = nu.size(1);
//
//    fftw_plan p3 = fftw_plan_dft_c2r_1d(2 * n_pts,
//            reinterpret_cast<fftw_complex*> (hfafreq),
//            x.data_ptr(), FFTW_ESTIMATE);
//    fftw_execute(p3);
//
//    t = ivg::Matrix(2 * n_pts, 1, 0.0);
//
//    return;
//}

void Tsa::fourierSpectrumNEQD(ivg::Matrix m2, ivg::Matrix sin_m2, ivg::Matrix &f, ivg::Matrix &px) {

    double mean_Interval = m2.diff().meanD();

    int ofac = 4; // oversampling factor
    double hifac = 1; // 
    //    double fhi = ??;// highest frequency that can be analyzed         //common way to chose fhi is to take fnyq that would be obtained if the N data points were evenly spaced over the same time interval
    //    double fnyq = ??; // Nyquist frequency 
    //    double hifac = fhi/fnyq; // high?? factor       Press et al. (1992)

    f = ivg::Matrix((pow(2 * mean_Interval, -1)) / (sin_m2.size(1) * ofac), (pow(2 * mean_Interval, -1)) / (sin_m2.size(1) * ofac), hifac * pow(2 * mean_Interval, -1), 1);


    sin_m2 -= sin_m2.meanD();
    px = ivg::Matrix(f.size(1), f.size(2), 100);
    ivg::Matrix e(m2.size(2), m2.size(1), 1);
    for (int k = 0; k < f.size(1); k++) {
        double wrun = f(k)*(2 * M_PI);

        double a = std::atan2(((m2 * (2 * wrun)).sin()).meanD(), ((m2 * (2 * wrun)).cos()).meanD());
        ivg::Matrix c = (m2 * wrun - a / 2).cos();

        double f1 = pow((sin_m2.transpose() * c)(0), 2);
        double f2 = (e * c.pow(2))(0);

        ivg::Matrix s = (m2 * wrun - a / 2).sin();

        double f3 = pow((sin_m2.transpose() * s)(0), 2);
        double f4 = (e * s.pow(2))(0);
        px(k) = 1 / (2 * pow(sin_m2.stdD(), 2)) * f1 / f2 + f3 / f4;

    }
    f *= (2 * M_PI);
    return;
}



// ---------------------------------------------------------------------------
void Tsa::calcTwoSampleAllanVarianceFrequency(const ivg::Matrix epo,const ivg::Matrix freq, ivg::Matrix &tau,ivg::Matrix &avar )
// two sample allan variance
// For geodetic time series use this function. 
// In case of clock data apply to frequency data, which is the .diff() operator on phase data.
//
// For more information, see:
// Fabian Czerwinski, Andrew C. Richardson, and Lene B. Oddershede,
// "Quantifying Noise in Optical Tweezers by Allan Variance,"
// Opt. Express 17, 13255-13269 (2009)
// http://dx.doi.org/10.1364/OE.17.013255
// ---------------------------------------------------------------------------
{
   //
   // determine minimum interval length
   //
   ivg::Matrix dt = epo.diff();
   //double dt_min = 5.0*dt.meanD(); // 5 values in an interval of mean length
   double dt_min = 3.0*dt.max();    // 3 times the max. time difference

   int nLevel = int( ( epo.max()-epo.min() )/dt_min )-1;

   avar.resize( nLevel,1,0.0 );
   tau.resize( nLevel,1,0.0 );

   //
   // calculate allan variance for different levels starting with the lowest
   // resolution
   //
   for( int i=2; i<=nLevel+1; ++i )
   {
      double tau_i = ( epo.max()-epo.min() )/double(i); // cluster length
      ivg::Matrix nCluster( i,1,0.0 );                       // obs. in clusters
      ivg::Matrix ave( i,1,0.0 );                            // cluster averages

      for( int j=0; j<epo.size(1); ++j )
      {
         // find interval of this observation
         double dt0 = epo(j) - epo(0);
         double idx1 = ( dt0/tau_i );
         int idx = int( dt0/tau_i );

         // with huge numbers, e.g., several days with 1 s resolution,
         // division results to integers although this is not right
         if( ( j == epo.size(1)-1 && fmod( dt0,tau_i ) < 1e-10 ) ||
             tau_i-fmod( dt0,tau_i ) < 1e-10 )
            idx--;

         // sum up values and number of observations in the corresponding
         // interval
         ave( idx ) += freq( j );
         nCluster( idx )++;
      }

      ave = ave.div_elem( nCluster );
      ivg::Matrix tmp = ave.diff();
      tmp = tmp^2.0;
      avar( i-2 ) = tmp.meanD()/2.0 ;
      tau( i-2 ) = tau_i;
   }
}
// ---------------------------------------------------------------------------
void Tsa::calcTwoSampleAllanVariancePhase(const ivg::Matrix t, const ivg::Matrix obs, ivg::Matrix &tau, ivg::Matrix &avar)
    // two sample allan variance
    // In case of clock data apply to phase data, which is the cumulative sum of frequency data.
    // Gives an identical estimate for the Allan Variance of geodetic time series data when using the cumulative sum (multiplied with stepsize) of the original geodetic observations. 
    //
    // For more information, see:
    // Estimating the Allan variance in the presence of long periods of missing data and outliers
    // Ilaria Sesia and Patrizia Tavella
    // Published 5 December 2008 • 2008 BIPM and IOP Publishing Ltd
    // Metrologia, Volume 45, Number 6 
    // http://dx.doi.org/10.1088/0026-1394/45/6/S19
    // ---------------------------------------------------------------------------
    {
        //
        // determine minimum interval length
        //
        double tau0 = t(2) - t(1);


        //   //double dt_min = 5.0*dt.meanD(); // 5 values in an interval of mean length
        //   double dt_min = 3.0*tau0;    // 3 times the max. time difference
        //   int nLevel = int( ( t.max()-t.min() )/dt_min )-1;
        //   int sz = nLevel;

        int nLevel = t.size(1);
        int sz;
        if (t.size(1) % 2 == 0) {
            sz = t.size(1) / 2 - 1;
        } else {
            sz = floor(t.size(1) / 2);
        }

        avar.resize(sz, 1, 0.0);
        tau.resize(sz, 1, 0.0);

        //
        // calculate allan variance for different levels starting with the lowest
        // resolution
        //

        for (int i = 1; i <= sz; ++i) {

            double taui = ((double) (i)) * tau0;
            ivg::Matrix sigma2((t.size(1) - 2 * i), 1, 0);
            for (int j = 1; j <= t.size(1) - 2 * (i); ++j) {
                sigma2(j - 1, 0) = pow((obs(j + 2 * i - 1) - 2 * obs(j + i - 1) + obs(j - 1)), 2);
            }

            tau(i - 1, 0) = taui;
            if (sigma2.size(1) == 1)
                avar(i - 1, 0) = sigma2(0) / (2 * pow(taui, 2)*(t.size(1) - 2 * (i)));
            else
                avar(i - 1, 0) = sigma2.sum_col()(0) / (2 * pow(taui, 2)*(t.size(1) - 2 * (i)));
        }
    }



} // namespace ivgat end