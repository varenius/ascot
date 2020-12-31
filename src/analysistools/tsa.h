#include "matrix.h"
#include "masterfile.h"
#include <complex>
#include <fftw3.h>


#ifndef TIMESERIESANALYSIS_H
#define TIMESERIESANALYSIS_H

namespace ivgat
{
    
enum stepMode {daily, as_reference};

class Tsa {
    
public:
    Tsa();

    /**
     *  \b Description: \n
     *        emp. crosscovariance function für equidistant data
     *        maximal lag is N/10
     *  \param [in] [ivg::Matrix] input first values-vector 
     *              [ivg::Matrix] input second values-vector
     *              [bool] bool to take all lags, not just until tenth of series (Nyquist and factor 5)
     *              [ivg::Matrix] lags vector [unit=integer lags]
     *  \return [ivg::Matrix] output values-vector
     */
//    ivg::Matrix crosscov(ivg::Matrix m1, ivg::Matrix m2, bool all_lags);
    ivg::Matrix crosscorr(ivg::Matrix m1, ivg::Matrix m2, bool all_lags, ivg::Matrix &lags);
    
    /**
     *  \b Description: \n
     *        emp. autocovariance function for equidistant data
     *        maximal lag is N/10
     *  \param [in] [ivg::Matrix] input values-vector 
     *  \return [ivg::Matrix] output values-vector
     */
    ivg::Matrix autocov(ivg::Matrix m1, bool all_lags);
    ivg::Matrix autocorr(ivg::Matrix m1, bool all_lags);
    
    /**
     *  \b Description: \n
     *        emp. autocovariance function für non-equidistant data with #classes
     *        maximal lag is N/10
     *  \param [in] [ivg::Matrix] input time-vector 
     *              [ivg::Matrix] input values-vector 
     *              [int] number of classes
     *              [ivg::Matrix] output time/lag-vector 
     *              [ivg::Matrix] output values-vector
     *  \return 
     */
    void autocovNEQD(ivg::Matrix t1 , ivg::Matrix m1,int numClasses,ivg::Matrix &tout , ivg::Matrix &mout);
    void autocorrNEQD(ivg::Matrix t1 , ivg::Matrix m1,int numClasses,ivg::Matrix &tout , ivg::Matrix &mout); 
    
    /**
     *  \b Description: \n
     *        emp. crosscovariance function für non-equidistant data with #classes.
     *        If second time_vector is 10 higher, max. correlation is at 10.
     *        maximal lag is N/10
     *  \param [in] [ivg::Matrix] input first time-vector 
     *              [ivg::Matrix] input first values-vector
     *              [ivg::Matrix] input second time-vector 
     *              [ivg::Matrix] input second values-vector 
     *              [int] number of classes
     *              [ivg::Matrix] output time/lag-vector 
     *              [ivg::Matrix] output values-vector
     *  \return 
     */
    void crosscovNEQD(ivg::Matrix t1 , ivg::Matrix m1,ivg::Matrix t2, ivg::Matrix m2,int numClasses,ivg::Matrix &tout , ivg::Matrix &mout);
    void crosscorrNEQD(ivg::Matrix t1 , ivg::Matrix m1,ivg::Matrix t2, ivg::Matrix m2,int numClasses,ivg::Matrix &tout , ivg::Matrix &mout);
    void crosscovNEQD2(const ivg::Matrix &t1 , const ivg::Matrix &m1,const ivg::Matrix &t2, const ivg::Matrix &m2,const int numClasses,ivg::Matrix &tout , ivg::Matrix &mout);
    void crosscorrNEQD2(const ivg::Matrix &t1 , const ivg::Matrix &m1,const ivg::Matrix &t2, const ivg::Matrix &m2,const int numClasses,ivg::Matrix &tout , ivg::Matrix &mout);
    
    void crosscovNEQD_Calc(const ivg::Matrix &t1 , const ivg::Matrix &m1,const ivg::Matrix &t2, const ivg::Matrix &m2,const int numClasses,ivg::Matrix &tout , ivg::Matrix &mout,int flag);
    
    
    /**
     *  \b Description: \n
     *        Moving Median of size window, i.e. floor(window/2) to each side
     *  \param [in] [ivg::Matrix] time-series to be filtered
     *              [int] window
     *  \return An instance of the class ivg::Matrix
     */
    ivg::Matrix movingMedian(ivg::Matrix m1, int window);

    /**
     *  \b Description: \n
     *        Moving Median of size aquivalent to double window (e.g. 7 hours) in time unit, i.e. window/2 to each side
     *  \param [in] [ivg::Matrix] input time-series x-values
     *               [ivg::Matrix] input time-series function values
     *              [double] window
     *               [ivg::Matrix] output time-series x values (i.e. new time axis)
     *               [ivg::Matrix] output time-series function values
     *  \return 
     */
    void movingMedian(ivg::Matrix t1, ivg::Matrix m1, double window, double stepsize, ivg::Matrix &t_neu, ivg::Matrix &fout);

    /**
     *
     * With "even_mjd" this method refers differences to whole mjd units/steps, e.g. one day from 0h to 24h (25 values)
     * If flag is not set to "even_mjd" (whole day) steps of pwlf could be 1 day from 0:30h, if data exists hourly at 0:30, 1:30, ...
     *
     *
     */
    void intervalMeanValues(ivg::Matrix t1, ivg::Matrix m1, double window, ivg::Matrix &t_neu, ivg::Matrix &out, stepMode flag);
    void intervalMeanBiases(ivg::Matrix t1, ivg::Matrix m1, ivg::Matrix t2, ivg::Matrix m2, double window, ivg::Matrix &t_neu, ivg::Matrix &out, stepMode flag);
    
     /**
     *  \b Description: \n
     *        Weighted Root Mean Squared Error
     *  \param [in] [ivg::Matrix] residuals
     *              [ivg::Matrix] standard deviations
     *  \return double WRMS-value
     */
    double WRMS(ivg::Matrix m1, ivg::Matrix std);
    
    /**
     *  \b Description: \n
     *        Root Mean Squared Error
     *  \param [in] [ivg::Matrix] residuals
     *  \return double RMS-value
     */
    double RMS(ivg::Matrix m1);

    /**
     *  \b Description: \n
     *        Fourier Amplitude Spectrum for equally-distant data
     *  \param [in] [ivg::Matrix] input time-series x-values
     *               [ivg::Matrix] input time-series function values
     *               [ivg::Matrix] frequency values (i.e. new x-axis)
     *               [ivg::Matrix] amplitude values
     *  \return 
     */
    void fourierAmplitudeSpectrumEQD(ivg::Matrix t, ivg::Matrix f, ivg::Matrix &nu, ivg::Matrix &A);
    
    /**
     *  \b Description: \n
     *        Fourier Amplitude Spectrum for equally-distant data
     *  \param [in] [ivg::Matrix] input time-series x-values
     *               [ivg::Matrix] input time-series function values
     *               [ivg::Matrix] frequency values (i.e. new x-axis)
     *               [ivg::Matrix] complex  spectrum
     *  \return 
     */
    void fourierSpectrumEQD(ivg::Matrix t, ivg::Matrix f, ivg::Matrix &nu, complex<double>*hfa_freq);
    
//    /**
//     *  \b Description: \n
//     *        Fourier Amplitude Spectrum for equally-distant data
//     *  \param [in] [ivg::Matrix] input time-series x-values
//     *               [ivg::Matrix] input time-series function values
//     *               [ivg::Matrix] frequency values (i.e. new x-axis)
//     *               [ivg::Matrix] complex  spectrum
//     *  \return 
//     */
//    void fourierSpectrumEQD2(ivg::Matrix t, ivg::Matrix f, ivg::Matrix &nu, fftw_complex *out);
    
//    /**
//     *  \b Description: \n
//     *        Fourier Amplitude Spectrum for equally-distant data
//     *  \param [in] [ivg::Matrix] frequency values 
//     *               [ivg::Matrix] complex  spectrum
//     *               [ivg::Matrix] time-values (i.e. new x-axis)
//     *               [ivg::Matrix] output time-series function values
//     *  \return 
//     */
//    void inverseFourierTransform(ivg::Matrix nu, complex<double>* hfa_freq, ivg::Matrix &t, ivg::Matrix &x);
    
    /**
     *  \b Description: \n
     *        Lomb-Scargle-Powerspectrum: : Fourier Amplitude Spectrum for non-equally-distant data
     *  \param [in] [ivg::Matrix] input time-series x-values
     *               [ivg::Matrix] input time-series function values
     *               [ivg::Matrix] frequency values (i.e. new x-axis)
     *               [ivg::Matrix] powerspectrum, i.e. amplitudes
     *  \return 
     */
    void fourierSpectrumNEQD(ivg::Matrix m2, ivg::Matrix sin_m2, ivg::Matrix &f, ivg::Matrix &px);
    
    
    /**
     *  \b Description: \n
     *        calculate Two-Sample-Allan-Variance
     *          Possible output plot is plot(tau,avar)
     *  \param [in] [ivg::Matrix] epo: time vector
     *               [ivg::Matrix] freq: frequencies (differencesin functional values)
     *               [ivg::Matrix] tau: time differences
     *               [ivg::Matrix] avar: allan variance values
     *               [ivg::Matrix] avarSig: allan variance sigma (stddev) values
     *  \return 
     */
    void calcTwoSampleAllanVarianceFrequency(const ivg::Matrix epo,const ivg::Matrix freq, ivg::Matrix &tau,ivg::Matrix &avar );
    void calcTwoSampleAllanVariancePhase(const ivg::Matrix t,const ivg::Matrix obs, ivg::Matrix &tau,ivg::Matrix &avar );
    
private:


};

} // namespace ivgat end

#endif /* TIMESERIESANALYSIS_H */

