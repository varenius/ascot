/* 
 * File:   AtmosphereParam.h
 * Author: armin
 *
 * Created on 13. Oktober 2015, 16:58
 */

#ifndef ATMOSPHEREPARAM_H
#define	ATMOSPHEREPARAM_H

#include"date.h"

namespace ltn{

class AtmosphereParam {
public:
    AtmosphereParam();
    virtual ~AtmosphereParam();

    AtmosphereParam(ivg::Date time, double temperature, double pressure, double humidity);
    
    //function
    void pressCor(double gamma, double g, double h, double h0);
    
    // Getter Setter
    void SetHumidity(double humidity);
    double GetHumidity() const;
    void SetPressure(double pressure);
    double GetPressure() const;
    void SetTemperature(double temperature);
    double GetTemperature() const;
    void SetTime(ivg::Date time);
    ivg::Date GetTime() const;

    
    
private:
    //Attributes
    ivg::Date time;
    double temperature;
    double pressure;
    double humidity;
    
    const double R = 287.058; //ideal gas constant

};

}

#endif	/* ATMOSPHEREPARAM_H */

