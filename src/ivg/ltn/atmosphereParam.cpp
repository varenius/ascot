/* 
 * File:   AtmosphereParam.cpp
 * Author: armin
 * 
 * Created on 13. Oktober 2015, 16:58
 */



#include "atmosphereParam.h"
#include "date.h"

namespace ltn{
  

AtmosphereParam::AtmosphereParam() {
    
}

AtmosphereParam::~AtmosphereParam() {
    
}

AtmosphereParam::AtmosphereParam(ivg::Date time, double temperature, double pressure, double humidity){
        this->time = time;
        this->temperature = temperature;
        this->pressure = pressure;
        this->humidity = humidity;
    }

void AtmosphereParam::pressCor(double gamma, double g, double h, double h0){
    
    double P0 = this->pressure;
    double T0 = this->temperature;

    this->pressure = P0 * pow( 1-gamma*(h-h0)/T0, g/(gamma*R) );
                
                
}

//Getter Setter
void AtmosphereParam::SetHumidity(double humidity) {
    this->humidity = humidity;
}

double AtmosphereParam::GetHumidity() const {
    return humidity;
}

void AtmosphereParam::SetPressure(double pressure) {
    this->pressure = pressure;
}

double AtmosphereParam::GetPressure() const {
    return pressure;
}

void AtmosphereParam::SetTemperature(double temperature) {
    this->temperature = temperature;
}

double AtmosphereParam::GetTemperature() const {
    return temperature;
}

void AtmosphereParam::SetTime(ivg::Date time) {
    this->time = time;
}

ivg::Date AtmosphereParam::GetTime() const {
    return time;
}

}

