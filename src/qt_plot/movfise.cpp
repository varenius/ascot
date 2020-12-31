/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   movfise.cpp
 * Author: schubert
 *
 * Created on 4. Mai 2016, 15:13
 */

#include <boost/algorithm/string/predicate.hpp>

#include "movfise.h"

Movfise::Movfise() {
    
    widget.setupUi(this);
    connect(widget.size_spinbox, SIGNAL(valueChanged(double)),this, SLOT(printUnitConversions(double)));
    connect(widget.stepsize_spinbox, SIGNAL(valueChanged(double)),this, SLOT(printUnitConversions(double)));
    connect(widget.calculatorPlainTextEdit, SIGNAL(textChanged()),this, SLOT(calculator()));
    widget.size_spinbox->setValue(widget.size_spinbox->value());
    widget.size_spinbox->value();
    widget.calculatorOutput->setTextInteractionFlags(Qt::TextSelectableByMouse);
}


Movfise::Movfise(QWidget *parent) : QDialog(parent){
    widget.setupUi(this);
    connect(widget.size_spinbox, SIGNAL(valueChanged(double)),this, SLOT(printUnitConversions(double)));
    connect(widget.stepsize_spinbox, SIGNAL(valueChanged(double)),this, SLOT(printUnitConversions(double)));
    connect(widget.calculatorPlainTextEdit, SIGNAL(textChanged()),this, SLOT(calculator()));
    widget.size_spinbox->value();
    widget.calculatorOutput->setTextInteractionFlags(Qt::TextSelectableByMouse);
}
//movfise::~Movfise() {
//}

void Movfise::printUnitConversions(double unit) {
    widget.textBrowser->clear();
    std::stringstream info_ss;
    double val = widget.size_spinbox->value();
//    info_ss << "Window Size: " << val << "sec, ";
//    info_ss << val / 60.0 << "min, ";
//    info_ss << val / 3600.0 << "hours, ";
//    info_ss << val / (24.0 * 3600.0) << "days, ";// << std::endl << std::endl;
//    info_ss << val / (24.0 * 3600.0*365.0) << "years" << std::endl<< std::endl;
    info_ss << "Window Size: " << val << "sec," << std::endl;
    info_ss << val / 60.0 << "min," << std::endl;
    info_ss << val / 3600.0 << "hours," << std::endl;
    info_ss << val / (24.0 * 3600.0) << "days," << std::endl;// << std::endl << std::endl;
    info_ss << val / (24.0 * 3600.0*365.0) << "years" << std::endl<< std::endl;
    double valStep = widget.stepsize_spinbox->value();
    if (valStep == 0) {
        info_ss << "Step Size: Original steps" << std::endl;
    } else {
        info_ss << "Step Size: " << valStep << "sec," << std::endl;
        info_ss << valStep / 60.0 << "min," << std::endl;
        info_ss << valStep / 3600.0 << "hours," << std::endl;
        info_ss << valStep / (24.0 * 3600.0) << "days," << std::endl;
        info_ss << valStep / (24.0 * 3600.0*365.0) << "years" << std::endl;
    }

    widget.textBrowser->setText(QString::fromStdString(info_ss.str()));
}

void Movfise::calculator() {
    widget.calculatorOutput->clear();
    std::string s = widget.calculatorPlainTextEdit->toPlainText().toStdString();

    double result = 1;

    std::string delimiter = "*";

    try {
        char buffer [50];
        if (s.find_first_not_of(".*0123456789") == std::string::npos) {
            if (!boost::starts_with(s, "*")) {
                if (!boost::ends_with(s, "*")) {
                    if (!boost::ends_with(s, ".")) {

                        size_t pos = 0;
                        std::string token;
                        while ((pos = s.find(delimiter)) != std::string::npos && (pos = s.find(delimiter)) != std::string::npos-1) {
                            token = s.substr(0, pos);
                            result *= std::stod(token);
                            s.erase(0, pos + delimiter.length());
                        }
                        if (s.compare("") != 0) {
                            result *= std::stod(s);

                            setlocale(LC_NUMERIC, "English");
                            sprintf(buffer, "= %f", result);
                        } else {
                            sprintf(buffer, "%s", "");
                        }
                    }
                }
            }





            //            std::stringstream info_ss;
            //            info_ss << "= " <<std::setprecision(10)<< result;
            //            std::cerr << "H2I";


        } else {
            sprintf(buffer, "%s", "Only mult. (*) allowed.");
        }
        widget.calculatorOutput->setText(QString::fromStdString(buffer));
    } catch (std::exception &ex) {
    }
}