/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   movfise.h
 * Author: schubert
 *
 * Created on 4. Mai 2016, 15:13
 */

#ifndef _MOVINGFILTERSETTINGS_H
#define _MOVINGFILTERSETTINGS_H

#include "ui_movfise.h"
#include <iostream>
#include <sstream> 
#include <string> 
#include <stdexcept>

class Movfise : public QDialog {
    Q_OBJECT
public:
    Movfise();
    Movfise(QWidget *parent);
//    virtual ~Movfise();

    Ui::Movfise widget;
private:
//signals:
//    void valueChanged(int val);
private slots:
     void printUnitConversions(double unit);
     void calculator();
};

#endif /* _MOVINGFILTERSETTINGS_H */
