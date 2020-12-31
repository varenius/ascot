#include <QApplication>
#include "matrix.h"
#include "figure.h"

int main(int argc, char **argv)
{
    QApplication a(argc, argv);

    //Generating Testdata
    //ivg::Matrix x1(100.0, 5, 500.0, 1);
    //ivg::Matrix y1 = x1;
    //y1.random();

    // FORMAT POSIBILITIES
    // -------------------
    // Qt::green , usw.
    // WIDTH = 2
    // Qt::SolidLine , DashLine , DotLine , DashDotLine , DashDotDotLine
    // QCPGraph::lsLine , lsNone , lsStepLeft  , lsStepRight  , lsStepCenter , lsImpulse
    // QCPScatterStyle::ssDisc , ssNone , ssDot , ssCross , ssPlus , ssCircle , ssSquare , ssDiamond , ssStar , usw.

    //Plot Data
    //ivg::Matrix values(100,7,1.0);
    //map<string,vector<int>> map1,map2,map3;
    //Statistics window1(values, map1, map2, map3);

    vector<double> x2 = {100.0, 300.0};
    vector<double> y2 = {0.5, 0.5};

   Figure window1;

    //Generate new y-values

    //Plot different Data
//    QPen pen2(Qt::red, 3, Qt::DashLine);
    window1.plot_data(x2,y2, { Qt::red, 4, Qt::DashDotLine, QCPGraph::lsLine, QCPScatterStyle::ssPlus, "A" });

//    QCustomPlot * plot = window1.get_plot();
//    plot->xAxis->setScaleType(QCPAxis::ScaleType::stLogarithmic);
//    plot->xAxis->setScaleLogBase(100.0);

//    Figure window2(x1,y1, { Qt::green, 2, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc });
//    window2.plot_data(x2,y2, { Qt::red, 4, Qt::DashDotLine, QCPGraph::lsLine, QCPScatterStyle::ssPlus, "A" });

    //Save as pdf
//    window1.save_pdf("/opt/bakkari/result.pdf");



//    plot->xAxis->setTickLabelType(QCPAxis::ltDateTime);
//    plot->xAxis->setDateTimeFormat("dd.MM.yyyy");

//    plot->replot();

    return a.exec();
}
