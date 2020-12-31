/* 
 * File:   plot.h
 * Author: iddink
 *
 * Created on 8. Juli 2015, 09:22
 */
#ifndef PLOT_H
#define	PLOT_H


#include <QWidget>
#include <QInputDialog>
#include <QMenu>
#include <QDate>
#include <QDateTime>
#include <QRubberBand>
#include <QPoint>
#include <QToolTip>
#include <QMessageBox>

#include "qcustomplot.h"
#include "matrix.h"

#include "exportdialog.h"
#include "graphdialog.h"
#include "projection.h"
#include "ui_plot.h"
#include "ui_axisrange.h"
#include "ui_highlight.h"
#include "ivg_const.h"
#include "ui_statistics_info.h"

#include "gsl/gsl_histogram.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_sf_legendre.h"


#include "ui_movfise.h"
#include <tsa.h>

struct QFormat
{
    // color for line AND marker (Qt::green)
    QColor color;
    // width of the line (2.0)
    qreal line_width;
    // style of the line (Qt::PenStyle::DashDotLine)
    Qt::PenStyle pen_style;
    // style of the line-connection (QCPGraph::LineStyle::lsImpulse)
    QCPGraph::LineStyle line_style;
    // style of the marker (QCPScatterStyle::ssDisc)
    QCPScatterStyle::ScatterShape marker;
    // string containing legend
    QString legend;
    // width of the line from the marker
    qreal marker_line_width;
    // marker size
    qreal marker_size;
};


class Plot  : public QWidget 
{
    Q_OBJECT
    
public:
        
    Plot(QWidget *parent = 0);
    
    Plot(ivg::Matrix x_ivg, ivg::Matrix y_ivg, QFormat format, QWidget *parent = 0);
    
    Plot(ivg::Matrix x_ivg, ivg::Matrix y_ivg, QWidget *parent = 0);
    
    Plot(const Plot& orig);
    
    virtual ~Plot();
    
    double histogram(ivg::Matrix y_ivg, int num_bins, QFormat format, int subplot_cnt=0, bool addToLegend = false,  bool inPercent = false);
    
    void boxplot(ivg::Matrix y_ivg, QFormat format, int subplot_cnt=0);
    
    void projection(QCPProjection projection, ivg::Matrix lam, ivg::Matrix phi, QFormat format, string title="Projection", vector<string> tooltips = vector<string>());
    
    void add_projected_colormap();
    
    QCPGraph* plot_data(ivg::Matrix x_ivg, ivg::Matrix y_ivg, int subplot_cnt=0, vector<string> tooltips = vector<string>());
    
    QCPGraph* plot_data(ivg::Matrix x_ivg, ivg::Matrix y_ivg, QFormat format, int subplot_cnt=0, vector<string> tooltips = vector<string>());
  
    QCPGraph*  plot_data(vector<double> x_std, vector<double> y_std, QFormat format, int subplot_cnt=0, vector<string> tooltips = vector<string>());
    
    void plot_mjd_series(ivg::Matrix mjd, ivg::Matrix y_ivg, QFormat format, int subplot_cnt=0, string time_format= std::string("dd-MMM-yy"), vector<string> tooltips = vector<string>());
    void plot_mjd_series_std(ivg::Matrix mjd, ivg::Matrix y_ivg, ivg::Matrix y_std, QFormat format, int subplot_cnt=0, string time_format= std::string("dd-MMM-yy"), vector<string> tooltips = vector<string>());
    void plot_mjd_series_hpd(ivg::Matrix mjd, ivg::Matrix y_ivg, ivg::Matrix hpdu, ivg::Matrix hpdl, QFormat format, int subplot_cnt=0, string time_format= std::string("dd-MMM-yy"), vector<string> tooltips = vector<string>());
    
    void plot_yerr(ivg::Matrix x_ivg, ivg::Matrix y_ivg, ivg::Matrix s_ivg, QFormat format, int subplot_cnt=0);
  
    void plot_yerr(vector<double> x_std, vector<double> y_std, vector<double> s_std, QFormat format, int subplot_cnt=0);
    
    // identify => e.g. {1}R1 means, look at row one of tooltips, is there "R1"?
    // identify => e.g. [>1]0.238 means, look at row one of tooltips, is the value bigger than 0.238 (works with <, >, =)
    // identify => e.g. {1}R1 \n {1}R 4  means, look at row one of tooltips, is there "R1" or "R4"? Need to be seperated with \n
    void highlight(string identify, QFormat format, bool remove_old, int subplot_cnt = 0);
    
    void imagesc(ivg::Matrix A, QCPColorGradient::GradientPreset color_gradient = QCPColorGradient::gpCold );
    
    int add_rowplot(string left_label, string bottom_label );
    
    void set_title(string title);
        
    QCustomPlot * get_plot(){return (_qcp); }
    
signals:
    
    // emit a signal in case a tooltip is showd by mouse over point
    void tooltipShowed(QString &tooltip);
    // emit a signal if a point is removed with CTRL + ALT
    void pointRemoved(QString &tooltip);
    
public slots:
    
    void selectionChanged();
    void contextMenuRequest(QPoint pos);
    
    void openAxisRangeDialog();
    void openAxisStyleDialog();
    void openExportDialog();
    void openHighlightDialog();
    
    void performGraphAnalysis();
        
    void mouseWheelTurned(QWheelEvent* event);
    void zoomBoxRequest(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent * event);
    void mouseReleaseEvent(QMouseEvent * event);
    void rescalePlot(QMouseEvent*);
    
    void titleDoubleClick(QMouseEvent* event, QCPPlotTitle* title);
    void axisLabelDoubleClick(QCPAxis *axis, QCPAxis::SelectablePart part);
    void legendDoubleClick(QCPLegend *legend, QCPAbstractLegendItem *item);
    
    void showPointToolTip(QMouseEvent *event);
    
    void removeSelectedGraph();
    
    void moveLegend();
    void removeLegend();
  
    void changeStyleSelectedGraph(); 
    void changeProjectionProperties();
    
    void switchErrorbars();
    void switchFreqPeriod(); // by Till
    void tickLabelsFreqDomain(); // by Till

    vector<string> tokenize_string( string &str, const char delimiter );
    
    map< QCPGraph *, vector<int> >* get_index_ptr(){ return &_index; };
    
private:
    
    // graphical user interface generated by QT Designer
    Ui::Plot _ui;
    
    // for all QCustomPlot interactions
    QCustomPlot *_qcp;
    
    // for further analysis plots, secondary plot
    Plot *_sec_plot;
    
    // adds a graph to a plot and returns the pointer to the added graph
    QCPGraph * add_graph(vector<double> x_std, vector<double> y_std, QCPAxis * keyAxis = 0, QCPAxis * valueAxis = 0);
    
    // adds an errorbar to an already existing graph
    void add_yerr(QCPGraph * graph, vector<double> x_std, vector<double> y_std, vector<double> s_std);
    
    // sets a specific format to a selected graph
    void setGraphStyle(QCPGraph * graph, QFormat format);
    
    QCPAxisRect * getSelectedAxisRect(QPoint &pos, QCPAxis::AxisType &axisType);
    
    // variables for zoom-box
    QRubberBand * _rubber_band;
    QPoint _click_origin;
    
    // increasing for each boxplot, used as key-element
    int _boxplot_cnt;
    // increasing for each rowplot
    int _subplot_cnt;
    // increasing color cnt
    int _color_cnt;
    
    // vector containing information about each plotted point
    // related to the x-position of a graph(QCPGraph *) in a subplot(AxisRect)
    map< QCPAxisRect* , map< QCPGraph *, map<double, string> > > _tooltips;
    
    // PROJECTION PART (all needed by plot.projection(.....))
    void _get_projected_xy(double &x, double &y, protype type, double lambda, double phi, double threshold);
    void _plot_projected_grid(QCPProjection projection);
    void _plot_coastline(protype type, QFormat format);
    bool _projected; // projection grid plotted
    QCPProjection _ref_pro; // saves latest projection configuration
    
    // map containing for each graph the original indices in ncfile
    map< QCPGraph *, vector<int> > _index;
    
};

#endif	/* PLOT_H */

