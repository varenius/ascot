#ifndef STATISTICS_H
#define STATISTICS_H

#include <QMainWindow>
#include <QTimer>
#include <QDialog>
#include <QInputDialog>
#include <QTabWidget>
#include <QFontDialog>
#include <QColorDialog>
#include <QMessageBox>
#include <QToolTip>
#include "qcustomplot.h"
#include "matrix.h"
#include "plot.h"
#include "tictoc.h"
#include "ivg_const.h"
#include "date.h"
#include "session.h"

#include <QStandardItemModel>
#include "ui_statistics.h"

struct AxisLayout
{
    string name;
    int data_idx;
    string unit;
    double fak;
    QCPAxis::LabelType label_type;
    string number_format;
    string time_format;
    
};

class Statistics : public QMainWindow
{
  Q_OBJECT
  
public:
  
  Statistics(ivg::Session *session_ptr, QWidget *parent = 0);
  explicit Statistics(QWidget *parent = 0);
  
  ~Statistics();
  
public slots:
    
  //Customized SLOTS
  void createFigure();
  void correctCBreak();
  void shiftResidUp(){ _shiftResid('U');};
  void shiftResidDown(){_shiftResid('D');};
  void showPointToolTip(QMouseEvent *event);
  void showBaselineInfo(QCPAbstractPlottable *plottable, QMouseEvent *event);
  void plotSiteSkyplot(QString &tooltip);
  void saveResid();
    
private:
    
  void _plot_sites();
  void _plot_sources();
  void _plot_baselines();
  void _init_residuals(ivg::Matrix &data, map< string,vector<int> > idx, residtype type);
  void _init_axis_layout();
  
  void _shiftResid(char direction);
    
  ivg::Session *_session;
  Ui::Statistics *ui;
  Plot *_last_dialog;
  map< delaytype , vector<Residual> > _residuals;
  vector<AxisLayout> _axis_layouts;
  map<string,int> _axis_assignment;
  QCPProjection _trf_projection,_crf_projection, _sky_projection;
  
  bool _bl_selected;
  
  // both vectors are needed for manual clock break correction
  // the i-th entry in _clock_break_epoch_mjd must correspond to the i-th 
  // entry in _clock_break_station
  std::vector<double> _clock_break_epoch_mjd;
  std::vector<std::string> _clock_break_station;
      
};

#endif // STATISTICS
