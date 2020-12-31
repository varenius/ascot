#ifndef STATIONDIALOG_H
#define STATIONDIALOG_H

#include <QGraphicsSimpleTextItem>
#include <iostream>
#include <QDialog>
#include <QGraphicsScene>
#include <QGraphicsEllipseItem>

#include <QPainter>
#include <QPrinter>
#include <QPrintEngine>

#include <tuple>

#include "transits.h"
#include "grid.h"
#include "geometry.h"
#include "ui_lp_sked_station.h"
#include "ivg_const.h"
#include "source.h"


namespace Ui {
class StationDialog;
}


class ObservationItem : public QGraphicsEllipseItem{

    private:
        int quasar_;
        QString name_;
        lps::Seconds start_; // in seconds after session start
        bool thisStation_;

    public:
        std::vector<ObservationItem*> copies;

       ObservationItem(int quasar, lps::Seconds start, QString name, const lps::Point & location, bool thisStation) :
           QGraphicsEllipseItem(location.x()-3,location.y()-3,6,6), quasar_(quasar), start_(start), name_(name), thisStation_(thisStation){
            this->setPen(QPen(Qt::red,2));
            this->setBrush(QBrush(Qt::white));
            setFlag(QGraphicsItem::ItemIsSelectable);
       }
       
       lps::Seconds start()const {return start_;}
       int quasar()const {return quasar_;}
       QString name() const {return name_;}
       bool thisStation() const{return thisStation_;};

       QVariant itemChange(GraphicsItemChange change, const QVariant &value)
       {
           if (change == QGraphicsItem::ItemSelectedChange)
           {
               if (value == true)
               {
                   std::cout << "selected: " << name_.toStdString() << " " <<  quasar_ << " observations: "<< copies.size()<< std::endl;
                   for(ObservationItem * item : copies){
                    item->setPen(QPen(Qt::blue,5));
                   }
               }
               else
               {
                   for(ObservationItem * item : copies){
                    item->setPen(QPen(Qt::red,2));
                   }
               }
           }

           return QGraphicsItem::itemChange(change, value);
       }
};

class TransitItem : public QGraphicsPathItem{

    private:
        QString name_;

    public:
        std::vector<TransitItem*> copies;

        TransitItem(const lps::Path & path, QPen pen, std::string souName):
        QGraphicsPathItem(path)
        {
            this->setPen(pen);
            this->setZValue(-2);
            this->setToolTip(QString::fromStdString(souName));
            setFlag(QGraphicsItem::ItemIsSelectable);
            this->name_ =QString::fromStdString(souName);
        }
        
       QString name() const {return name_;}

       QVariant itemChange(GraphicsItemChange change, const QVariant &value)
       {
           if (change == QGraphicsPathItem::ItemSelectedChange)
           {
               if (value == true)
               {
                   std::cout << "selected transit of source: " << name_.toStdString() << std::endl;
                   for(TransitItem * item : copies){
                    item->setPen(QPen(Qt::blue,5));
                   }
               }
               else
               {
                   for(TransitItem * item : copies){
                    item->setPen(QPen(Qt::red,2));
                   }
               }
           }

           return QGraphicsPathItem::itemChange(change, value);
       }
};

class StationDialog : public QDialog
{
    Q_OBJECT

public:
    explicit StationDialog(ivg::Session* session, lps::TemporalGrid* tg, QWidget *parent = nullptr);
    virtual ~StationDialog();

    QGraphicsScene scene;

    void setupStarMap();
    
    void show_tree(){
        ui->groupBox->setChecked(true);
        ui->groupBox_2->setChecked(false);
    }
    
    void addElevationMask(const std::vector<lps::Position> &path, QPen pen);
    void addTraverse(QPen pen);
    
    QGraphicsPathItem* addTransit(const std::vector<lps::Position> &path, QPen pen,  std::string souName);
    
    QGraphicsPathItem* addSunTransit(const std::vector<lps::Position> & path, QPen pen);

    std::vector<QGraphicsPathItem *> addTransits(int station, const lps::Transits& transits, QPen pen);

    QGraphicsSimpleTextItem* addLabel(double azimuth, double elevation, QString text, int alignment = 0);
    void drawPoint(double azimuth, double elevation);
    void drawPointCart(double x, double y);
   
    QGraphicsPathItem* addCell(const lps::Path &rect);
    
    QGraphicsPathItem* addCell(const lps::Path& rect, int level, int temporal_idx, bool valid);
    
    ObservationItem *addObservation(ivg::Source& source, double azimuth, double elevation, lps::Seconds start, std::vector<unsigned> temporal_idx, bool thisStation = true);


    QGraphicsRectItem *drawRect(const lps::Rect &rect);

    void clear(){
        
        for(QGraphicsItem * item : removable){
            scene.removeItem(item);
        }
        removable.clear();
        obs_info.clear();
        traverse_info.clear();
        tree_info.clear();

        tmpObsCon.clear();
        tmpObsCon.resize(n_temporal_grids);
        
    }
    
    void print_pdf(std::string filePath);
    
    void set_default_path(std::string path){   ui->lineEdit->setText(QString::fromStdString(path)); };

private slots:
    void on_groupBox_toggled(bool checked);
    void on_groupBox_2_toggled(bool checked);
    void on_groupBox_3_toggled(bool checked);

    void on_checkBox_toggled(bool checked);
    void on_checkBox_2_toggled(bool checked);
    void on_checkBox_3_toggled(bool checked);
    
    void on_checkBox_4_toggled(bool checked);
    
    void on_checkBox_5_toggled(bool checked);
    
    void on_checkBox_6_toggled(bool checked);
    
    void on_checkBox_7_toggled(bool checked);
    
    void save_pdf(){print_pdf(ui->lineEdit->text().toStdString());};
    
    void refresh_map();
    
    void setCellVisibility();
    void setObservationVisibility();
    void setTransitVisibility();
    void setTraverseVisibility();
    void updateIntervalLabelText(int t);

public:    
    std::vector<QGraphicsItem*> removable;
    Ui::StationDialog *ui;
    QGraphicsItem* tree;
    QGraphicsItem* map;
    QGraphicsItem* transits;
    QGraphicsItem* observations;
    QGraphicsItem* elevationMask;
    QGraphicsItem* sun;
    QGraphicsItem* pole;
    QGraphicsItem* obsConnection;
       
    ivg::Session* session;
    int n_temporal_grids;
    
    const lps::TemporalGrid* tg;
    
    // tree level, temporal idx, valid
    std::map<QGraphicsPathItem*, std::tuple<int,int,bool>> tree_info;
    
    std::map<ObservationItem*, std::vector<unsigned> > obs_info;
    std::map<QGraphicsPathItem*, int> transit_info;
    std::map<QGraphicsPathItem*, int> traverse_info;
    
    std::vector<lps::Path> tmpObsCon;

};

#endif // STATIONDIALOG_H
