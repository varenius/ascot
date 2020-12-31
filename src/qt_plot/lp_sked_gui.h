#ifndef GUI_H
#define GUI_H

#include <QObject>
#include <QTime>
#include <QGraphicsRectItem>
#include <QColor>
#include <QApplication>
#include <QThread>

#include "lp_sked_worker.h"
#include "lp_sked_station.h"

#include "session.h"
#include "milp.h"
#include "grid.h"



namespace lps{

class GUI : public QObject
{
    Q_OBJECT
    QThread workerThread;
    
public:
    explicit GUI(ivg::Session* session, QObject* parent = nullptr);

    ~GUI(){
        //std::cerr << "~GUI" <<  " " << stationViews.size() << std::endl;
        
        workerThread.quit();

//        for(int i = 0; i < stationViews.size(); ++i){
//            delete stationViews[i];
//        }
            
        //delete stationViews[1];
        
        for( StationDialog* sd : stationViews ){
            //delete sd;
            
      
            sd->tree_info.clear();
            sd->obs_info.clear();
            sd->transit_info.clear();
            sd->traverse_info.clear();
            sd->tmpObsCon.clear();
            
            for(QGraphicsItem* item : sd->scene.items()){
                sd->scene.removeItem(item);
            }

            delete sd->tree;
            delete sd->map;
            delete sd->transits;
            delete sd->observations;
            delete sd->elevationMask;
            delete sd->sun;
            delete sd->pole;
            delete sd->obsConnection;

            delete sd->ui;
                        
        }
        
    }
    
    int run(int argc, char *argv[]);

    void initGUI();
    
    void showSession();
    
    std::vector<StationDialog*> stationViews;
    
    std::map<unsigned, std::vector<ivg::Analysis_station*> > twinMap;
       
    ivg::Session* session;

    lps::TemporalGrid tg;
signals:

public slots:
    void foundSolution(std::vector<lps::StationActivity> activity);
    void selectedCell(int level, int temporal_idx, lps::Path rect, int sta_idx, bool valid);
    void print_skyplots();
};

} //namespace

#endif // GUI_H
