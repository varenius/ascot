#ifndef WORKER_H
#define WORKER_H

#include <QObject>
#include <gurobi_c++.h>
#include <QThread>
#include <fstream>

#include "session.h"
#include "milp.h"
#include "geometry.h"
#include "grid.h"

namespace lps{

class Worker : public QObject
{
    Q_OBJECT
public:
    explicit Worker(QObject *parent = 0);

    Worker(){}
    Worker(ivg::Session* session):session(session){};

    ivg::Session* session;

signals:
    void foundSolution(std::vector<lps::StationActivity> activity);
    void selectedCell(int level, int temporal_idx, lps::Path rect, int station, bool valid);
public slots:
   void process();
};

} //namespace
#endif // WORKER_H


