#include <./qwidget.h>

#include "lp_sked_station.h"

StationDialog::StationDialog(ivg::Session* session, lps::TemporalGrid* tg, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::StationDialog)
{
    setMouseTracking(true);
    
    ui->setupUi(this);
    ui->graphicsView->setScene(&scene);
    ui->graphicsView->setRenderHint(QPainter::Antialiasing);
    
    tree = new QGraphicsRectItem();
    scene.addItem(tree);
    tree->setVisible(false);
    
    map = new QGraphicsRectItem();
    scene.addItem(map);
    map->setVisible(true);
    
    transits = new QGraphicsRectItem();
    scene.addItem(transits);
    transits->setVisible(true);
    
    observations = new QGraphicsRectItem();
    scene.addItem(observations);
    observations->setVisible(true);
    
    elevationMask = new QGraphicsRectItem();
    scene.addItem(elevationMask);
    elevationMask->setVisible(false);
    
    sun = new QGraphicsRectItem();
    scene.addItem(sun);
    sun->setVisible(false);
    
    pole = new QGraphicsRectItem();
    scene.addItem(pole);
    pole->setVisible(true);
    
    obsConnection = new QGraphicsRectItem();
    scene.addItem(obsConnection);
    obsConnection->setVisible(false);
   
    
    setupStarMap();
   
    
    this->tg = tg;
    n_temporal_grids = tg->get_number_of_intervals();

       
    connect(ui->pushButton, SIGNAL(clicked()), this, SLOT(save_pdf()));
    connect(ui->pushButton_2, SIGNAL(clicked()), this, SLOT(refresh_map()));
    connect(ui->pushButton_3, SIGNAL(clicked()), this, SLOT(setCellVisibility()));
    
    ui->horizontalSlider->setMaximum(n_temporal_grids-1);
    
    connect(ui->horizontalSlider, SIGNAL(valueChanged(int)),
            ui->lcdNumber, SLOT( display(int) ));
        
    connect(ui->horizontalSlider, SIGNAL(valueChanged(int)),
            this, SLOT( setObservationVisibility() ));
    connect(ui->horizontalSlider, SIGNAL(valueChanged(int)),
            this, SLOT( setCellVisibility() ));
    connect(ui->horizontalSlider, SIGNAL(valueChanged(int)),
            this, SLOT( setTransitVisibility() ));
    connect(ui->horizontalSlider, SIGNAL(valueChanged(int)),
            this, SLOT( setTraverseVisibility() ));
    
    connect(ui->horizontalSlider, SIGNAL(valueChanged(int)),
            this, SLOT(  updateIntervalLabelText(int) ));
    
 
    this->session = session;
    
    tmpObsCon.resize(n_temporal_grids);
    
    
//    this->setStyleSheet("QGroupBox { "
//                        "border: 1px solid gray; "
//                        "border-radius: 5px;" 
//                        "}");
}



StationDialog::~StationDialog() 
{

    std::cerr << ">>>>>>>>>>>>>> ~StationDialog( ) " << std::endl; 

    clear();
 
    delete tree;
    delete map;
    delete transits;
    delete observations;
    delete elevationMask;
    delete sun;
    delete pole;
    delete obsConnection;

    delete ui;        
}

double scale =2.0;

void StationDialog::refresh_map(){
    for(QGraphicsItem * item : map->childItems()){
            scene.removeItem(item);
    }
    setupStarMap();
}

void StationDialog::setupStarMap(){
    QTransform transform;
    transform.scale(scale,scale);
    map->setZValue(-4);
    
    //elevation
    if( ui->spinBox->value() > 0){
        double inc = 180.0/ui->spinBox->value();
        for(double el=0; el <= 180.0; el+= inc ){
           QGraphicsEllipseItem * item = scene.addEllipse(-el/2*scale,-el/2*scale,el*scale,el*scale, QPen(Qt::black,.5));
           item->setZValue(-4);
           item->setParentItem(map);
        }
    }

    // azimuth
    if( ui->spinBox_2->value() > 0){
        double inc = 360.0/ui->spinBox_2->value();
        for(double az=0; az < 359.999; az+= inc ){
            QGraphicsItem * item = scene.addLine(0,0,std::sin(az*ivg::d2rad)*90*scale,-std::cos(az*ivg::d2rad)*90*scale,QPen(Qt::black,.5));
            item->setZValue(-4);
            item->setParentItem(map);
         
            QGraphicsSimpleTextItem* l1 = addLabel(az, 0, QString::number(round(az)), floor(az/90)  );
            l1->setParentItem(map);
            
        }
    }
}

void StationDialog::setCellVisibility(){
    for(auto& a : tree_info){
        bool inRange = std::get<0>(a.second) >= ui->spinBox_3->value() && std::get<0>(a.second) <= ui->spinBox_4->value() && ui->horizontalSlider->value() == std::get<1>(a.second);
            
        if(ui->checkBox_8->isChecked()){
            a.first->setVisible( std::get<2>(a.second) &&  inRange );
        } else {
            a.first->setVisible( inRange );
        }
        
    }
}

void StationDialog::setObservationVisibility(){
    int t = ui->horizontalSlider->value();
    
    std::map< unsigned, std::vector< ObservationItem*>> twin_info;
    
    for(auto& a : obs_info){
        if(ui->groupBox_3->isChecked()){
            bool vis = std::find(a.second.begin(), a.second.end(), ui->horizontalSlider->value()) != a.second.end();
            a.first->setVisible( vis );
            if( vis && ui->checkBox_2->isChecked() ){
                unsigned id = ceil(a.first->start()/10.0);
                twin_info[id].push_back(a.first);
            }
        } else {
            a.first->setVisible( true );
            if(  ui->checkBox_2->isChecked() ){
                unsigned id = ceil(a.first->start()/10.0);
                twin_info[id].push_back(a.first);
            }
        }
    }
    
    if( ui->checkBox_2->isChecked() ){
        int i = 0;
        for(auto& a : twin_info){
            for(ObservationItem* item: a.second){
                item->setPen(QPen(QColor(ivg::get_color_value(i).c_str()), 4));
            }
            ++i;
        }
    } else {
        for(auto& a : obs_info){
            if(a.first->thisStation()){
                a.first->setPen(QPen(Qt::red, 2));
            } else {
                a.first->setPen(QPen(Qt::magenta, 2));
            }
        }
    }
    
}

void StationDialog::setTransitVisibility(){
    int t = ui->horizontalSlider->value();
    for(auto& a : transit_info){
        if(ui->groupBox_3->isChecked()){
            a.first->setVisible( a.second == t );
        } else {
            a.first->setVisible( true );
        }
    }
}

void StationDialog::setTraverseVisibility(){
    int t = ui->horizontalSlider->value();
    for(auto& a : traverse_info){
        if(ui->groupBox_3->isChecked()){
            a.first->setVisible( a.second == t );
        } else {
            a.first->setVisible( true );
        }
    }
}

void StationDialog::updateIntervalLabelText(int t)
{
    std::pair<lps::Seconds, lps::Seconds> intervalBorder =  tg->getIntervalBorder( t );
    ui->label_5->setText(QString::number(intervalBorder.first) + " - " + QString::number(intervalBorder.second));
}


QGraphicsSimpleTextItem* StationDialog::addLabel(double azimuth, double elevation,  QString text, int alignment){
    double x = std::sin(azimuth*ivg::d2rad)*(90-elevation)*scale;
    double y = -std::cos(azimuth*ivg::d2rad)*(90-elevation)*scale;
    QGraphicsSimpleTextItem* item = scene.addSimpleText(text);
    
    item->setFont(QFont(QApplication::font().family(), 13));
    QRectF br = item->sceneBoundingRect();
    switch (alignment){
        case 0:{
            item->setPos(  x, y- br.height() );
            break;
        }
        case 1:{
            item->setPos(  x, y );
            break;
        }
        case 2:{
            item->setPos(  x- br.width(), y );
            break;
        }
        case 3:{
            item->setPos(  x- br.width(), y- br.height() );
            break;
        }
            
    }

    return item;
}

void StationDialog::drawPoint(double azimuth, double elevation)
{
    double x = std::sin(azimuth*ivg::d2rad)*(90-elevation)*scale;
    double y = -std::cos(azimuth*ivg::d2rad)*(90-elevation)*scale;
    QGraphicsItem * item = scene.addEllipse(x-4, y-4, 8, 8, QPen(Qt::black,2), QBrush(Qt::black));
    item->setParentItem(pole);
}

void StationDialog::drawPointCart(double x, double y)
{
    scene.addEllipse(x*scale-2, y*scale-2, 4, 4, QPen(Qt::green,1));
}

QGraphicsRectItem* StationDialog::drawRect(const lps::Rect & rect){
    QTransform transform;
    transform.scale(scale,scale);
    return scene.addRect(transform.mapRect(rect));
}

QGraphicsPathItem *StationDialog::addCell(const lps::Path& rect)
{
    QTransform transform;
    transform.scale(scale,scale);
    QGraphicsPathItem * item= new QGraphicsPathItem(transform.map(rect),tree);
    removable.push_back(item);
    return item;
}

QGraphicsPathItem *StationDialog::addCell(const lps::Path& rect, int level, int temporal_idx, bool valid)
{
    QGraphicsPathItem* item = addCell(rect);
    tree_info[item] = std::make_tuple(level, temporal_idx, valid);
    if( ui->groupBox_3->isChecked()  )
        item->setVisible(temporal_idx == ui->horizontalSlider->value());
    
    return item;
}

ObservationItem *StationDialog::addObservation(ivg::Source& source, double azimuth, double elevation, lps::Seconds start, std::vector<unsigned> temporal_idx, bool thisStation)
{
    double x = std::sin(azimuth*ivg::d2rad)*(90-elevation)*scale;
    double y = -std::cos(azimuth*ivg::d2rad)*(90-elevation)*scale;
    
    QString scr_name = source.get_name(ivg::srcname::ivs).c_str();
    ObservationItem* item = new ObservationItem(source.get_idx(), start, scr_name, lps::Point(x,y), thisStation );
    scene.addItem(item);

    removable.push_back(item);
    item->setParentItem(observations);
            
    obs_info[item] = temporal_idx;
    
    if( ui->groupBox_3->isChecked()  )
        item->setVisible( std::find(temporal_idx.begin(), temporal_idx.end(), ui->horizontalSlider->value()) != temporal_idx.end() );
    
    if(thisStation){
        for(unsigned& tidx : temporal_idx){
            if( tmpObsCon[tidx].elementCount() == 0 ){
                tmpObsCon[tidx].moveTo(x, y);
            } else {
                tmpObsCon[tidx].lineTo(x, y);
            }
        }
    }
    
    
    return item;
}

lps::Path transformPath(const std::vector<lps::Position> & path){
    lps::Path result;
    for(size_t i=0; i < path.size(); ++i){
        const lps::Position & p = path[i];
        double x = std::sin(p.azimuth()*ivg::d2rad)*(90-p.elevation())*scale;
        double y = -std::cos(p.azimuth()*ivg::d2rad)*(90-p.elevation())*scale;
        if(i==0){
            result.moveTo(x,y);
        }else{
            result.lineTo(x,y);
        }
    }
    return result;
}

void StationDialog::addElevationMask(const std::vector<lps::Position> &path, QPen pen){
    lps::Path p = transformPath(path);
    QGraphicsPathItem* item= scene.addPath(p,pen);
    item->setZValue(-2);
    item->setParentItem(elevationMask);

}

void StationDialog::addTraverse(QPen pen){
    for(int i = 0; i < tmpObsCon.size(); ++i){
        QGraphicsPathItem* item= scene.addPath(tmpObsCon[i],pen);
        item->setZValue(-2); 
        item->setParentItem(obsConnection);
        traverse_info[item] = i;
        removable.push_back(item);
    }
}

std::vector<QGraphicsPathItem*> StationDialog::addTransits(int station, const lps::Transits &transits, QPen pen)
{
    std::vector<QGraphicsPathItem*> result;
    
    lps::Seconds temporal_grid_resolution = ((int)(*session->get_setup())["SKED"]["temporal_resolution"]) * 60;

    for(const lps::Transit & transit : transits.get()){
        if(transit.get_sta_idx() == station){
            for(int t = 0; t < n_temporal_grids;  ++t){
                ivg::Date start = session->getStart();
                ivg::Date end = session->getStart();
                
                std::pair<lps::Seconds, lps::Seconds> intervalBorder =  tg->getIntervalBorder( t );
   
                start.add_secs(intervalBorder.first);
                end.add_secs(intervalBorder.second);
                
                std::vector<lps::Position> tmp = transit.get_path_part(start, end);
                if(tmp.size() > 0){
                    QGraphicsPathItem* item = addTransit(tmp,pen,transit.get_sou_name());
                    result.push_back(item);
                    transit_info[item] = t;
                }
            }
        }
    }
    return result;
}

QGraphicsPathItem* StationDialog::addTransit(const std::vector<lps::Position> & path, QPen pen, std::string souName){
    lps::Path p = transformPath(path);
    
    TransitItem* item =  new TransitItem(p, pen, souName);
    scene.addItem(item);
    
    item->setZValue(-2);
    item->setParentItem(transits);
        
    return item;

}

QGraphicsPathItem* StationDialog::addSunTransit(const std::vector<lps::Position> & path, QPen pen){
    lps::Path p = transformPath(path);
    QGraphicsPathItem* item= scene.addPath(p,pen);
    item->setZValue(-2);
    item->setParentItem(sun);

    return item;

}


void StationDialog::on_groupBox_toggled(bool checked)
{
    setCellVisibility();
    tree->setVisible(checked);
}

void StationDialog::on_groupBox_2_toggled(bool checked)
{
    map->setVisible(checked);
}

void StationDialog::on_groupBox_3_toggled(bool checked)
{
    this->setObservationVisibility();
    this->setTransitVisibility();
    this->setCellVisibility();
}

void StationDialog::on_checkBox_toggled(bool checked){
    setTraverseVisibility();
    obsConnection->setVisible(checked);
}

void StationDialog::on_checkBox_2_toggled(bool checked){
    setObservationVisibility();
}

void StationDialog::on_checkBox_3_toggled(bool checked)
{
    transits->setVisible(checked);
}

void StationDialog::on_checkBox_4_toggled(bool checked)
{
    setObservationVisibility();
    observations->setVisible(checked);
    ui->checkBox_2->setCheckable(checked);
    ui->checkBox_2->setEnabled(checked);
}

void StationDialog::on_checkBox_5_toggled(bool checked)
{
    elevationMask->setVisible(checked);
}

void StationDialog::on_checkBox_6_toggled(bool checked)
{
    sun->setVisible(checked);
}

void StationDialog::on_checkBox_7_toggled(bool checked)
{
    pole->setVisible(checked);
}

void StationDialog::print_pdf(std::string filePath){
    
    log<DETAIL>("*** printing ") % filePath;
    
    const QString pdfCreator="";
    const QString pdfTitle=this->windowTitle();
    
    QPrinter printer(QPrinter::ScreenResolution);
    printer.setOutputFileName( QString::fromStdString(filePath) );
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setFullPage(true);
    printer.setColorMode(QPrinter::Color);
    printer.printEngine()->setProperty(QPrintEngine::PPK_Creator, pdfCreator);
    printer.printEngine()->setProperty(QPrintEngine::PPK_DocumentName, pdfTitle);

    lps::Rect boundingRect = lps::Rect( lps::Point(-90,-90), lps::Point(0,0));
    printer.setPaperSize(boundingRect.size(), QPrinter::DevicePixel);
    
    QPainter painter(&printer);
    painter.setRenderHint(QPainter::Antialiasing);
    scene.render(&painter);  
}