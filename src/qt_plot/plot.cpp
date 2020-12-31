#include <QtCore/qdatetime.h>

#include "plot.h"
#include "movfise.h"


// ...........................................................................
Plot::Plot(QWidget *parent) : QWidget(parent)
// ...........................................................................
{

    // setup design of plot-widget
    _ui.setupUi(this);

    // used for the zoombox (middle mouse button)
    _rubber_band = new QRubberBand(QRubberBand::Rectangle, this);

    // to increase the key for each new boxplot
    _boxplot_cnt = 0;

    // in order to increase for each used color
    _color_cnt = 0;

    // set projected to false
    _projected = false;

    // get pointer to QCustomPlot for further operations
    _qcp = _ui.customPlot;
    _qcp->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom));

    //Set labels and title
    _qcp->plotLayout()->insertRow(0);
    _qcp->plotLayout()->addElement(0, 0, new QCPPlotTitle(_qcp, " "));
    _qcp->xAxis->setLabel("x");
    _qcp->yAxis->setLabel("y");

    //Set legend and show it
    _qcp->legend->setVisible(true);
    QFont legendFont = font();
    legendFont.setPointSize(10);
    _qcp->legend->setFont(legendFont);
    _qcp->legend->setSelectedFont(legendFont);
    _qcp->legend->setSelectableParts(QCPLegend::spItems); // legend box shall not be selectable, only legend items

    // set some more important settings
    _qcp->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP:: iSelectAxes | QCP::iSelectLegend | QCP::iSelectPlottables | QCP::iSelectItems);
    _qcp->setInteraction(QCP::iMultiSelect, true);
    _qcp->setContextMenuPolicy(Qt::CustomContextMenu);

    //connect slot that ties axis selections
    connect(_qcp, SIGNAL(selectionChangedByUser()), this, SLOT(selectionChanged()));
    //setup context menu popup
    connect(_qcp, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequest(QPoint)));
    //setup scroll zoom actions
    connect(_qcp, SIGNAL(mouseWheel(QWheelEvent*)), this, SLOT(mouseWheelTurned(QWheelEvent*)));
    //setup zoombox action
    connect(_qcp, SIGNAL(mousePress(QMouseEvent*)), this, SLOT(zoomBoxRequest(QMouseEvent*)));
    connect(_qcp, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMoveEvent(QMouseEvent*)));
    connect(_qcp, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseReleaseEvent(QMouseEvent*)));
    //setup doubleclick event to see all data
    connect(_qcp, SIGNAL(mouseDoubleClick(QMouseEvent*)), this, SLOT(rescalePlot(QMouseEvent*)));

    //setup label/legend/title edit by double click
    connect(_qcp, SIGNAL(titleDoubleClick(QMouseEvent*,QCPPlotTitle*)), this, SLOT(titleDoubleClick(QMouseEvent*,QCPPlotTitle*)));
    connect(_qcp, SIGNAL(axisDoubleClick(QCPAxis*,QCPAxis::SelectablePart,QMouseEvent*)), this, SLOT(axisLabelDoubleClick(QCPAxis*,QCPAxis::SelectablePart)));
    connect(_qcp, SIGNAL(legendDoubleClick(QCPLegend*,QCPAbstractLegendItem*,QMouseEvent*)), this, SLOT(legendDoubleClick(QCPLegend*,QCPAbstractLegendItem*)));

    // add trash layers where old projected grids and coastlines are stored
    _qcp->addLayer("trash_grid",_qcp->layer("grid"),QCustomPlot::LayerInsertMode::limAbove);
    _qcp->layer("trash_grid")->setVisible(false);

    _qcp->addLayer("trash_coastline",_qcp->layer("trash_grid"),QCustomPlot::LayerInsertMode::limAbove);
    _qcp->layer("trash_coastline")->setVisible(false);

}
Plot::Plot(ivg::Matrix x_ivg, ivg::Matrix y_ivg, QFormat format, QWidget *parent) : Plot(parent)
// ...........................................................................
{
    plot_data(x_ivg, y_ivg, format);
}
// ...........................................................................
Plot::Plot(ivg::Matrix x_ivg, ivg::Matrix y_ivg, QWidget *parent) : Plot(parent)
// ...........................................................................
{

    QFormat tmp = { QColor(Qt::red), 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, " ", 1.0, 5.0 };
    plot_data(x_ivg, y_ivg, tmp);
}
// ...........................................................................
Plot::Plot(const Plot& orig) {
}
// ...........................................................................
Plot::~Plot() {
}
// ...........................................................................
double Plot::histogram(ivg::Matrix y_ivg, int num_bins, QFormat format, int subplot_cnt, bool addToLegend, bool inPercent)
// ...........................................................................
{
    // pop up new plot
    (*this).show();

    // generate new barplot for the histogramm
    QCPBars *histogram = new QCPBars(_qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom), _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atLeft));
    _qcp->addPlottable(histogram);

    // initialize histrogram based on min,max,num_bins
    gsl_histogram * h = gsl_histogram_alloc(num_bins);
    gsl_histogram_set_ranges_uniform (h, y_ivg.min(), y_ivg.max());

    // fill histrogram with data from y_ivg matrix, element by element
    for(auto &y: y_ivg)
        gsl_histogram_increment (h, y);
 
    // compute histrogram-plot data
    QVector<double> ticks, data;
    double lower,upper,amount_max=0;
    for(int i=0; i<num_bins; i++)
    {
        gsl_histogram_get_range (h,i,&lower,&upper);

        ticks << (double) (upper+lower)/2.0;
        double amount = gsl_histogram_get(h,i);

        if(amount>amount_max)
            amount_max=amount;

        data << amount;
    }

    gsl_histogram_free (h);
    
    if( inPercent ) {
        unsigned int n_data_points = y_ivg.length();
        for(double& d: data){
            d = d/n_data_points*100.0;
        }
    }
    
    // appearance settings
    histogram->setPen(QPen(format.color, format.line_width, format.pen_style));
    QColor tmp = format.color;
    tmp.setAlpha(50);
    histogram->setBrush(tmp);
    
    if(!format.legend.isEmpty())
        histogram->setName(format.legend);
    
    if(!addToLegend)
        histogram->removeFromLegend();
        
    // get a nice appearance
    histogram->setWidth(upper-lower);
    histogram->setData(ticks, data );
    histogram->rescaleAxes();

    _qcp->replot();

    // return bin-width for calculating gauss-curve
    return upper-lower;
}
// ...........................................................................
void Plot::boxplot(ivg::Matrix y_ivg, QFormat format, int subplot_cnt)
// ...........................................................................
{
    // increase boxplot counter for correct key (increasing with each boxplot)
    _boxplot_cnt++;

    // get amount of elements
    int n = y_ivg.length();

    // copy content from matrix to array
    double data[n];
    copy(y_ivg.begin() ,y_ivg.end(), data);

    // sort data in order to be able to get median
    gsl_sort(data, 1, n);

    // determine min, max, median, lower and upper quantile of dataset
    double min = y_ivg.min();
    double max = y_ivg.max();
    double median = gsl_stats_median_from_sorted_data (data,  1, n);
    double upperq = gsl_stats_quantile_from_sorted_data (data, 1, n, 0.75);
    double lowerq = gsl_stats_quantile_from_sorted_data (data, 1, n, 0.25);

    // pop up new plot
    (*this).show();

    // generate new boxplot
    QCPStatisticalBox *boxplot = new QCPStatisticalBox(_qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom), _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atLeft));
    _qcp->addPlottable(boxplot);

    // get nice appearance
    boxplot->setPen(QPen(format.color, format.line_width, format.pen_style));
    boxplot->setBrush(QBrush(format.color,Qt::Dense6Pattern));
    boxplot->setKey(_boxplot_cnt);
    boxplot->setMinimum(min);
    boxplot->setLowerQuartile(lowerq);
    boxplot->setMedian(median);
    boxplot->setUpperQuartile(upperq);
    boxplot->setMaximum(max);

    // remove corresponding legend item
    _qcp->legend->removeAt(_qcp->legend->elementCount()-1);
    _qcp->rescaleAxes();
    _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom)->setRange(-1,_boxplot_cnt+2);
    _qcp->replot();

}
// ...........................................................................
void Plot::projection(QCPProjection projection, ivg::Matrix lam, ivg::Matrix phi, QFormat format, string title, vector<string> tooltips)
// ...........................................................................
{
    // threshold for mollweide projection iteration loop
    #define THRESHOLD 1e-3
    #define TEXTSIZE 28

    if(projection.arrow_value == 0.0)
    {
        if(_ref_pro.plot_coast) // in case of TRF
            projection.arrow_value = 1.0; // 1[cm] default
        else // in case of CRF
            projection.arrow_value = 1000.0; // 1000[microas] default
    }
      
    // if we got tooltips for the positions, add click event
    if( tooltips.size() == lam.length() )
        connect(_qcp, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(showPointToolTip(QMouseEvent*)));

    if(phi.max() <= M_PI_2 && phi.min() >= -M_PI_2) // lam.max() <= M_PI && lam.min() >= -M_PI && 
    {
        // only plot grid and coastline if it's the first time
        if(_projected == false)
        {
            _ref_pro = projection;
            _plot_projected_grid(_ref_pro);
            // only plot coastline if desired
            if(_ref_pro.plot_coast)
                _plot_coastline(_ref_pro.type, { projection.coast_color, projection.coast_width, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, " " });

            _projected = true;
        }

        map<double, string> sorted_tooltips;
        // check correct size of the matrices
        if(lam.cols() <= 2 && phi.cols() <= 2 && lam.rows() == phi.rows())
        {
            double x,y,end_x,end_y;
            QVector<double> x_vec,y_vec;
            for(int i=0; i<lam.rows(); i++)
            {
                _get_projected_xy(x, y, _ref_pro.type, lam(i,0),  phi(i,0), 1e-3);
                x_vec.push_back(x);
                y_vec.push_back(y);

                // we need to know the position-relation between the data and the tooltips
                // for that reason we use the key-value x as identifier
                if( tooltips.size() == lam.length() )
                    sorted_tooltips[x] = tooltips.at(i);

                // if we got two columns for matrices lam and phi, use the second column for the arrow
                if(lam.cols() == 2 && phi.cols() == 2 )
                {
                    double dx = lam(i,1);
                    double dy = phi(i,1);

                    std::complex<double> mycomplex_dxdy (dx,dy);
                    double xy_dist = sqrt(std::norm(mycomplex_dxdy));

                    double arc = atan2(dx,dy);
                    double t = fmod(arc,(2*M_PI));
                    
                    cerr << "d_ra/d_phi [microas]: "<< lam(i,1)*ivg::rad2mas*1000 << "|" << phi(i,1)*ivg::rad2mas*1000  << " | dist [mas]: " << xy_dist*ivg::rad2mas*1000 << endl;
                    
                   // add the arrow:
                   QCPItemLine *arrow = new QCPItemLine(_qcp);
                   _qcp->addItem(arrow);
                   arrow->setPen(QPen(format.color, format.line_width, Qt::PenStyle::SolidLine));
                   arrow->start->setCoords(x,y);

                   double new_end_x = x + xy_dist * 1e08 *_ref_pro.arrow_scale * sin(t);
                   double new_end_y = y + xy_dist * 1e08 *_ref_pro.arrow_scale * cos(t);

                    //                   cerr << "x: " << x << "|" << new_end_x << endl;
                    //                   cerr << "y: " << y << "|" << new_end_y << endl;

                   arrow->end->setCoords(new_end_x,new_end_y);
                    arrow->setHead(QCPLineEnding(_ref_pro.arrow_end, _ref_pro.arrow_width, _ref_pro.arrow_length, false));

                    // last iteration we add a reference arrow
                   if(i == lam.rows()-1)
                   {
                        // plot the arrow
                        QCPItemLine *ref_arrow = new QCPItemLine(_qcp);
                        _qcp->addItem(ref_arrow);
                        ref_arrow->setPen(QPen(Qt::black, format.line_width, Qt::PenStyle::SolidLine));
                        ref_arrow->start->setCoords(-2.8,1.0);
                        
                        // depending on CRF or TRF the arrow is based on different units
                        double dy_ref;
                        if(_ref_pro.plot_coast)
                            dy_ref = (_ref_pro.arrow_value) * 1.56789571971343e-09; // from [cm] to [rad] on equator (40074km)
                        else
                            dy_ref = (_ref_pro.arrow_value/1000.0) * ivg::mas2rad; // from [microas] to [rad]
                        
                        ref_arrow->end->setCoords(-2.8, 1.0  + dy_ref * 1e08 *_ref_pro.arrow_scale);
                        ref_arrow->setHead(QCPLineEnding(_ref_pro.arrow_end, _ref_pro.arrow_width, _ref_pro.arrow_length, false));

                        // plot the text
                        QCPItemText *textLabel = new QCPItemText(_qcp);
                        _qcp->addItem(textLabel);
                        textLabel->setSelectable(false);
                        textLabel->setPositionAlignment(Qt::AlignLeft | Qt::AlignTop);
                        textLabel->position->setType(QCPItemPosition::ptPlotCoords);
                        textLabel->setRotation(-90.0);
                        textLabel->position->setCoords(-2.7, 1.0); // place position at center/top of axis rect
                        stringstream ss;
                        ss << _ref_pro.arrow_value;
                        QString unit;
                        if(_ref_pro.plot_coast)
                            unit = QString::fromUtf8("[cm]");
                        else
                            unit = QString::fromUtf8("[\u00B5as]");
                        
                        QString text = QString::fromStdString(ss.str()) + unit;
                        textLabel->setText(text);
                        textLabel->setFont(QFont(font().family(), TEXTSIZE, QFont::Weight::Normal)); // make font a bit larger
                   }
                   
                   
                }
            }



            // all points in lam and phi are set up as ONE graph
            QCPGraph *graph = _qcp->addGraph(_qcp->xAxis, _qcp->yAxis);
            graph->setData(x_vec, y_vec);
            graph->setName(format.legend);
            setGraphStyle(graph, format);

            // save tooltips for each added graph
            _tooltips[_qcp->axisRect(0)][graph] = sorted_tooltips;

            // show the results
            set_title(title);
            (*this).show();
            _qcp->replot();
}
        else
            throw runtime_error( "void Plot::projection(protype type, ivg::Matrix lam, ivg::Matrix phi, QFormat format, string title): wrong matrix size of lambda and/or phi");
    }
    else
        throw runtime_error( "void Plot::projection(protype type, ivg::Matrix lam, ivg::Matrix phi, QFormat format, string title): lambda and/or phi in a not expected range");
}
// ...........................................................................
void Plot::_plot_coastline(protype type, QFormat format)
// ...........................................................................
{
    // create new layer where the coastline is plotted on
    _qcp->addLayer("pro_coastline",_qcp->layer("trash_coastline"), QCustomPlot::LayerInsertMode::limAbove ); 

    // hard implemented path to the coastlinedata (segments containing latitude and longitude)  
    stringstream path_ss;
    if(const char* env_ascot_dir = std::getenv("ASCOT_DIR"))
        path_ss << env_ascot_dir << "/src/qt_plot/world.gmt";

    string path = path_ss.str();

    ifstream inStream;
    string line;

    double x,y,t;
    QVector<double> x_vec,y_vec,t_vec;

    // reads the file and each segment is ploted as individual curve
    inStream.open(path.c_str(), ios::in);
    if( !inStream.is_open() )
            throw runtime_error( "void Plot::_plot_coastline(QFormat format) in plot.cpp: Failed to open file: "+path);
    else
    {
        while( getline( inStream, line ) )
        {
            // collect data as long it's the same segment
            if(line.substr(0,1) != ">")
            {
                ivg::Matrix lam_phi( 1, 2 );
                istringstream lam_phi_line( line );
                lam_phi_line >> lam_phi( 0 ) >> lam_phi( 1 );

                _get_projected_xy(x, y, type, lam_phi(0)*(M_PI/180),  lam_phi(1)*(M_PI/180), 1e-3);

                t_vec.push_back(t);
                x_vec.push_back(x);
                y_vec.push_back(y);
            }
            // plot the data as a new curve if segment is finished
            else
            {
                QCPCurve *newCurve = new QCPCurve(_qcp->xAxis, _qcp->yAxis);
                _qcp->addPlottable(newCurve);

                newCurve->setPen(QPen(format.color, format.line_width, format.pen_style));
                newCurve->setData(t_vec,x_vec,y_vec);
                newCurve->setSelectable(false);
                newCurve->setAntialiased(true);
                newCurve->setLayer("pro_coastline");
                _qcp->legend->removeAt(_qcp->legend->elementCount()-1);

                t++;

                // clear vectors to be prepared for the new segment
                t_vec.clear();
                x_vec.clear();
                y_vec.clear();
            }
        }
        _qcp->replot();
    }
}
// ...........................................................................
void Plot::_plot_projected_grid(QCPProjection projection)
// ...........................................................................
{
    #define TEXTSIZE 15

    // create new layer where the grid is plotted on
    _qcp->addLayer("pro_grid",_qcp->layer("trash_grid"), QCustomPlot::LayerInsertMode::limAbove ); 

    // pen definition for the grid (color, width, type)
    QPen pen(projection.grid_color, projection.grid_width, Qt::SolidLine);

    double phi_zero = 90.0;
    double phi_start = -90.0;
    double phi_inc = 1.0;
    double lam_inc = 360.0;
    // in case of polarstereo-plot the loops have to run slightly different
    if( projection.type == protype::polarstereo )
    {
        phi_zero = 0.0;
        phi_start = 0.0;
        phi_inc = 45.0;
        lam_inc = 5.0;
    }

    // plotting a grid based on predefined steps/spaces
    // plotting latitude-lines depending on selected projection
    int run=0;
    double x,y;
    for(double lambda= -180.0; lambda <= 180.0; lambda += projection.lam_space)
    {
        QVector<double> x_vec,y_vec,t_vec;
        for(double phi= phi_start; phi <= 90.0; phi += phi_inc)
        {
            _get_projected_xy(x, y, projection.type, lambda*(M_PI/180.0),  phi*(M_PI/180.0), 1e-3);

            t_vec.push_back(run);
            x_vec.push_back(x);
            y_vec.push_back(y);
        };

        QCPCurve *newCurve = new QCPCurve(_qcp->xAxis, _qcp->yAxis);
        _qcp->addPlottable(newCurve);

        newCurve->setData(t_vec,x_vec,y_vec);
        newCurve->setSelectable(false);
        newCurve->setPen(pen);
        newCurve->setAntialiased(true);

        newCurve->setLayer("pro_grid");
        _qcp->legend->removeAt(_qcp->legend->elementCount()-1);
        run++;
    }

    // plotting longitude-lines depending on selected projection
    for(double phi= phi_start; phi <= 90.0; phi += projection.phi_space)
    {
        QVector<double> x_vec,y_vec,t_vec;
        for(double lambda= -180.0 - projection.grid_extension; lambda <= 180.0 + projection.grid_extension; lambda += lam_inc + (2*projection.grid_extension))
        {
            _get_projected_xy(x, y,projection.type, lambda*(M_PI/180.0),  phi*(M_PI/180.0), 1e-3);

            t_vec.push_back(run);
            x_vec.push_back(x);
            y_vec.push_back(y);
        };

        QCPCurve *newCurve = new QCPCurve(_qcp->xAxis, _qcp->yAxis);
        _qcp->addPlottable(newCurve);

        newCurve->setData(t_vec,x_vec,y_vec);
        newCurve->setSelectable(false);
        newCurve->setPen(pen);
        newCurve->setAntialiased(true);

        newCurve->setLayer("pro_grid");
        _qcp->legend->removeAt(_qcp->legend->elementCount()-1);

        run++;
    }

    // set text (12h and 0degree) if desired
    if(projection.plot_text == true)
    {
        _get_projected_xy(x, y,projection.type, 0,  phi_zero*(M_PI/180.0), 1e-3);

        QCPItemText *ra_text = new QCPItemText(_qcp);
        _qcp->addItem(ra_text);
        ra_text->position->setType(QCPItemPosition::PositionType::ptPlotCoords);
        ra_text->position->setCoords(x,y);
        ra_text->setPositionAlignment(Qt::AlignBottom|Qt::AlignHCenter);
        ra_text->setText("12h");
        if( projection.type == protype::polarstereo )
            ra_text->setText(QString::fromUtf8("0\u00B0"));
        ra_text->setFont(QFont(font().family(), TEXTSIZE));
        ra_text->setLayer("pro_grid");


        _get_projected_xy(x, y,projection.type, -180*(M_PI/180.0), 0, 1e-3);

        QCPItemText *dec_text1 = new QCPItemText(_qcp);
        _qcp->addItem(dec_text1);
        dec_text1->position->setType(QCPItemPosition::PositionType::ptPlotCoords);
        dec_text1->position->setCoords(x,y);
        dec_text1->setPositionAlignment(Qt::AlignBottom|Qt::AlignRight);
        dec_text1->setText(QString::fromUtf8("0\u00B0"));
        if( projection.type == protype::polarstereo )
            dec_text1->setText("");
        dec_text1->setFont(QFont(font().family(), TEXTSIZE));
        dec_text1->setLayer("pro_grid");
    }

    // get nice appearance for a each projection type
    _qcp->rescaleAxes();
    if( projection.type == protype::polarstereo )
    {
        _qcp->xAxis->setRange(-1.8,1.8);
        _qcp->yAxis->setRange(-1.8,1.8);
    }
    else
    {
        _qcp->xAxis->setRange(-3.1,3.1);
        _qcp->yAxis->setRange(-1.7,1.9);
    }
    _qcp->xAxis->grid()->setVisible(false);
    _qcp->yAxis->grid()->setVisible(false);
    _qcp->xAxis->setVisible(false);
    _qcp->yAxis->setVisible(false);
    _qcp->legend->setVisible(false);

}
// ...........................................................................
void Plot::_get_projected_xy(double &x, double &y, protype type, double lambda, double phi, double threshold)
// ...........................................................................
{
    // depending on the selected projection, the transformation between lambda/phi and x/y is performed
    if(type == protype::mollweide)
    {
        // declaration
        double theta, d_theta_s, theta_s, diff;

        // first guess
        theta_s = 2.0*asin((2.0*phi)/M_PI);

        diff = 999;
        while(abs(diff) > threshold)
        {
            theta = 0.5*theta_s;

            d_theta_s = -((theta_s + sin(theta_s) - M_PI * sin(phi))/(1.0+cos(theta_s)));
            theta_s = theta_s + d_theta_s;

            diff = (2.0*theta + sin(2.0*theta)) -  (M_PI * sin(phi));
        }

        x = (2.0*sqrt(2.0)*lambda*cos(phi))/M_PI;
        y = sqrt(2.0) * sin(phi);
    }
    else if(type == protype::naturalearth)
    {
        double phi2 = phi * phi;
        double phi4 = phi2 * phi2;
        x = lambda * (0.8707 + phi2 * (-0.131979 + phi2 * (-0.013791 + phi4 * phi2 * (0.003971 + phi2 * -0.001529))));
        y = phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 + -0.005916 * phi4)));

    }
    else if(type == protype::mercator)
    {
        x = lambda;
        y = 0.5 * log((1+sin(phi))/(1-sin(phi)));
    }
    else if(type == protype::polarstereo)
    {      
        x = ((M_PI/2.0)-phi)*sin(lambda);
        y = ((M_PI/2.0)-phi)*cos(lambda);
    }
}
// ...........................................................................
void Plot::add_projected_colormap()
// ...........................................................................
{
    _qcp->addLayer("pro_colormap",_qcp->layer("trash_grid"), QCustomPlot::LayerInsertMode::limBelow ); 

    QCPColorMap *colorMap = new QCPColorMap(_qcp->xAxis, _qcp->yAxis);
    _qcp->addPlottable(colorMap);

    colorMap->data()->setSize(100, 100);
    colorMap->data()->setRange(QCPRange(-3.1, 3.1), QCPRange(-1.7, 1.7));

    //    for (int x=0; x<100; ++x)
    //        for (int y=0; y<100; ++y)
    //                colorMap->data()->setCell(x, y, 1.0);


    double z;
    for(double lambda= -180.0; lambda <= 180.0; lambda += 30.0)
    {
        double x,y;
        for(double phi= -90.0; phi <= 90.0; phi += 20.0)
        {
            _get_projected_xy(x, y, protype::naturalearth, lambda*(M_PI/180.0),  phi*(M_PI/180.0), 1e-3);
            z = drand48();

            //            gsl_sf_legendre_sphPlm

            colorMap->data()->setData(x, y, z);

        }
    }


    //    for (int x=0; x<500; ++x)
    //    {
    //        for (int y=0; y<500; ++y)
    //        {
    //            if(x == y)
    //                 colorMap->data()->setData(x, y, 2);
    //            else
    //                colorMap->data()->setCell(x, y, qCos(x/10.0)+qSin(y/10.0));
    //        }
    //    }

    // add a color scale:
    QCPColorScale *colorScale = new QCPColorScale( _qcp );
    _qcp->plotLayout()->addElement(1, 1, colorScale); // add it to the right of the main axis rect
    colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
    colorMap->setColorScale(colorScale); // associate the color map with the color scale
    colorScale->axis()->setLabel("matrix elements");

    colorMap->setGradient(QCPColorGradient::gpCold);
    colorMap->rescaleDataRange(true);
    colorMap->setInterpolate(true);

    colorMap->setLayer("pro_colormap");

    _qcp->rescaleAxes();
    _qcp->replot();

}
// ...........................................................................
QCPGraph* Plot::plot_data(ivg::Matrix x_ivg, ivg::Matrix y_ivg, int subplot_cnt, vector<string> tooltips)
// ...........................................................................
{
    // easiest way to plot data without defining any format in the function-call

    QFormat tmp = { QColor(ivg::color_values.at(_color_cnt).c_str()), 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, " ", 1.0 , 5.0 }; 
    _color_cnt++;

    return plot_data(x_ivg, y_ivg, tmp, subplot_cnt, tooltips);
}
// ...........................................................................
QCPGraph* Plot::plot_data(ivg::Matrix x_ivg, ivg::Matrix y_ivg, QFormat format, int subplot_cnt, vector<string> tooltips)
// ...........................................................................
{
    return plot_data(vector<double>(x_ivg.begin(),x_ivg.end()),vector<double> (y_ivg.begin(),y_ivg.end()),format, subplot_cnt, tooltips);
}
// ...........................................................................
QCPGraph* Plot::plot_data(vector<double> x_std, vector<double> y_std, QFormat format, int subplot_cnt, vector<string> tooltips)
// ...........................................................................
{
    // if we got tooltips, add them and the signal
    map<double, string> sorted_tooltips;
    if(tooltips.size() == x_std.size())
    {
        for(int i=0; i<tooltips.size(); i++)
            sorted_tooltips[x_std.at(i)] = tooltips.at(i);

    }
    // in case of no tooltips, create tooltip based on x/y position
    else if(tooltips.size() == 0)
    {
        for(int i=0; i<x_std.size(); i++)
        {
            stringstream ss;
            ss << "x: " << setprecision(10) << scientific <<  x_std.at(i) << endl;
            ss << "y: " << setprecision(10) << scientific << y_std.at(i);
            sorted_tooltips[x_std.at(i)] = ss.str();
        }
    }
    else
        cerr << "!!! No tooltips available because unequal length of tooltips and datapoints: " << tooltips.size() << " vs " << x_std.size() << endl;

    QCPGraph * added_graph = add_graph(x_std, y_std, _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom), _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atLeft));

    _tooltips[_qcp->axisRect(subplot_cnt)][added_graph] = sorted_tooltips;
    connect(_qcp, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(showPointToolTip(QMouseEvent*)));

    setGraphStyle(added_graph, format);
    added_graph->setErrorType(QCPGraph::ErrorType::etNone);
    _qcp->replot();
    
    return added_graph;
}
// ...........................................................................
void Plot::plot_yerr(ivg::Matrix x_ivg, ivg::Matrix y_ivg, ivg::Matrix s_ivg, QFormat format, int subplot_cnt)
// ...........................................................................
{
    plot_yerr(vector<double>(x_ivg.begin(),x_ivg.end()),
             vector<double> (y_ivg.begin(),y_ivg.end()),
             vector<double> (s_ivg.begin(),s_ivg.end()),format, subplot_cnt);
}
// ...........................................................................
void Plot::plot_yerr(vector<double> x_std, vector<double> y_std,
        vector<double> s_std, QFormat format, int subplot_cnt)
// ...........................................................................
{
    QCPGraph * added_graph = add_graph(x_std, y_std, _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom), _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atLeft));
    add_yerr(added_graph, x_std, y_std, s_std);
    setGraphStyle(added_graph, format);
    _qcp->replot();
}
// ...........................................................................
void Plot::add_yerr(QCPGraph * graph, vector<double> x_std, vector<double> y_std,
        vector<double> s_std)
// ...........................................................................
{
    (*this).show();

    QVector<double> x = QVector<double>::fromStdVector(x_std);
    QVector<double> y = QVector<double>::fromStdVector(y_std);
    QVector<double> s = QVector<double>::fromStdVector(s_std);

    graph->setErrorType(QCPGraph::ErrorType::etValue);
    graph->setDataValueError(x, y, s);
}
// ...........................................................................
void Plot::plot_mjd_series(ivg::Matrix mjd, ivg::Matrix y_ivg, QFormat format, int subplot_cnt, string time_format, vector<string> tooltips)
// ...........................................................................
{
    _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom)->setTickLabelType(QCPAxis::ltDateTime);
    _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom)->setDateTimeFormat(QString::fromStdString(time_format));

    vector<double> time_qcustomplot;
    for(int i=0; i<mjd.length(); i++)
        time_qcustomplot.push_back((mjd(i) - 40587.0 ) * 86400.0);

    plot_data(time_qcustomplot,vector<double> (y_ivg.begin(),y_ivg.end()),format, subplot_cnt, tooltips);

}
// ...........................................................................
void Plot::plot_mjd_series_std(ivg::Matrix mjd, ivg::Matrix y_ivg, ivg::Matrix y_std, QFormat format, int subplot_cnt, string time_format, vector<string> tooltips)
// ...........................................................................
{
    _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom)->setTickLabelType(QCPAxis::ltDateTime);
    _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom)->setDateTimeFormat(QString::fromStdString(time_format));

    vector<double> time_qcustomplot;
    for(int i=0; i<mjd.length(); i++)
        time_qcustomplot.push_back((mjd(i) - 40587.0 ) * 86400.0);

    plot_data(time_qcustomplot,vector<double> (y_ivg.begin(),y_ivg.end()),format, subplot_cnt, tooltips);

    QVector<double> x = QVector<double>::fromStdVector(time_qcustomplot);
    QVector<double> y = QVector<double>::fromStdVector(vector<double> (y_ivg.begin(),y_ivg.end()));
    QVector<double> s = QVector<double>::fromStdVector(vector<double> (y_std.begin(),y_std.end()));

    _qcp->graph()->setDataValueError(x, y, s);
    _qcp->graph()->setErrorType(QCPGraph::ErrorType::etValue);
    _qcp->graph()->rescaleAxes(true,false);
    _qcp->replot();
}
// ...........................................................................
void Plot::plot_mjd_series_hpd(ivg::Matrix mjd, ivg::Matrix y_ivg, ivg::Matrix hpdu, ivg::Matrix hpdl, QFormat format, int subplot_cnt, string time_format, vector<string> tooltips)
// ...........................................................................
{
    _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom)->setTickLabelType(QCPAxis::ltDateTime);
    _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom)->setDateTimeFormat(QString::fromStdString(time_format));
        
    vector<double> time_qcustomplot;
    for(int i=0; i<mjd.length(); i++)
        time_qcustomplot.push_back((mjd(i) - 40587.0 ) * 86400.0);
    
    plot_data(time_qcustomplot,vector<double> (y_ivg.begin(),y_ivg.end()),format, subplot_cnt, tooltips);
    
    QVector<double> x = QVector<double>::fromStdVector(time_qcustomplot);
    QVector<double> y = QVector<double>::fromStdVector(vector<double> (y_ivg.begin(),y_ivg.end()));
    QVector<double> su = QVector<double>::fromStdVector(vector<double> (hpdu.begin(),hpdu.end()));
    QVector<double> sl = QVector<double>::fromStdVector(vector<double> (hpdl.begin(),hpdl.end()));

    _qcp->graph()->setDataValueError(x, y, su, sl);
    _qcp->graph()->setErrorType(QCPGraph::ErrorType::etValue);
    _qcp->graph()->rescaleAxes(true,false);
    _qcp->replot();
}
// ...........................................................................
void Plot::highlight(string identify, QFormat format, bool remove_old, int subplot_cnt)
// ...........................................................................
{
    vector<double> new_x,new_y,new_s;
    map<double, string> new_sorted_tooltips;

    QList<QCPGraph*> graph_list;

    if(_qcp->selectedGraphs().size()>0)
        graph_list = _qcp->selectedGraphs();
    else
        graph_list = _qcp->axisRect(subplot_cnt)->graphs();


    for(auto *graph: graph_list)
    {
        map<double, string> sorted_tooltips = _tooltips[_qcp->axisRect(subplot_cnt)][graph];

        // in case of seperated tokens to identify
        vector<string> check_for_tokens = tokenize_string( identify, '\n' );
        for(auto &tt: sorted_tooltips)
        {
            // get strings looking for and looking in
            string check_in = tt.second;
            string check_for;

            bool highlight = false;
            for(int i=0; i<check_for_tokens.size();i++)
            {
                identify = check_for_tokens.at(i);

                // tokenize the tooltips by delimiter new-line \n            
                vector<string> tokens = tokenize_string( tt.second, '\n' );

                // in case of {1}, only check a specific row of the tooltips, seperated by \n
                if(identify.substr(0,1) == "{" && identify.substr(2,1) == "}")
                {
                    int row = stoi(identify.substr(1,1));

                    // requested row must fit to amount of tokens
                    if(tokens.size() >= row)
                    {
                        check_in = tokens.at(row-1);
                        check_for = identify.substr(3);

                        // if mulitple string e.g "{5}Wz$Al$65"
                        if(check_for.find("$") != string::npos)
                        {
                            vector<string> and_tokens = tokenize_string( check_for, '$' );
                            for(int t=0; t<and_tokens.size(); t++)
                            {
                                highlight = check_in.find(and_tokens.at(t)) != string::npos;
                                if(highlight == false)
                                    break;
                            }
                        }
                        else
                        {
                            // check if we need to highlight the data-point
                            highlight = check_in.find(check_for) != string::npos;
                        }
                    }
                }
                // in case of [>1], only check a specific row of the tooltips, seperated by \n
                else if(identify.substr(0,1) == "[" && identify.substr(3,1) == "]")
                {
                    int row = stoi(identify.substr(2,1));

                    // requested row must fit to amount of tokens
                    if(tokens.size() >= row)
                    {
                        double value = stod(identify.substr(4));
                        double double_token = stod(tokens.at(row-1));
                        // check if we need to highlight the data-point
                        if(identify.substr(1,1) == ">")
                            highlight = double_token > value;
                        else if(identify.substr(1,1) == "<")
                            highlight = double_token < value;
                        else if(identify.substr(1,1) == "=") 
                            highlight = double_token == value;
                    }
                }

                if(highlight == true)
                    break;
            }
            // in case of highlighting, store data-point
            if(highlight == true)
            {
                new_x.push_back(graph->data()->find(tt.first)->key);
                new_y.push_back(graph->data()->find(tt.first)->value);
                new_s.push_back(graph->data()->find(tt.first)->valueErrorMinus);
                new_sorted_tooltips[graph->data()->find(tt.first)->key] = tt.second;
                // depending on selection, remove old data-point or not
                if(remove_old)
                    graph->data()->remove(tt.first);
            }
        }
    }

    // add neq graph
    QVector<double> x = QVector<double>::fromStdVector(new_x);
    QVector<double> y = QVector<double>::fromStdVector(new_y);
    QVector<double> s = QVector<double>::fromStdVector(new_s);
    
    QCPGraph *new_graph = _qcp->addGraph(_qcp->axisRect(subplot_cnt)->axis(QCPAxis::atBottom), _qcp->axisRect(subplot_cnt)->axis(QCPAxis::atLeft));
    new_graph->setErrorType(QCPGraph::ErrorType::etValue);
    new_graph->setDataValueError(x, y, s);

    _tooltips[_qcp->axisRect(subplot_cnt)][new_graph] = new_sorted_tooltips;
    connect(_qcp, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(showPointToolTip(QMouseEvent*)));

    setGraphStyle(new_graph, format);
    _qcp->replot();
}
// ...........................................................................
void Plot::imagesc(ivg::Matrix A, QCPColorGradient::GradientPreset color_gradient)
// ...........................................................................
{
    A = A.transpose();

    _qcp->yAxis->setRangeReversed( true );

    QCPColorMap *colorMap = new QCPColorMap(_qcp->xAxis, _qcp->yAxis);

    _qcp->addPlottable(colorMap);
    colorMap->data()->setSize(A.rows(), A.cols());
    colorMap->data()->setRange(QCPRange(0, A.rows()-1), QCPRange(0, A.cols()-1));
    for (int x=0; x<A.rows(); ++x)
       for (int y=0; y<A.cols(); ++y)
          colorMap->data()->setCell(x, y, A(x,y));

    // add a color scale:
    QCPColorScale *colorScale = new QCPColorScale( _qcp );
    _qcp->plotLayout()->addElement(1, 1, colorScale); // add it to the right of the main axis rect
    colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
    colorMap->setColorScale(colorScale); // associate the color map with the color scale
    colorMap->setInterpolate(false);
    colorScale->axis()->setLabel("matrix elements");

    colorMap->setGradient(color_gradient);
    colorMap->rescaleDataRange(true);
    _qcp->legend->setVisible(false);
    _qcp->rescaleAxes();

    //    _qcp->xAxis->setScaleRatio(_qcp->yAxis,1.0);
    _qcp->replot();


    (*this).show();
}
// ...........................................................................
int Plot::add_rowplot( string left_label, string bottom_label )
// ...........................................................................
{
    _subplot_cnt++;

    QCPAxisRect *wideAxisRect = new QCPAxisRect(_qcp);
    wideAxisRect->axis(QCPAxis::AxisType::atLeft)->setLabel(QString::fromStdString(left_label));
    wideAxisRect->axis(QCPAxis::AxisType::atBottom)->setLabel(QString::fromStdString(bottom_label));

    _qcp->plotLayout()->addElement(_qcp->axisRectCount()+1, 0, wideAxisRect);
    _qcp->replot();

    return _qcp->axisRectCount()-1;
}
// ...........................................................................
void Plot::set_title(string title)
// ...........................................................................
{
    _qcp->plotLayout()->removeAt(0);
    _qcp->plotLayout()->addElement(0, 0, new QCPPlotTitle(_qcp, QString::fromStdString(title)));
    _qcp->replot();
}
// ...........................................................................
void Plot::setGraphStyle(QCPGraph * graph, QFormat format)
// ...........................................................................
{
    // defines properties of the LINE
    QPen pen_line(format.color, format.line_width, format.pen_style);
    graph->setPen(pen_line);
    graph->setLineStyle(format.line_style);

    // defines properties of the MARKER
    QCPScatterStyle myScatter;
    myScatter.setShape(format.marker);
    myScatter.setSize(format.marker_size);

    QPen pen_marker(format.color, format.marker_line_width, Qt::PenStyle::SolidLine);
    myScatter.setPen(pen_marker);

    graph->setScatterStyle(myScatter);

    graph->setName(format.legend);

}
// ...........................................................................
QCPGraph * Plot::add_graph(vector<double> x_std, vector<double> y_std, QCPAxis * keyAxis, QCPAxis * valueAxis)
// ...........................................................................
{
    (*this).show();

    QVector<double> x = QVector<double>::fromStdVector(x_std);
    QVector<double> y = QVector<double>::fromStdVector(y_std);

    QCPGraph *graph = _qcp->addGraph(keyAxis, valueAxis);
    graph->setData(x, y);
    graph->rescaleAxes();

    return graph;
}

// ...........................................................................
void Plot::selectionChanged()
// ...........................................................................
{
    for(int i=0; i< _qcp->axisRectCount(); i++)
    {
        QCPAxis * x_axis = _qcp->axisRect(i)->axis(QCPAxis::AxisType::atBottom,0);
        QCPAxis * y_axis = _qcp->axisRect(i)->axis(QCPAxis::AxisType::atLeft,0);

        if (x_axis->selectedParts().testFlag(QCPAxis::spAxis) || x_axis->selectedParts().testFlag(QCPAxis::spTickLabels))
            x_axis->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
        else if (y_axis->selectedParts().testFlag(QCPAxis::spAxis) || y_axis->selectedParts().testFlag(QCPAxis::spTickLabels))
            y_axis->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
    }

    // synchronize selection of graphs with selection of corresponding legend items:
    for (int i=0; i< _qcp->graphCount(); ++i)
    {
        QCPGraph *graph = _qcp->graph(i);
        QCPPlottableLegendItem *item = _qcp->legend->itemWithPlottable(graph);
      if (item->selected() || graph->selected())
      {
            item->setSelected(true);
            graph->setSelected(true);
        }
    }
}
// ...........................................................................
void Plot::contextMenuRequest(QPoint pos)
// ...........................................................................
{
    QMenu *menu = new QMenu(this);
    menu->setAttribute(Qt::WA_DeleteOnClose);

  if (_qcp->legend->selectTest(pos, false) >= 0) // context menu on legend requested
  {
    if (_qcp->selectedGraphs().size() > 0){
            menu->addAction("Change style of selected graph", this, SLOT(changeStyleSelectedGraph()));
            menu->addAction("Remove selected graph", this, SLOT(removeSelectedGraph()));
            menu->addAction("Open new window of graph", this, SLOT(performGraphAnalysis()))->setData("SIMPLE");
            menu->addAction("Toggle errorbars", this, SLOT(switchErrorbars()));

            QMenu *submenu = new QMenu(this);
            submenu->setAttribute(Qt::WA_DeleteOnClose);
            submenu->setTitle("Time Series Analysis Tools");

            menu->addMenu(submenu);
            submenu->addAction("Show Statictics", this, SLOT(performGraphAnalysis()))->setData("STATS");
            submenu->addAction("Perform FFT", this, SLOT(performGraphAnalysis()))->setData("FFT");
            submenu->addAction("Perform Moving Median dt", this, SLOT(performGraphAnalysis()))->setData("MOVMEDIANTIME");
            //submenu->addAction("Perform Moving Median #val", this, SLOT(performGraphAnalysis()))->setData("MOVMEDIANVAL");
    }
    else
    {
        menu->addAction("Move to top left", this, SLOT(moveLegend()))->setData((int)(Qt::AlignTop|Qt::AlignLeft));
        menu->addAction("Move to top center", this, SLOT(moveLegend()))->setData((int)(Qt::AlignTop|Qt::AlignHCenter));
        menu->addAction("Move to top right", this, SLOT(moveLegend()))->setData((int)(Qt::AlignTop|Qt::AlignRight));
        menu->addAction("Move to bottom right", this, SLOT(moveLegend()))->setData((int)(Qt::AlignBottom|Qt::AlignRight));
        menu->addAction("Move to bottom left", this, SLOT(moveLegend()))->setData((int)(Qt::AlignBottom|Qt::AlignLeft));
            menu->addAction("Toggle legend", this, SLOT(removeLegend()));
        }

  }
  else  // general context menu on graphs requested
    {
    if (_qcp->selectedGraphs().size() > 0){
            menu->addAction("Change style of selected graph", this, SLOT(changeStyleSelectedGraph()));
            menu->addAction("Remove selected graph", this, SLOT(removeSelectedGraph()));
            menu->addAction("Open new window of graph", this, SLOT(performGraphAnalysis()))->setData("SIMPLE");
            menu->addAction("Toggle errorbars", this, SLOT(switchErrorbars()));

            QMenu *submenu = new QMenu(this);
            submenu->setAttribute(Qt::WA_DeleteOnClose);
            submenu->setTitle("Time Series Analysis Tools");

            menu->addMenu(submenu);
            submenu->addAction("Show Statictics", this, SLOT(performGraphAnalysis()))->setData("STATS");
            submenu->addAction("Perform FFT", this, SLOT(performGraphAnalysis()))->setData("FFT");
            submenu->addAction("Perform Moving Median dt", this, SLOT(performGraphAnalysis()))->setData("MOVMEDIANTIME");
            //submenu->addAction("Perform Moving Median #val", this, SLOT(performGraphAnalysis()))->setData("MOVMEDIANVAL");
            
        } else {
            menu->addAction("Export plot", this, SLOT(openExportDialog()));
            menu->addAction("Highlight data", this, SLOT(openHighlightDialog()));
            menu->addAction("Toggle errorbars", this, SLOT(switchErrorbars()));
            menu->addAction("Toggle legend", this, SLOT(removeLegend()));
       if(_projected == true)
                menu->addAction("Projection properties", this, SLOT(changeProjectionProperties()));
        }
    }

    QCPAxis::AxisType axisType;
    QCPAxisRect * axisRect = getSelectedAxisRect(pos,axisType);
    if( axisType == QCPAxis::AxisType::atLeft || axisType == QCPAxis::AxisType::atBottom )
    {
        menu->addAction("Change Ticks", this, SLOT(openAxisRangeDialog()))->setData(pos);
        menu->addAction("Change Style", this, SLOT(openAxisStyleDialog()))->setData(pos);
    }
    if( axisType == QCPAxis::AxisType::atBottom )
    {
        QMenu *submenu = new QMenu(this);
        submenu->setAttribute(Qt::WA_DeleteOnClose);
        submenu->setTitle("Axis Domain");
        menu->addMenu(submenu);
        
        submenu->addAction("Switch frequency/period", this, SLOT(switchFreqPeriod()));
        submenu->addAction("Set frequency tick labels [Hz]", this, SLOT(tickLabelsFreqDomain()))->setData("freq_sec");
        submenu->addAction("Set frequency tick labels [per min]", this, SLOT(tickLabelsFreqDomain()))->setData("freq_min");
        submenu->addAction("Set frequency tick labels [per hour]", this, SLOT(tickLabelsFreqDomain()))->setData("freq_hour");
        submenu->addAction("Set frequency tick labels [per day]", this, SLOT(tickLabelsFreqDomain()))->setData("freq_day");
        submenu->addAction("Set frequency tick labels [per year]", this, SLOT(tickLabelsFreqDomain()))->setData("freq_year");
        submenu->addAction("Set frequency tick labels [per decade]", this, SLOT(tickLabelsFreqDomain()))->setData("freq_dec");
        submenu->addAction("Set period tick labels [secs]", this, SLOT(tickLabelsFreqDomain()))->setData("period_sec");
        submenu->addAction("Set period tick labels [mins]", this, SLOT(tickLabelsFreqDomain()))->setData("period_min");
        submenu->addAction("Set period tick labels [hours]", this, SLOT(tickLabelsFreqDomain()))->setData("period_hour");
        submenu->addAction("Set period tick labels [days]", this, SLOT(tickLabelsFreqDomain()))->setData("period_day");
        submenu->addAction("Set period tick labels [years]", this, SLOT(tickLabelsFreqDomain()))->setData("period_year");
        submenu->addAction("Set period tick labels [decades]", this, SLOT(tickLabelsFreqDomain()))->setData("period_dec");
    }

    menu->popup(_qcp->mapToGlobal(pos));

}
// ...........................................................................
void Plot::changeProjectionProperties()
// ...........................................................................
{
    // initialize new dialog containing changable settings
    Projection *projection_dialog = new Projection(this);

    // in order to know if there is a porjectiontype change
    protype old_projection_type = _ref_pro.type;

    // set pointer which will be changed after dialog is accepted
    projection_dialog->setProjection(&_ref_pro);

    // show dialog
    projection_dialog->show();

    // only if user pushed OK, apply changes
    if (projection_dialog->exec() == QDialog::Accepted) 
    {
        // remove old grid and generate new grid based on changed _ref_pro
        _qcp->removeLayer(_qcp->layer("pro_grid"));
        _plot_projected_grid(_ref_pro);

        // only plot coastline if we really need to do it!
//        if(old_projection_type != _ref_pro.type)
//        {
            _qcp->removeLayer(_qcp->layer("pro_coastline"));
            _plot_coastline(_ref_pro.type, { Qt::black, _ref_pro.grid_width, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssNone, " " });
//        }
        _qcp->replot();
    }
}
// ...........................................................................
void Plot::changeStyleSelectedGraph()
// ...........................................................................
{
    if (_qcp->selectedGraphs().size() == 1)
    {
        GraphDialog *graph_dialog = new GraphDialog(this);

        //set selected graph as graph to edit
        graph_dialog->setGraph(_qcp->selectedGraphs().first());
        //show dialog
        graph_dialog->show();

        //only if user pushed OK, apply changes
        if (graph_dialog->exec() == QDialog::Accepted) 
        {
            graph_dialog->updateGraph();
            _qcp->replot();
        }
    }
}
// ...........................................................................
void Plot::moveLegend()
// ...........................................................................
{
    // make sure this slot is really called by a context menu action, so it carries the data we need
  if (QAction* contextAction = qobject_cast<QAction*>(sender())) 
  {
        bool ok;
        int dataInt = contextAction->data().toInt(&ok);
    if (ok)
    {
            _qcp->axisRect()->insetLayout()->setInsetAlignment(0, (Qt::Alignment)dataInt);
            _qcp->replot();
        }
    }
}
// ...........................................................................
void Plot::removeLegend()
// ...........................................................................
{
    if(_qcp->legend->visible())
        _qcp->legend->setVisible(false);
    else
        _qcp->legend->setVisible(true);

    _qcp->replot();
}
// ...........................................................................
void Plot::removeSelectedGraph()
// ...........................................................................
{
  if (_qcp->selectedGraphs().size() > 0)
  {
        _qcp->removeGraph(_qcp->selectedGraphs().first());
        _qcp->replot();
        
        _index.erase(_qcp->selectedGraphs().first());
    }
}
// ...........................................................................
void Plot::switchErrorbars()
// ...........................................................................
{
  if (_qcp->selectedGraphs().size() == 1)
  {
    if(_qcp->selectedGraphs().first()->errorType() == QCPGraph::ErrorType::etNone)
            _qcp->selectedGraphs().first()->setErrorType(QCPGraph::ErrorType::etValue);
    else if(_qcp->selectedGraphs().first()->errorType() == QCPGraph::ErrorType::etValue)
            _qcp->selectedGraphs().first()->setErrorType(QCPGraph::ErrorType::etNone);
    }
  else
  {
        // if errorbars from all graphs should ne switched
      for(int i=0; i<_qcp->graphCount(); i++)
      {
        if(_qcp->graph(i)->errorType() == QCPGraph::ErrorType::etNone)
                _qcp->graph(i)->setErrorType(QCPGraph::ErrorType::etValue);
        else if(_qcp->graph(i)->errorType() == QCPGraph::ErrorType::etValue)
                _qcp->graph(i)->setErrorType(QCPGraph::ErrorType::etNone);
        }
    }

    _qcp->replot();
}
// ...........................................................................

void Plot::switchFreqPeriod()
// ...........................................................................
{
    // get data [mjd, value, std] of all selected graphs
    ivg::Matrix mjd, secs, value, std;
    vector<string> tooltips;

//    std::cerr << "Switched to " << "???" << std::endl;  //TODO: status

    QColor new_color;
    QList<QCPGraph*> graph_list = _qcp->xAxis->graphs();
    double maxX = 0.0;
    double maxY = 0.0;
    double minY = 1.0/0.0;
    for (auto *graph : graph_list) {
        // get axis of current plot
        QCPAxis *key_axis = graph->keyAxis();
        QCPAxis *value_axis = graph->valueAxis();

        // all data points of selected graph
        QCPDataMap * data = graph->data();
        QVector<double> x_reziprok(data->count(), 0);
        QVector<double> y_neu(data->count(), 0);
        double xMax = 0;
        double yMax = 0;
        int k = 0;
        // Invert keys
        for (auto &point : (*data)) {
            x_reziprok.insert(data->count() - k - 1, 1.0 / point.key);
            y_neu.insert(data->count() - k - 1, point.value);

            if (1.0 / point.key > xMax && 1.0 / point.key < 1.0/0.0)
                xMax = 1.0 / point.key;
            if (point.value > yMax && 1.0 / point.key < 1.0/0.0)
                yMax = point.value;

            k++;
        }
        graph->setData(x_reziprok, y_neu);

        double yMin = std::numeric_limits<double>::max();

        foreach(double p, y_neu) {
            yMin = qMin(yMin, p);
        }

        if (yMax > maxY)
            maxY = yMax;
        if (yMin < minY)
            minY = yMin;
        if (xMax > maxX)
            maxX = xMax;

    }

    _qcp->xAxis->setRange(0, maxX);
    _qcp->yAxis->setRange(minY, maxY);

    _qcp->replot();
}
// ...........................................................................
void Plot::tickLabelsFreqDomain()
// ...........................................................................
{
    QAction* contextAction = qobject_cast<QAction*>(sender());
    QString perform_type = (QString) contextAction->data().toString();
    QList<QCPGraph*> graph_list = _qcp->xAxis->graphs();
    double maxX = 0.0;
    double minX = 1.0/0.0;
    for (auto *graph : graph_list) {
        QCPAxis *key_axis = graph->keyAxis();
        QCPAxis *value_axis = graph->valueAxis();
        
        // all data points of selected graph
        QCPDataMap * data = graph->data();
        double xMax = 0;
        double xMin = 1.0/0.0;
        for (auto &point : (*data)) {
            if (point.key >= key_axis->range().lower && point.key <= key_axis->range().upper) {
                if (point.key > xMax && point.key < 1.0 / 0.0)
                    xMax = point.key;
                if (point.key < xMin && point.key < 1.0 / 0.0)
                    xMin = point.key;
            }
        }
        if (xMax > maxX)
            maxX = xMax;
        if (xMin < minX)
            minX = xMin;
    }
    
    QVector<double> vec;
    QVector<QString> strVec;
    double unit;
    string unit_str;
    if(perform_type=="freq_sec") {
        unit = 1.0;
        unit_str = "Hz";
    }else if (perform_type=="freq_min") {
        unit = 1.0/60.0;
        unit_str = "pm";
    }else if(perform_type=="freq_hour") {
        unit = 1.0/3600.0;
        unit_str = "ph";
     } else if(perform_type=="freq_day") {
        unit = 1.0/(24.0*3600.0);
        unit_str = "pd";
     } else if(perform_type=="freq_year") {
        unit = 1.0/(365.0*24.0*3600.0);
        unit_str = "py";
    } else if(perform_type=="freq_dec") {
        unit = 1.0/(10.0*365.0*24.0*3600.0);
        unit_str = "pdec";
    }else if (perform_type=="period_min") {
        unit = 60.0;
        unit_str = "min";
    }else if(perform_type=="period_hour") {
        unit = 3600.0;
        unit_str = "hr";
     } else if(perform_type=="period_day") {
        unit = 24.0*3600.0;
        unit_str = "days";
     } else if(perform_type=="period_year") {
        unit = 365.0*24.0*3600.0;
        unit_str = "yr";
    } else if(perform_type=="period_dec") {
        unit = 10.0*365.0*24.0*3600.0;
        unit_str = "dec";
    }else {
        unit = 1.0;
        unit_str = "sec";
    }
    int numTicks = 6;
    
    double pw = pow(10,floor(log10(maxX/(unit*numTicks))));
    double schritt = round(maxX/(unit*numTicks*pw)) * pw*unit;
    if(schritt==0){
        schritt = maxX / numTicks;
    }
    
    int k = 0;
    for (double d = minX; d <= maxX; d = d + schritt) {
        vec.push_back(d);
        char stemp[100] = "";
        if(schritt/unit>1) {
            snprintf(stemp,100,"%.0f",d/unit);
        }else if (schritt/unit>0.1){
            snprintf(stemp,100,"%.1f",d/unit);
        }else if (schritt/unit>0.01){
            snprintf(stemp,100,"%.2f",d/unit);
        }else {
            snprintf(stemp,100,"%.2e",d/unit);
        }
        string s=stemp;
        strVec.push_back(QString::fromStdString(s+unit_str));
        k++;
    }
    _qcp->xAxis->setAutoTicks(false);
    _qcp->xAxis->setAutoTickLabels(false);
    _qcp->xAxis->setTickVector(vec);
    _qcp->xAxis->setTickVectorLabels(strVec);

    _qcp->replot();
    _qcp->xAxis->setAutoTicks(true);
    _qcp->xAxis->setAutoTickLabels(true);
}
// ...........................................................................
void Plot::performGraphAnalysis()
// ...........................................................................
{

  if (_qcp->selectedGraphs().size() > 0)
  {
        QAction* contextAction = qobject_cast<QAction*>(sender());
        QString perform_type = (QString) contextAction->data().toString();
      
        _sec_plot = new Plot();
        
        // get data [mjd, value, std] of all selected graphs
        ivg::Matrix mjd, secs, value, std;
        vector<string> tooltips;

        QColor new_color;
        QList<QCPGraph*> graph_list = _qcp->selectedGraphs();
        for(auto *graph: graph_list)
        {
            // get axis of current plot
            QCPAxis *key_axis = graph->keyAxis();
            QCPAxis *value_axis = graph->valueAxis();

            // all data points of selected graph
            QCPDataMap * data = graph->data();
            new_color = graph->pen().color();

            // reference mjd to calculate sec since 1970 to mjd
            double mjd_1970 = 40587.0; // generated with ivg::Date( 1970, 1.0 ).get_double_mjd();

            map<double, string> sorted_tooltips = _tooltips[key_axis->axisRect()][graph];

            // fill new data matrices only of the visibile area in current plot
            int row=0;
            for(auto &point: (*data))
            {
                if(point.key >= key_axis->range().lower && point.key <= key_axis->range().upper
                    && point.value >= value_axis->range().lower && point.value <= value_axis->range().upper)
                {         
                    secs.append_rows(point.key);
                    mjd.append_rows((point.key / 86400.0) + mjd_1970);
                    value.append_rows(point.value);
                    std.append_rows(point.valueErrorMinus);
                    tooltips.push_back(sorted_tooltips[point.key]);
                }
            }
        }

        // initialize time series analyzer
        ivgat::Tsa tsa;
        if(perform_type == "SIMPLE")
        {
            // plot data in a secondary new plot
            QFormat format1 = {new_color, 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, "X", 1.0, 3.0};
            _sec_plot->plot_mjd_series_std(mjd, value, std, format1, 0, "dd-MMM-yy", tooltips);
            _sec_plot->_qcp->legend->setVisible(false);

            string tileString = "WRMS = " + std::to_string(tsa.WRMS(value, std)) + ", RMS = " + std::to_string(tsa.RMS(value));
            _sec_plot->set_title(tileString);
            _sec_plot->show();
        }
        else if(perform_type == "FFT")
        {
            // plot data in a secondary new plot
            QFormat format1 = {new_color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, "X", 1.0, 3.0};

            ivg::Matrix AA,nu2;
            tsa.fourierAmplitudeSpectrumEQD(secs, value, nu2, AA);

            _sec_plot->plot_data(nu2.transpose(), AA, format1);
            _sec_plot->set_title("Spektrum");
            _sec_plot->get_plot()->yAxis->setLabel("Amp(nu)");
            _sec_plot->get_plot()->xAxis->setLabel("nu");
            _sec_plot->get_plot()->replot();
            _sec_plot->_qcp->legend->setVisible(true);
            _sec_plot->show();
            
        }
        else if(perform_type == "STATS")
        {
            Ui::Statistics_info ui;
            QDialog *highlight = new QDialog();
            ui.setupUi(highlight);
            highlight->show();
            ui.textBrowser->clear();
            
            stringstream info_ss;
            info_ss << "Statistics summary:\n" << endl;
            info_ss <<    "Mean: " << value.meanD() << endl;
            info_ss <<    "Standard deviation: " << value.stdD() << endl;
            info_ss <<    "Median: " << value.median() << endl;
            info_ss <<    "RMS: " << tsa.RMS(value) << endl;
            info_ss <<    "WRMS: " << tsa.WRMS(value,std) << endl;
    
            ui.textBrowser->setText(QString::fromStdString(info_ss.str()));
        }
        else if(perform_type == "MOVMEDIANTIME")
        {
            // plot data in a secondary new plot
            QFormat format1 = {new_color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, "median", 1.0, 3.0};

             
////            Ui::Movfise ui;
//            Movfise ui;
////            ui.widget
////            ui.widget
//            QDialog *moving = new QDialog();
//            ui.setupUi(moving);
//            moving->show();
            
            Movfise *export_dialog = new Movfise(this);
            export_dialog->show();

        if (export_dialog->exec()== QDialog::Accepted)
            {
                Ui::Movfise *ui = &export_dialog->widget;
//                ui->
                double window = ui->size_spinbox->value();
                double stepsize = ui->stepsize_spinbox->value();
                bool chked = ui->same_checkbox->isChecked();

                ivg::Matrix AA,nu2;
                tsa.movingMedian(secs, value, window, stepsize, nu2, AA);

                string titleString = "MovMedian win=" + std::to_string(window);
                QString qstring = QString::fromStdString(titleString);
                if (chked){
                    // plot data in a secondary new plot
                    format1 = {new_color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, qstring, 1.0, 3.0};
                    this->plot_data(nu2, AA, format1);
                }else {
                    // plot data in a secondary new plot
                    format1 = {new_color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, qstring, 1.0, 3.0};
                    _sec_plot->plot_data(nu2, AA, format1);
                    _sec_plot->set_title(titleString);
                    _sec_plot->get_plot()->replot();
                    _sec_plot->_qcp->legend->setVisible(true);
                    _sec_plot->show();
                }
            }
        }
        else if(perform_type == "MOVMEDIANVAL")
        {
             // plot data in a secondary new plot
            //  QFormat format1 = { QColor(ivg::color_values.at(_color_cnt).c_str()), 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, "", 1.0 , 3.0  };
            QFormat format1 = {new_color, 2.0, Qt::SolidLine, QCPGraph::lsLine, QCPScatterStyle::ssDisc, "median", 1.0, 3.0};

            ivg::Matrix AA,nu2;

            int window = 20;
            AA = tsa.movingMedian(value, window);
            nu2 = mjd.get_sub((window - 1) / 2 - 1, 0, mjd.size(1)-(window - 1) / 2 - 2, 0); 
            string tileString = "Moving Median: " + std::to_string(window); 

    //        double rms = tsa.RMS(AA);
    //        tileString += " values, RMS = " + std::to_string(rms);

            _sec_plot->plot_data(nu2, AA, format1);
            _sec_plot->set_title(tileString);
            _sec_plot->get_plot()->replot();
            _sec_plot->_qcp->legend->setVisible(true);
            _sec_plot->show();
        }
    }
}
// ...........................................................................
void Plot::openExportDialog()
// ...........................................................................
{
    ExportDialog *export_dialog = new ExportDialog(this);
    export_dialog->show();

    if (export_dialog->exec() == QDialog::Accepted) {
        Ui::Export *exp = &export_dialog->ui;

        bool ok = false;
        if( exp->RB_pdf->isChecked() )
            ok = _qcp->savePdf( exp->LE_pathfile->text(), exp->CB_cospen, exp->SB_width->value(), exp->SB_height->value() , QString(), exp->LE_pdftitle->text() ) ;
        else if(exp->RB_png->isChecked() )
            ok = _qcp->savePng( exp->LE_pathfile->text(), exp->SB_width->value(), exp->SB_height->value(), exp->SB_scale->value(), exp->SB_quality->value() ) ;
        else if(exp->RB_bmp->isChecked() )
            ok = _qcp->saveBmp( exp->LE_pathfile->text(), exp->SB_width->value(), exp->SB_height->value(), exp->SB_scale->value()) ;
        else if(exp->RB_jpg->isChecked() )
            ok = _qcp->saveJpg( exp->LE_pathfile->text(), exp->SB_width->value(), exp->SB_height->value(), exp->SB_scale->value(), exp->SB_quality->value()) ;

        QMessageBox msgBox;
        if(ok == false)
        {
            msgBox.setIcon(QMessageBox::Icon::Warning);
            msgBox.setText("Export not successful. Please try again.");
        }
        else
        {
            msgBox.setIcon(QMessageBox::Icon::Information);
            msgBox.setText("Plot exported to "+exp->LE_pathfile->text());
        }
        msgBox.exec();
    }
}
// ...........................................................................
void Plot::openHighlightDialog()
// ...........................................................................
{
    // create and show new highlight dialog
    Ui::Highlight ui;
    QDialog *highlight = new QDialog();
    ui.setupUi(highlight);
    highlight->show();

    // in case of pushed OK-button
    if (highlight->exec() == QDialog::Accepted) {

        // now get the desired graphstyle for the new highlighted datapoints
        GraphDialog *graph_dialog = new GraphDialog(this);
        graph_dialog->setGraph(_qcp->graph());
        graph_dialog->show();

        //only if user pushed OK, apply changes
        if (graph_dialog->exec() == QDialog::Accepted) 
        {
            string identify;
            bool remove = ui.remove_checkbox->isChecked();
            // in case of highlighting base on string comparison
            if(ui.type1_radiobutton->isChecked())
            {
                identify = ui.identify_lineedit->text().toStdString();
                // in case of specific row != 0, only look in a specific row from tooltip
                int row = ui.row1_spinbox->value();
               if( row != 0)
               {
                    stringstream ss;
                    ss << "{" << row << "}" << identify;
                    identify = ss.str();
                }
            }
            // in case of highlighting based on doube-value comparison
            else if(ui.type2_radiobutton->isChecked())
            {
                double value = ui.value_spinbox->value();
                int row = ui.row2_spinbox->value();
               if( row != 0)
               {
                    stringstream ss;
                    ss << "[>" << row << "]" << value;
                    identify = ss.str();
                }
            }

            // if checked, highlight data in all subplots
            int subplot_end = 0;
            if(ui.allsubs_checkbox->isChecked())
                subplot_end = _subplot_cnt;

            for(int subplot=0; subplot <= subplot_end; subplot++ )
            {
                 (*this).highlight( identify, { QColor(Qt::green), 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssCircle, QString::fromStdString(identify), 1.0, 6.0 }, remove, subplot);
                graph_dialog->updateGraph(_qcp->graph());
            }

            _qcp->replot();
        }
    }
}
// ...........................................................................
void Plot::openAxisRangeDialog()
// ...........................................................................
{
    QAction* contextAction = qobject_cast<QAction*>(sender());

    QPoint pos = (QPoint) contextAction->data().toPoint();

    QCPAxis::AxisType axisType;
    QCPAxisRect * axisRect = getSelectedAxisRect(pos,axisType);
    QCPRange tmp_range = axisRect->axis(axisType,0)->range();

    Ui::Range ui_tmp;
    QDialog *axisRange = new QDialog();
    ui_tmp.setupUi(axisRange);

    // if tick label is not a value but a string e.g. 14-JAN-07
    if(axisRect->axis(axisType,0)->tickLabelType() == QCPAxis::LabelType::ltDateTime)
    {
        ui_tmp.LE_max->setText(QDateTime::fromTime_t((uint)tmp_range.upper).toString("dd-MM-yyyy")); // axisRect->axis(axisType,0)->dateTimeFormat()
        ui_tmp.LE_min->setText(QDateTime::fromTime_t((uint)tmp_range.lower).toString("dd-MM-yyyy"));
    }
    else
    {
        ui_tmp.LE_max->setText(QString::number(tmp_range.upper));
        ui_tmp.LE_min->setText(QString::number(tmp_range.lower));
    }
    axisRange->show();

    if (axisRange->exec() == QDialog::Accepted)
    {

        double max,min;
        if(axisRect->axis(axisType,0)->tickLabelType() == QCPAxis::LabelType::ltDateTime)
        {
            max = (double)QDateTime(QDate::fromString(ui_tmp.LE_max->text(),  "dd-MM-yyyy")).toTime_t();
            min = (double)QDateTime(QDate::fromString(ui_tmp.LE_min->text(),  "dd-MM-yyyy")).toTime_t();

        }
        else
        {
            max = ui_tmp.LE_max->text().toDouble();
            min = ui_tmp.LE_min->text().toDouble();
        }

        axisRect->axis(axisType,0)->setRange(min,max);

        axisRect->axis(axisType,0)->setTickLabelRotation(ui_tmp.LE_rot->text().toDouble());

        _qcp->replot();
    }
}
// ...........................................................................
void Plot::openAxisStyleDialog()
// ...........................................................................
{
    bool ok;
    QFont font = QFontDialog::getFont( &ok, QFont("Helvetica [Cronyx]", 10), this);
    if (ok) {

        QAction* contextAction = qobject_cast<QAction*>(sender());
        QPoint pos = (QPoint) contextAction->data().toPoint();

        for(int i=0; i< _qcp->axisRectCount(); i++)
        {
            QCPAxisRect * axisRect = _qcp->axisRect(i);
            if(axisRect->selectTest(pos, false) >=0)
            {
                axisRect->axis(QCPAxis::AxisType::atLeft,0)->setSelectedLabelFont(font);
                axisRect->axis(QCPAxis::AxisType::atLeft,0)->setSelectedTickLabelFont(font);
                axisRect->axis(QCPAxis::AxisType::atLeft,0)->setLabelFont(font);
                axisRect->axis(QCPAxis::AxisType::atLeft,0)->setTickLabelFont(font);
                axisRect->axis(QCPAxis::AxisType::atBottom,0)->setSelectedLabelFont(font);
                axisRect->axis(QCPAxis::AxisType::atBottom,0)->setSelectedTickLabelFont(font);
                axisRect->axis(QCPAxis::AxisType::atBottom,0)->setLabelFont(font);
                axisRect->axis(QCPAxis::AxisType::atBottom,0)->setTickLabelFont(font);
            }
        }
        _qcp->replot();
    }
}
// ...........................................................................
void Plot::zoomBoxRequest(QMouseEvent* event)
// ...........................................................................
{
    if (event->button() == Qt::MiddleButton )
    {
        _click_origin = event->pos();
        _rubber_band->setGeometry(QRect(_click_origin, QSize()));
        _rubber_band->show();
    }
}
// ...........................................................................
void Plot::mouseMoveEvent(QMouseEvent * event)
// ...........................................................................
{
    if (_rubber_band->isVisible())
        _rubber_band->setGeometry(QRect(_click_origin, event->pos()).normalized());
}
// ...........................................................................
void Plot::mouseReleaseEvent(QMouseEvent * event)
// ...........................................................................
{
    if (_rubber_band->isVisible())
    {
        if (event->button() == Qt::MiddleButton)
        {
            // loop over all axisRects (axex that are rectamgular) probably x and y Axis ?
            for(int i=0; i< _qcp->axisRectCount(); i++)
            {
                QCPAxisRect * axisRect = _qcp->axisRect(i);
                if(axisRect->selectTest(_click_origin, false) >=0)
                {
                    const QRect & zoomRect = _rubber_band->geometry();
                    int xp1, yp1, xp2, yp2;
                    zoomRect.getCoords(&xp1, &yp1, &xp2, &yp2);
                    auto x1 = axisRect->axis(QCPAxis::AxisType::atBottom,0)->pixelToCoord(xp1); // x_min
                    auto x2 = axisRect->axis(QCPAxis::AxisType::atBottom,0)->pixelToCoord(xp2); // x_max
                    auto y1 = axisRect->axis(QCPAxis::AxisType::atLeft,0)->pixelToCoord(yp1); // y_min
                    auto y2 = axisRect->axis(QCPAxis::AxisType::atLeft,0)->pixelToCoord(yp2); // y_max

                    // depending on pressed buttons, remove data, create new graph or add data points to existing graph
                    // STRG+ALT = removes data in rubberband
                    // STRG+SHIFT = add new graph for data in rubberband
                    // STRG+SHIFT+ALT = add data in rubberband to selcted graph
                    if(event->modifiers().testFlag(Qt::ControlModifier) == true || event->modifiers().testFlag(Qt::AltModifier) == true)
                    {
                        // get data [mjd, value, std] of all selected graphs
                        ivg::Matrix key, value, std;
                        vector<string> tooltips;
                        vector<int> orig_idx;
			
                        QList<QCPGraph*> graph_list;
                        if(_qcp->selectedGraphs().size() > 0 )
                            graph_list = _qcp->selectedGraphs();
                        else
                            graph_list = axisRect->graphs();

                        if(event->modifiers().testFlag(Qt::ControlModifier) == true && event->modifiers().testFlag(Qt::AltModifier) == true &&
                                event->modifiers().testFlag(Qt::ShiftModifier) == true)
                            graph_list = axisRect->graphs();
			
                        for(auto *graph: graph_list)
                        {
                            // get axis of current plot
                            QCPAxis *key_axis = graph->keyAxis();
                            QCPAxis *value_axis = graph->valueAxis();

                            // all data points of graph
                            QCPDataMap * data = graph->data();
                            // remove all points within rubberband
                            int count = 0;
                            for(auto &point: (*data))
                            {
                                if(point.key >= x1 && point.key <= x2 && point.value >= y2 && point.value <= y1)
                                {
                                    map<double, string> sorted_tooltips = _tooltips[key_axis->axisRect()][graph];
                                    
                                    // --- remove selected points if ALT and CTRL is pressed ---
                                    if(event->modifiers().testFlag(Qt::ControlModifier) == true && event->modifiers().testFlag(Qt::AltModifier) == true
                                            && event->modifiers().testFlag(Qt::ShiftModifier) == false)
                                    {
                                        data->remove(point.key);
                                        if(!sorted_tooltips[point.key].empty())
                                        {
                                            QString tt = QString::fromStdString(sorted_tooltips[point.key]);
                                            emit pointRemoved(tt);
                                            sorted_tooltips.erase(point.key);
                                        }
                                    }
                                    // --- create new graph if CTRL and SHIFT is pressed ---
                                    else if(event->modifiers().testFlag(Qt::ControlModifier) == true && event->modifiers().testFlag(Qt::ShiftModifier) == true)
                                    {
				    
                                        // fill new data matrices only of the visibile area in current plot
                                        key.append_rows(point.key);
				
                                        value.append_rows(point.value);
					
                                        std.append_rows(point.valueErrorMinus);
					
                                        tooltips.push_back(sorted_tooltips[point.key]);
					if  (_index[graph].size()>0)
					  orig_idx.push_back(_index[graph].at(count));
					
                                        
                                    }
                                }
                                count++;
                            }
                        }
		
                        QVector<double> x = QVector<double>::fromStdVector(vector<double> (key.begin(),key.end()));
                        QVector<double> y = QVector<double>::fromStdVector(vector<double> (value.begin(),value.end()));
                        QVector<double> s = QVector<double>::fromStdVector(vector<double> (std.begin(),std.end()));
			
                        if(event->modifiers().testFlag(Qt::ControlModifier) == true && event->modifiers().testFlag(Qt::AltModifier) == true &&
                                event->modifiers().testFlag(Qt::ShiftModifier) == true)
                        {
                            if(_qcp->selectedGraphs().size() == 1)
                                _qcp->selectedGraphs().first()->addData(x,y);
                        }
                        // --- create new graph if CTRL and SHIFT is pressed ---
                        else if(event->modifiers().testFlag(Qt::ControlModifier) == true && event->modifiers().testFlag(Qt::ShiftModifier) == true)
                        {
			  
                            QCPGraph *new_graph = _qcp->addGraph(axisRect->axis(QCPAxis::AxisType::atBottom,0), axisRect->axis(QCPAxis::AxisType::atLeft,0));
                            _index[new_graph] = orig_idx;
                            _color_cnt++;
			    
                            QFormat format = { QColor(ivg::color_values.at(_color_cnt).c_str()), 1.0, Qt::SolidLine, QCPGraph::lsNone, QCPScatterStyle::ssDisc, "X", 1.0 , 3.0  };
                            setGraphStyle(new_graph, format);
                            new_graph->setErrorType(QCPGraph::ErrorType::etNone);
                            _qcp->graph()->setDataValueError(x, y, s);
			    
                        }
                        _qcp->replot();
		        
                    }
                    else
                    {
                        axisRect->axis(QCPAxis::AxisType::atBottom,0)->setRange(x1, x2);
                        axisRect->axis(QCPAxis::AxisType::atLeft,0)->setRange(y1, y2);
                    }

                    _rubber_band->hide();
		    
                    _qcp->replot();
                }
            }
        }
    }
}
// ...........................................................................
void Plot::rescalePlot(QMouseEvent * event)
// ...........................................................................
{
    if (event->button() == Qt::MiddleButton )
    {
        // depens on current plot which scales need to be used
        // in case of projecteion plot
        if(_projected)
        {
            if( _ref_pro.type == protype::polarstereo )
            {
                _qcp->xAxis->setRange(-1.8,1.8);
                _qcp->yAxis->setRange(-1.8,1.8);
            }
        else
            {
                _qcp->xAxis->setRange(-3.1,3.1);
                _qcp->yAxis->setRange(-1.7,1.9);
            }
        }
        // in case of regular plot
        else
            _qcp->rescaleAxes();

        _qcp->replot();
    }
}
// ...........................................................................
void Plot::mouseWheelTurned(QWheelEvent* event)
// ...........................................................................
{
    for(int i=0; i< _qcp->axisRectCount(); i++)
    {
        QCPAxisRect * axisRect = _qcp->axisRect(i);
        if(axisRect->selectTest(event->pos(), false) >=0)
        {
            QCPAxis * x_axis = axisRect->axis(QCPAxis::AxisType::atBottom,0);
            QCPAxis * y_axis = axisRect->axis(QCPAxis::AxisType::atLeft,0);

            if(x_axis->selectedParts().testFlag(QCPAxis::spAxis))
                axisRect->setRangeZoom(x_axis->orientation());
            else if(y_axis->selectedParts().testFlag(QCPAxis::spAxis))
                axisRect->setRangeZoom(y_axis->orientation());
            else
                axisRect->setRangeZoom(Qt::Horizontal | Qt::Vertical);
        }
    }
}
// ...........................................................................
void Plot::titleDoubleClick(QMouseEvent* event, QCPPlotTitle* title)
// ...........................................................................
{
    Q_UNUSED(event)

            bool ok;
    QFont font = QFontDialog::getFont( &ok, QFont("Helvetica [Cronyx]", 10), this);
    if (ok) {

        title->setFont(font);

        QString newTitle = QInputDialog::getText(this, " ", "New plot title:", QLineEdit::Normal, title->text(), &ok);
        if (ok)
        {
            title->setText(newTitle);
            _qcp->replot();
        }
    }


}
// ...........................................................................
void Plot::axisLabelDoubleClick(QCPAxis *axis, QCPAxis::SelectablePart part)
// ...........................................................................
{
    // Set an axis label by double clicking on it
    if (part == QCPAxis::spAxisLabel) // only react when the actual axis label is clicked, not tick label or axis backbone
    {
        bool ok;
        QString newLabel = QInputDialog::getText(this, " ", "New axis label:", QLineEdit::Normal, axis->label(), &ok);
    if (ok)
    {
            axis->setLabel(newLabel);
            _qcp->replot();
        }
    }
}
// ...........................................................................
void Plot::legendDoubleClick(QCPLegend *legend, QCPAbstractLegendItem *item)
// ...........................................................................
{
    // Rename a graph by double clicking on its legend item
    Q_UNUSED(legend)
    if (item) // only react if item was clicked (user could have clicked on border padding of legend where there is no item, then item is 0)
    {
        QCPPlottableLegendItem *plItem = qobject_cast<QCPPlottableLegendItem*>(item);
        bool ok;
    QFont font = QFontDialog::getFont( &ok, QFont("Helvetica [Cronyx]", 10), this);
        if (ok) {

            plItem->setFont(font);
            plItem->setSelectedFont(font);

            QString newName = QInputDialog::getText(this, " ", "New graph name:", QLineEdit::Normal, plItem->plottable()->name(), &ok);
        if (ok)
        {
                plItem->plottable()->setName(newName);
                _qcp->replot();
            }
        }
    }
}
// ...........................................................................
void Plot::showPointToolTip(QMouseEvent *event)
// ...........................................................................
{
    // get position of event
    QPoint pos =  event->pos();

    QCPAxis::AxisType axisType;
    QCPAxisRect * axisRect = getSelectedAxisRect(pos,axisType);

    // only in case if event is over an axisRect
    if(axisRect != NULL)
    {
        // get coordinates of event-position
        //        double x = axisRect->axis(QCPAxis::atBottom,0)->pixelToCoord(event->pos().x());
        //        double y = axisRect->axis(QCPAxis::atLeft,0)->pixelToCoord(event->pos().y());

        double x_mp = event->pos().x(); // mp = mouse pixel
        double y_mp = event->pos().y(); // mp = mouse pixel

        QList<QCPGraph*> selected_graphs = _qcp->selectedGraphs();

        // get all data plotted in the figure
        for(auto const &graph : selected_graphs ) 
        { 
            if( graph->selectTest(event->pos(), false) < _qcp->selectionTolerance() )
            {
                map<double, string> tooltips = _tooltips[axisRect][graph];

                // search for neartes point by euclidean distance in the pixel system
                double dist = 1e5;
                double final_key;
                for(QMap<double,QCPData>::iterator data = graph->data()->begin(); data != graph->data()->end(); ++data)
                {
                    double y_pp = axisRect->axis(QCPAxis::atLeft,0)->coordToPixel(data.value().value);
                    double x_pp = axisRect->axis(QCPAxis::atBottom,0)->coordToPixel(data.value().key);

                    double new_dist = sqrt( pow(y_pp - y_mp, 2.0 ) + pow(x_pp - x_mp, 2.0 ) );

                    if( new_dist < dist)
                    {
                        dist = new_dist;
                        final_key = data.value().key;
                    }
                }
                // show tooltip of nearest point
                QString tt = QString::fromStdString(tooltips[final_key]);
                QToolTip::showText(QCursor::pos(), tt);
                emit tooltipShowed(tt);
            }
        }
    }

}
// ...........................................................................
QCPAxisRect * Plot::getSelectedAxisRect(QPoint &pos, QCPAxis::AxisType &axisType)
// ...........................................................................
{
    for(int i=0; i< _qcp->axisRectCount(); i++)
    {
        QCPAxisRect * axisRect = _qcp->axisRect(i);
        if(axisRect->selectTest(pos, false) >=0)
        {
            if(axisRect->axis(QCPAxis::AxisType::atLeft,0)->selectTest(pos, false) >=0)
                axisType = QCPAxis::AxisType::atLeft;
            else if(axisRect->axis(QCPAxis::AxisType::atBottom,0)->selectTest(pos, false) >=0)
                axisType = QCPAxis::AxisType::atBottom;

            return axisRect;
        }
    }

    return NULL;
}

// ...........................................................................
vector<string> Plot::tokenize_string( string &str, const char delimiter )
// ...........................................................................
{

    string token;
    vector<string> tokens;
    stringstream tt_stream(str);
    while (getline(tt_stream, token, delimiter))
        tokens.push_back(token);

    return tokens;
}

