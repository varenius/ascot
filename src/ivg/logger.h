#ifndef FORMATTED_LOG_T_H
#define	FORMATTED_LOG_T_H

/**
*
* @brief Logger class
* @author AI - bakkari developer team
* @date 2015-04-02
* @version 0.1
*
*/
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <QtCore/QString>
#include <QtWidgets/QLineEdit>

using namespace std;

enum loglevel {
    NOTHING,
    WARNING,
    RESULT,
    INFO,
    DETAIL,
    ALL,
    SAVE,
    QT
};


//global variable, set in main
extern loglevel g_verbose;

namespace ivg
{
     
class Logger
// ===========================================================================
{
    public:
        
    // save-information added to the line if loglevel SAVE is called
    static string add_to_saveline;
    // QLineEdit emits signal if log<XXX> is used and qt_output is defined ( != NULL )
    static QLineEdit *qt_output;
           
    // assigns each level a specific color
    static array<string,8> level_color;
    
    static std::string get_color(string s){ return _color.at(s);};
        
    Logger( loglevel level, string msg ) : fmt(msg), level(level) {}
    Logger( loglevel level, string path, string msg ) : path(path), fmt(msg), level(level) {}
    
    ~Logger() {
        if ( level <= g_verbose && g_verbose != 5 && qt_output == NULL ){
            std::cout << level_color[level] <<  fmt <<  _color["white"] <<endl;	 
        // in case of g_verbose also the loglevel will be shown
        } else if ( level <= g_verbose && g_verbose == 5 ){
	    std::cout << level_color[level] << "<LEVEL" << level << "> " << fmt << _color["white"] <<endl;
        // in case of qt output for log-browser
        } else if ( level <= g_verbose && qt_output != NULL && g_verbose == 7  ){
            qt_output->setText(QString::number((int)level)+QString::fromStdString(fmt));
        }     
        // in case of loglevel SAVE (e.g. for VASCC2016) use path to save a file
        if(level == 6 && !path.empty())
        {
            ofstream outstream ( path.c_str(), ios::out | ios::app ); 
            outstream << fmt;
            outstream << add_to_saveline;
            outstream.close();
        }
    }      
    
    template <typename T> Logger& operator %(T value) {
        ostringstream ss;
        ss << value;
        fmt += ss.str();
        return *this;
    }    
    
    protected:
        loglevel level;
        string fmt;
        string path;
    
    private:
        static map<string,string> _color;
};   
}

//easy use of logger
template <loglevel level>
ivg::Logger log(string msg) {
    return ivg::Logger( level, msg );
}
  
//easy use of save-logger (loglevel = 6)
template <loglevel level>
ivg::Logger log(string path, string msg) {
    return ivg::Logger( level, path, msg );
}

#endif	/* LOGGER_H */

