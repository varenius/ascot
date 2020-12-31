# include "logger.h"
string ivg::Logger::add_to_saveline = "";
QLineEdit* ivg::Logger::qt_output = NULL;

// ANSI escape codes
std::map<string,string> ivg::Logger::_color = map<string,string>({
                                                                 {"boldred",   "\x1b[31;1m"},
                                                                 {"red",       "\x1b[31m"},
                                                                 {"boldwhite", "\x1b[37;1m"},
                                                                 {"white",     "\x1b[37;0m"},
                                                                 {"boldgreen", "\x1b[32;1m"},
                                                                 {"green",     "\x1b[32m"},
                                                                 {"blue",      "\x1b[34m"},
                                                                 {"boldblue",  "\x1b[34;1m"}
                                                                });

 // color for each log level
std::array<string,8> ivg::Logger::level_color = std::array<string,8>({ _color["white"],
                                                                       _color["boldred"],
                                                                       _color["boldblue"],
                                                                       _color["white"],
                                                                       _color["white"],
                                                                       _color["white"],
                                                                       _color["white"],
                                                                       _color["white"]
                                                                       }); // white
                                                             
