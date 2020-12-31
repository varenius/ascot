#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <vector>
#include <sstream>
#include <string>

namespace lps{
    
using Seconds = double;

const Seconds DurationDay = 24*60*60;

const double mjd_J2000 = 51544.0; // J2000 in mjd

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> tokens;
    while (std::getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}


template <class Iter>
class range {
    Iter b;
    Iter e;
public:

    range(Iter b, Iter e) : b(b), e(e) {}

    Iter begin() { return b; }
    Iter end() { return e; }
    Iter cbegin() const {return b;}
    Iter cend()   const {return e;}
};

template <class Container>
range<typename Container::iterator>
make_range(Container& c, std::size_t b, std::size_t e) {
    return range<typename Container::iterator> (c.begin()+b, c.begin()+e);
}

template <class Container>
range<typename Container::const_iterator>
make_const_range(Container& c, std::size_t b, std::size_t e){
    return range<typename Container::const_iterator> (c.cbegin()+b, c.cbegin()+e);
}

}

#endif // DEFINITIONS_H
