#ifndef TICTOC_H
#define TICTOC_H

#include <sys/time.h>
#include <iostream>
#include <iomanip>

using namespace std;

class tictoc
{
    public:
        void tic()
        {
            gettimeofday( &_start, 0 );
        }

        double toc()
        {
            gettimeofday( &_stop, 0 );
            // calculate duration
            // seconds
            _duration = _stop.tv_sec - _start.tv_sec;
            // microseconds
            _duration += ( _stop.tv_usec - _start.tv_usec ) / 1e6;
            return _duration;
        }

        void cerr_toc(string info="here")
        {
            cerr << "Absolute time required until " << info << ": " << toc() <<
                 " seconds." << endl;
        }

        void cerr_diff(string info="here")
        {
            double last = _duration;
            cerr << "Relative time required until " << info << ": " << fixed << std::setprecision(8) << toc()-last <<
                 " seconds." << endl;
        }

        double duration()
        {
            return _duration;
        }
        void reset()
        {
            _start.tv_sec = 0;
            _start.tv_usec = 0;

            _stop.tv_sec = 0;
            _stop.tv_usec = 0;

            _duration = 0.0;
        }

    private:
        timeval _start;
        timeval _stop;
        // duration:
        double _duration;
};

#endif // TICTOC_H
