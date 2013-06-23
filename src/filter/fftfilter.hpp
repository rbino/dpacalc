/*
Copyright (C) 2013 Riccardo Binetti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
*/
#pragma once
#include "dpacalc.h"
#include "base.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/info_parser.hpp"
#include "boost/foreach.hpp"
#include "unsupported/Eigen/FFT"
#include <iostream>
#include <fstream>
#define TUKEY_ALPHA 0.5

enum windowShape {
    RECT,
    HAMMING,
    HANN,
    TUKEY
};

typedef struct {
    windowShape shape;
    double freq1;
    double freq2;
} filterParam;

using namespace Eigen;
using namespace boost::property_tree;

namespace Filters
{
    class fftfilter : public base
    {
        public:
            fftfilter ( TCLAP::CmdLine& cmd, shared_ptr<SamplesInput::base> _input) :
                base ( cmd, _input ),
                filterConfArg ( "c", "filter-conf", "Filter configuration filename", false, "", "path")
                { cmd.add(filterConfArg);
                };
            virtual void init();
            void applyFilter(shared_ptr<TraceWithData>& tracewd);

        protected:
            TCLAP::ValueArg<std::string> filterConfArg;
            ifstream config;
            vector<filterParam> filterParamVect;
            void tokenize(const string& str, vector<string>& tokens, const string& delimiters);
            vector<TraceValueType> filter;
            void initializeToZero(vector<TraceValueType>& filt, unsigned long long length);
            unsigned long long nextPow2(unsigned long long n){
                return pow(2, ceil(log2(n)));
            }
            unsigned long long fftLength;
            double fNyq;
            void generateWindows(vector<TraceValueType>& filt, vector<filterParam>& parameters);
            void debugPrint(Trace& trace, string filename);
    };
}
