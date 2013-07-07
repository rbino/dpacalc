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



using namespace Eigen;
using namespace boost::property_tree;

namespace Filters
{
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
        double tukeyAlpha;
    } filterParam;

    enum padType {
        ZERO,
        MEAN,
        HOLD
    };

    class fftfilter : public base
    {
        public:
            fftfilter ( TCLAP::CmdLine& cmd, shared_ptr<SamplesInput::base> _input) :
                base ( cmd, _input ),
#if defined(CONFIG_FILTER_OUTPUT_DISK)
                filterOutFileConfArg("t", "filter-output", "Filter output file", false, "", "path"),
#endif
                filterConfArg( "c", "filter-conf", "Filter configuration filename", false, "", "path")
                { cmd.add(filterConfArg);
#if defined(CONFIG_FILTER_OUTPUT_DISK)
                  cmd.add(filterOutFileConfArg);
#endif
                };
            virtual void init();
            void applyFilter(shared_ptr<TraceWithData>& tracewd);
            void initFilterOutput();
            void writeFilteredTrace(shared_ptr<TraceWithData> tracewd, unsigned int id);
            void* getFilteredPointer(unsigned int& newsize);

        protected:
#if defined(CONFIG_FILTER_OUTPUT_DISK)
            TCLAP::ValueArg<std::string> filterOutFileConfArg;
#endif
            TCLAP::ValueArg<std::string> filterConfArg;
            ifstream config;
            vector<filterParam> filterParamVect;
            shared_ptr<Trace> filter;
            mutex writeMutex;
            padType padtype;
            TraceValueType padding;
            void initializeToZero(shared_ptr<Trace>& filt, unsigned long long length);
            unsigned long long nextPow2(unsigned long long n){
                return pow(2, ceil(log2(n)));
            }
            unsigned long long fftLength;
            double fNyq;
            void generateWindows(shared_ptr<Trace>& filt, vector<filterParam>& parameters);
            void debugPrint(shared_ptr<Trace>& trace, string filename);
            void debugPrint(shared_ptr<ComplexTrace>& trace, string filename);
            void combineFilter(unsigned long pos, TraceValueType windowValue);
            TraceValueType maxBin;
            unsigned long long getSampleOffset ( unsigned long long trace, unsigned long long samplenum ) {
                //trace and samplenum are zero-based
                return sizeof ( struct fileheaders ) + trace * ( sizeof(TraceValueType) * SamplesPerTrace + DATA_SIZE_BYTE ) + sizeof(TraceValueType) * samplenum;
            }
            unsigned long long getDataOffset ( unsigned long long trace ) {
                //trace and samplenum are zero-based
                return sizeof ( struct fileheaders ) + trace * ( sizeof(TraceValueType) * SamplesPerTrace + DATA_SIZE_BYTE ) + sizeof(TraceValueType) * SamplesPerTrace;
            }
            unsigned long long getBufferDimension () {
                //trace and samplenum are zero-based
                return sizeof ( struct fileheaders ) + NumTraces * ( sizeof(TraceValueType) * SamplesPerTrace + DATA_SIZE_BYTE );
            }
#if defined(CONFIG_FILTER_OUTPUT_DISK)
            ofstream outFile;
            int outFd;
#endif
            void* outBuffer;

    };
}
