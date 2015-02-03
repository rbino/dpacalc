/*
Copyright (C) 2014 Riccardo Binetti

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
#include "math.h"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/info_parser.hpp"
#include "boost/foreach.hpp"
#include "unsupported/Eigen/FFT"
#include <iostream>
#include <fstream>
#include <utility>


using namespace Eigen;
using namespace boost::property_tree;

namespace FilterFind
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

    bool filterParamCompare(const filterParam &firstElem, const filterParam &secondElem);

    enum padType {
        ZERO,
        MEAN,
        HOLD
    };

    class fftfilter_nary : public base
    {
        public:
             /**
             * @brief fftfilter_nary Creates an instance of the fftfilter_nary class
             * @param cmd The TCLAP::CmdLine which is used to add command line arguments
             */
            fftfilter_nary ( TCLAP::CmdLine& cmd ) :
                base ( cmd ),
                filterConfArg( "c", "filter-conf", "Filter search configuration filename", true, "", "path")
                { cmd.add(filterConfArg);
                };
            virtual void init(unsigned long long _samplespertrace, unsigned long long _numtraces);
            /**
             * @brief applyFilter Applies the member filter to the trace (coefficient wise multiplication)
             * @param tracewd shared_ptr to the TraceWithData you want to apply the filter to (only the Trace is affected)
             */
            void applyFilter(shared_ptr<TraceWithData> &tracewd);
            /**
             * @brief initFilterOutput Initializes the file or buffer (depending on the config) in which the filtered input is written
             */
            void initFilterOutput();
            /**
             * @brief writeFilteredTrace Writes the filtered trace to the output file or buffer (depending on the config)
             * @param tracewd The trace you want to write to the output
             * @param id The progressive id of the trace (to calculate the offset)
             */
            void writeFilteredTrace(shared_ptr<TraceWithData> tracewd, unsigned int id);
            /**
             * @brief getFilteredPointer Returns the void* pointer to the mmaped file or the buffer (depending on the config) which points to the filtered input and sets the new size
             * @param newsize This variable (passed by reference) is set to the new size of the input
             * @return the void* pointer to the mmaped file or the buffer (depending on the config) which points to the filtered input
             */
            void* getFilteredPointer(unsigned int& newsize);

            bool hasFinished();

            void setBaseline(unsigned long steps);

            FilterBand beginStep();

            virtual void endStep(unsigned int successTraces);

            virtual bool isLastStep();

            virtual void end();

        protected:

            TCLAP::ValueArg<std::string> filterConfArg;
            ifstream config;
            ofstream configGoodOut;
            ofstream configGoodUglyOut;
            deque<filterParam> toBeProcessed;
            deque<filterParam> alreadyProcessed;
            deque<filterParam> goodAndUglyFilters;
            deque<filterParam> goodFilters;
            filterParam currentFilter;
            mutex writeMutex;
            padType padtype;
            char padconf;
            unsigned int currentStep;
            windowShape filterShape;
            float filterTukeyAlpha;
            unsigned int filterSteps;
            unsigned int filterDivisions;
            float overlap;
            shared_ptr<Trace> filter;
            unsigned long baseline;
            TraceValueType padding;
            unsigned int maxBadDepth;
            /**
             * @brief nextPow2 Calculates the next power of 2 of n
             * @param n The number which we want to increase to the next power of 2
             * @return the next power of 2 of n
             */
            unsigned long long nextPow2(unsigned long long n){
                return pow(2, ceil(log2(n)));
            }
            unsigned long long fftLength;
            double fNyq;
            /**
             * @brief generateWindows Generates filter Windows and puts it in the filter
             * @param filt shared_ptr to the filter to be generated
             * @param parameters Vector of filterParam for the generation of the windows
             */
            void generateWindows(shared_ptr<Trace>& filt, filterParam &windowParam);
            /**
             * @brief debugPrint Utility function to print the values of a Trace (filter or trace) in a gnuplot friendly format
             * @param trace shared_ptr to the trace you want to print
             * @param filename Output filename
             */
            void debugPrint(shared_ptr<Trace>& trace, string filename);
            /**
             * @brief debugPrint Utility function to print the absolute values of a ComplexTrace (fft) in a gnuplot friendly format
             * @param trace shared_ptr to the ComplexTrace you want to print
             * @param filename Output filename
             */
            void debugPrint(shared_ptr<ComplexTrace>& trace, string filename);
            /**
             * @brief combineFilter Combines (depending on the config) the filter bin in position pos with the calculated windowValue
             * @param pos The index of the filter bin which you want to combine with the window value
             * @param windowValue The calculated window value
             */
            void combineFilter(unsigned long pos, TraceValueType windowValue);
            TraceValueType maxBin;
            /**
             * @brief getTraceOffset Returns the trace offset in the output buffer (or file)
             * @param trace The progressive id of the trace
             * @return The offset corresponding to the beginning of the trace with the passed id
             */
            unsigned long long getTraceOffset ( unsigned long long trace ) {
                //trace and samplenum are zero-based
                return sizeof ( struct fileheaders ) + trace * ( sizeof(TraceValueType) * SamplesPerTrace + DATA_SIZE_BYTE );
            }
            /**
             * @brief getDataOffset Returns the data offset in the output buffer (or file)
             * @param trace The progressive id of the trace
             * @return The offset corresponding to the beginning of the data with the passed id
             */
            unsigned long long getDataOffset ( unsigned long long trace ) {
                //trace and samplenum are zero-based
                return sizeof ( struct fileheaders ) + trace * ( sizeof(TraceValueType) * SamplesPerTrace + DATA_SIZE_BYTE ) + sizeof(TraceValueType) * SamplesPerTrace;
            }
            /**
             * @brief getBufferDimension Calculates the buffer dimension to allocate the output buffer (or file)
             * @return The buffer (or file) dimension
             */
            unsigned long long getBufferDimension () {
                //trace and samplenum are zero-based
                return sizeof ( struct fileheaders ) + NumTraces * ( sizeof(TraceValueType) * SamplesPerTrace + DATA_SIZE_BYTE );
            }

            void writeConfig(ofstream& outfile, deque<filterParam> params, char padding, double fsampling);

            void* outBuffer;

    };
}
