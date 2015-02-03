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
#include "fftfilter.hpp"
#include <sys/mman.h>
#include <fcntl.h>

void Filters::fftfilter::init(unsigned long long _samplespertrace, unsigned long long _numtraces){
    SamplesPerTrace = _samplespertrace;
    NumTraces = _numtraces;
    config.open(filterConfArg.getValue().c_str());
    if (!config.is_open()){
        if(!filterConfArg.isSet()){
            cerr << "If you filter (-i), you must specify a filter configuration file (-c filename)" << endl;
        } else {
            cerr << "Cannot open filter configuration " << filterConfArg.getValue() << endl;
        }
        exit(3);
    }
#if defined(CONFIG_FILTER_OUTPUT_DISK)
    if (!filterOutFileConfArg.isSet()){
        cerr << "Filter output file not set. Specify it with -t" << endl;
        exit(3);
    }
#endif
    try
    {
        ptree pt;
        info_parser::read_info(config, pt);
        try
        {
            SamplingFrequency = pt.get<float>("samplingfrequency");
            fNyq = SamplingFrequency/2;
            char pad = pt.get<char>("paddingtype");
            switch(pad){
                case 'z':
                    padtype = ZERO;
                    padding = 0;
                    break;
                case 'm':
                    padtype = MEAN;
                    break;
                case 'h':
                    padtype = HOLD;
                    break;
                default:
                    cerr << "paddingtype must be z (zero), m (mean) or h (hold last value)" << endl;
                    exit(3);
            }

            ptree filters (pt.get_child("filters"));
            ptree::const_iterator itFilt;
            for ( itFilt = filters.begin(); itFilt != filters.end(); ++itFilt ){
                filterParam param;
                char shape = itFilt->second.get<char>("window");
                switch (shape){
                    case 'r':
                        param.shape = RECT;
                        break;
                    case 'h':
                        param.shape = HAMMING;
                        break;
                    case 'H':
                        param.shape = HANN;
                        break;
                    case 't':
                        param.shape = TUKEY;
                        param.tukeyAlpha = itFilt->second.get<double>("alpha", 0.5);
                        if (param.tukeyAlpha <= 0){
                            cerr << "alpha must be > 0" << endl;
                            exit(3);
                        }
                        break;
                    default:
                        cerr << "Filter configuration error: Window shape must be r, h, H or t" << endl;
                        exit ( 3 );
                }
                string type = itFilt->second.get<string>("type");
                if (type.compare("bp") == 0) {
                    param.freq1 = itFilt->second.get<double>("frequencies.low");
                    param.freq2 = itFilt->second.get<double>("frequencies.high");
                    if (param.freq1 > param.freq2){
                        cerr << "Low frequency must be lower than high frequency" << endl;
                        exit (3);
                    }
                } else if (type.compare("lp") == 0) {
                    param.freq2 = itFilt->second.get<double>("frequencies.high");
                    param.freq1 = -param.freq2;
                } else if (type.compare("hp") == 0) {
                    param.freq1 = itFilt->second.get<double>("frequencies.low");
                    param.freq2 = SamplingFrequency - param.freq1;
                } else {
                    cerr << "Filter configuration error: Filter type must be bp, lp or hp" << endl;
                    exit (3);
                }
                filterParamVect.push_back(param);
            }
        } catch( ptree_error e) {
            cerr << "Filter configuration error. Check the config file" << endl;
            exit (3);
        }
    } catch (info_parser::info_parser_error e) {
        cerr << "Cannot parse filter configuration" << endl;
        exit ( 3 );
    }
    maxBin = 1;
    /* FFT padded to the next power of 2 to avoid boundary effects and to speed up computation.
     * SamplesPerTrace+1 takes care of traces whose length is already a power of 2 */
    fftLength = nextPow2(SamplesPerTrace+1);
    filter.reset(new Trace(fftLength));
    *filter = Trace::Zero(fftLength);
    generateWindows(filter, filterParamVect);
}

void Filters::fftfilter::applyFilter(shared_ptr<TraceWithData>& tracewd){
    FFT<TraceValueType> fft;
    shared_ptr<ComplexTrace> freqvec (new ComplexTrace(fftLength));
    unsigned long padLength = fftLength - SamplesPerTrace;
    if (padtype == MEAN){
        padding = (tracewd->trace)->mean();
    } else if (padtype == HOLD){
        padding = (*tracewd->trace)(SamplesPerTrace-1);
    }
    Trace padTrace = Trace::Constant(padLength, padding);
    Trace paddedTrace(fftLength);
    paddedTrace << *(tracewd->trace),padTrace;
    *(tracewd->trace) = paddedTrace;
    fft.fwd(*freqvec, *(tracewd->trace));
    if (freqvec->size() == filter->size()){
        (*freqvec) = freqvec->cwiseProduct(*filter);
    } else {
        cerr << "FFT length not equal to filter length" << endl;
        exit(3);
    }
    fft.inv(*(tracewd->trace), *freqvec);
    /* Back to the original size, remove padding */
    (tracewd->trace)->conservativeResize(SamplesPerTrace);
}

void Filters::fftfilter::debugPrint(shared_ptr<Trace>& trace, string filename){
    ofstream debug(filename);
    if (debug.is_open()){
        for (int i=0; i < trace->size(); i++){
            debug << i << "\t" << (*trace) (i) << endl;
        }
    }
    debug.close();
}

void Filters::fftfilter::debugPrint(shared_ptr<ComplexTrace>& trace, string filename){
    ofstream debug(filename);
    if (debug.is_open()){
        for (int i=0; i < trace->size(); i++){
            std::complex<TraceValueType> tmp = (*trace)(i);
            debug << i << "\t" << abs(tmp) << endl;
        }
    }
    debug.close();
}

void Filters::fftfilter::generateWindows(shared_ptr<Trace>& filt, vector<filterParam>& parameters){
    for (vector<filterParam>::iterator windowParam = parameters.begin(); windowParam != parameters.end(); ++windowParam){
        int nBins = (int) ceil((windowParam->freq2 - windowParam->freq1)/fNyq * fftLength/2);
        Trace window(nBins);
        if (nBins == 1){
            window(1) = 1;
        } else {
            switch (windowParam->shape){
                case RECT:
                    for (int n=0; n < nBins; n++){
                        window(n) = 1;
                    }
                    break;
                case HAMMING:
                    for (int n=0; n < nBins; n++){
                        window(n) = 0.54 - 0.46 * cos( (2*M_PI*n) / ( nBins-1 ) );
                    }
                    break;
                case HANN:
                    for (int n=0; n < nBins; n++){
                        window(n) = 0.5 * (1 - cos(2*M_PI*n/(nBins-1)));
                    }
                    break;
                case TUKEY:
                    for (int n=0; n < nBins; n++)
                        if (n <= windowParam->tukeyAlpha*(nBins-1)/2){
                            window(n) = 0.5 * ( 1 + cos(M_PI * ((2*n/(windowParam->tukeyAlpha*(nBins-1)))-1)));
                        } else if (n <= (nBins-1)*(1-(windowParam->tukeyAlpha/2))) {
                            window(n) = 1;
                        }   else {
                            window(n) = 0.5 * ( 1 + cos(M_PI * ((2*n/(windowParam->tukeyAlpha*(nBins-1)))- (2/windowParam->tukeyAlpha) + 1)));
                        }
                    break;
            }
        }
        int startBin = (int) ceil((windowParam->freq1 * fftLength/2) / fNyq);
        for (int n=0; n < nBins; n++){
            combineFilter((startBin+n)%fftLength, window(n));   //Modulo fftlength so the bins are placed circularly in the right position
        }
    }
#if defined(CONFIG_FILTER_COMBINE_NORMALIZE)
    *filter/=maxBin;
#endif
    unsigned long nyquistBin = filter->size() / 2;
    filter->conservativeResize(nyquistBin);
    Trace bilateralFilter(nyquistBin*2);
    bilateralFilter << (*filter),filter->reverse();
    (*filter) = bilateralFilter;
}

void Filters::fftfilter::combineFilter(unsigned long pos, TraceValueType windowValue){
    (*filter) (pos) += windowValue;
#if defined(CONFIG_FILTER_COMBINE_NORMALIZE)
    if ((*filter)(pos) > maxBin){
        /* Check the highest value in the filter to normalize */
        maxBin = (*filter)(pos);
    }
#elif defined(CONFIG_FILTER_COMBINE_CLAMP)
    if ((*filter) (pos) > 1){
        (*filter) (pos) = 1;
    }
#endif
}

void Filters::fftfilter::initFilterOutput(){
#if defined(CONFIG_FILTER_OUTPUT_DISK)

    outFd = open (filterOutFileConfArg.getValue().c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRWXU | S_IRGRP | S_IROTH);
    if (outFd == -1){
        cerr << "Could not open filter output file" << endl;
        exit(3);
    }
    int seekRes;
    int writeRes;
    seekRes = lseek(outFd, getBufferDimension() - 1, SEEK_SET);
    writeRes = write(outFd, "", 1);
    if (writeRes == -1 || seekRes == -1){
        cerr << "Error inflating the file" << endl;
        exit(3);
    }
    outBuffer = mmap(NULL, getBufferDimension(), PROT_READ | PROT_WRITE, MAP_SHARED, outFd, 0);
#elif defined(CONFIG_FILTER_OUTPUT_RAM)
    outBuffer = malloc(getBufferDimension());
    if (outBuffer ==  NULL){
        cout << "Out of memory" << endl;
        exit(3);
    }
#endif
    struct fileheaders header;
    header.numtraces = NumTraces;
    header.numsamples_per_trace = SamplesPerTrace;
#if defined(CONFIG_TRACETYPE_FLOAT)
    header.datatype = 'f';
#elif defined(CONFIG_TRACETYPE_DOUBLE)
    header.datatype = 'd';
#endif
    header.knowndatalength = DATA_SIZE_BYTE;
    memcpy(outBuffer, &header, sizeof(struct fileheaders));
}

void Filters::fftfilter::writeFilteredTrace(shared_ptr<TraceWithData> tracewd, unsigned int id){
    writeMutex.lock();
    memcpy((char*)outBuffer + getTraceOffset(id), (tracewd->trace)->data(), (tracewd->trace)->size() * sizeof(TraceValueType) );
    BitsetToBuffer<KEY_SIZE_BYTE>((*tracewd->data), (char*) outBuffer + getDataOffset(id));
    writeMutex.unlock();
}

void* Filters::fftfilter::getFilteredPointer(unsigned int& newsize){
    newsize = getBufferDimension();
    return outBuffer;
}

