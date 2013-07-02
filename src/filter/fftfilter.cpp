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

void Filters::fftfilter::init(){
    filter.reset(new Trace(SamplesPerTrace));
    config.open(filterConfArg.getValue().c_str());
    try
    {
        ptree pt;
        info_parser::read_info(config, pt);
        try
        {
            SamplingFrequency = pt.get<float>("samplingfrequency");
            fNyq = SamplingFrequency/2;
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
        cerr << "Cannot open filter configuration file" << endl;
        exit ( 3 );
    }
    maxBin = 1;
    fftLength = nextPow2(input->SamplesPerTrace);
    initializeToZero(filter, fftLength);
    generateWindows(filter, filterParamVect);
    debugPrint(*filter, "/home/rbino/dpaoutput/filterDebug");
}

void Filters::fftfilter::initializeToZero(shared_ptr<Trace>& trace, unsigned long long length){
    for (unsigned long long i=0; i < length; i++){
        trace->push_back(0);
    }
}

void Filters::fftfilter::applyFilter(shared_ptr<TraceWithData>& tracewd){
    FFT<TraceValueType> fft;
    shared_ptr<ComplexTrace> freqvec (new ComplexTrace(fftLength));
    unsigned long zeroPadLength = fftLength - input->SamplesPerTrace;
    shared_ptr<Trace> zeroPad (new Trace(zeroPadLength));
    initializeToZero(zeroPad, zeroPadLength);
    debugPrint(tracewd->trace, "/home/rbino/dpaoutput/preFftTrace");
    tracewd->trace.insert(tracewd->trace.end(), zeroPad->begin(), zeroPad->end());
    fft.fwd(*freqvec, tracewd->trace);
    if (freqvec->size() == filter->size()){
        for(unsigned int i=0; i<freqvec->size(); i++){
            (*freqvec) [i] = (*freqvec) [i] * (*filter) [i];
        }
    } else {
        cerr << "FFT length not equal to filter length" << endl;
        exit(3);
    }
    fft.inv(tracewd->trace, *freqvec);
    tracewd->trace.resize(input->SamplesPerTrace);
    debugPrint(tracewd->trace, "/home/rbino/dpaoutput/postFftTrace");
}

void Filters::fftfilter::debugPrint(Trace& trace, string filename){
    int counter = 0;
    ofstream debug(filename);
    if (debug.is_open()){
        BOOST_FOREACH(TraceValueType sample, trace)
        {
         debug << counter << "\t" << sample << endl;
         counter++;
        }
    }
}

void Filters::fftfilter::generateWindows(shared_ptr<Trace>& filt, vector<filterParam>& parameters){
    for (vector<filterParam>::iterator windowParam = parameters.begin(); windowParam != parameters.end(); ++windowParam){
        int nBins = (int) ceil((windowParam->freq2 - windowParam->freq1)/fNyq * fftLength/2);
        vector<TraceValueType> window;
        if (nBins == 1){
            window.push_back(1);
        } else {
            switch (windowParam->shape){
                case RECT:
                    for (int n=0; n < nBins; n++){
                        window.push_back(1);
                    }
                    break;
                case HAMMING:
                    for (int n=0; n < nBins; n++){
                        window.push_back( 0.54 - 0.46 * cos( (2*M_PI*n) / ( nBins-1 ) ) );
                    }
                    break;
                case HANN:
                    for (int n=0; n < nBins; n++){
                        window.push_back( 0.5 * (1 - cos(2*M_PI*n/(nBins-1))));
                    }
                    break;
                case TUKEY:
                    for (int n=0; n < nBins; n++)
                        if (n <= windowParam->tukeyAlpha*(nBins-1)/2){
                            window.push_back( 0.5 * ( 1 + cos(M_PI * ((2*n/(windowParam->tukeyAlpha*(nBins-1)))-1))));
                        } else if (n <= (nBins-1)*(1-(windowParam->tukeyAlpha/2))) {
                            window.push_back(1);
                        }   else {
                            window.push_back(0.5 * ( 1 + cos(M_PI * ((2*n/(windowParam->tukeyAlpha*(nBins-1)))- (2/windowParam->tukeyAlpha) + 1))));
                        }
                    break;
            }
        }
        int startBin = (int) ceil((windowParam->freq1 * fftLength/2) / fNyq);
        for (int n=0; n < nBins; n++){
            combineFilter(startBin+n%fftLength, window[n]);
        }

    }
#if defined(CONFIG_FILTER_COMBINE_NORMALIZE)
    for (unsigned long i=0; i<filter.size(); i++){
        filter[i] = filter[i] / maxBin;
    }
#endif
    std::size_t const nyquistBin = filter->size() / 2;
    std::vector<TraceValueType> lower_half(filter->begin(), filter->begin() + nyquistBin + 1);
    std::vector<TraceValueType> upper_half(lower_half);
    upper_half.pop_back();
    reverse(upper_half.begin(), upper_half.end());
    upper_half.pop_back();
    lower_half.insert(lower_half.end(), upper_half.begin(), upper_half.end());
    *filter = lower_half;
}

void Filters::fftfilter::combineFilter(unsigned long pos, TraceValueType windowValue){
    (*filter) [pos] += windowValue;
#if defined(CONFIG_FILTER_COMBINE_NORMALIZE)
    if (filter[pos] > maxBin){
        maxBin = filter[pos];
    }
#elif defined(CONFIG_FILTER_COMBINE_CLAMP)
    if ((*filter) [pos] > 1){
        (*filter) [pos] = 1;
    }
#endif
}

