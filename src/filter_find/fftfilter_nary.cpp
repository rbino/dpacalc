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
#include "fftfilter_nary.hpp"
#include <sys/mman.h>
#include <fcntl.h>

void FilterFind::fftfilter_nary::init(){
    SamplesPerTrace = input->SamplesPerTrace;
    NumTraces = input->RealNumTraces;
    config.open(filterConfArg.getValue().c_str());
    if (!config.is_open()){
        cerr << "Cannot open filter configuration " << filterConfArg.getValue() << endl;
        exit(3);
    }
    try
    {
        ptree pt;
        info_parser::read_info(config, pt);
        try
        {
            SamplingFrequency = pt.get<float>("samplingfrequency");
            fNyq = SamplingFrequency/2;
            padconf = pt.get<char>("paddingtype");
            switch(padconf){
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
            char shape = pt.get<char>("window");
            switch (shape){
                case 'r':
                    filterShape = RECT;
                    break;
                case 'h':
                    filterShape = HAMMING;
                    break;
                case 'H':
                    filterShape = HANN;
                    break;
                case 't':
                    filterShape = TUKEY;
                    filterTukeyAlpha = pt.get<double>("alpha", 0.5);
                    if (filterTukeyAlpha <= 0){
                        cerr << "alpha must be > 0" << endl;
                        exit(3);
                    }
                    break;
                default:
                    cerr << "Filter configuration error: Window shape must be r, h, H or t" << endl;
                    exit ( 3 );
            }
            filterSteps = pt.get<unsigned int>("steps");
            filterDivisions = pt.get<unsigned int>("divisions");
            overlap = pt.get<float>("overlap");
            maxBadDepth = pt.get<unsigned int>("maxbaddepth");
            configGoodOut.open(pt.get<string>("goodout"));
            if (!configGoodOut.is_open()){
                cerr << "Cannot open good filter configuration output " << endl;
                exit(3);
            }
            configGoodUglyOut.open(pt.get<string>("gooduglyout"));
            if (!configGoodUglyOut.is_open()){
                cerr << "Cannot open ugly+good filter configuration output " << endl;
                exit(3);
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
    filterParam init;
    init.freq1 = 0;
    init.freq2 = fNyq;
    init.shape = filterShape;
    init.tukeyAlpha = filterTukeyAlpha;
    alreadyProcessed.push_back(init);

}

void FilterFind::fftfilter_nary::applyFilter(shared_ptr<TraceWithData> &tracewd){
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

void FilterFind::fftfilter_nary::debugPrint(shared_ptr<Trace>& trace, string filename){
    ofstream debug(filename);
    if (debug.is_open()){
        for (int i=0; i < trace->size(); i++){
            debug << i << "\t" << (*trace) (i) << endl;
        }
    }
    debug.close();
}

void FilterFind::fftfilter_nary::debugPrint(shared_ptr<ComplexTrace>& trace, string filename){
    ofstream debug(filename);
    if (debug.is_open()){
        for (int i=0; i < trace->size(); i++){
            std::complex<TraceValueType> tmp = (*trace)(i);
            debug << i << "\t" << abs(tmp) << endl;
        }
    }
    debug.close();
}

void FilterFind::fftfilter_nary::generateWindows(shared_ptr<Trace>& filt, filterParam& windowParam){
    filt.reset(new Trace(fftLength));
    *filt = Trace::Zero(fftLength);
    int nBins = (int) ceil((windowParam.freq2 - windowParam.freq1)/fNyq * fftLength/2);
    Trace window(nBins);
    if (nBins == 1){
        window(1) = 1;
    } else {
        switch (windowParam.shape){
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
                    if (n <= windowParam.tukeyAlpha*(nBins-1)/2){
                        window(n) = 0.5 * ( 1 + cos(M_PI * ((2*n/(windowParam.tukeyAlpha*(nBins-1)))-1)));
                    } else if (n <= (nBins-1)*(1-(windowParam.tukeyAlpha/2))) {
                        window(n) = 1;
                    }   else {
                        window(n) = 0.5 * ( 1 + cos(M_PI * ((2*n/(windowParam.tukeyAlpha*(nBins-1)))- (2/windowParam.tukeyAlpha) + 1)));
                    }
                break;
        }
    }
    int startBin = (int) ceil((windowParam.freq1 * fftLength/2) / fNyq);
    for (int n=0; n < nBins; n++){
        combineFilter((startBin+n)%fftLength, window(n));   //Modulo fftlength so the bins are placed circularly in the right position
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

void FilterFind::fftfilter_nary::combineFilter(unsigned long pos, TraceValueType windowValue){
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

void FilterFind::fftfilter_nary::initFilterOutput(){
    outBuffer = malloc(getBufferDimension());
    if (outBuffer ==  NULL){
        cout << "Out of memory" << endl;
        exit(3);
    }
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

void FilterFind::fftfilter_nary::writeFilteredTrace(shared_ptr<TraceWithData> tracewd, unsigned int id){
    writeMutex.lock();
    memcpy((char*)outBuffer + getTraceOffset(id), (tracewd->trace)->data(), (tracewd->trace)->size() * sizeof(TraceValueType) );
    BitsetToBuffer<KEY_SIZE_BYTE>((*tracewd->data), (char*) outBuffer + getDataOffset(id));
    writeMutex.unlock();
}

void* FilterFind::fftfilter_nary::getFilteredPointer(unsigned int& newsize){
    newsize = getBufferDimension();
    return outBuffer;
}

bool FilterFind::fftfilter_nary::hasFinished(){
    if (currentStep >= filterSteps && toBeProcessed.empty()){
        return true;
    } else {
        return false;
    }
}

void FilterFind::fftfilter_nary::setBaseline(unsigned long steps){
    baseline = steps;
}

FilterBand FilterFind::fftfilter_nary::beginStep(){

    if (toBeProcessed.empty()){
        currentStep++;
        while (!alreadyProcessed.empty()){
            filterParam iter = alreadyProcessed.front();
            alreadyProcessed.pop_front();
            int winLen = iter.freq2 - iter.freq1;
            int curLo = iter.freq1;
            int subWinLen = winLen / ( (filterDivisions) - (filterDivisions-1) * overlap);
            for (unsigned int i = 0; i < filterDivisions; i++){
                filterParam filt;
                filt.shape = filterShape;
                filt.tukeyAlpha = filterTukeyAlpha;
                filt.freq1 = floor(curLo);
                filt.freq2 = ceil(curLo + subWinLen);
                curLo = curLo + subWinLen - (subWinLen * overlap);
                toBeProcessed.push_back(filt);
            }
        }
    }
    currentFilter = toBeProcessed.front();
    toBeProcessed.pop_front();
    generateWindows(filter,currentFilter);
    FilterBand curBand;
    curBand.first = currentFilter.freq1;
    curBand.second = currentFilter.freq2;
    return curBand;
}

void FilterFind::fftfilter_nary::endStep(unsigned int successTraces){
    if (successTraces < baseline){
        if (currentStep < filterSteps){
            alreadyProcessed.push_back(currentFilter);
        } else {
            goodFilters.push_back(currentFilter);
            goodAndUglyFilters.push_back(currentFilter);
        }
    } else if (successTraces < NumTraces) {
        if (currentStep < filterSteps) {
            alreadyProcessed.push_back(currentFilter);
        } else {
            goodAndUglyFilters.push_back(currentFilter);
        }
    } else {
        if (currentStep <= maxBadDepth) {
            alreadyProcessed.push_back(currentFilter);
        }
    }
}

bool FilterFind::fftfilter_nary::isLastStep(){
    if (currentStep == filterSteps){
        return true;
    } else {
        return false;
    }
}

bool FilterFind::filterParamCompare(const filterParam& firstElem, filterParam& secondElem) {
    return firstElem.freq1 < secondElem.freq1;
}


void FilterFind::fftfilter_nary::writeConfig(ofstream& outfile, deque<filterParam> params, char padding, double fsampling){
    ptree pt;
    int count = 0;
    filterParam cur;
    string curFilt;
    pt.put("paddingtype", padding);
    pt.put("samplingfrequency", fsampling);
    sort(params.begin(),params.end(), filterParamCompare);
    while (!params.empty()){
        cur = params.front();
        params.pop_front();
        while ((cur.freq2 > params.front().freq1) && !params.empty()){
            cur.freq2 = params.front().freq2;
            params.pop_front();
        }
        curFilt = "f" + to_string(count);
        pt.put("filters."+ curFilt + ".type", "bp");
        pt.put("filters."+ curFilt + ".frequencies.low", cur.freq1);
        pt.put("filters."+ curFilt + ".frequencies.high", cur.freq2);
        switch(cur.shape){
            case RECT:
                pt.put("filters."+ curFilt + ".window", "r");
                break;
            case HAMMING:
                pt.put("filters."+ curFilt + ".window", "h");
                break;
            case HANN:
                pt.put("filters."+ curFilt + ".window", "H");
                break;
            case TUKEY:
                pt.put("filters."+ curFilt + ".window", "t");
                pt.put("filters."+ curFilt + ".alpha", cur.tukeyAlpha);
                break;
        }
        count++;
    }
    write_info(outfile,pt);
}


void FilterFind::fftfilter_nary::end(){
    if (goodFilters.empty()){
        cout << "No good filters. Writing empty filter configuration." << endl;
    } else {
        writeConfig(configGoodOut, goodFilters, padconf, SamplingFrequency);
    }
    if (goodAndUglyFilters.empty()){
        cout << "No good and ugly filters. Writing empty filter configuration." << endl;
    } else {
        writeConfig(configGoodUglyOut, goodAndUglyFilters, padconf, SamplingFrequency);
    }
    configGoodOut.close();
    configGoodUglyOut.close();
}
