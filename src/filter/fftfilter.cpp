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
    config.open(filterConfArg.getValue().c_str());
    if (config.is_open()){
        config >> SamplingFrequency;
        fNyq = SamplingFrequency / 2;
        string filterParamStr;
        while (config >> filterParamStr){
            vector<string> tokens;
            filterParam param;
            tokenize(filterParamStr, tokens, ":");
            switch (tokens[0][0]){
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
                    break;
                default:
                    cerr << "Filter configuration error: Window shape must be r, h, H or t" << endl;
                    exit ( 3 );
            }
            float tmp;
            tmp = ::atof(tokens[1].c_str());
            param.freq1 = tmp;
            if (tokens[2].compare("LP") == 0){
                param.freq2 = tmp;
                param.freq1 = - param.freq2;
            } else if (tokens[2].compare("HP") == 0){
                param.freq2 = SamplingFrequency - param.freq1;
            } else {
                tmp = ::atof(tokens[2].c_str());
                param.freq2 = tmp;
            }
            filterParamVect.push_back(param);
       }
    } else {
        cerr << "Cannot open filter configuration file" << endl;
        exit ( 3 );
    }
    fftLength = nextPow2(input->SamplesPerTrace);
    initializeFilter(filter, fftLength);
    generateWindows(filter, filterParamVect);

}

void Filters::fftfilter::tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void Filters::fftfilter::initializeFilter(vector<TraceValueType>& filt, unsigned long long length){
    for (unsigned long long i=0; i < length; i++){
        filt.push_back(0);
    }
}

void Filters::fftfilter::generateWindows(vector<TraceValueType>& filt, vector<filterParam>& parameters){
    for (vector<filterParam>::iterator windowParam = parameters.begin(); windowParam != parameters.end(); ++windowParam){
        int nBins = (int) ceil((windowParam->freq2 - windowParam->freq1)/fNyq * fftLength/2);
        vector<TraceValueType> window;
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
                    if (n <= TUKEY_ALPHA*(nBins-1)/2){
                        window.push_back( 0.5 * ( 1 + cos(M_PI * ((2*n/(TUKEY_ALPHA*(nBins-1)))-1))));
                    } else if (n <= (nBins-1)*(1-(TUKEY_ALPHA/2))) {
                        window.push_back(1);
                    }   else {
                        window.push_back(0.5 * ( 1 + cos(M_PI * ((2*n/(TUKEY_ALPHA*(nBins-1)))- (2/TUKEY_ALPHA) + 1))));
                    }
                break;
        }
        int startBin = (int) ceil((windowParam->freq1 * fftLength/2) / fNyq);
        for (int n=0; n < nBins; n++){
            filter[(startBin+n)%fftLength] += window[n]; //TODO: combine
        }

    }

    std::size_t const nyquistBin = filter.size() / 2;
    std::vector<TraceValueType> lower_half(filter.begin(), filter.begin() + nyquistBin + 1);
    std::vector<TraceValueType> upper_half(lower_half);
    upper_half.pop_back();
    reverse(upper_half.begin(), upper_half.end());
    upper_half.pop_back();
    lower_half.insert(lower_half.end(), upper_half.begin(), upper_half.end());
    filter = lower_half;

}

