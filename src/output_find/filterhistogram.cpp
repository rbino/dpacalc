/*
Copyright (C) 2014	Riccardo Binetti

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
#include "dpacalc.h"
#include "filterhistogram.hpp"
#include <mutex>
#include <bits/stl_map.h>
using namespace std;
void OutputFind::filterhistogram::init()
{
	dataoutp.open ( dataNameArg.getValue() );
	if ( !dataoutp.is_open() ) {
		cerr << "Please provide a valid data output filename" << endl;
		exit ( 1 );
    }
}

bool OutputFind::bandWithTracesSort(BandWithTraces band1, BandWithTraces band2){
    return band1.band.first < band2.band.first;
}

void OutputFind::filterhistogram::writeBand ( FilterBand band, unsigned int ntraces )
{
    BandWithTraces res;
    res.band = band;
    res.ntraces = ntraces;
    results.push_back(res);
}

void OutputFind::filterhistogram::end()
{
    sort(results.begin(), results.end(), bandWithTracesSort);
    unsigned int prevlo = results.front().band.first + 10;
    float prevhi = results.front().band.second;
    for (vector<BandWithTraces>::iterator it = results.begin(); it != results.end(); ++it){
        if (prevhi < it->band.first) {
           dataoutp << std::fixed << prevhi << "\t" << "0" << endl;
        }
        if (floor(it->band.first != prevlo)){
           dataoutp << std:: fixed << it->band.first << "\t" << it->ntraces << endl;
        }
        prevhi = it->band.second;
        prevlo = floor(it->band.first);
    }
    dataoutp << results.back().band.second << "\t" << results.back().ntraces << endl;
    OutputFind::base::end();
    dataoutp.close();
}

