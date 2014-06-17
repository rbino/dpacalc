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
#include "output_fftfilterconf.hpp"
#include <mutex>
#include <bits/stl_map.h>
using namespace std;
void OutputFilter::output_fftfilterconf::init()
{
	dataoutp.open ( dataNameArg.getValue() );
	if ( !dataoutp.is_open() ) {
		cerr << "Please provide a valid data output filename" << endl;
		exit ( 1 );
	}
	ofstream scriptoutp ( scriptNameArg.getValue() );
	if ( !scriptoutp.is_open() ) {
		cerr << "Please provide a valid script output filename" << endl;
		exit ( 1 );
	}
    confdataoutp.open( confidenceDataNameArg.getValue()  );
    if ( !confdataoutp.is_open() ) {
        cerr << "Please provide a valid confidence data output filename" << endl;
        exit ( 1 );
    }
    ofstream confidencescriptoutp ( confidenceScriptNameArg.getValue() );
    if ( !confidencescriptoutp.is_open() ) {
        cerr << "Please provide a valid confidence script output filename" << endl;
        exit ( 1 );
    }
    if (alphaArg.getValue() <= 0 || alphaArg.getValue() >= 1){
        cerr << "Alpha must be 0 < alpha < 1" << endl;
        exit ( 1 );
    }
    scriptoutp << "set term png size 3000,1500 crop;" << endl;
    scriptoutp << "set output \"output_prog.png\";" << endl;
	scriptoutp << "set autoscale;" << endl;
	scriptoutp << "set xtic auto;" << endl;
	scriptoutp << "set ytic auto;" << endl;
	scriptoutp << "set key outside right;" << endl;
    scriptoutp << "set title \"dpacalc_prog graphical output\";" << endl;
    scriptoutp << "set xlabel \"Number of traces\";" << endl;
    scriptoutp << "set ylabel \"Max Pearson coefficient\";" << endl << endl << endl << endl << endl;
	scriptoutp << "plot ";
	for ( unsigned long long k = 0; k < KEYNUM; k++ ) {
		scriptoutp  << " \"" << dataNameArg.getValue() << "\" u 1:" << k + 2 << " t \"" << keygen->getKeyAsString ( k ) << "\" with lines";
		if ( k != KEYNUM - 1 ) {
			scriptoutp  << ",";
		}
	}
	scriptoutp.close();

    confidencescriptoutp << "set term pngcairo dashed size 3000,1500 crop;" << endl;
    confidencescriptoutp << "set output \"output_prog_conf.png\";" << endl;
    confidencescriptoutp << "set autoscale;" << endl;
    confidencescriptoutp << "set xtic auto;" << endl;
    confidencescriptoutp << "set ytic auto;" << endl;
    confidencescriptoutp << "set style line 1 lt 2 lw 2 pt 1 linecolor rgb \"green\";" << endl;
    confidencescriptoutp << "set style line 2 lt 1 lw 2 pt 1 linecolor rgb \"green\";" << endl;
    confidencescriptoutp << "set style line 3 lt 2 lw 2 pt 1 linecolor rgb \"red\";" << endl;
    confidencescriptoutp << "set style line 4 lt 1 lw 2 pt 1 linecolor rgb \"red\";" << endl;
    confidencescriptoutp << "set key outside right;" << endl;
    confidencescriptoutp << "set title \"dpacalc_prog graphical output with confidence interval\";" << endl;
    confidencescriptoutp << "set xlabel \"Number of traces\";" << endl;
    confidencescriptoutp << "set ylabel \"Max Pearson coefficient and confidence interval with alpha=" << alphaArg.getValue() << "\";" << endl << endl << endl << endl << endl;
    confidencescriptoutp << "plot ";
    for ( unsigned long long k = 0; k < 6; k++ ) {
        confidencescriptoutp  << " \"" << confidenceDataNameArg.getValue() << "\" u 1:" << k + 2 << " t \"" << keygen->getKeyAsString ( k ) << "\" with lines ls ";
        /* Ugly as hell, but it's to have the line styles right */
        switch (k) {
        case 0:
        case 2:
            confidencescriptoutp << "1,";
            break;
        case 1:
            confidencescriptoutp << "2,";
            break;
        case 3:
            confidencescriptoutp << "3,";
            break;
        case 4:
            confidencescriptoutp << "4,";
            break;
        case 5:
            confidencescriptoutp << "3";
            break;
        }
    }
    confidencescriptoutp.close();

}

void OutputFilter::output_fftfilterconf::WriteBatch ( unsigned long long id, shared_ptr< StatisticIndexMatrix >& s )
{
    Eigen::Matrix < StatisticValueType, 1, KEYNUM > maxPearson;
    ( *s ) = s->cwiseAbs();
    maxPearson = s->colwise().maxCoeff();
    checkMutex.lock();
       for(unsigned long long key = 0; key < KEYNUM; key++){
           //TODO: Maybe transform the two vectors in a 2xKEYNUM Matrix, compute the colwise max and assign it to bestPearson
           if (bestPearson (key) < maxPearson (key)){
               bestPearson (key) = maxPearson (key);
           }
       }
    checkMutex.unlock();
}

void OutputFilter::output_fftfilterconf::endTraceBlock()
{
    dataoutp << currentTraces;
    for (long long key = 0; key < bestPearson.cols(); key++){
        dataoutp << "\t" << bestPearson (key);
    }
    dataoutp << endl;
    // Little trick to find max and second max: we find max, than we set it to 0 (we don't need it anymore anyway)
    unsigned long long maxRow, maxCol;
    StatisticValueType best = bestPearson.maxCoeff(&maxRow, &maxCol);
    bestPearson (maxCol) = 0;
    ConfidencePair bestConf = getConfidence(best, currentTraces, alphaArg.getValue());

    StatisticValueType secondBest = bestPearson.maxCoeff();
    ConfidencePair secondBestConf = getConfidence(secondBest, currentTraces, alphaArg.getValue());

    confdataoutp << currentTraces << "\t" << bestConf.first << "\t" << best << "\t" << bestConf.second << "\t" << secondBestConf.first << "\t" << secondBest << "\t" << secondBestConf.second << endl;
    bestPearson = Trace::Zero(KEYNUM);
}

void OutputFilter::output_fftfilterconf::end()
{
    OutputFilter::base::end();
    confdataoutp.close();
	dataoutp.close();
}

ConfidencePair OutputFilter::output_fftfilterconf::getConfidence(StatisticValueType r, unsigned long n, StatisticValueType alpha)
{
    //TODO: if n < 4 I think there will be problems

    // Confidence level
    StatisticValueType gamma = 1-alpha;
    // Fisher transform r to map correctly the domains
    StatisticValueType Fr = atanh(r);
    StatisticValueType z_score = sqrt(n-3)*Fr;
    // obtain Gaussian quantile of transformed z_score
    StatisticValueType half_width=sqrt(2)*boost::math::erf_inv(2*((1+gamma)/2)-1);
    // compute confidence interval in Fisher domain
    StatisticValueType inf=z_score-half_width;
    StatisticValueType sup=z_score+half_width;
    // invert Fisher map to get the domain straight again
    inf=tanh(inf/sqrt(n-3));
    sup=tanh(sup/sqrt(n-3));
    return ConfidencePair(inf,sup);
}

