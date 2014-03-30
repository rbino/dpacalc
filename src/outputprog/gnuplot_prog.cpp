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
#include "gnuplot_prog.hpp"
#include <mutex>
#include <bits/stl_map.h>
using namespace std;
void OutputProg::gnuplot_prog::init()
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
}

void OutputProg::gnuplot_prog::WriteBatch ( unsigned long long id, shared_ptr< StatisticIndexMatrix >& s )
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

void OutputProg::gnuplot_prog::endTraceBlock()
{
    dataoutp << currentTraces;
    for (long long key = 0; key < bestPearson.cols(); key++){
        dataoutp << "\t" << bestPearson (key);
    }
    dataoutp << endl;
    bestPearson = Trace::Zero(KEYNUM);
}

void OutputProg::gnuplot_prog::end()
{
	OutputProg::base::end();
	dataoutp.close();
}

