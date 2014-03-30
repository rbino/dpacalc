/*
Copyright (C) 2012	Massimo Maggi

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
void Output::gnuplot_prog::init()
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
	scriptoutp << "set output \"output.png\";" << endl;
	scriptoutp << "set autoscale;" << endl;
	scriptoutp << "set xtic auto;" << endl;
	scriptoutp << "set ytic auto;" << endl;
	scriptoutp << "set key outside right;" << endl;
	scriptoutp << "set title \"DPAcalc graphical output\";" << endl;
	scriptoutp << "set xlabel \"Time\";" << endl;
	scriptoutp << "set ylabel \"Pearson coefficient\";" << endl << endl << endl << endl << endl;
	scriptoutp << "plot ";
	for ( unsigned long long k = 0; k < KEYNUM; k++ ) {
		scriptoutp  << " \"" << dataNameArg.getValue() << "\" u 1:" << k + 2 << " t \"" << keygen->getKeyAsString ( k ) << "\" with lines";
		if ( k != KEYNUM - 1 ) {
			scriptoutp  << ",";
		}
	}
	scriptoutp.close();
}
//this queue system is not designed to enhance output performance, but to ensure output correctness.
// (batches with higher ID cannot write to the output before previous batches)
void Output::gnuplot_prog::WriteBatch ( unsigned long long id, shared_ptr< StatisticIndexMatrix >& s )
{
	queueMutex.lock();
	//if(id==2)id=3;else if(id==3)id=2; // Ugly, destroys final output, but clearly demonstrates the problem.
	if ( id == doneId + 1 ) {
		this->RealWriteBatch ( id, s );
		doneId = id;
		this->emptyQueue();
	} else {
		wqueue.insert ( pair<unsigned long long, shared_ptr< StatisticIndexMatrix > > ( id, shared_ptr<StatisticIndexMatrix> ( s ) ) );
	}
	queueMutex.unlock();
}
void Output::gnuplot_prog::emptyQueue()
{
	long id = doneId;
	while ( wqueue.count ( ++id ) ) {
		this->RealWriteBatch ( id, wqueue[id] );
		doneId = id;
		wqueue.erase ( id );
	}
}

void Output::gnuplot_prog::RealWriteBatch ( unsigned long long id, shared_ptr< StatisticIndexMatrix >& s )
{
	for ( long long r = 0; r < s->rows(); r++ ) {
		dataoutp << ( id * BATCH_SIZE ) + r ;
		for ( long long c = 0; c < s->cols(); c++ ) {
			dataoutp << "\t" << ( *s ) ( r, c );
		}
		dataoutp << endl;
	}
}
void Output::gnuplot_prog::end()
{
	Output::base::end();
	dataoutp.close();
}

