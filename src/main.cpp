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
#include "statisticaltest/base.hpp"
#include <iostream>
#include <fstream>
#include "includes.h"
#include "main.hpp"
#include <sys/time.h>
#define VERSION "1.0beta"
using namespace Eigen;
using namespace std;

void DPA::doRun() //BBQ-style. This method can be started multiple times in different threads. Elaborates BATCH_SIZE time samples, for each possible key value.
{
	unsigned long long myid;
	shared_ptr<TracesMatrix> traces;
	shared_ptr<StatisticIndexMatrix> sm;
	int num = input->read ( &myid, &traces );
	sm.reset ( new StatisticIndexMatrix ( num, KEYNUM ) );
	stat->generate ( sm, traces, num );
	outp->WriteBatch ( myid, sm );
}

void DPA::prefetch()
{
	while ( input->CurrentSample < input->SamplesPerTrace ) {
		input->populateQueue();
	}
}

int DPA::main ( int argc, char** argv )
{
	TCLAP::CmdLine cmd ( "DPA calc", ' ', VERSION );
	exec = shared_ptr<ExecMethod::base> ( new ExecMethod::EXECCLASS ( cmd ) );
	input = shared_ptr<SamplesInput::base> ( new SamplesInput::INPUTCLASS ( cmd ) );
    filter = shared_ptr<Filters::base> ( new Filters::FILTERCLASS ( cmd ) );
	keygen = shared_ptr<KeyGenerators::base> ( new KeyGenerators::KEYGENCLASS ( cmd ) );
	interm = shared_ptr<GenerateIntermediateValues::base> ( new GenerateIntermediateValues::GENINTERMCLASS ( cmd, keygen ) );
	genpm = shared_ptr<GeneratePowerModel::base> ( new GeneratePowerModel::GENPOWERMODELCLASS ( cmd ) );
	stat = shared_ptr<Statistic::base> ( new Statistic::STATISTICCLASS ( cmd ) );
	outp = shared_ptr<Output::base> ( new Output::OUTPUTCLASS ( cmd, keygen ) );
	this->ShowCompileTimeOptions();
	try {
		cmd.parse ( argc, argv );
	} catch ( TCLAP::ArgException& e ) {
		cerr << "Error " << e.error() << " in command line argument " << e.argId() << std::endl;
		return 1;
	}
	input->init();
	timeval start, end;
	gettimeofday ( &start, NULL );
	keygen->init();
	interm->init();
	genpm->init();
	outp->init();
	numbatches = ( input->SamplesPerTrace / BATCH_SIZE ) + ( ( ( input->SamplesPerTrace % BATCH_SIZE ) == 0 ) ? 0 : 1 );
	cout << "Reading known data..." << endl;
	data = input->readData();
	cout << "Done. Calculating intermediate values.....[single threaded]" << endl;
	intval.reset (  new IntermediateValueMatrix ( input->NumTraces, KEYNUM ) );
	interm->generate ( data, intval );
	data.reset(); // I don't need that data anymore.
	cout << "Done. Calculating power model.....[single threaded]" << endl;
	pm.reset ( new PowerModelMatrix ( input->NumTraces, KEYNUM ) );
	genpm->generate ( intval, pm );
	intval.reset(); // I don't need that data anymore, let's free some memory!
	cout << "Done. Initializing statistic test [single threaded]:" << endl;
	// StatisticIndexMatrix size should be a multiple of BATCH_SIZE
	unsigned long sz = input->SamplesPerTrace;
	if ( sz % BATCH_SIZE > 0 ) { sz += ( BATCH_SIZE - ( sz % BATCH_SIZE ) ) ; }
	stat->init ( pm );
	cout << "Done. Starting statistic test pass 1 [multithreaded]" << endl;
	exec->RunAndWait ( numbatches );
	gettimeofday ( &end, NULL );
	outp->end();
	cout << "Elaboration took " << timevaldiff ( &start, &end ) << " milliseconds." << endl;
	return 0;
}
void DPA::ShowCompileTimeOptions()
{
	cout << "DPAcalc was compiled with : " << endl;
	cout << "Batch size : " << BATCH_SIZE << endl;
	cout << "Number of bit of the key to guess : " << KEY_HYP_BIT << endl;
	cout << "Size of known data : " << DATA_SIZE_BIT << " bit " << endl;
	cout << "Size of key : " << KEY_SIZE_BIT << " bit " << endl;
	cout << endl;
	cout << "Name of the class that reads input file: " << INPUTCLASS_STR << endl;
	cout << "Name of the class that generates intermediate values: " << GENINTERMCLASS_STR << endl;
	cout << "Name of the class that generates power model: " << GENPOWERMODELCLASS_STR << endl;
	cout << "Name of the class that calculates statistic data: " << STATISTICCLASS_STR << endl;
	cout << "Name of the class that manages parallelization: " << EXECCLASS_STR << endl;
	cout << "Name of the class that writes output: " << OUTPUTCLASS_STR << endl;
	cout << endl;
}

int main ( int argc, char** argv )
{
	DPA* me = DPA::instance();
	return me->main ( argc, argv );
}
