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
#include "dpacalc.hpp"
#include <sys/time.h>
#define VERSION "1.0beta"
using namespace Eigen;
using namespace std;

int DPA::main ( int argc, char** argv )
{
    auto filtFunc = [&] {
        shared_ptr<TraceWithData> trace(new TraceWithData);
        traceMutex.lock();
        input->readTraceWithData(trace, curTrace);
        int localTrace = curTrace++;
        traceMutex.unlock();
        filter->applyFilter(trace);
        filter->writeFilteredTrace(trace, localTrace);
    };
    auto prefetchFunc = [&] {
        while ( input->CurrentSample < input->SamplesPerTrace ) {
            input->populateQueue();
        }
    };
    auto runFunc = [&] {
        unsigned long long myid;
        shared_ptr<TracesMatrix> traces;
        shared_ptr<StatisticIndexMatrix> sm;
        int num = input->read ( &myid, &traces );
        sm.reset ( new StatisticIndexMatrix ( num, KEYNUM ) );
        stat->generate ( sm, traces, num );
        outp->WriteBatch ( myid, sm );
    };
	TCLAP::CmdLine cmd ( "DPA calc", ' ', VERSION );
	exec = shared_ptr<ExecMethod::base> ( new ExecMethod::EXECCLASS ( cmd ) );
	input = shared_ptr<SamplesInput::base> ( new SamplesInput::INPUTCLASS ( cmd ) );
    filter = shared_ptr<Filters::base> ( new Filters::FILTERCLASS ( cmd ) );
	keygen = shared_ptr<KeyGenerators::base> ( new KeyGenerators::KEYGENCLASS ( cmd ) );
	interm = shared_ptr<GenerateIntermediateValues::base> ( new GenerateIntermediateValues::GENINTERMCLASS ( cmd, keygen ) );
	genpm = shared_ptr<GeneratePowerModel::base> ( new GeneratePowerModel::GENPOWERMODELCLASS ( cmd ) );
	stat = shared_ptr<Statistic::base> ( new Statistic::STATISTICCLASS ( cmd ) );
	outp = shared_ptr<Output::base> ( new Output::OUTPUTCLASS ( cmd, keygen ) );
    TCLAP::SwitchArg filterSwitch("i", "filter-input", "If set, the input is filtered. You must provide a configuration file with -c");
    cmd.add(filterSwitch);
	this->ShowCompileTimeOptions();
	try {
		cmd.parse ( argc, argv );
	} catch ( TCLAP::ArgException& e ) {
		cerr << "Error " << e.error() << " in command line argument " << e.argId() << std::endl;
		return 1;
	}
	input->init();
    timeval start, endfilter, end;
	gettimeofday ( &start, NULL );
    if (filterSwitch.isSet()){
        filter->init(input->SamplesPerTrace, input->NumTraces);
    } else {
        filter.reset();
    }
	keygen->init();
	interm->init();
	genpm->init();
	outp->init();
    if (filterSwitch.isSet()){
        cout << "Filtering..." << endl;
        filter->initFilterOutput();
        curTrace = 0;
        exec->RunAndWait(input->NumTraces, filtFunc, NULL);
        unsigned int newsize;
        void* newpointer = filter->getFilteredPointer(newsize);
        input->changeFileOffset(newpointer, newsize);
        gettimeofday ( &endfilter, NULL );
        cout << "Filtering took " << timevaldiff ( &start, &endfilter ) << " milliseconds." << endl;
        cout << "Done. ";
    }
	numbatches = ( input->SamplesPerTrace / BATCH_SIZE ) + ( ( ( input->SamplesPerTrace % BATCH_SIZE ) == 0 ) ? 0 : 1 );
	cout << "Reading known data..." << endl;
	data = input->readData();
	cout << "Done. Calculating intermediate values.....[single threaded]" << endl;
	intval.reset (  new IntermediateValueMatrix ( input->NumTraces, KEYNUM ) );
	interm->generate ( data, intval );
    data.reset(); // tracewd->I don't need that data anymore.
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
    exec->RunAndWait ( numbatches,  runFunc, prefetchFunc);
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
#if defined(CONFIG_FILTER_OUTPUT_DISK)
    cout << "Filter output on disk" << endl;
#elif defined(CONFIG_FILTER_OUTPUT_RAM)
    cout << "Filter output on RAM" << endl;
#endif
#if defined(CONFIG_FILTER_COMBINE_NOTHING)
    cout << "Filter combining: do nothing" << endl;
#elif defined(CONFIG_FILTER_COMBINE_NORMALIZE)
    cout << "Filter combining: normalize" << endl;
#elif defined(CONFIG_FILTER_COMBINE_CLAMP)
    cout<< "Filter combining: clamp" << endl;
#endif


	cout << endl;
	cout << "Name of the class that reads input file: " << INPUTCLASS_STR << endl;
    cout << "Name of the class that filters the data: " << FILTERCLASS_STR << endl;
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
