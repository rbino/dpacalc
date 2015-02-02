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
#include "dpacalc.h"
#include "statisticaltest/base.hpp"
#include <iostream>
#include <fstream>
#include "includes.h"
#include "dpacalc_findfilter.hpp"
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
        stat->progressiveGenerate ( sm, traces, num, myid );
        verif->batchBest(myid, sm);
    };
	TCLAP::CmdLine cmd ( "DPA calc", ' ', VERSION );
	exec = shared_ptr<ExecMethod::base> ( new ExecMethod::EXECCLASS ( cmd ) );
    input = shared_ptr<SamplesInputFind::base> ( new SamplesInputFind::INPUTFINDCLASS ( cmd ) );
    filter = shared_ptr<FilterFind::base> ( new FilterFind::FILTERFINDCLASS ( cmd, input ) );
	keygen = shared_ptr<KeyGenerators::base> ( new KeyGenerators::KEYGENCLASS ( cmd ) );
	interm = shared_ptr<GenerateIntermediateValuesProg::base> ( new GenerateIntermediateValuesProg::GENINTERMPROGCLASS ( cmd, keygen ) );
	genpm = shared_ptr<GeneratePowerModelProg::base> ( new GeneratePowerModelProg::GENPOWERMODELPROGCLASS ( cmd ) );
    stat = shared_ptr<StatisticProg::base> ( new StatisticProg::STATISTICPROGCLASS ( cmd ) );
    verif = shared_ptr<VerifyAttack::base> (new VerifyAttack::VERIFICATIONCLASS ( cmd ) );
    outp = shared_ptr<OutputFind::base> ( new OutputFind::OUTPUTFINDCLASS ( cmd ) );
    TCLAP::ValueArg<unsigned int> traceJump("t", "trace-step", "How many traces are added at every progressive round", true, 0, "int");
    cmd.add(traceJump);
	this->ShowCompileTimeOptions();
	try {
		cmd.parse ( argc, argv );
	} catch ( TCLAP::ArgException& e ) {
		cerr << "Error " << e.error() << " in command line argument " << e.argId() << std::endl;
		return 1;
    }
    input->init();
    timeval start, startbatch, end, endbatch, endfilter;
	gettimeofday ( &start, NULL );
    filter->init();
    keygen->init();
    interm->init();
    genpm->init();
	outp->init();

    input->NumTraces = 0;
    unsigned int step = traceJump.getValue();
    bool success = false;
    cout << "dpacalc_prog: calculating baseline" << endl;
    while(input->NumTraces < input->RealNumTraces && !success){
        gettimeofday ( &startbatch, NULL );
        input->reinit();
        // changeNumTraces clamps the number of traces to RealNumTraces
        input->increaseNumTraces(step);
        verif->currentTraces = input->NumTraces;
        numbatches = ( input->SamplesPerTrace / BATCH_SIZE ) + ( ( ( input->SamplesPerTrace % BATCH_SIZE ) == 0 ) ? 0 : 1 );
        cout << "dpacalc_prog: now processing " << input->NumTraces << "/" << input->RealNumTraces << " unfiltered traces..." << endl;
        //cout << "Reading known data..." << endl;
        data = input->readProgressiveData(step);
        //cout << "Done. Calculating intermediate values.....[single threaded]" << endl;
        interm->progressiveGenerate ( data, intval, step );
        //cout << "Done. Calculating power model.....[single threaded]" << endl;
        genpm->progressiveGenerate ( intval, pm, step );
        //cout << "Done. Initializing statistic test [single threaded]:" << endl;
        // StatisticIndexMatrix size should be a multiple of BATCH_SIZE
        unsigned long sz = input->SamplesPerTrace;
        if ( sz % BATCH_SIZE > 0 ) { sz += ( BATCH_SIZE - ( sz % BATCH_SIZE ) ) ; }
        stat->init ( pm, step, numbatches );
        //cout << "Done. Starting statistic test pass 1 [multithreaded]" << endl;
        exec->RunAndWait ( numbatches,  runFunc, prefetchFunc);
        //cout << " Done!" << endl;
        success = verif->verify(stat);
        gettimeofday ( &endbatch, NULL );
        if (success){
            successfulTraces = input->NumTraces;
            cout << "Attack successful!" << endl;
        } else {
            cout << "Attack failed" << endl;
        }
        cout << "Batch elaboration of " << input->NumTraces << " traces took " << timevaldiff ( &startbatch, &endbatch ) << " milliseconds." << endl;
    }
    cout << "Baseline found: " << successfulTraces << " traces" << endl << endl;
    cout << "dpacalc_prog: beginning to search for the best filter" << endl;
    filter->setBaseline(successfulTraces);
    while (!filter->hasFinished()){
        success = false;
        input->NumTraces = input->RealNumTraces;
        pm.reset();
        intval.reset();
        stat->reset();
        FilterBand curBand;
        curTrace = 0;
        curBand = filter->beginStep();
        cout << "Trying with band " << std::fixed << curBand.first << "-" << curBand.second << " Hz" << endl;
        filter->initFilterOutput();
        unsigned int newsize;
        void* newpointer = filter->getFilteredPointer(newsize);
        input->changeFileOffset(newpointer, newsize);
        //cout << "Done. ";

        input->NumTraces = 0;

        while(input->NumTraces < input->RealNumTraces && !success){
            gettimeofday ( &startbatch, NULL );
            input->reinit();
            // changeNumTraces clamps the number of traces to RealNumTraces
            input->increaseNumTraces(step);
            verif->currentTraces = input->NumTraces;
            numbatches = ( input->SamplesPerTrace / BATCH_SIZE ) + ( ( ( input->SamplesPerTrace % BATCH_SIZE ) == 0 ) ? 0 : 1 );
            cout << "dpacalc_prog: now processing " << input->NumTraces << "/" << input->RealNumTraces << " traces..." << endl;
            //cout << "Reading known data..." << endl;
            //cout << "Filtering..." << endl;
            exec->RunAndWait(step, filtFunc, NULL);
            gettimeofday ( &endfilter, NULL );
            cout << "Filtering took " << timevaldiff ( &startbatch, &endfilter ) << " milliseconds." << endl;
            data = input->readProgressiveData(step);
            //cout << "Done. Calculating intermediate values.....[single threaded]" << endl;
            interm->progressiveGenerate ( data, intval, step );
            //cout << "Done. Calculating power model.....[single threaded]" << endl;
            genpm->progressiveGenerate ( intval, pm, step );
            //cout << "Done. Initializing statistic test [single threaded]:" << endl;
            // StatisticIndexMatrix size should be a multiple of BATCH_SIZE
            unsigned long sz = input->SamplesPerTrace;
            if ( sz % BATCH_SIZE > 0 ) { sz += ( BATCH_SIZE - ( sz % BATCH_SIZE ) ) ; }
            stat->init ( pm, step, numbatches );
            //cout << "Done. Starting statistic test pass 1 [multithreaded]" << endl;
            exec->RunAndWait ( numbatches,  runFunc, prefetchFunc);
            //cout << " Done!" << endl;
            success = verif->verify(stat);
            gettimeofday ( &endbatch, NULL );
            if (success){
                cout << "Attack successful!" << endl;
            } else {
                cout << "Attack failed" << endl;
            }
            cout << "Batch elaboration of " << input->NumTraces << " traces took " << timevaldiff ( &startbatch, &endbatch ) << " milliseconds." << endl;
        }

        if (filter->isLastStep()){
            outp->writeBand(curBand, input->NumTraces % input->RealNumTraces );
        }
        filter->endStep(input->NumTraces);
        input->resetFileOffset();
    }
    filter->end();
    gettimeofday ( &end, NULL );
	outp->end();
    cout << "Total elaboration took " << timevaldiff ( &start, &end ) << " milliseconds." << endl;
	return 0;
}
void DPA::ShowCompileTimeOptions()
{
    cout << "DPAcalc_prog was compiled with : " << endl;
	cout << "Batch size : " << BATCH_SIZE << endl;
	cout << "Number of bit of the key to guess : " << KEY_HYP_BIT << endl;
	cout << "Size of known data : " << DATA_SIZE_BIT << " bit " << endl;
	cout << "Size of key : " << KEY_SIZE_BIT << " bit " << endl;
    cout << "Filter output on RAM" << endl;

	cout << endl;
	cout << "Name of the class that reads input file: " << INPUTCLASS_STR << endl;
    cout << "Name of the class that filters the data: " << FILTERCLASS_STR << endl;
	cout << "Name of the class that generates intermediate values: " << GENINTERMCLASS_STR << endl;
	cout << "Name of the class that generates power model: " << GENPOWERMODELCLASS_STR << endl;
    cout << "Name of the class that calculates statistic data: " << STATISTICPROGCLASS_STR << endl;
	cout << "Name of the class that manages parallelization: " << EXECCLASS_STR << endl;
    cout << "Name of the class that writes progressive output: " << OUTPUTPROGCLASS_STR << endl;
	cout << endl;
}

int main ( int argc, char** argv )
{
	DPA* me = DPA::instance();
	return me->main ( argc, argv );
}
