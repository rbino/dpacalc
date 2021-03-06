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
#include "bin1_prog.hpp"
#include <algorithm>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

long long unsigned int SamplesInputProg::bin1_prog::read ( long long unsigned int* id, shared_ptr< TracesMatrix >* traces )
{
	queueelement qe;
	while ( qe.size == -1 ) {
		queuemutex.lock();
		if ( !readytraces.empty() ) {
			qe = readytraces.front();
			readytraces.pop();
			queuemutex.unlock();
		} else {
			queuemutex.unlock();
			populateQueue();
		}
	}
	*id = qe.id;
	( *traces ) = qe.traces;
	return qe.size;
}
void SamplesInputProg::bin1_prog::init()
{
	struct fileheaders header;
	long long TotalFileSize;
    offsetUnmap = true;
	SamplesInputProg::base::init();
	inputfd = open ( nameArg.getValue().c_str(), O_RDONLY );
	if ( inputfd == -1 ) {
		cerr << "Cannot open " << nameArg.getValue() << endl;
		exit ( 1 );
	}
	if ( ::read ( inputfd, ( void* ) &header, sizeof ( header ) ) != sizeof ( header ) ) {
		cerr << "Invalid header size" << endl;
		exit ( 1 );
	}
	SamplesPerTrace = header.numsamples_per_trace;
	NumTraces = header.numtraces;
    RealNumTraces = header.numtraces;
	switch ( header.datatype ) {
		case 'b':
			samplesize = 1;
			sampletype = 'b';
			break;
		case 'c':
			samplesize = 2;
			sampletype = 'c';
			break;
		case 'f':
			samplesize = 4;
			sampletype = 'f';
			break;
		case 'd':
			samplesize = 8;
			sampletype = 'd';
			break;
		default:
			cerr << "Sample data type should be b,c,f, or d. " << endl;
			exit ( 1 );
	}
	TotalFileSize = ( long long ) sizeof ( header ) + ( long long ) header.numtraces * ( ( long long ) header.knowndatalength + samplesize * ( long long ) header.numsamples_per_trace );
	RealFileSize = lseek ( inputfd, 0, SEEK_END );
	cout << "File size is " << RealFileSize << endl << " Header size: " << sizeof ( header ) << endl << header.numtraces << " traces, single trace size: " << ( int ) header.knowndatalength << " (known data) + " << header.numsamples_per_trace << "*" << samplesize << " bytes." << endl;
	if ( RealFileSize != TotalFileSize ) {
		cerr << "File size should be " << TotalFileSize << " but it is " << RealFileSize << endl << " Header size: " << sizeof ( header ) << endl << header.numtraces << " traces, single trace size: " << ( int ) header.knowndatalength << " (known data) + " << header.numsamples_per_trace << "*" << samplesize << " bytes." << endl;
		exit ( 2 );
	}
	if ( header.knowndatalength != DATA_SIZE_BYTE ) {
		cerr << "The length of each known data should be " << DATA_SIZE_BYTE << " but it is " << ( int ) header.knowndatalength << ". Change data file or recompile this program with appropriate options." << endl;
		exit ( 2 );
	}
	fileoffset = mmap ( NULL, RealFileSize, PROT_READ, MAP_SHARED, inputfd, 0 );
	if ( fileoffset == MAP_FAILED ) {
		cerr << "Cannot memory map input file. Cannot continue" << endl;
		exit ( 3 );
	}
    originalfileoffset = fileoffset;
    originalsamplesize = samplesize;
    OriginalRealFileSize = RealFileSize;
    originalsampletype = sampletype;
	if ( mlockArg.getValue() ) {
		cout << "mlock-ing" << endl;
		mlock ( fileoffset, RealFileSize );
		cout << "mlock-ed" << endl;
	}
}

void SamplesInputProg::bin1_prog::populateQueue()
{
	unsigned long long cur_trace;
	unsigned long long mysample;
	queueelement qe = queueelement();
	if ( readytraces.size() >= INPUT_QUEUE_SIZE ) {
		return;
	}
	input_mutex.lock();
	if ( CurrentSample >= SamplesPerTrace ) {
		input_mutex.unlock();
		return ;
	}
	unsigned long long num = min<unsigned long long> ( BATCH_SIZE, SamplesPerTrace - CurrentSample );
	mysample = CurrentSample;
	CurrentSample += num;
	++CurrentId;
	qe.size = num;
	qe.id = CurrentId;
    qe.traces = shared_ptr<TracesMatrix> ( new TracesMatrix ( CurrentStep, BATCH_SIZE ) );
	//  cout << "I'm going to allocate a " << NumTraces << " * " << BATCH_SIZE << " * " << sizeof(TraceValueType) << " = " << (NumTraces*BATCH_SIZE*sizeof(TraceValueType)/1024) << " kb matrix"<<endl;
    for ( cur_trace = NumTraces - CurrentStep; cur_trace < NumTraces; cur_trace++ ) {
		switch ( sampletype ) {
			case 'b':
				readSamples<uint8_t> ( qe.traces, cur_trace, mysample, num );
				break;
			case 'c':
				readSamples<uint16_t> ( qe.traces, cur_trace, mysample, num );
				break;
			case 'f':
				readSamples<float> ( qe.traces, cur_trace, mysample, num );
				break;
			case 'd':
				readSamples<double> ( qe.traces, cur_trace, mysample, num );
				break;
		}
	}
	queuemutex.lock();
	readytraces.push ( qe );
	queuemutex.unlock();
	input_mutex.unlock();
	return ;
}

void SamplesInputProg::bin1_prog::changeFileOffset(void *newOffset, long long newSize){
    fileoffset = newOffset;
    RealFileSize = newSize;
    samplesize = sizeof(TraceValueType);
    if (samplesize == sizeof(float)){
        sampletype = 'f';
    }
    else if (samplesize == sizeof(double)){
        sampletype = 'd';
    }
    offsetUnmap = false;
}

template <class T>void SamplesInputProg::bin1_prog::readSamples ( shared_ptr<TracesMatrix>& traces, unsigned long curtrace, unsigned long startingsample, unsigned long numsamples )
{
	T* buffer;
	//File is big enough, checked right after open.
	buffer = ( T* ) ( ( char* ) fileoffset + getSampleOffset ( curtrace, startingsample ) );
	for ( unsigned long i = 0; i < numsamples; i++ ) {
        ( *traces ) ( curtrace - NumTraces + CurrentStep, i ) = buffer[i];
	}
}

void SamplesInputProg::bin1_prog::readTraceWithData(shared_ptr<TraceWithData>& tracewd, unsigned long id){
    switch ( originalsampletype ) {
        case 'b':
            readTraceWithDataImplem<uint8_t> (tracewd, id );
            break;
        case 'c':
            readTraceWithDataImplem<uint16_t> (tracewd, id );
            break;
        case 'f':
            readTraceWithDataImplem<float> (tracewd, id );
            break;
        case 'd':
            readTraceWithDataImplem<double> (tracewd, id );
            break;
    }
}

template <class T> void SamplesInputProg::bin1_prog::readTraceWithDataImplem(shared_ptr<TraceWithData>& tracewd, unsigned long id){
    char* traceData;
    T* buffer;
    (tracewd -> trace).reset(new Trace(SamplesPerTrace));
    (tracewd -> data).reset(new DataValueType());
    //File is big enough, checked right after open.
    traceData = ( char* ) originalfileoffset +  getOriginalDataOffset ( id );
    BufferToBitset<DATA_SIZE_BYTE> (traceData, *(tracewd->data));
    buffer = ( T* ) ( ( char* ) originalfileoffset + getOriginalSampleOffset ( id, 0 ) );
    for ( unsigned long i = 0; i < SamplesPerTrace; i++ ) {
        (*tracewd -> trace) (i) = (TraceValueType) buffer[i];
    }
}

std::shared_ptr< DataMatrix > SamplesInputProg::bin1_prog::readProgressiveData(unsigned int step)
{
	char* buffer;
	if ( data.get() != NULL ) {
        input_mutex.lock();
        for ( unsigned long cur_trace = NumTraces - step; cur_trace < NumTraces; cur_trace++ ) {
            buffer = ( char* ) fileoffset +  getDataOffset ( cur_trace );
            BufferToBitset<DATA_SIZE_BYTE> ( buffer, ( *data ) [cur_trace - NumTraces + step] );
        }
        input_mutex.unlock();
		return shared_ptr<DataMatrix> ( data );
	}
	input_mutex.lock();
    data.reset ( new DataMatrix ( step ) );
    for ( unsigned long cur_trace = 0; cur_trace < step; cur_trace++ ) {
		buffer = ( char* ) fileoffset +  getDataOffset ( cur_trace );
		BufferToBitset<DATA_SIZE_BYTE> ( buffer, ( *data ) [cur_trace] );
	}
	input_mutex.unlock();
	return shared_ptr<DataMatrix> ( data );
}
SamplesInputProg::bin1_prog::~bin1_prog()
{
	if ( mlockArg.getValue() ) {
		munlock ( fileoffset, RealFileSize );
	}
    if (offsetUnmap){
        munmap ( fileoffset, RealFileSize );
    } else {
        free(fileoffset);
    }
	close ( inputfd );
}

void SamplesInputProg::bin1_prog::increaseNumTraces(unsigned int& step){
    if (NumTraces + step >= RealNumTraces){
        step = RealNumTraces - NumTraces;
    }
    CurrentStep = step;
    NumTraces += step;
}

void SamplesInputProg::bin1_prog::reinit(){
    CurrentSample = 0;
    CurrentId = -1;
}
