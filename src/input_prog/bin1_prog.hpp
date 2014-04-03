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
#pragma once
#include "dpacalc.h"
#include <iostream>
#include <queue>
#include <utility>
#include "base.hpp"
using namespace Eigen;
using namespace std;
namespace SamplesInputProg
{
/*
#pragma pack(push)

#pragma pack(1)
	struct fileheaders {
		uint32_t numtraces;
		uint32_t numsamples_per_trace;
		char datatype;
		uint8_t knowndatalength;
	};

#pragma pack(pop)
*/
	struct queueelement {
		unsigned long long id;
		long long size;
		shared_ptr<TracesMatrix> traces;
		queueelement() {
			id = 0;
			size = -1;
		}
	};
	class bin1_prog : public base
	{
		public:
			virtual unsigned long long read ( unsigned long long* id, shared_ptr<TracesMatrix>* traces );
			bin1_prog ( TCLAP::CmdLine& cmd ) :
				base ( cmd ),
				nameArg ( "f", "filename", "Input file name", true, "", "path" ),
				mlockArg ( "m", "mlock", "mlock entire input file in ram? Use only if you know what you're doing.", false ),
				queuemutex(), readytraces() {
				cmd.add ( nameArg );
				cmd.add ( mlockArg );
			};
			virtual void init();
			~bin1_prog();
            shared_ptr<DataMatrix> readProgressiveData(unsigned int step);
            /**
             * @brief readTraceWithData Reads a trace with the associated data
             * @param tracewd shared_ptr the the trace which is filled by the function
             * @param id The progressive id of the trace
             */
            void readTraceWithData(shared_ptr<TraceWithData>& tracewd, unsigned long id);
            /**
             * @brief changeFileOffset Changes the fileoffset member variable (a pointer that points to the input traces)
             * @param newOffset The new pointer, which will point to the filtered traces
             * @param newSize The new size of the buffer (or mmaped file)
             */
            void changeFileOffset(void* newOffset, long long newSize);
            /**
             * @brief changeNumTraces Changes the number of traces to be processed
             * @param newOffset The new number of traces
             */
            void changeNumTraces(unsigned long long newNum);
            void reinit();
		protected:
			TCLAP::ValueArg<std::string> nameArg;
			TCLAP::SwitchArg mlockArg;
			std::mutex input_mutex;
			int inputfd;
			template <class T> void readSamples ( shared_ptr<TracesMatrix>& traces, unsigned long curtrace, unsigned long startingsample, unsigned long numsamples );
			char sampletype;
			int samplesize;
			void* fileoffset;
            bool offsetUnmap;
			shared_ptr<DataMatrix> data;
			void populateQueue();
			mutex queuemutex;
			long long RealFileSize;
			queue<queueelement> readytraces;
            unsigned long long getSampleOffset ( unsigned long long trace, unsigned long long samplenum ) {
                //trace and samplenum are zero-based
                return sizeof ( struct fileheaders ) + trace * ( samplesize * SamplesPerTrace + DATA_SIZE_BYTE ) + samplesize * samplenum;
            }
            unsigned long long getDataOffset ( unsigned long long trace ) {
                //trace and samplenum are zero-based
                return sizeof ( struct fileheaders ) + trace * ( samplesize * SamplesPerTrace + DATA_SIZE_BYTE ) + samplesize * SamplesPerTrace;
            }
            template <class T> void readTraceWithDataImplem(shared_ptr<TraceWithData>& tracewd, unsigned long id);
	};


}
