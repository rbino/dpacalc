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
#include <mutex>

using namespace Eigen;
using namespace std;
namespace SamplesInput
{
	class base
	{
		public:
			unsigned long long SamplesPerTrace; //from metadata, dimension N of matrix T
			unsigned long long NumTraces; //dimension N of matrix T
			unsigned long long CurrentSample;  //dimension N of matrix T
			unsigned long long CurrentTrace;
			base ( TCLAP::CmdLine& cmd ) {};
			virtual void init() {
				CurrentSample = 0;
				CurrentTrace = 0;
				CurrentId = -1;
				NumTraces = 0;
				SamplesPerTrace = 0;
			};
			virtual unsigned long long read ( unsigned long long* id, shared_ptr<TracesMatrix>* traces ) = 0;
			virtual shared_ptr<DataMatrix> readData() = 0;
			virtual void populateQueue() {};
			virtual void readTraceWithData(shared_ptr<TraceWithData>& tracewd, unsigned long id) = 0;
			virtual void changeFileOffset(void* newOffset, long long newSize) = 0;

		protected:
			shared_ptr<DataMatrix> data;
			unsigned long long CurrentId;
	};

}
