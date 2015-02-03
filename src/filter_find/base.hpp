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
#pragma once
#include "dpacalc.h"
#include "input_find/base.hpp"
#include <utility>

using namespace std;
using namespace Eigen;
namespace FilterFind
{
	class base
	{
		public:
            base ( TCLAP::CmdLine& cmd ) {};
            virtual void init(unsigned long long _samplespertrace, unsigned long long _numtraces) {};
						virtual void applyFilter(shared_ptr<TraceWithData>& tracewd) = 0;
						virtual void* getFilteredPointer(unsigned int& newsize) = 0;
						virtual void initFilterOutput() = 0;
            virtual void writeFilteredTrace(shared_ptr<TraceWithData> tracewd, unsigned int id) = 0;
            virtual bool hasFinished() = 0;
            virtual void setBaseline(unsigned long steps) = 0;
            virtual FilterBand beginStep() = 0;
            virtual void endStep(unsigned int successTraces) = 0;
            virtual bool isLastStep() = 0;
            virtual void end() = 0;

		protected:
			unsigned long long SamplesPerTrace; //from metadata, dimension N of matrix T
			unsigned long long NumTraces; //dimension N of matrix T
            double SamplingFrequency;
	};
}
