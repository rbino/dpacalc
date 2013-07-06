/*
Copyright (C) 2013	Riccardo Binetti

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
#include "input/base.hpp"

using namespace std;
using namespace Eigen;
namespace Filters
{
	class base
	{
		public:
			base ( TCLAP::CmdLine& cmd, shared_ptr<SamplesInput::base> _input ) : input ( _input ) {};
            virtual void init() {};
						virtual void applyFilter(shared_ptr<TraceWithData>& tracewd) = 0;
						virtual void* getFilteredPointer(unsigned int& newsize) = 0;
						virtual void initFilterOutput() = 0;
            virtual void writeFilteredTrace(shared_ptr<TraceWithData> tracewd, unsigned int id) = 0;

          //  virtual void applyFilter() = 0;
		protected:
			shared_ptr<SamplesInput::base> input;
			unsigned long long SamplesPerTrace; //from metadata, dimension N of matrix T
			unsigned long long NumTraces; //dimension N of matrix T
            double SamplingFrequency;
	};
}
