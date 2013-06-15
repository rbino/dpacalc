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

using namespace Eigen;
namespace Filters
{
	class base
	{
		public:
			unsigned long long SamplesPerTrace; //from metadata, dimension N of matrix T
			unsigned long long NumTraces; //dimension N of matrix T
            unsigned long long SamplingFrequency;
			base ( TCLAP::CmdLine& cmd ) {};
          //  virtual void applyFilter() = 0;
	};
}
