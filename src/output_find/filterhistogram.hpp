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
#include "base.hpp"
#include "math.h"
#include <boost/math/special_functions/erf.hpp>
#include <fstream>
#include <map>
#include <mutex>
using namespace std;



namespace OutputFind
{

    class filterhistogram: public base
	{
		public:
            filterhistogram ( TCLAP::CmdLine& cmd) : base ( cmd),
                dataNameArg ( "o", "output", "Gnuplot data (filter bands and number of traces for successful attack) file name", true, "", "path" )
            {
                cmd.add ( dataNameArg );
                prevHi = 0;
			};
			virtual void init();
			virtual void end();
            virtual void writeBand ( FilterBand band, unsigned int ntraces );
		protected:
			TCLAP::ValueArg<std::string> dataNameArg;
			std::ofstream dataoutp;
            std::vector<BandWithTraces> results;
            float prevHi;

	};

}
