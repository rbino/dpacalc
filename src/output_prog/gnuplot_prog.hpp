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
namespace OutputProg
{
    class gnuplot_prog: public base
	{
		public:
            gnuplot_prog ( TCLAP::CmdLine& cmd, shared_ptr<KeyGenerators::base> _keygen ) : base ( cmd, _keygen ),
                dataNameArg ( "o", "output", "Gnuplot data (all keys with correlation coefficient) file name", true, "", "path" ),
                scriptNameArg ( "s", "script-output", "Gnuplot script output (all keys with correlation coefficient) file name", true, "", "path" ),
                confidenceDataNameArg ( "r", "confidence-output", "Gnuplot data (best and second best keys with confidence interval) file name", true, "", "path" ),
                confidenceScriptNameArg ( "g", "confidence-script-output", "Gnuplot script output (best and second best keys with confidence interval) file name", true, "", "path" ),
                alphaArg ( "a", "alpha", "The alpha to compute the (1 - alpha) confidence interval", true, 0.05, "0-1" ),
                bestPearson() {
                bestPearson = Trace::Zero(KEYNUM);
				cmd.add ( dataNameArg );
				cmd.add ( scriptNameArg );
                cmd.add ( confidenceDataNameArg );
                cmd.add ( confidenceScriptNameArg );
                cmd.add ( alphaArg );
                currentTraces = 0;
			};
			virtual void init();
			virtual void end();
			virtual void WriteBatch ( unsigned long long id, shared_ptr<StatisticIndexMatrix>& s );
			virtual void endTraceBlock();
		protected:
			TCLAP::ValueArg<std::string> dataNameArg;
			TCLAP::ValueArg<std::string> scriptNameArg;
            TCLAP::ValueArg<std::string> confidenceDataNameArg;
            TCLAP::ValueArg<std::string> confidenceScriptNameArg;
            TCLAP::ValueArg<StatisticValueType> alphaArg;
			std::ofstream dataoutp;
            std::ofstream confdataoutp;
            Eigen::Matrix <StatisticValueType, 1, KEYNUM > bestPearson;
            std::mutex checkMutex;
            ConfidencePair getConfidence(StatisticValueType r, unsigned long n, StatisticValueType alpha);

	};

}
