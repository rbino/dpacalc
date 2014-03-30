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
#include "base.hpp"
#include <fstream>
#include <map>
#include <mutex>
using namespace std;
namespace Output
{
    class gnuplot_prog: public base
	{
		public:
            gnuplot_prog ( TCLAP::CmdLine& cmd, shared_ptr<KeyGenerators::base> _keygen ) : base ( cmd, _keygen ),
				dataNameArg ( "o", "output", "Gnuplot data file name", true, "", "path" ),
				scriptNameArg ( "s", "script-output", "Gnuplot script output file name", true, "", "path" ), wqueue(), queueMutex() {
				doneId = -1;
				cmd.add ( dataNameArg );
				cmd.add ( scriptNameArg );
			};
			virtual void init();
			virtual void end();
			virtual void WriteBatch ( unsigned long long id, shared_ptr<StatisticIndexMatrix>& s );
		protected:
			TCLAP::ValueArg<std::string> dataNameArg;
			TCLAP::ValueArg<std::string> scriptNameArg;
			std::ofstream dataoutp;
			unsigned long long doneId;
			std::map<unsigned long long, shared_ptr<StatisticIndexMatrix>> wqueue;
			std::mutex queueMutex;
			virtual void RealWriteBatch ( unsigned long long id, shared_ptr<StatisticIndexMatrix>& s );
			virtual void emptyQueue();
	};

}
