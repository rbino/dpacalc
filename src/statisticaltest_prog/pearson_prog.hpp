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
#include <map>
using namespace Eigen;
using namespace std;
namespace StatisticProg
{
	class pearson_prog: public base
	{
		public:
			virtual void generate ( shared_ptr<StatisticIndexMatrix>& stat, shared_ptr<TracesMatrix>& traces, unsigned long numvalid );
			pearson_prog ( TCLAP::CmdLine& cmd ) : base ( cmd ) {};
			virtual void init ( shared_ptr<PowerModelMatrix>& _pm );
		protected:
			Eigen::Matrix<TraceValueType, Dynamic, Dynamic> pmexpect;
			Eigen::Matrix<TraceValueType, 1, Dynamic> pmexpect_bykey;
	};
}

