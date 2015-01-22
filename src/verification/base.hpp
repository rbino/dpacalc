/*
Copyright (C) 2014 Riccardo Binetti

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
#include "statisticaltest_prog/base.hpp"

namespace VerifyAttack
{
	class base
	{
		public:
			base ( TCLAP::CmdLine& cmd ) {};
			virtual void init() {};
      virtual void batchBest(unsigned long long id, shared_ptr<StatisticIndexMatrix>& s ) {};
			virtual bool verify(shared_ptr<StatisticProg::base> stat) = 0;
			unsigned long long currentTraces;

		protected:
	
	};

}
