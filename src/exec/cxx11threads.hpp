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
#include <functional>

namespace ExecMethod
{
	class cxx11threads: public base
	{
		public:

            cxx11threads ( TCLAP::CmdLine& cmd ) : base ( cmd ), procArg ( "p", "processors", "Number of processors to use (0=autodetect)", false, 0, "int" ) {
				cmd.add ( procArg );
			};
            /**
             * @brief RunAndWait Runs f1 once and f2 numberoftimes times in different threads
             * @param numberoftimes The number of times to run f2
             * @param f1 The function to be executed once in a single thread. If it's NULL, nothing it's executed in its thread.
             * @param f2 The function to be executed numberoftimes times in a thread per CPU core. If it's NULL, nothing it's executed in its threads.
             */
            virtual void RunAndWait ( unsigned long numberoftimes, std::function<void()> f1, std::function<void()>  f2 );
		protected:
			TCLAP::ValueArg<int> procArg;
	};
}
