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
#include <mutex>
#include "input_find/base.hpp"
#include "filter_find/base.hpp"
#include "keygen/base.hpp"
#include "gen_intermediate_prog/base.hpp"
#include "gen_powermodel_prog/base.hpp"
#include "statisticaltest_prog/base.hpp"
#include "exec/base.hpp"
#include "verification/base.hpp"
#include "output_find/base.hpp"

using namespace Eigen;
using namespace std;

class DPA
{
	public:

		int main ( int argc, char** argv );
		static DPA* instance() {
			static DPA theInstance;
			return &theInstance;
        }
	protected:
		shared_ptr<DataMatrix> data;
		shared_ptr<IntermediateValueMatrix> intval;
		shared_ptr<PowerModelMatrix> pm;
		unsigned long numbatches;
		shared_ptr<ExecMethod::base> exec;
        shared_ptr<SamplesInputFind::base> input;
        shared_ptr<FilterFind::base> filter;
		shared_ptr<KeyGenerators::base> keygen;
		shared_ptr<GenerateIntermediateValuesProg::base> interm;
		shared_ptr<GeneratePowerModelProg::base> genpm;
        shared_ptr<StatisticProg::base> stat;
        shared_ptr<VerifyAttack::base> verif;
        shared_ptr<OutputFind::base> outp;
		virtual void ShowCompileTimeOptions();
        unsigned long curTrace;
        unsigned long successfulTraces;
        mutex traceMutex;
	private:
		DPA() {}
};
