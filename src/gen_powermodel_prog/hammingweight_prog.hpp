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
#include "base.hpp"
using namespace Eigen;
using namespace std;
namespace GeneratePowerModelProg
{
	class hammingweight_prog : public base
	{
		public:
            virtual void progressiveGenerate ( shared_ptr<IntermediateValueMatrix>& intval, shared_ptr<PowerModelMatrix>& pm, unsigned int step );
			hammingweight_prog ( TCLAP::CmdLine& cmd ) : base ( cmd ) {}
        protected:
            virtual void fill ( shared_ptr<IntermediateValueMatrix>& intval, shared_ptr<PowerModelMatrix>& pm, unsigned long startTrace, unsigned int step );
	};
}

