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
#include "base.hpp"
using namespace Eigen;
using namespace std;
namespace GenerateIntermediateValuesProg
{
	class aes128round1_prog : public base
	{
		public:

            virtual void progressiveGenerate ( shared_ptr<DataMatrix>& knowndata, shared_ptr<IntermediateValueMatrix>& intval, unsigned int step );
			virtual void init();
			aes128round1_prog ( TCLAP::CmdLine& cmd, shared_ptr<KeyGenerators::base> _keygen ) : base ( cmd, _keygen ),
				whichsboxArg ( "b", "sbox", "From which SBOX output should I start to correlate?", false, 0, "0-15" ),
				sboxnumArg ( "v", "sboxnum", "How many consecutive SBOXes should I consider?", false, 1, "1-8" ) {
				cmd.add ( whichsboxArg );
				cmd.add ( sboxnumArg );
			}
		protected:
            virtual void fill ( shared_ptr<DataMatrix>& knowndata, shared_ptr<IntermediateValueMatrix>& intval, unsigned long startTrace, unsigned int step );
			TCLAP::ValueArg<int> whichsboxArg;
			TCLAP::ValueArg<int> sboxnumArg;
	};
}

