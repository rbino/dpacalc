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
#include "dpacalc.h"
#include "hammingweight_prog.hpp"
void GeneratePowerModelProg::hammingweight_prog::fill ( shared_ptr<IntermediateValueMatrix>& intval, shared_ptr<PowerModelMatrix>& pm, unsigned long startTrace, unsigned int step )
{
	std::bitset<sizeof ( IntermediateValueType ) *8> tmp;
	IntermediateValueType value;
    for ( unsigned long d = startTrace; d < startTrace + step; d++ ) {
		for ( unsigned long long k = 0; k < KEYNUM; k++ ) {
			value = ( *intval ) ( d, k );
			BufferToBitset<sizeof ( IntermediateValueType ) > ( ( char* ) &value, tmp );
			( *pm ) ( d, k ) = tmp.count();
		}
	}
}

void GeneratePowerModelProg::hammingweight_prog::progressiveGenerate ( shared_ptr<IntermediateValueMatrix>& intval, shared_ptr<PowerModelMatrix>& pm, unsigned int step )
{
    if (pm.get() == NULL){
        pm.reset(new PowerModelMatrix ( step, KEYNUM ));
    }
    fill(intval, pm, 0, step);
    return;
    /*
    unsigned long curTraces = pm->rows();
    pm->conservativeResize(curTraces+step,NoChange);
    fill(intval, pm, curTraces, step);
    */
}
