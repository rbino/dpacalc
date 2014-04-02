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
#include "dpacalc.h"
#include "aes128round1_prog.hpp"
#include "aes.h"
void GenerateIntermediateValuesProg::aes128round1_prog::init()
{
	if ( whichsboxArg.getValue() + sboxnumArg.getValue() > AES_STATE_BYTES_NO || whichsboxArg.getValue() < 0 || sboxnumArg.getValue() < 1 ) {
		cerr << "Maximum number of SBOXes exceeded." << endl;
		exit ( 1 );
	}
}

void GenerateIntermediateValuesProg::aes128round1_prog::generate ( shared_ptr<DataMatrix>& knowndata, shared_ptr<IntermediateValueMatrix>& intval )
{
	std::bitset<KEY_SIZE_BIT> key;
	uint8_t fullaeskey[AES_STATE_BYTES_NO];
	uint8_t fullaesdata[AES_STATE_BYTES_NO];
	memset ( fullaeskey, 0, AES_STATE_BYTES_NO );
	uint8_t UnrolledRoundKey[AES_128_ROUND_KEY_NO][AES_STATE_BYTES_NO];
	uint8_t* dataptr = &fullaesdata[whichsboxArg.getValue()];
	for ( KeyIndexType keyidx = 0; keyidx < KEYNUM ; ++keyidx ) {
		key = keygen->getKeyFromIndex ( keyidx );
		BitsetToBuffer<KEY_SIZE_BYTE> ( key, ( char* ) &fullaeskey );
		KeyExpansion ( fullaeskey, UnrolledRoundKey );
		for ( unsigned long trcidx = 0; trcidx < knowndata->size(); trcidx++ ) {
			BitsetToBuffer<DATA_SIZE_BYTE> ( ( *knowndata ) [trcidx], ( char* ) &fullaesdata );
			AddRoundKey ( fullaesdata, UnrolledRoundKey[0] );
			( *intval ) ( trcidx, keyidx ) = 0;
			switch ( sboxnumArg.getValue() ) {
				case 8:
					( *intval ) ( trcidx, keyidx ) += ( ( IntermediateValueType ) SBOX[ ( * ( dataptr + 7 ) )] ) << 56;
				case 7:
					( *intval ) ( trcidx, keyidx ) += ( ( IntermediateValueType ) SBOX[ ( * ( dataptr + 6 ) )] ) << 48;
				case 6:
					( *intval ) ( trcidx, keyidx ) += ( ( IntermediateValueType ) SBOX[ ( * ( dataptr + 5 ) )] ) << 40;
				case 5:
					( *intval ) ( trcidx, keyidx ) += ( ( IntermediateValueType ) SBOX[ ( * ( dataptr + 4 ) )] ) << 32;
				case 4:
					( *intval ) ( trcidx, keyidx ) += ( ( IntermediateValueType ) SBOX[ ( * ( dataptr + 3 ) )] ) << 24;
				case 3:
					( *intval ) ( trcidx, keyidx ) += ( ( IntermediateValueType ) SBOX[ ( * ( dataptr + 2 ) )] ) << 16;
				case 2:
					( *intval ) ( trcidx, keyidx ) += ( ( IntermediateValueType ) SBOX[ ( * ( dataptr + 1 ) )] ) << 8;
				case 1:
					( *intval ) ( trcidx, keyidx ) += ( ( IntermediateValueType ) SBOX[ ( * ( dataptr ) )] );
			}
		}
	}
}

