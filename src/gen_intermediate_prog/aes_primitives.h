#pragma once
#include <stdint.h>
#include "aes_constants.h"
#include "aes_helpers.h"


void MixColumns ( uint8_t InValue[] )
{
	ColumnStructure_t TempInColumn;
	ColumnStructure_t TempOutColumn;
	uint8_t   i;
	for ( i = 0; i < 4; i++ ) {
		TempInColumn.E1 = InValue[i * 4];
		TempInColumn.E2 = InValue[i * 4 + 1];
		TempInColumn.E3 = InValue[i * 4 + 2];
		TempInColumn.E4 = InValue[i * 4 + 3];
		TempOutColumn = OneMixColumn ( TempInColumn );
		InValue[i * 4] = TempOutColumn.E1;
		InValue[i * 4 + 1] = TempOutColumn.E2;
		InValue[i * 4 + 2] = TempOutColumn.E3;
		InValue[i * 4 + 3] = TempOutColumn.E4;
	}
}


void AddRoundKey ( uint8_t InValue[], uint8_t RoundKey[] )
{
	uint8_t i;
	for ( i = 0; i < 16; i++ )
	{ InValue[i] ^= RoundKey[i]; }
}

void SubBytes ( uint8_t InValue[] )
{
	uint8_t i;
	for ( i = 0; i < 16; i++ )
	{ InValue[i] = SBOX[InValue[i]]; }
}

void ShiftRows ( uint8_t InValue[] )
{
	uint8_t TempValue[16];
	uint8_t i;
	for ( i = 0; i < 16; i++ )
	{ TempValue[i] = InValue[i]; }
	for ( i = 0; i < 16; i++ )
	{ InValue[i] = TempValue[ ( i + ( ( i % 4 ) * 4 ) ) % 16]; }
}

void InShiftRows ( uint8_t InValue[] )
{
	uint8_t TempValue[16];
	uint8_t i;
	for ( i = 0; i < 16; i++ )
	{ TempValue[i] = InValue[i]; }
	for ( i = 0; i < 16; i++ )
	{ InValue[i] = TempValue[ ( i + 16 - ( ( i % 4 ) * 4 ) ) % 16]; }
}

void InSubBytes ( uint8_t InValue[] )
{
	uint8_t i;
	for ( i = 0; i < 16; i++ )
	{ InValue[i] = InSBOX[InValue[i]]; }
}

void InMixColumns ( uint8_t InValue[] )
{
	ColumnStructure_t TempInColumn;
	ColumnStructure_t TempOutColumn;
	uint8_t   i;
	for ( i = 0; i < 4; i++ ) {
		TempInColumn.E1 = InValue[i * 4];
		TempInColumn.E2 = InValue[i * 4 + 1];
		TempInColumn.E3 = InValue[i * 4 + 2];
		TempInColumn.E4 = InValue[i * 4 + 3];
		TempOutColumn = OneInMixColumn ( TempInColumn );
		InValue[i * 4] = TempOutColumn.E1;
		InValue[i * 4 + 1] = TempOutColumn.E2;
		InValue[i * 4 + 2] = TempOutColumn.E3;
		InValue[i * 4 + 3] = TempOutColumn.E4;
	}
}
#pragma GCC diagnostic warning "-Warray-bounds"
void KeyExpansion ( uint8_t Key[], uint8_t RoundKey[11][16] )
{
	uint8_t i;
	uint8_t j;
	uint8_t Temp;
	uint8_t RConCount;
	RConCount = 1;
	for ( i = 0; i < 16; i++ )
	{ RoundKey[0][i] = Key[i]; }
	//Every time i change I'm in a new set of 16 bytes
	for ( i = 1; i < 11; i++ ) {
		for ( j = 0; j < 4; j++ ) {
			//Use RoundKey as a temp variable
			//It is not important if 12 + j*4 exceeds the bound of 16
			//it stays anyway in the proper area of the variable RoundKey
			RoundKey[i][j * 4] = RoundKey[i - 1][12 + j * 4];
			RoundKey[i][j * 4 + 1] = RoundKey[i - 1][12 + j * 4 + 1];
			RoundKey[i][j * 4 + 2] = RoundKey[i - 1][12 + j * 4 + 2];
			RoundKey[i][j * 4 + 3] = RoundKey[i - 1][12 + j * 4 + 3];
			//Group of four byte in position multiple of 4
			//Need to apply SubWord and RotateWord
			if ( j == 0 ) {
				Temp = RoundKey[i][j * 4];
				RoundKey[i][j * 4] = SBOX[RoundKey[i][j * 4 + 1]] ^ RCon[RConCount];
				RoundKey[i][j * 4 + 1] = SBOX[RoundKey[i][j * 4 + 2]];
				RoundKey[i][j * 4 + 2] = SBOX[RoundKey[i][j * 4 + 3]];
				RoundKey[i][j * 4 + 3] = SBOX[Temp];
				RConCount++;
			}
			RoundKey[i][j * 4] ^= RoundKey[i - 1][j * 4];
			RoundKey[i][j * 4 + 1] ^= RoundKey[i - 1][j * 4 + 1];
			RoundKey[i][j * 4 + 2] ^= RoundKey[i - 1][j * 4 + 2];
			RoundKey[i][j * 4 + 3] ^= RoundKey[i - 1][j * 4 + 3];
		}
	}
}
