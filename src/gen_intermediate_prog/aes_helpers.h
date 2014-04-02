/*
This file includes helper functions and the full polynomial multiplication table
*/

#include <stdint.h>
#include <stdio.h>
#include "aes_structs.h"

uint8_t Mult2[256];
uint8_t Mult3[256];
uint8_t Multe[256];
uint8_t Mult9[256];
uint8_t Multd[256];
uint8_t Multb[256];



void PrintVectorFile ( uint8_t* DaStampare, FILE* fp )
{
	int i;
	for ( i = 0; i < 16; i++ )
	{ fprintf ( fp, "%02X ", DaStampare[i] ); }
	printf ( "\n" );
}

int checkAllZero ( uint8_t diff[16] )
{
	int i = 0;
	for ( i = 0; i < 16; i++ )
		if ( diff[i] != 0 )
		{ return 0; }
	return 1;
}


void printStateToFile ( unsigned char state[16], FILE* fp )
{
	fprintf ( fp, "|%02X %02X %02X %02X|\n", state[0], state[4], state[8], state[12] );
	fprintf ( fp, "|%02X %02X %02X %02X|\n", state[1], state[5], state[9], state[13] );
	fprintf ( fp, "|%02X %02X %02X %02X|\n", state[2], state[6], state[10], state[14] );
	fprintf ( fp, "|%02X %02X %02X %02X|\n", state[3], state[7], state[11], state[15] );
	fprintf ( fp, "------------------------------\n" );
}

uint8_t returnDecimalFromString ( char hex )
{
	if ( hex == 'a' )
	{ return 0; }
	else if ( hex == '1' )
	{ return 1; }
	else if ( hex == '2' )
	{ return 2; }
	else if ( hex == '3' )
	{ return 3; }
	else if ( hex == '4' )
	{ return 4; }
	else if ( hex == '5' )
	{ return 5; }
	else if ( hex == '6' )
	{ return 6; }
	else if ( hex == '7' )
	{ return 7; }
	else if ( hex == '8' )
	{ return 8; }
	else if ( hex == '9' )
	{ return 9; }
	else if ( hex == 'A' )
	{ return 10; }
	else if ( hex == 'B' )
	{ return 11; }
	else if ( hex == 'C' )
	{ return 12; }
	else if ( hex == 'D' )
	{ return 13; }
	else if ( hex == 'E' )
	{ return 14; }
	else if ( hex == 'F' )
	{ return 15; }
	else if ( hex == 'a' )
	{ return 10; }
	else if ( hex == 'b' )
	{ return 11; }
	else if ( hex == 'c' )
	{ return 12; }
	else if ( hex == 'd' )
	{ return 13; }
	else if ( hex == 'e' )
	{ return 14; }
	else if ( hex == 'f' )
	{ return 15; }
	return 0;
}

void PrintVector ( uint8_t* DaStampare )
{
	int i;
	for ( i = 0; i < 16; i++ )
	{ printf ( "%02X ", DaStampare[i] ); }
	printf ( "\n" );
}

ColumnStructure_t OneMixColumn ( ColumnStructure_t InValue )
{
	ColumnStructure_t OutValue;
	OutValue.E1 = Mult2[InValue.E1] ^ Mult3[InValue.E2] ^ InValue.E3 ^ InValue.E4;
	OutValue.E2 = InValue.E1 ^ Mult2[InValue.E2] ^ Mult3[InValue.E3] ^ InValue.E4;
	OutValue.E3 = InValue.E1 ^ InValue.E2 ^ Mult2[InValue.E3] ^ Mult3[InValue.E4];
	OutValue.E4 = Mult3[InValue.E1] ^ InValue.E2 ^ InValue.E3 ^ Mult2[InValue.E4];
	return ( OutValue );
}

ColumnStructure_t OneInMixColumn ( ColumnStructure_t InValue )
{
	ColumnStructure_t OutValue;
	OutValue.E1 = Multe[InValue.E1] ^ Multb[InValue.E2] ^ Multd[InValue.E3] ^ Mult9[InValue.E4];
	OutValue.E2 = Mult9[InValue.E1] ^ Multe[InValue.E2] ^ Multb[InValue.E3] ^ Multd[InValue.E4];
	OutValue.E3 = Multd[InValue.E1] ^ Mult9[InValue.E2] ^ Multe[InValue.E3] ^ Multb[InValue.E4];
	OutValue.E4 = Multb[InValue.E1] ^ Multd[InValue.E2] ^ Mult9[InValue.E3] ^ Multe[InValue.E4];
	return ( OutValue );
}

uint8_t PolyMultiply ( uint8_t FirstPoly, uint8_t SecondPoly )
{
	uint8_t p = 0;
	uint8_t Counter;
	uint8_t Hi_bit_set;
	for ( Counter = 0; Counter < 8; Counter++ ) {
		if ( ( SecondPoly & 1 ) == 1 )
		{ p ^= FirstPoly; }
		Hi_bit_set = ( FirstPoly & 0x80 );
		FirstPoly <<= 1;
		if ( Hi_bit_set == 0x80 )
		{ FirstPoly ^= 0x1b; }
		SecondPoly >>= 1;
	}
	return p;
}


void FillHelperData()
{
	uint16_t i;
	for ( i = 0; i < 256; i++ ) {
		Mult2[i] = PolyMultiply ( 2, i );
		Mult3[i] = PolyMultiply ( 3, i );
		Multe[i] = PolyMultiply ( 0xe, i );
		Mult9[i] = PolyMultiply ( 9, i );
		Multd[i] = PolyMultiply ( 0xd, i );
		Multb[i] = PolyMultiply ( 0xb, i );
	}
}

