#include "aes_primitives.h"

/**
* Runs a full AES 128 encryption, requires an unrolled key
*/

void AES_128 ( uint8_t InValue[], uint8_t UnrolledRoundKey[AES_128_ROUND_KEY_NO][AES_STATE_BYTES_NO], uint8_t OutValue[] )
{
	uint8_t i;
	for ( i = 0; i < 16; i++ )
	{ OutValue[i] = InValue[i]; }
	AddRoundKey ( OutValue, UnrolledRoundKey[0] );
	for ( i = 1; i <= 10; i++ ) {
		SubBytes ( OutValue );
		ShiftRows ( OutValue );
		if ( i != 10 )
		{ MixColumns ( OutValue ); }
		AddRoundKey ( OutValue, UnrolledRoundKey[i] );
	}
}

/**
* Runs a full AES 192 encryption, requires an unrolled key
*/

void AES_192 ( uint8_t InValue[], uint8_t UnrolledRoundKey[AES_128_ROUND_KEY_NO][AES_STATE_BYTES_NO], uint8_t OutValue[] )
{
	uint8_t i;
	for ( i = 0; i < 16; i++ )
	{ OutValue[i] = InValue[i]; }
	AddRoundKey ( OutValue, UnrolledRoundKey[0] );
	for ( i = 1; i <= 12; i++ ) {
		SubBytes ( OutValue );
		ShiftRows ( OutValue );
		if ( i != 12 )
		{ MixColumns ( OutValue ); }
		AddRoundKey ( OutValue, UnrolledRoundKey[i] );
	}
}

/**
* Runs a full AES 256 encryption, requires an unrolled key
*/

void AES_256 ( uint8_t InValue[], uint8_t UnrolledRoundKey[AES_256_ROUND_KEY_NO][AES_STATE_BYTES_NO], uint8_t OutValue[] )
{
	uint8_t i;
	for ( i = 0; i < 16; i++ )
	{ OutValue[i] = InValue[i]; }
	AddRoundKey ( OutValue, UnrolledRoundKey[0] );
	for ( i = 1; i <= 14; i++ ) {
		SubBytes ( OutValue );
		ShiftRows ( OutValue );
		if ( i != 14 )
		{ MixColumns ( OutValue ); }
		AddRoundKey ( OutValue, UnrolledRoundKey[i] );
	}
}
