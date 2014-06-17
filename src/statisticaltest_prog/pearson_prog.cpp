/*
Copyright (C) 2014	Riccardo Binetti

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
#include "pearson_prog.hpp"

void StatisticProg::pearson_prog::init ( shared_ptr< PowerModelMatrix >& _pm )
{
	StatisticProg::base::init ( _pm );
	Eigen::Matrix<TraceValueType, 1, Dynamic> pmaverage = Matrix<TraceValueType, 1, Dynamic> (1,KEYNUM);
	pmexpect = Matrix<TraceValueType, Dynamic, Dynamic> ( pm->rows(), KEYNUM );
	pmexpect_bykey = Matrix<TraceValueType, Dynamic, Dynamic> ( 1, KEYNUM );
	for ( long d = 0; d < KEYNUM; d++ ) {
		pmaverage ( 0, d ) = pm->col ( d ).array().sum() / pm->col ( d ).array().count();
	}
	for ( long d = 0; d < pm->cols(); d++ ) {
		pmexpect.col ( d ) = pm->col ( d ).array() - pmaverage ( 0, d );
	}
	for ( long keyh = 0; keyh < KEYNUM; keyh++ ) {
		pmexpect_bykey ( 0, keyh ) = pmexpect.col ( keyh ).squaredNorm();
	}
}

void StatisticProg::pearson_prog::generate ( shared_ptr<StatisticIndexMatrix>& stat, shared_ptr< TracesMatrix >& traces, long unsigned int numvalid )
{
	/*PSEUDOCODE:
	 For each time instant:
	 Calculating averages avg(t(j))
	 For each time and key hypotesis:
	stat(j,l) = dividend/divisor
	dividend = sum(i){t(i,j)-mediat(j))*pmexpect(i,l)}
	divisor = sqrt(sum(i)^2{t(i,j)-*mediat(j)} * sum(i){pmexpect_squared(i,l)} )
	Real code calculates all dividens and divisors together for each row, and then divides entire rows,
	to enable Eigen to better optimize those operations.
	*/
	assert ( numvalid <= BATCH_SIZE );
	TraceValueType tavg;
	auto dividendi = shared_ptr<Matrix<StatisticValueType, 1, Dynamic> > ( new Matrix<StatisticValueType, 1, Dynamic>(1,KEYNUM) );
	auto divisori = shared_ptr<Matrix<StatisticValueType, 1, Dynamic> > ( new Matrix<StatisticValueType, 1, Dynamic>(1,KEYNUM) );
	for ( unsigned long long time = 0; time < numvalid; time++ ) {
		tavg = traces->col ( time ).array().sum() / traces->col ( time ).array().count();
		( *divisori ) = ( ( traces->col ( time ).array() - tavg ).matrix().squaredNorm() * pmexpect_bykey ).array().sqrt();
		( *dividendi ).noalias() = ( traces->col ( time ).array() - tavg ).matrix().transpose() * pmexpect;
		stat->row ( time ) = ( *dividendi ).array() / ( *divisori ).array();
	}
}

