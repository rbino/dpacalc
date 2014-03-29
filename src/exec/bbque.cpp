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
#include "bbque.hpp"
#include "dpacalc.hpp"
void prefetchthread();
void threadfunction();

void ExecMethod::bbque::RunAndWait ( unsigned long numberoftimes )
{
	std::thread prefetcht = thread ( prefetchthread );
	RTLIB_Services_t* rtlib;
	RTLIB_Init ( "DPAcalc", &rtlib );
	assert ( rtlib );
	shared_ptr<ExecMethod::BbqueClass> bc = shared_ptr<ExecMethod::BbqueClass> ( new ExecMethod::BbqueClass ( "DPAcalc", "DPAcalcRecipe", rtlib, numberoftimes ) );
	assert ( bc->isRegistered() );
	bc->Start();
	bc->WaitCompletion();
	prefetcht.join();
	//BbqueClass destruction happens here thanks to shared_ptr
}

ExecMethod::BbqueClass::BbqueClass ( const string& name, const string& recipe, RTLIB_Services_t* rtlib, unsigned long _tot ) : BbqueEXC ( name, recipe, rtlib ), mtx()
{
	tot = _tot;
	num = 0;
}
RTLIB_ExitCode_t ExecMethod::BbqueClass::onConfigure ( uint8_t awm_id )
{
	numthreads = awm_id + 1;
	return RTLIB_OK;
}
RTLIB_ExitCode_t ExecMethod::BbqueClass::onMonitor()
{
	return RTLIB_OK;
}
RTLIB_ExitCode_t ExecMethod::BbqueClass::onRun()
{
	unsigned int oldnum;
	mtx.lock();
	if ( num >= tot ) {
		mtx.unlock();
		return RTLIB_EXC_WORKLOAD_NONE;
	}
	oldnum = num;
	vector<thread> thrs;
	for ( ; num < oldnum + numthreads; num++ ) {
		if ( num >= tot ) {
			break;
		}
		thrs.push_back ( thread ( threadfunction ) );
	}
	mtx.unlock();
	for ( auto t = thrs.begin(); t != thrs.end(); ++t ) {
		t->join();
	}
	return RTLIB_OK;
}
RTLIB_ExitCode_t ExecMethod::BbqueClass::onSetup()
{
	return RTLIB_OK;
}

void prefetchthread()
{
	DPA::instance()->prefetch();
}

void threadfunction()
{
	DPA::instance()->doRun();
}
