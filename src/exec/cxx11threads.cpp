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
#include "cxx11threads.hpp"
#include <thread>
#include <unistd.h>
#include "main.hpp"
#include <functional>
void prefetchthread();
void threadfunction();
void filterfunction();
unsigned long batchmax = 0;
unsigned long batchcur = 0;
static std::function<void()> fun1 = NULL;
static std::function<void()> fun2 = NULL;
mutex cur_mutex;
mutex filter_mutex;
void ExecMethod::cxx11threads::RunAndWait ( unsigned long numberoftimes, std::function<void()> f1, std::function<void()>  f2 )
{
	batchmax = numberoftimes;
	batchcur = 0;
    fun1 = f1;
    fun2 = f2;
	int numCPU = procArg.getValue();
	if ( numCPU == 0 ) {
		numCPU = sysconf ( _SC_NPROCESSORS_ONLN );
	}
	vector<thread> thrs;
	thrs.push_back ( thread ( prefetchthread ) );
	for ( int i = 0; i < numCPU; i++ ) {
		thrs.push_back ( thread ( threadfunction ) );
	}
	for ( auto t = thrs.begin(); t != thrs.end(); ++t ) {
		t->join();
	}
}

void prefetchthread()
{
    if (fun2) {
        fun2();
    }
}

void threadfunction()
{
	for ( ;; ) {
		cur_mutex.lock();
		if ( batchcur < batchmax ) {
			++batchcur;
		} else {
			cur_mutex.unlock();
			break;
		}
		cur_mutex.unlock();
        if (fun1){
            fun1();
        }
	}
}
