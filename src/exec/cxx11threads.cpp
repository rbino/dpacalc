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
void prefetchthread();
void threadfunction();
void filterfunction();
unsigned long batchmax = 0;
unsigned long batchcur = 0;
unsigned long filtercur = 0;
unsigned long filtermax = 0;
mutex cur_mutex;
mutex filter_mutex;
void ExecMethod::cxx11threads::RunAndWait ( unsigned long numberoftimes )
{
	batchmax = numberoftimes;
	batchcur = 0;
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

void ExecMethod::cxx11threads::RunFilter (unsigned long numberoftimes){
    filtercur = 0;
    filtermax = numberoftimes;
    int numCPU = procArg.getValue();
    if ( numCPU == 0 ) {
        numCPU = sysconf ( _SC_NPROCESSORS_ONLN );
    }
    vector<thread> thrs;
    for ( int i = 0; i < numCPU; i++ ) {
        thrs.push_back ( thread ( filterfunction ) );
    }
    for ( auto t = thrs.begin(); t != thrs.end(); ++t ) {
        t->join();
    }
}

void prefetchthread()
{
	DPA::instance()->prefetch();
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
		DPA::instance()->doRun();
	}
}

void filterfunction()
{
    for ( ;; ) {
        filter_mutex.lock();
        if (filtercur < filtermax) {
            ++filtercur;
        } else {
            filter_mutex.unlock();
            break;
        }
        filter_mutex.unlock();
        DPA::instance()->doFilter();
    }
}
