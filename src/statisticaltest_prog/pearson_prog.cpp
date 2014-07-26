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

void StatisticProg::pearson_prog::init ( shared_ptr< PowerModelMatrix >& _pm, unsigned int step, unsigned long nbatch )
{
    StatisticProg::base::init ( _pm, step, nbatch);

    curStep = step;
    curNtraces += step;

    if (first){
        first = false;
        sum_hyp.reset(new Matrix<TraceValueType, 1, Dynamic> (1,KEYNUM));
        sum_hypSquared.reset(new Matrix<TraceValueType, 1, Dynamic> (1,KEYNUM));
        for ( long d = 0; d < KEYNUM; d++ ){
            (*sum_hyp) ( 0, d ) = pm->col (d).topRows(step).array().sum();
            (*sum_hypSquared) ( 0, d ) = pm->col (d).topRows(step).squaredNorm();
        }
        sum_traces.resize(nbatch);
        sum_tracesSquared.resize(nbatch);
        sum_tracehyp.resize(nbatch);
    } else {
        for ( long d = 0; d < KEYNUM; d++ ){
            (*sum_hyp) ( 0, d ) += pm->col (d).topRows(step).array().sum();
            (*sum_hypSquared) ( 0, d ) += pm->col (d).topRows(step).squaredNorm();
        }
    }
}

void StatisticProg::pearson_prog::progressiveGenerate ( shared_ptr<StatisticIndexMatrix>& stat, shared_ptr< TracesMatrix >& traces, long unsigned int numvalid, unsigned long long id )
{
    assert ( numvalid <= BATCH_SIZE );
    if(sum_traces[id].get() == NULL){
        sum_traces[id].reset(new Matrix<StatisticValueType, Dynamic, Dynamic> (1,numvalid));
        sum_tracesSquared[id].reset(new Matrix<StatisticValueType, Dynamic, Dynamic> (1,numvalid));
        sum_tracehyp[id].reset(new Matrix<StatisticValueType, Dynamic, Dynamic> (numvalid,KEYNUM));
        for ( unsigned long long time = 0; time < numvalid; time++ ) {
            (*sum_traces[id]) (0, time) = traces->col( time ).array().sum();
            (*sum_tracesSquared[id]) (0, time) = traces->col( time ).squaredNorm();
        }
        (*sum_tracehyp[id]) = traces->transpose().topRows(numvalid) * pm->topRows(curStep);
    } else {
        for ( unsigned long long time = 0; time < numvalid; time++ ) {
            (*sum_traces[id]) (0, time) += traces->col( time ).array().sum();
            (*sum_tracesSquared[id]) (0, time) += traces->col( time ).squaredNorm();
        }
        (*sum_tracehyp[id]) += traces->transpose().topRows(numvalid) * pm->topRows(curStep);
    }

    auto mean_traces = shared_ptr<Matrix<StatisticValueType, 1, Dynamic> > ( new Matrix<StatisticValueType, 1, Dynamic>(1,numvalid) );
    auto mean_hyp = shared_ptr<Matrix<StatisticValueType, 1, Dynamic> > ( new Matrix<StatisticValueType, 1, Dynamic>(1,KEYNUM) );
    (*mean_traces) = (*sum_traces[id]) / curNtraces;
    (*mean_hyp) = (*sum_hyp) / curNtraces;
    auto variance_traces = shared_ptr<Matrix<StatisticValueType, 1, Dynamic> > ( new Matrix<StatisticValueType, 1, Dynamic>(1,numvalid) );
    auto variance_hyp = shared_ptr<Matrix<StatisticValueType, 1, Dynamic> > ( new Matrix<StatisticValueType, 1, Dynamic>(1,KEYNUM) );
    (*variance_traces) = (*sum_tracesSquared[id]) + (mean_traces->array().square().matrix() * curNtraces)  - (2  * mean_traces->cwiseProduct(*sum_traces[id]));
    (*variance_hyp)= (*sum_hypSquared) + (mean_hyp->array().square().matrix() * curNtraces)  - (2  * mean_hyp->cwiseProduct(*sum_hyp));
    for (unsigned int i=0; i < KEYNUM; i++){
        if ((*variance_hyp) (0, i) == 0){
            (*variance_hyp) (0, i) = numeric_limits<StatisticValueType>::infinity();
        }
    }
    for (unsigned int i=0; i < numvalid; i++){
        if ((*variance_traces) (0, i) == 0){
            (*variance_traces) (0, i) = numeric_limits<StatisticValueType>::infinity();
        }
    }
    auto joint_variance = shared_ptr<Matrix<StatisticValueType, Dynamic, Dynamic> > ( new Matrix<StatisticValueType, Dynamic, Dynamic>(numvalid,KEYNUM) );
    (*joint_variance) = (variance_traces->transpose() * (*variance_hyp)).array().sqrt();
    (*stat) = ((*sum_tracehyp[id]) + ((mean_traces->transpose() * (*mean_hyp)) * curNtraces) - (mean_traces->transpose() * (*sum_hyp)) - (sum_traces[id]->transpose()*(*mean_hyp))).cwiseQuotient(*joint_variance);
}

