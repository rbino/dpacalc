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
#pragma once
#include <bitset>
#include <memory>
#include <iostream>
#include <Eigen/Dense>
#include <tclap/CmdLine.h>
#include "config.h"
#include "utils.hpp"

typedef std::bitset<DATA_SIZE_BIT> DataValueType;



#define DATA_SIZE_BYTE (DATA_SIZE_BIT >> 3)
#define KEY_SIZE_BYTE (KEY_SIZE_BIT >> 3)
#define KEYNUM (1 <<KEY_HYP_BIT)

typedef std::bitset<KEY_SIZE_BIT> KeyDataType;
typedef unsigned long KeyIndexType;
typedef unsigned long long IntermediateValueType;
typedef TraceValueType StatisticValueType;

typedef Eigen::Matrix<TraceValueType, Eigen::Dynamic, BATCH_SIZE> TracesMatrix;
typedef std::vector< DataValueType > DataMatrix;
typedef Eigen::Matrix<TraceValueType, 1, Eigen::Dynamic> Trace;
typedef Eigen::Matrix<std::complex<TraceValueType>, 1, Eigen::Dynamic> ComplexTrace;
typedef Eigen::Matrix<IntermediateValueType, Eigen::Dynamic, Eigen::Dynamic> IntermediateValueMatrix;
typedef Eigen::Matrix<TraceValueType, Eigen::Dynamic, Eigen::Dynamic> PowerModelMatrix;
typedef Eigen::Matrix<StatisticValueType, Eigen::Dynamic, Eigen::Dynamic> StatisticIndexMatrix;
typedef std::pair<StatisticValueType,StatisticValueType> ConfidencePair;
typedef std::pair<float,float> FilterBand;
typedef struct {
    FilterBand band;
    unsigned int ntraces;
} BandWithTraces;
typedef struct {
    std::shared_ptr<DataValueType> data;
    std::shared_ptr<Trace> trace;
} TraceWithData;
#pragma pack(push)
#pragma pack(1)
struct fileheaders {
    uint32_t numtraces;
    uint32_t numsamples_per_trace;
    char datatype;
    uint8_t knowndatalength;
};
#pragma pack(pop)



