/*
 * Copyright 1993-2014 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */



////////////////////////////////////////////////////////////////////////////////
// Shortcut definition
////////////////////////////////////////////////////////////////////////////////
typedef unsigned int uint;



///////////////////////////////////////////////////////////////////////////////
// Sort result validation routines
////////////////////////////////////////////////////////////////////////////////
//Sorted keys array validation (check for integrity and proper order)
extern "C" uint validateSortedKeys(
    uint *resKey,
    uint *srcKey,
    uint batchSize,
    uint arrayLength,
    uint numValues,
    uint dir
);

extern "C" int validateValues(
    uint *resKey,
    uint *resVal,
    uint *srcKey,
    uint batchSize,
    uint arrayLength
);



////////////////////////////////////////////////////////////////////////////////
// CUDA sorting networks
////////////////////////////////////////////////////////////////////////////////

extern "C" uint bitonicSortShared(
	float *d_DstKey,
	uint *d_DstVal,
	float *d_SrcKey,
	uint *d_SrcVal,
	uint batchSize,
	uint arrayLength,
	uint dir
);

extern "C" uint bitonicSort(
	float *d_DstKey,
    uint *d_DstVal,
	float *d_SrcKey,
    uint *d_SrcVal,
    uint batchSize,
    uint arrayLength,
    uint dir
);

extern "C" void oddEvenMergeSort(
	float *d_DstKey,
    uint *d_DstVal,
	float *d_SrcKey,
    uint *d_SrcVal,
    uint batchSize,
    uint arrayLength,
    uint dir
);
