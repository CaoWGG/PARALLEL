#include <stdio.h>
#include <algorithm>
#include <math.h>
#include "CMU418intrin.h"
#include "logger.h"
using namespace std;


void absSerial(float* values, float* output, int N) {
    for (int i=0; i<N; i++) {
	float x = values[i];
	if (x < 0) {
	    output[i] = -x;
	} else {
	    output[i] = x;
	}
    }
}

// implementation of absolute value using 15418 instrinsics
void absVector(float* values, float* output, int N) {
    __cmu418_vec_float x;
    __cmu418_vec_float result;
    __cmu418_vec_float zero = _cmu418_vset_float(0.f);
    __cmu418_mask maskAll, maskIsNegative, maskIsNotNegative;

    //  Note: Take a careful look at this loop indexing.  This example
    //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
    //  Why is that the case?
    for (int i=0; i<N; i+=VECTOR_WIDTH) {

	// All ones
	maskAll = _cmu418_init_ones();

	// All zeros
	maskIsNegative = _cmu418_init_ones(0);

	// Load vector of values from contiguous memory addresses
	_cmu418_vload_float(x, values+i, maskAll);               // x = values[i];

	// Set mask according to predicate
	_cmu418_vlt_float(maskIsNegative, x, zero, maskAll);     // if (x < 0) {

	// Execute instruction using mask ("if" clause)
	_cmu418_vsub_float(result, zero, x, maskIsNegative);      //   output[i] = -x;

	// Inverse maskIsNegative to generate "else" mask
	maskIsNotNegative = _cmu418_mask_not(maskIsNegative);     // } else {

	// Execute instruction ("else" clause)
	_cmu418_vload_float(result, values+i, maskIsNotNegative); //   output[i] = x; }

	// Write results back to memory
	_cmu418_vstore_float(output+i, result, maskAll);
    }
}

// Accepts an array of values and an array of exponents
// For each element, compute values[i]^exponents[i] and clamp value to
// 4.18.  Store result in outputs.
// Uses iterative squaring, so that total iterations is proportional
// to the log_2 of the exponent
void clampedExpSerial(float* values, int* exponents, float* output, int N) {
    for (int i=0; i<N; i++) {
	float x = values[i];
	float result = 1.f;
	int y = exponents[i];
	float xpower = x;
	while (y > 0) {
    	    if (y & 0x1) {
		result *= xpower;
		if (result > 4.18f) {
		    result = 4.18f;
		    break;
		}
            }
	    xpower = xpower * xpower;
	    y >>= 1;
	}
	output[i] = result;
    }
}

void clampedExpVector(float* values, int* exponents, float* output, int N) {
    // TODO: Implement your vectorized version of clampedExpSerial here
    //  ...
	__cmu418_vec_float x, xpower, res;
	__cmu418_vec_int y, e, isOneInt;
	__cmu418_mask isContinue, dataMask, tmp;
	
	// const value
	__cmu418_mask oneMask=_cmu418_init_ones(), zeroMask=_cmu418_init_ones(0);
	__cmu418_vec_int oneInt=_cmu418_vset_int(1), zeroInt=_cmu418_vset_int(0);
	__cmu418_vec_float vec418=_cmu418_vset_float(4.18f);


	for (int i=0;i<N;i+=VECTOR_WIDTH){
		dataMask=_cmu418_init_ones(N-i<VECTOR_WIDTH?N-i:VECTOR_WIDTH);
		isContinue=dataMask;
		_cmu418_vload_float(x,values+i,dataMask);
		_cmu418_vset_float(res,1.0f,dataMask);
		_cmu418_vload_int(y,exponents+i,dataMask);
		_cmu418_vload_float(xpower,values+i,dataMask);

		while (_cmu418_cntbits(isContinue)>0){
			// calculation
			_cmu418_vbitand_int(isOneInt,y,oneInt,dataMask);
			_cmu418_veq_int(tmp,isOneInt,oneInt,dataMask);
			_cmu418_vmult_float(res,res,xpower,tmp);

			_cmu418_vgt_float(tmp, res, vec418, dataMask);
			_cmu418_vset_float(res,4.18f,tmp);
			_cmu418_vmult_float(xpower,xpower,xpower,dataMask);
			_cmu418_vshiftright_int(y,y,oneInt,dataMask);

			tmp = _cmu418_mask_not(tmp);
			isContinue = _cmu418_mask_and(tmp,isContinue);
			_cmu418_vgt_int(tmp,y,zeroInt,isContinue);
			isContinue = _cmu418_mask_and(tmp,isContinue);
		}
		_cmu418_vstore_float(output+i, res, dataMask);
	}
}


float arraySumSerial(float* values, int N) {
    float sum = 0;
    for (int i=0; i<N; i++) {
	sum += values[i];
    }

    return sum;
}

// Assume N % VECTOR_WIDTH == 0
// Assume VECTOR_WIDTH is a power of 2
float arraySumVector(float* values, int N) {
    // TODO: Implement your vectorized version here
    // ...
	float ans=0.f;
	float tmp[VECTOR_WIDTH];
	__cmu418_vec_float x, y;
	__cmu418_mask dataMask, oneHot, ones1, ones2;

	for (int i=0;i<N;i+=VECTOR_WIDTH){
		int endId=N-i<VECTOR_WIDTH?N-i:VECTOR_WIDTH;
		dataMask=_cmu418_init_ones(endId);
		_cmu418_vload_float(x,values+i,dataMask);
		
		while (endId>1){
			if (endId%2==1){
				// add 0 to pos endId; endId+=1
				ones1=_cmu418_init_ones(endId+1);
				ones2=_cmu418_init_ones(endId);
				oneHot=_cmu418_mask_xor(ones1, ones2);
				_cmu418_vset_float(x,0.f,oneHot);
				endId+=1;
			}
			_cmu418_hadd_float(y,x);
			_cmu418_interleave_float(x,y);
			endId=endId/2;
		}
		ones1=_cmu418_init_ones(1);
		_cmu418_vstore_float(tmp,x,ones1);
		ans+=tmp[0];
	}

    return ans;
}
