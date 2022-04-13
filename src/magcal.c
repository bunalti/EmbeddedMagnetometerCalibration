// Copyright (c) 2014, Freescale Semiconductor, Inc.
// All rights reserved.
// vim: set ts=4:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Freescale Semiconductor, Inc. nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL FREESCALE SEMICONDUCTOR, INC. BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// This file contains magnetic calibration functions.  It is STRONGLY RECOMMENDED
// that the casual developer NOT TOUCH THIS FILE.  The mathematics behind this file
// is extremely complex, and it will be very easy (almost inevitable) that you screw
// it up.
//
// Haha - This file has been edited!  Please do not blame or pester NXP (formerly
//        Freescale) about the "almost inevitable" issues!

#include "magcal.h"

#define FXOS8700_UTPERCOUNT  0.1
#define DEFAULTB 50.0			// default geomagnetic field (uT)
#define X 0                         // vector components
#define Y 1
#define Z 2
#define ONETHIRD 0.33333333        // one third
#define ONESIXTH 0.166666667       // one sixth
#define MINMEASUREMENTS4CAL 40      // minimum number of measurements for 4 element calibration
#define MINMEASUREMENTS7CAL 100     // minimum number of measurements for 7 element calibration
#define MINMEASUREMENTS10CAL 150    // minimum number of measurements for 10 element calibration
#define MINBFITUT 22.0             // minimum geomagnetic field B (uT) for valid calibration
#define MAXBFITUT 67.0             // maximum geomagnetic field B (uT) for valid calibration
#define FITERRORAGINGSECS 7200.0   // 2 hours: time for fit error to increase (age) by e=2.718



// run the magnetic calibration
int MagCal_Run(void)
{
	int i, j;			// loop counters
	int isolver;		// magnetic solver used
	int count=0;


	// count number of data points
	for (i=0; i < MAGBUFFSIZE; i++) {
		if (magcal.valid[i]) count++;
	}


  ESP_LOGI("MAGCAL", "Count: %ld", count);
	if (count < MINMEASUREMENTS4CAL) return 0;

	if (magcal.ValidMagCal) {
		// age the existing fit error to avoid one good calibration locking out future updates
		magcal.FitErrorAge *= 1.02f;
	}

	// is enough data collected
	if (count < MINMEASUREMENTS7CAL) {
		isolver = 4;
    ESP_LOGI("MAGCAL", "Solver 4");
		fUpdateCalibration4INV(&magcal); // 4 element matrix inversion calibration
		if (magcal.trFitErrorpc < 12.0f) magcal.trFitErrorpc = 12.0f;
	} else if (count < MINMEASUREMENTS10CAL) {
		isolver = 7;
    ESP_LOGI("MAGCAL", "Solver 7");
		fUpdateCalibration7EIG(&magcal); // 7 element eigenpair calibration
		if (magcal.trFitErrorpc < 7.5f) magcal.trFitErrorpc = 7.5f;
	} else {
		isolver = 10;
    ESP_LOGI("MAGCAL", "Solver 10");
		fUpdateCalibration10EIG(&magcal); // 10 element eigenpair calibration
	}

	// the trial geomagnetic field must be in range (earth is 22uT to 67uT)
	if ((magcal.trB >= MINBFITUT) && (magcal.trB <= MAXBFITUT))	{
		// always accept the calibration if
		//  1: no previous calibration exists
		//  2: the calibration fit is reduced or
		//  3: an improved solver was used giving a good trial calibration (4% or under)
		if ((magcal.ValidMagCal == 0) ||
				(magcal.trFitErrorpc <= magcal.FitErrorAge) ||
				((isolver > magcal.ValidMagCal) && (magcal.trFitErrorpc <= 4.0F))) {
			// accept the new calibration solution
			//printf("new magnetic cal, B=%.2f uT\n", magcal.trB);
			magcal.ValidMagCal = isolver;
			magcal.FitError = magcal.trFitErrorpc;
			if (magcal.trFitErrorpc > 2.0f) {
				magcal.FitErrorAge = magcal.trFitErrorpc;
			} else {
				magcal.FitErrorAge = 2.0f;
			}
			magcal.B = magcal.trB;
			magcal.FourBsq = 4.0F * magcal.trB * magcal.trB;
			for (i = X; i <= Z; i++) {
				magcal.V[i] = magcal.trV[i];
				for (j = X; j <= Z; j++) {
					magcal.invW[i][j] = magcal.trinvW[i][j];
				}
			}
			return 1; // indicates new calibration applied
		}
	}
	return 0;
}



// 4 element calibration using 4x4 matrix inverse
void fUpdateCalibration4INV(MagCalibration_t *MagCal)
{
	float fBp2;					// fBp[X]^2+fBp[Y]^2+fBp[Z]^2
	float fSumBp4;				// sum of fBp2
	float fscaling;				// set to FUTPERCOUNT * FMATRIXSCALING
	float fE;					// error function = r^T.r
	int16_t iOffset[3];			// offset to remove large DC hard iron bias in matrix
	int16_t iCount;				// number of measurements counted
	int i, j, k;				// loop counters

	// working arrays for 4x4 matrix inversion
	float *pfRows[4];
	int8_t iColInd[4];
	int8_t iRowInd[4];
	int8_t iPivot[4];

	// compute fscaling to reduce multiplications later
	fscaling = FXOS8700_UTPERCOUNT / DEFAULTB;

	// the trial inverse soft iron matrix invW always equals
	// the identity matrix for 4 element calibration
	f3x3matrixAeqI(MagCal->trinvW);

	// zero fSumBp4=Y^T.Y, vecB=X^T.Y (4x1) and on and above
	// diagonal elements of matA=X^T*X (4x4)
	fSumBp4 = 0.0F;
	for (i = 0; i < 4; i++) {
		MagCal->vecB[i] = 0.0F;
		for (j = i; j < 4; j++) {
			MagCal->matA[i][j] = 0.0F;
		}
	}

	// the offsets are guaranteed to be set from the first element but to avoid compiler error
	iOffset[X] = iOffset[Y] = iOffset[Z] = 0;

	// use from MINEQUATIONS up to MAXEQUATIONS entries from magnetic buffer to compute matrices
	iCount = 0;
	for (j = 0; j < MAGBUFFSIZE; j++) {
		if (MagCal->valid[j]) {
			// use first valid magnetic buffer entry as estimate (in counts) for offset
			if (iCount == 0) {
				for (k = X; k <= Z; k++) {
					iOffset[k] = MagCal->BpFast[k][j];
				}
			}

			// store scaled and offset fBp[XYZ] in vecA[0-2] and fBp[XYZ]^2 in vecA[3-5]
			for (k = X; k <= Z; k++) {
				MagCal->vecA[k] = (float)((int32_t)MagCal->BpFast[k][j]
					- (int32_t)iOffset[k]) * fscaling;
				MagCal->vecA[k + 3] = MagCal->vecA[k] * MagCal->vecA[k];
			}

			// calculate fBp2 = Bp[X]^2 + Bp[Y]^2 + Bp[Z]^2 (scaled uT^2)
			fBp2 = MagCal->vecA[3] + MagCal->vecA[4] + MagCal->vecA[5];

			// accumulate fBp^4 over all measurements into fSumBp4=Y^T.Y
			fSumBp4 += fBp2 * fBp2;

			// now we have fBp2, accumulate vecB[0-2] = X^T.Y =sum(Bp2.Bp[XYZ])
			for (k = X; k <= Z; k++) {
				MagCal->vecB[k] += MagCal->vecA[k] * fBp2;
			}

			//accumulate vecB[3] = X^T.Y =sum(fBp2)
			MagCal->vecB[3] += fBp2;

			// accumulate on and above-diagonal terms of matA = X^T.X ignoring matA[3][3]
			MagCal->matA[0][0] += MagCal->vecA[X + 3];
			MagCal->matA[0][1] += MagCal->vecA[X] * MagCal->vecA[Y];
			MagCal->matA[0][2] += MagCal->vecA[X] * MagCal->vecA[Z];
			MagCal->matA[0][3] += MagCal->vecA[X];
			MagCal->matA[1][1] += MagCal->vecA[Y + 3];
			MagCal->matA[1][2] += MagCal->vecA[Y] * MagCal->vecA[Z];
			MagCal->matA[1][3] += MagCal->vecA[Y];
			MagCal->matA[2][2] += MagCal->vecA[Z + 3];
			MagCal->matA[2][3] += MagCal->vecA[Z];

			// increment the counter for next iteration
			iCount++;
		}
	}

	// set the last element of the measurement matrix to the number of buffer elements used
	MagCal->matA[3][3] = (float) iCount;

	// store the number of measurements accumulated
	MagCal->MagBufferCount = iCount;

	// use above diagonal elements of symmetric matA to set both matB and matA to X^T.X
	for (i = 0; i < 4; i++) {
		for (j = i; j < 4; j++) {
			MagCal->matB[i][j] = MagCal->matB[j][i]
				= MagCal->matA[j][i] = MagCal->matA[i][j];
		}
	}

	// calculate in situ inverse of matB = inv(X^T.X) (4x4) while matA still holds X^T.X
	for (i = 0; i < 4; i++) {
		pfRows[i] = MagCal->matB[i];
	}
	fmatrixAeqInvA(pfRows, iColInd, iRowInd, iPivot, 4);

	// calculate vecA = solution beta (4x1) = inv(X^T.X).X^T.Y = matB * vecB
	for (i = 0; i < 4; i++) {
		MagCal->vecA[i] = 0.0F;
		for (k = 0; k < 4; k++) {
			MagCal->vecA[i] += MagCal->matB[i][k] * MagCal->vecB[k];
		}
	}

	// calculate P = r^T.r = Y^T.Y - 2 * beta^T.(X^T.Y) + beta^T.(X^T.X).beta
	// = fSumBp4 - 2 * vecA^T.vecB + vecA^T.matA.vecA
	// first set P = Y^T.Y - 2 * beta^T.(X^T.Y) = SumBp4 - 2 * vecA^T.vecB
	fE = 0.0F;
	for (i = 0; i < 4; i++) {
		fE += MagCal->vecA[i] * MagCal->vecB[i];
	}
	fE = fSumBp4 - 2.0F * fE;

	// set vecB = (X^T.X).beta = matA.vecA
	for (i = 0; i < 4; i++) {
		MagCal->vecB[i] = 0.0F;
		for (k = 0; k < 4; k++) {
			MagCal->vecB[i] += MagCal->matA[i][k] * MagCal->vecA[k];
		}
	}

	// complete calculation of P by adding beta^T.(X^T.X).beta = vecA^T * vecB
	for (i = 0; i < 4; i++) {
		fE += MagCal->vecB[i] * MagCal->vecA[i];
	}

	// compute the hard iron vector (in uT but offset and scaled by FMATRIXSCALING)
	for (k = X; k <= Z; k++) {
		MagCal->trV[k] = 0.5F * MagCal->vecA[k];
	}

	// compute the scaled geomagnetic field strength B (in uT but scaled by FMATRIXSCALING)
	MagCal->trB = sqrtf(MagCal->vecA[3] + MagCal->trV[X] * MagCal->trV[X] +
			MagCal->trV[Y] * MagCal->trV[Y] + MagCal->trV[Z] * MagCal->trV[Z]);

	// calculate the trial fit error (percent) normalized to number of measurements
	// and scaled geomagnetic field strength
	MagCal->trFitErrorpc = sqrtf(fE / (float) MagCal->MagBufferCount) * 100.0F /
			(2.0F * MagCal->trB * MagCal->trB);

	// correct the hard iron estimate for FMATRIXSCALING and the offsets applied (result in uT)
	for (k = X; k <= Z; k++) {
		MagCal->trV[k] = MagCal->trV[k] * DEFAULTB
			+ (float)iOffset[k] * FXOS8700_UTPERCOUNT;
	}

	// correct the geomagnetic field strength B to correct scaling (result in uT)
	MagCal->trB *= DEFAULTB;
}










// 7 element calibration using direct eigen-decomposition
void fUpdateCalibration7EIG(MagCalibration_t *MagCal)
{
	float det;					// matrix determinant
	float fscaling;				// set to FUTPERCOUNT * FMATRIXSCALING
	float ftmp;					// scratch variable
	int16_t iOffset[3];			// offset to remove large DC hard iron bias
	int16_t iCount;				// number of measurements counted
	int i, j, k, m, n;			// loop counters


  
	// compute fscaling to reduce multiplications later
	fscaling = FXOS8700_UTPERCOUNT / DEFAULTB;

	// the offsets are guaranteed to be set from the first element but to avoid compiler error
	iOffset[X] = iOffset[Y] = iOffset[Z] = 0;

	// zero the on and above diagonal elements of the 7x7 symmetric measurement matrix matA
	for (m = 0; m < 7; m++) {
		for (n = m; n < 7; n++) {
			MagCal->matA[m][n] = 0.0F;
		}
	}

	// place from MINEQUATIONS to MAXEQUATIONS entries into product matrix matA
	iCount = 0;
	for (j = 0; j < MAGBUFFSIZE; j++) {
		if (MagCal->valid[j]) {
			// use first valid magnetic buffer entry as offset estimate (bit counts)
			if (iCount == 0) {
				for (k = X; k <= Z; k++) {
					iOffset[k] = MagCal->BpFast[k][j];
				}
			}

			// apply the offset and scaling and store in vecA
			for (k = X; k <= Z; k++) {
				MagCal->vecA[k + 3] = (float)((int32_t)MagCal->BpFast[k][j]
					- (int32_t)iOffset[k]) * fscaling;
				MagCal->vecA[k] = MagCal->vecA[k + 3] * MagCal->vecA[k + 3];
			}

			// accumulate the on-and above-diagonal terms of
			// MagCal->matA=Sigma{vecA^T * vecA}
			// with the exception of matA[6][6] which will sum to the number
			// of measurements and remembering that vecA[6] equals 1.0F
			// update the right hand column [6] of matA except for matA[6][6]
			for (m = 0; m < 6; m++) {
				MagCal->matA[m][6] += MagCal->vecA[m];
			}
			// update the on and above diagonal terms except for right hand column 6
			for (m = 0; m < 6; m++) {
				for (n = m; n < 6; n++) {
					MagCal->matA[m][n] += MagCal->vecA[m] * MagCal->vecA[n];
				}
			}

			// increment the measurement counter for the next iteration
			iCount++;
		}
	}

	// finally set the last element matA[6][6] to the number of measurements
	MagCal->matA[6][6] = (float) iCount;

	// store the number of measurements accumulated
	MagCal->MagBufferCount = iCount;

	// copy the above diagonal elements of matA to below the diagonal
	for (m = 1; m < 7; m++) {
		for (n = 0; n < m; n++) {
			MagCal->matA[m][n] = MagCal->matA[n][m];
		}
	}

	// set tmpA7x1 to the unsorted eigenvalues and matB to the unsorted eigenvectors of matA
	eigencompute(MagCal->matA, MagCal->vecA, MagCal->matB, 7);

	// find the smallest eigenvalue
	j = 0;
	for (i = 1; i < 7; i++) {
		if (MagCal->vecA[i] < MagCal->vecA[j]) {
			j = i;
		}
	}

	// set ellipsoid matrix A to the solution vector with smallest eigenvalue,
	// compute its determinant and the hard iron offset (scaled and offset)
	f3x3matrixAeqScalar(MagCal->A, 0.0F);
	det = 1.0F;
	for (k = X; k <= Z; k++) {
		MagCal->A[k][k] = MagCal->matB[k][j];
		det *= MagCal->A[k][k];
		MagCal->trV[k] = -0.5F * MagCal->matB[k + 3][j] / MagCal->A[k][k];
	}

	// negate A if it has negative determinant
	if (det < 0.0F) {
		f3x3matrixAeqMinusA(MagCal->A);
		MagCal->matB[6][j] = -MagCal->matB[6][j];
		det = -det;
	}

	// set ftmp to the square of the trial geomagnetic field strength B
	// (counts times FMATRIXSCALING)
	ftmp = -MagCal->matB[6][j];
	for (k = X; k <= Z; k++) {
		ftmp += MagCal->A[k][k] * MagCal->trV[k] * MagCal->trV[k];
	}

	// calculate the trial normalized fit error as a percentage
	MagCal->trFitErrorpc = 50.0F *
		sqrtf(fabs(MagCal->vecA[j]) / (float) MagCal->MagBufferCount) / fabs(ftmp);

	// normalize the ellipsoid matrix A to unit determinant
	f3x3matrixAeqAxScalar(MagCal->A, powf(det, -(ONETHIRD)));

	// convert the geomagnetic field strength B into uT for normalized
	// soft iron matrix A and normalize
	MagCal->trB = sqrtf(fabs(ftmp)) * DEFAULTB * powf(det, -(ONESIXTH));

	// compute trial invW from the square root of A also with normalized
	// determinant and hard iron offset in uT
	f3x3matrixAeqI(MagCal->trinvW);
	for (k = X; k <= Z; k++) {
		MagCal->trinvW[k][k] = sqrtf(fabs(MagCal->A[k][k]));
		MagCal->trV[k] = MagCal->trV[k] * DEFAULTB + (float)iOffset[k] * FXOS8700_UTPERCOUNT;
	}
}





// 10 element calibration using direct eigen-decomposition
void fUpdateCalibration10EIG(MagCalibration_t *MagCal)
{
	float det;					// matrix determinant
	float fscaling;				// set to FUTPERCOUNT * FMATRIXSCALING
	float ftmp;					// scratch variable
	int16_t iOffset[3];			// offset to remove large DC hard iron bias in matrix	
	int16_t iCount;				// number of measurements counted
	int i, j, k, m, n;			// loop counters

	// compute fscaling to reduce multiplications later
	fscaling = FXOS8700_UTPERCOUNT / DEFAULTB;

	// the offsets are guaranteed to be set from the first element but to avoid compiler error
	iOffset[X] = iOffset[Y] = iOffset[Z] = 0;

	// zero the on and above diagonal elements of the 10x10 symmetric measurement matrix matA
	for (m = 0; m < 10; m++) {
		for (n = m; n < 10; n++) {
			MagCal->matA[m][n] = 0.0F;
		}
	}

	// sum between MINEQUATIONS to MAXEQUATIONS entries into the 10x10 product matrix matA
	iCount = 0;
	for (j = 0; j < MAGBUFFSIZE; j++) {
		if (MagCal->valid[j]) {
			// use first valid magnetic buffer entry as estimate for offset
			// to help solution (bit counts)
			if (iCount == 0) {
				for (k = X; k <= Z; k++) {
					iOffset[k] = MagCal->BpFast[k][j];
				}
			}

			// apply the fixed offset and scaling and enter into vecA[6-8]
			for (k = X; k <= Z; k++) {
				MagCal->vecA[k + 6] = (float)((int32_t)MagCal->BpFast[k][j]
					- (int32_t)iOffset[k]) * fscaling;
			}

			// compute measurement vector elements vecA[0-5] from vecA[6-8]
			MagCal->vecA[0] = MagCal->vecA[6] * MagCal->vecA[6];
			MagCal->vecA[1] = 2.0F * MagCal->vecA[6] * MagCal->vecA[7];
			MagCal->vecA[2] = 2.0F * MagCal->vecA[6] * MagCal->vecA[8];
			MagCal->vecA[3] = MagCal->vecA[7] * MagCal->vecA[7];
			MagCal->vecA[4] = 2.0F * MagCal->vecA[7] * MagCal->vecA[8];
			MagCal->vecA[5] = MagCal->vecA[8] * MagCal->vecA[8];

			// accumulate the on-and above-diagonal terms of matA=Sigma{vecA^T * vecA}
			// with the exception of matA[9][9] which equals the number of measurements
			// update the right hand column [9] of matA[0-8][9] ignoring matA[9][9]
			for (m = 0; m < 9; m++) {
				MagCal->matA[m][9] += MagCal->vecA[m];
			}
			// update the on and above diagonal terms of matA ignoring right hand column 9
			for (m = 0; m < 9; m++) {
				for (n = m; n < 9; n++) {
					MagCal->matA[m][n] += MagCal->vecA[m] * MagCal->vecA[n];
				}
			}

			// increment the measurement counter for the next iteration
			iCount++;
		}
	}

	// set the last element matA[9][9] to the number of measurements
	MagCal->matA[9][9] = (float) iCount;

	// store the number of measurements accumulated
	MagCal->MagBufferCount = iCount;

	// copy the above diagonal elements of symmetric product matrix matA to below the diagonal
	for (m = 1; m < 10; m++) {
		for (n = 0; n < m; n++) {
			MagCal->matA[m][n] = MagCal->matA[n][m];
		}
	}

	// set MagCal->vecA to the unsorted eigenvalues and matB to the unsorted
	// normalized eigenvectors of matA
	eigencompute(MagCal->matA, MagCal->vecA, MagCal->matB, 10);

	// set ellipsoid matrix A from elements of the solution vector column j with
	// smallest eigenvalue
	j = 0;
	for (i = 1; i < 10; i++) {
		if (MagCal->vecA[i] < MagCal->vecA[j]) {
			j = i;
		}
	}
	MagCal->A[0][0] = MagCal->matB[0][j];
	MagCal->A[0][1] = MagCal->A[1][0] = MagCal->matB[1][j];
	MagCal->A[0][2] = MagCal->A[2][0] = MagCal->matB[2][j];
	MagCal->A[1][1] = MagCal->matB[3][j];
	MagCal->A[1][2] = MagCal->A[2][1] = MagCal->matB[4][j];
	MagCal->A[2][2] = MagCal->matB[5][j];

	// negate entire solution if A has negative determinant
	det = f3x3matrixDetA(MagCal->A);
	if (det < 0.0F) {
		f3x3matrixAeqMinusA(MagCal->A);
		MagCal->matB[6][j] = -MagCal->matB[6][j];
		MagCal->matB[7][j] = -MagCal->matB[7][j];
		MagCal->matB[8][j] = -MagCal->matB[8][j];
		MagCal->matB[9][j] = -MagCal->matB[9][j];
		det = -det;
	}

	// compute the inverse of the ellipsoid matrix
	f3x3matrixAeqInvSymB(MagCal->invA, MagCal->A);

	// compute the trial hard iron vector in offset bit counts times FMATRIXSCALING
	for (k = X; k <= Z; k++) {
		MagCal->trV[k] = 0.0F;
		for (m = X; m <= Z; m++) {
			MagCal->trV[k] += MagCal->invA[k][m] * MagCal->matB[m + 6][j];
		}
		MagCal->trV[k] *= -0.5F;
	}

	// compute the trial geomagnetic field strength B in bit counts times FMATRIXSCALING
	MagCal->trB = sqrtf(fabs(MagCal->A[0][0] * MagCal->trV[X] * MagCal->trV[X] +
			2.0F * MagCal->A[0][1] * MagCal->trV[X] * MagCal->trV[Y] +
			2.0F * MagCal->A[0][2] * MagCal->trV[X] * MagCal->trV[Z] +
			MagCal->A[1][1] * MagCal->trV[Y] * MagCal->trV[Y] +
			2.0F * MagCal->A[1][2] * MagCal->trV[Y] * MagCal->trV[Z] +
			MagCal->A[2][2] * MagCal->trV[Z] * MagCal->trV[Z] - MagCal->matB[9][j]));

	// calculate the trial normalized fit error as a percentage
	MagCal->trFitErrorpc = 50.0F * sqrtf(
		fabs(MagCal->vecA[j]) / (float) MagCal->MagBufferCount) /
		(MagCal->trB * MagCal->trB);

	// correct for the measurement matrix offset and scaling and
	// get the computed hard iron offset in uT
	for (k = X; k <= Z; k++) {
		MagCal->trV[k] = MagCal->trV[k] * DEFAULTB + (float)iOffset[k] * FXOS8700_UTPERCOUNT;
	}

	// convert the trial geomagnetic field strength B into uT for
	// un-normalized soft iron matrix A
	MagCal->trB *= DEFAULTB;

	// normalize the ellipsoid matrix A to unit determinant and
	// correct B by root of this multiplicative factor
	f3x3matrixAeqAxScalar(MagCal->A, powf(det, -(ONETHIRD)));
	MagCal->trB *= powf(det, -(ONESIXTH));

	// compute trial invW from the square root of fA (both with normalized determinant)
	// set vecA to the unsorted eigenvalues and matB to the unsorted eigenvectors of matA
	// where matA holds the 3x3 matrix fA in its top left elements
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			MagCal->matA[i][j] = MagCal->A[i][j];
		}
	}
	eigencompute(MagCal->matA, MagCal->vecA, MagCal->matB, 3);

	// set MagCal->matB to be eigenvectors . diag(sqrt(sqrt(eigenvalues))) =
	//   matB . diag(sqrt(sqrt(vecA))
	for (j = 0; j < 3; j++) { // loop over columns j
		ftmp = sqrtf(sqrtf(fabs(MagCal->vecA[j])));
		for (i = 0; i < 3; i++) { // loop over rows i
			MagCal->matB[i][j] *= ftmp;
		}
	}

	// set trinvW to eigenvectors * diag(sqrt(eigenvalues)) * eigenvectors^T =
	//   matB * matB^T = sqrt(fA) (guaranteed symmetric)
	// loop over rows
	for (i = 0; i < 3; i++) {
		// loop over on and above diagonal columns
		for (j = i; j < 3; j++) {
			MagCal->trinvW[i][j] = 0.0F;
			// accumulate the matrix product
			for (k = 0; k < 3; k++) {
				MagCal->trinvW[i][j] += MagCal->matB[i][k] * MagCal->matB[j][k];
			}
			// copy to below diagonal element
			MagCal->trinvW[j][i] = MagCal->trinvW[i][j];
		}
	}
}

void raw_data_reset(void)
{

  memset(&magcal, 0, sizeof(magcal));
  magcal.V[2] = 80.0f;  // initial guess
  magcal.invW[0][0] = 1.0f;
  magcal.invW[1][1] = 1.0f;
  magcal.invW[2][2] = 1.0f;
  magcal.FitError = 100.0f;
  magcal.FitErrorAge = 100.0f;
  magcal.B = 50.0f;
}

void apply_calibration(int16_t rawx, int16_t rawy, int16_t rawz, Point_t *out)
{
  float x, y, z;

  x = ((float)rawx * UT_PER_COUNT) - magcal.V[0];
  y = ((float)rawy * UT_PER_COUNT) - magcal.V[1];
  z = ((float)rawz * UT_PER_COUNT) - magcal.V[2];
  out->x = x * magcal.invW[0][0] + y * magcal.invW[0][1] + z * magcal.invW[0][2];
  out->y = x * magcal.invW[1][0] + y * magcal.invW[1][1] + z * magcal.invW[1][2];
  out->z = x * magcal.invW[2][0] + y * magcal.invW[2][1] + z * magcal.invW[2][2];
}

// function sets the 3x3 matrix A to the identity matrix
void f3x3matrixAeqI(float A[][3])
{
  float *pAij;  // pointer to A[i][j]
  int8_t i, j;  // loop counters

  for (i = 0; i < 3; i++) {
    // set pAij to &A[i][j=0]
    pAij = A[i];
    for (j = 0; j < 3; j++) {
      *(pAij++) = 0.0F;
    }
    A[i][i] = 1.0F;
  }
}

// function sets the matrix A to the identity matrix
void fmatrixAeqI(float *A[], int16_t rc)
{
  // rc = rows and columns in A

  float *pAij;  // pointer to A[i][j]
  int8_t i, j;  // loop counters

  for (i = 0; i < rc; i++) {
    // set pAij to &A[i][j=0]
    pAij = A[i];
    for (j = 0; j < rc; j++) {
      *(pAij++) = 0.0F;
    }
    A[i][i] = 1.0F;
  }
}

// function sets every entry in the 3x3 matrix A to a constant scalar
void f3x3matrixAeqScalar(float A[][3], float Scalar)
{
  float *pAij;  // pointer to A[i][j]
  int8_t i, j;  // counters

  for (i = 0; i < 3; i++) {
    // set pAij to &A[i][j=0]
    pAij = A[i];
    for (j = 0; j < 3; j++) {
      *(pAij++) = Scalar;
    }
  }
}

// function multiplies all elements of 3x3 matrix A by the specified scalar
void f3x3matrixAeqAxScalar(float A[][3], float Scalar)
{
  float *pAij;  // pointer to A[i][j]
  int8_t i, j;  // loop counters

  for (i = 0; i < 3; i++) {
    // set pAij to &A[i][j=0]
    pAij = A[i];
    for (j = 0; j < 3; j++) {
      *(pAij++) *= Scalar;
    }
  }
}

// function negates all elements of 3x3 matrix A
void f3x3matrixAeqMinusA(float A[][3])
{
  float *pAij;  // pointer to A[i][j]
  int8_t i, j;  // loop counters

  for (i = 0; i < 3; i++) {
    // set pAij to &A[i][j=0]
    pAij = A[i];
    for (j = 0; j < 3; j++) {
      *pAij = -*pAij;
      pAij++;
    }
  }
}

// function directly calculates the symmetric inverse of a symmetric 3x3 matrix
// only the on and above diagonal terms in B are used and need to be specified
void f3x3matrixAeqInvSymB(float A[][3], float B[][3])
{
  float fB11B22mB12B12; // B[1][1] * B[2][2] - B[1][2] * B[1][2]
  float fB12B02mB01B22; // B[1][2] * B[0][2] - B[0][1] * B[2][2]
  float fB01B12mB11B02; // B[0][1] * B[1][2] - B[1][1] * B[0][2]
  float ftmp;       // determinant and then reciprocal

  // calculate useful products
  fB11B22mB12B12 = B[1][1] * B[2][2] - B[1][2] * B[1][2];
  fB12B02mB01B22 = B[1][2] * B[0][2] - B[0][1] * B[2][2];
  fB01B12mB11B02 = B[0][1] * B[1][2] - B[1][1] * B[0][2];

  // set ftmp to the determinant of the input matrix B
  ftmp = B[0][0] * fB11B22mB12B12 + B[0][1] * fB12B02mB01B22 + B[0][2] * fB01B12mB11B02;

  // set A to the inverse of B for any determinant except zero
  if (ftmp != 0.0F) {
    ftmp = 1.0F / ftmp;
    A[0][0] = fB11B22mB12B12 * ftmp;
    A[1][0] = A[0][1] = fB12B02mB01B22 * ftmp;
    A[2][0] = A[0][2] = fB01B12mB11B02 * ftmp;
    A[1][1] = (B[0][0] * B[2][2] - B[0][2] * B[0][2]) * ftmp;
    A[2][1] = A[1][2] = (B[0][2] * B[0][1] - B[0][0] * B[1][2]) * ftmp;
    A[2][2] = (B[0][0] * B[1][1] - B[0][1] * B[0][1]) * ftmp;
  } else {
    // provide the identity matrix if the determinant is zero
    f3x3matrixAeqI(A);
  }
}

// function calculates the determinant of a 3x3 matrix
float f3x3matrixDetA(float A[][3])
{
  return (A[X][X] * (A[Y][Y] * A[Z][Z] - A[Y][Z] * A[Z][Y]) +
      A[X][Y] * (A[Y][Z] * A[Z][X] - A[Y][X] * A[Z][Z]) +
      A[X][Z] * (A[Y][X] * A[Z][Y] - A[Y][Y] * A[Z][X]));
}

// function computes all eigenvalues and eigenvectors of a real symmetric matrix A[0..n-1][0..n-1]
// stored in the top left of a 10x10 array A[10][10]
// A[][] is changed on output.
// eigval[0..n-1] returns the eigenvalues of A[][].
// eigvec[0..n-1][0..n-1] returns the normalized eigenvectors of A[][]
// the eigenvectors are not sorted by value
void eigencompute(float A[][10], float eigval[], float eigvec[][10], int8_t n)
{
  // maximum number of iterations to achieve convergence: in practice 6 is typical
#define NITERATIONS 15

  // various trig functions of the jacobi rotation angle phi
  float cot2phi, tanhalfphi, tanphi, sinphi, cosphi;
  // scratch variable to prevent over-writing during rotations
  float ftmp;
  // residue from remaining non-zero above diagonal terms
  float residue;
  // matrix row and column indices
  int8_t ir, ic;
  // general loop counter
  int8_t j;
  // timeout ctr for number of passes of the algorithm
  int8_t ctr;

  // initialize eigenvectors matrix and eigenvalues array
  for (ir = 0; ir < n; ir++) {
    // loop over all columns
    for (ic = 0; ic < n; ic++) {
      // set on diagonal and off-diagonal elements to zero
      eigvec[ir][ic] = 0.0F;
    }

    // correct the diagonal elements to 1.0
    eigvec[ir][ir] = 1.0F;

    // initialize the array of eigenvalues to the diagonal elements of m
    eigval[ir] = A[ir][ir];
  }

  // initialize the counter and loop until converged or NITERATIONS reached
  ctr = 0;
  do {
    // compute the absolute value of the above diagonal elements as exit criterion
    residue = 0.0F;
    // loop over rows excluding last row
    for (ir = 0; ir < n - 1; ir++) {
      // loop over above diagonal columns
      for (ic = ir + 1; ic < n; ic++) {
        // accumulate the residual off diagonal terms which are being driven to zero
        residue += fabs(A[ir][ic]);
      }
    }

    // check if we still have work to do
    if (residue > 0.0F) {
      // loop over all rows with the exception of the last row (since only rotating above diagonal elements)
      for (ir = 0; ir < n - 1; ir++) {
        // loop over columns ic (where ic is always greater than ir since above diagonal)
        for (ic = ir + 1; ic < n; ic++) {
          // only continue with this element if the element is non-zero
          if (fabs(A[ir][ic]) > 0.0F) {
            // calculate cot(2*phi) where phi is the Jacobi rotation angle
            cot2phi = 0.5F * (eigval[ic] - eigval[ir]) / (A[ir][ic]);

            // calculate tan(phi) correcting sign to ensure the smaller solution is used
            tanphi = 1.0F / (fabs(cot2phi) + sqrtf(1.0F + cot2phi * cot2phi));
            if (cot2phi < 0.0F) {
              tanphi = -tanphi;
            }

            // calculate the sine and cosine of the Jacobi rotation angle phi
            cosphi = 1.0F / sqrtf(1.0F + tanphi * tanphi);
            sinphi = tanphi * cosphi;

            // calculate tan(phi/2)
            tanhalfphi = sinphi / (1.0F + cosphi);

            // set tmp = tan(phi) times current matrix element used in update of leading diagonal elements
            ftmp = tanphi * A[ir][ic];

            // apply the jacobi rotation to diagonal elements [ir][ir] and [ic][ic] stored in the eigenvalue array
            // eigval[ir] = eigval[ir] - tan(phi) *  A[ir][ic]
            eigval[ir] -= ftmp;
            // eigval[ic] = eigval[ic] + tan(phi) * A[ir][ic]
            eigval[ic] += ftmp;

            // by definition, applying the jacobi rotation on element ir, ic results in 0.0
            A[ir][ic] = 0.0F;

            // apply the jacobi rotation to all elements of the eigenvector matrix
            for (j = 0; j < n; j++) {
              // store eigvec[j][ir]
              ftmp = eigvec[j][ir];
              // eigvec[j][ir] = eigvec[j][ir] - sin(phi) * (eigvec[j][ic] + tan(phi/2) * eigvec[j][ir])
              eigvec[j][ir] = ftmp - sinphi * (eigvec[j][ic] + tanhalfphi * ftmp);
              // eigvec[j][ic] = eigvec[j][ic] + sin(phi) * (eigvec[j][ir] - tan(phi/2) * eigvec[j][ic])
              eigvec[j][ic] = eigvec[j][ic] + sinphi * (ftmp - tanhalfphi * eigvec[j][ic]);
            }

            // apply the jacobi rotation only to those elements of matrix m that can change
            for (j = 0; j <= ir - 1; j++) {
              // store A[j][ir]
              ftmp = A[j][ir];
              // A[j][ir] = A[j][ir] - sin(phi) * (A[j][ic] + tan(phi/2) * A[j][ir])
              A[j][ir] = ftmp - sinphi * (A[j][ic] + tanhalfphi * ftmp);
              // A[j][ic] = A[j][ic] + sin(phi) * (A[j][ir] - tan(phi/2) * A[j][ic])
              A[j][ic] = A[j][ic] + sinphi * (ftmp - tanhalfphi * A[j][ic]);
            }
            for (j = ir + 1; j <= ic - 1; j++) {
              // store A[ir][j]
              ftmp = A[ir][j];
              // A[ir][j] = A[ir][j] - sin(phi) * (A[j][ic] + tan(phi/2) * A[ir][j])
              A[ir][j] = ftmp - sinphi * (A[j][ic] + tanhalfphi * ftmp);
              // A[j][ic] = A[j][ic] + sin(phi) * (A[ir][j] - tan(phi/2) * A[j][ic])
              A[j][ic] = A[j][ic] + sinphi * (ftmp - tanhalfphi * A[j][ic]);
            }
            for (j = ic + 1; j < n; j++) {
              // store A[ir][j]
              ftmp = A[ir][j];
              // A[ir][j] = A[ir][j] - sin(phi) * (A[ic][j] + tan(phi/2) * A[ir][j])
              A[ir][j] = ftmp - sinphi * (A[ic][j] + tanhalfphi * ftmp);
              // A[ic][j] = A[ic][j] + sin(phi) * (A[ir][j] - tan(phi/2) * A[ic][j])
              A[ic][j] = A[ic][j] + sinphi * (ftmp - tanhalfphi * A[ic][j]);
            }
          }   // end of test for matrix element already zero
        }   // end of loop over columns
      }   // end of loop over rows
    }  // end of test for non-zero residue
  } while ((residue > 0.0F) && (ctr++ < NITERATIONS)); // end of main loop
}

// function uses Gauss-Jordan elimination to compute the inverse of matrix A in situ
// on exit, A is replaced with its inverse
void fmatrixAeqInvA(float *A[], int8_t iColInd[], int8_t iRowInd[], int8_t iPivot[], int8_t isize)
{
  float largest;          // largest element used for pivoting
  float scaling;          // scaling factor in pivoting
  float recippiv;         // reciprocal of pivot element
  float ftmp;           // temporary variable used in swaps
  int8_t i, j, k, l, m;     // index counters
  int8_t iPivotRow, iPivotCol;  // row and column of pivot element

  // to avoid compiler warnings
  iPivotRow = iPivotCol = 0;

  // initialize the pivot array to 0
  for (j = 0; j < isize; j++) {
    iPivot[j] = 0;
  }

  // main loop i over the dimensions of the square matrix A
  for (i = 0; i < isize; i++) {
    // zero the largest element found for pivoting
    largest = 0.0F;
    // loop over candidate rows j
    for (j = 0; j < isize; j++) {
      // check if row j has been previously pivoted
      if (iPivot[j] != 1) {
        // loop over candidate columns k
        for (k = 0; k < isize; k++) {
          // check if column k has previously been pivoted
          if (iPivot[k] == 0) {
            // check if the pivot element is the largest found so far
            if (fabs(A[j][k]) >= largest) {
              // and store this location as the current best candidate for pivoting
              iPivotRow = j;
              iPivotCol = k;
              largest = (float) fabs(A[iPivotRow][iPivotCol]);
            }
          } else if (iPivot[k] > 1) {
            // zero determinant situation: exit with identity matrix
            fmatrixAeqI(A, isize);
            return;
          }
        }
      }
    }
    // increment the entry in iPivot to denote it has been selected for pivoting
    iPivot[iPivotCol]++;

    // check the pivot rows iPivotRow and iPivotCol are not the same before swapping
    if (iPivotRow != iPivotCol) {
      // loop over columns l
      for (l = 0; l < isize; l++) {
        // and swap all elements of rows iPivotRow and iPivotCol
        ftmp = A[iPivotRow][l];
        A[iPivotRow][l] = A[iPivotCol][l];
        A[iPivotCol][l] = ftmp;
      }
    }

    // record that on the i-th iteration rows iPivotRow and iPivotCol were swapped
    iRowInd[i] = iPivotRow;
    iColInd[i] = iPivotCol;

    // check for zero on-diagonal element (singular matrix) and return with identity matrix if detected
    if (A[iPivotCol][iPivotCol] == 0.0F) {
      // zero determinant situation: exit with identity matrix
      fmatrixAeqI(A, isize);
      return;
    }

    // calculate the reciprocal of the pivot element knowing it's non-zero
    recippiv = 1.0F / A[iPivotCol][iPivotCol];
    // by definition, the diagonal element normalizes to 1
    A[iPivotCol][iPivotCol] = 1.0F;
    // multiply all of row iPivotCol by the reciprocal of the pivot element including the diagonal element
    // the diagonal element A[iPivotCol][iPivotCol] now has value equal to the reciprocal of its previous value
    for (l = 0; l < isize; l++) {
      A[iPivotCol][l] *= recippiv;
    }
    // loop over all rows m of A
    for (m = 0; m < isize; m++) {
      if (m != iPivotCol) {
        // scaling factor for this row m is in column iPivotCol
        scaling = A[m][iPivotCol];
        // zero this element
        A[m][iPivotCol] = 0.0F;
        // loop over all columns l of A and perform elimination
        for (l = 0; l < isize; l++) {
          A[m][l] -= A[iPivotCol][l] * scaling;
        }
      }
    }
  } // end of loop i over the matrix dimensions

  // finally, loop in inverse order to apply the missing column swaps
  for (l = isize - 1; l >= 0; l--) {
    // set i and j to the two columns to be swapped
    i = iRowInd[l];
    j = iColInd[l];

    // check that the two columns i and j to be swapped are not the same
    if (i != j) {
      // loop over all rows k to swap columns i and j of A
      for (k = 0; k < isize; k++) {
        ftmp = A[k][i];
        A[k][i] = A[k][j];
        A[k][j] = ftmp;
      }
    }
  }
}

// function re-orthonormalizes a 3x3 rotation matrix
void fmatrixAeqRenormRotA(float A[][3])
{
  float ftmp;         // scratch variable

  // normalize the X column of the low pass filtered orientation matrix
  ftmp = sqrtf(A[X][X] * A[X][X] + A[Y][X] * A[Y][X] + A[Z][X] * A[Z][X]);
  if (ftmp > CORRUPTMATRIX) {
    // normalize the x column vector
    ftmp = 1.0F / ftmp;
    A[X][X] *= ftmp;
    A[Y][X] *= ftmp;
    A[Z][X] *= ftmp;
  } else {
    // set x column vector to {1, 0, 0}
    A[X][X] = 1.0F;
    A[Y][X] = A[Z][X] = 0.0F;
  }

  // force the y column vector to be orthogonal to x using y = y-(x.y)x
  ftmp = A[X][X] * A[X][Y] + A[Y][X] * A[Y][Y] + A[Z][X] * A[Z][Y];
  A[X][Y] -= ftmp * A[X][X];
  A[Y][Y] -= ftmp * A[Y][X];
  A[Z][Y] -= ftmp * A[Z][X];

  // normalize the y column vector
  ftmp = sqrtf(A[X][Y] * A[X][Y] + A[Y][Y] * A[Y][Y] + A[Z][Y] * A[Z][Y]);
  if (ftmp > CORRUPTMATRIX) {
    // normalize the y column vector
    ftmp = 1.0F / ftmp;
    A[X][Y] *= ftmp;
    A[Y][Y] *= ftmp;
    A[Z][Y] *= ftmp;
  } else {
    // set y column vector to {0, 1, 0}
    A[Y][Y] = 1.0F;
    A[X][Y] = A[Z][Y] = 0.0F;
  }

  // finally set the z column vector to x vector cross y vector (automatically normalized)
  A[X][Z] = A[Y][X] * A[Z][Y] - A[Z][X] * A[Y][Y];
  A[Y][Z] = A[Z][X] * A[X][Y] - A[X][X] * A[Z][Y];
  A[Z][Z] = A[X][X] * A[Y][Y] - A[Y][X] * A[X][Y];
}
