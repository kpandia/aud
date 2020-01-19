/*
    This is a header file for BatchProcessWaveform.c

    Copyright (C) 2002-2016 Speech and Music Technology Lab,
    Indian Institute of Technology Madras
    
    Contributed by Hema A Murthy <hema@cse.iitm.ac.in>

    This file is part of KWS.

    KWS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    KWS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with KWS.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifndef BATCH_PROCESS_WAVEFORM_H
#define BATCH_PROCESS_WAVEFORM_H
F_VECTOR *FrameComputeWaveform(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLinearCepstrumMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLinearCepstrumRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLinearCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLinearDeltaCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLinearDeltaDeltaCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLinearAugmentedCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeResidualGDelay(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeHilbertEnvelopeResidualGDelay(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeGDelay(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResGdCepstrumMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResGdCepstrumRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResGdCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResGdDeltaCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResGdDeltaDeltaCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResGdAugmentedCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeResidualModGDelay(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeResidualModGDelaySmooth(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMinGDelay(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLPGDelay(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGDelayLP(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGDelayLPSmooth(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeSourceModGDelay(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGDelayPitch(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMinGDelayPitch(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeCepstrumPitch(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeSourceLPModGDelay(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGDelay(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGDelayMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMinGDelayMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGDelaySmooth(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResModGdCepstrumMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResModGdCepstrumRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResModGdCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResModGdDeltaCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResModGdDeltaDeltaCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelResModGdAugmentedCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLPCepstrumMeanFFT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);F_VECTOR *FrameComputeLPCepstrumRawFFT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLPCepstrumFFT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaLPCepstrumFFT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaLPCepstrumFFT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedLPCepstrumFFT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFFTCepstrumMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);F_VECTOR *FrameComputeFFTCepstrumRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFFTCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaFFTCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaFFTCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedFFTCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFFTSpectrumMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);F_VECTOR *FrameComputeFFTCepstrumRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFFTSpectrumRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFFTSpectrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaFFTSpectrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaFFTSpectrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedFFTSpectrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLPCepstrumMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);F_VECTOR *FrameComputeLPCepstrumRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLPCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaLPCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaLPCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedLPCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFilterbankEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFilterbankLogEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaFilterbankLogEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaFilterbankLogEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedFilterbankLogEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelCepstrumMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelCepstrumRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelDeltaCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelDeltaDeltaCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelAugmentedCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelRootCepstrumMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelRootCepstrumRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelRootCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelRootDeltaCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelRootDeltaDeltaCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelRootAugmentedCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdLogSmthCepstrumNcNRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdLogSmthCepstrumNcNMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdLogSmthCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaModGdLogSmthCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaModGdLogSmthCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedModGdLogSmthCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdCepstrumNcNRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdCepstrumNcNMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaModGdCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaModGdCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedModGdCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdCepstrumDCTRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdCepstrumDCTMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdCepstrumDCT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaModGdCepstrumDCT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaModGdCepstrumDCT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedModGdCepstrumDCT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdCepstrumLPDCTRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdCepstrumLPDCTMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdCepstrumLPDCT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaModGdCepstrumLPDCT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaModGdCepstrumLPDCT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedModGdCepstrumLPDCT(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeProductGdCepstrumNcNRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeProductGdCepstrumNcNMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeProductGdCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaProductGdCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaProductGdCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedProductGdCepstrumNcN(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMinGdCepstrumRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMinGdCepstrumMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMinGdCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaMinGdCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaMinGdCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeAugmentedMinGdCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFilterbankLogEnergyMean(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFilterbankLogEnergyRaw(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFilterbankLogEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelSlope(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelDeltaSlope(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelDeltaDeltaSlope(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMelAugmentedSlope(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeZeroCrossing(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeSpectralFlatness(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeSpectralFlux(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeSignificantChange(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLogEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaLogEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeDeltaDeltaLogEnergy(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFFTSpectrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLPSpectrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLPResAutoCorrSource(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeStdCepstrumSource(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeMinGdSource(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeLPModGdSource(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeModGdSource(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFlatSpectrumLogCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeFlatSpectrumRootCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeSmthSpectrumLogCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
F_VECTOR *FrameComputeSmthSpectrumRootCepstrum(ASDF *asdf, int frameIndex, F_VECTOR *fvect);
#endif







