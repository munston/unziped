{-# LANGUAGE ScopedTypeVariables #-}

module Fourier
    ( -- Types
      FourierCoeffs
    , FourierState(..)
    , FourierConfig(..)
      -- Core functions
    , initializeFourierState
    , initializeFourierStateWithLength
    , dftVandermonde
    , evaluateFourier
    , evaluateFourierSequence
    , evaluateFourierState
      -- Utility
    , nextPowerOf2
    ) where

import Data.Complex
import Numeric.FFT (fft, ifft)

-- Types
type FourierCoeffs = [Complex Double]

-- Store both coefficients and the FFT length used
data FourierState = FourierState
    { fsCoeffs :: FourierCoeffs    -- FFT coefficients
    , fsOrigLen :: Int              -- Original data length
    , fsPaddedLen :: Int            -- Padded length (power of 2)
    } deriving (Show, Eq)

data FourierConfig = FourierConfig
    { fourierOrder :: Int      -- Number of frequencies to keep
    , fourierDeltaT :: Double  -- Time step
    } deriving (Show, Eq)

-- ============================================================================
-- FFT Initialization
-- ============================================================================

-- Find next power of 2 >= n
nextPowerOf2 :: Int -> Int
nextPowerOf2 n = head [2^k | k <- [0..], 2^k >= n]

-- Initialize Fourier state from real data using FFT
-- Returns FFT coefficients with proper normalization tracking
initializeFourierState :: [Double] -> FourierState
initializeFourierState realData =
    let origLen = length realData
        paddedLen = nextPowerOf2 origLen
        -- Pad with zeros to power of 2
        paddedData = realData ++ replicate (paddedLen - origLen) 0
        complexData = map (:+ 0) paddedData
        -- Compute FFT (returns paddedLen coefficients)
        fullFFT = fft complexData
    in FourierState
        { fsCoeffs = fullFFT
        , fsOrigLen = origLen
        , fsPaddedLen = paddedLen
        }

-- Initialize with explicit target padded length (for matching reduced/full periods)
initializeFourierStateWithLength :: Int -> [Double] -> FourierState
initializeFourierStateWithLength targetPaddedLen realData =
    let origLen = length realData
        -- Pad to target length (must be power of 2)
        paddedData = realData ++ replicate (targetPaddedLen - origLen) 0
        complexData = map (:+ 0) paddedData
        -- Compute FFT
        fullFFT = fft complexData
    in FourierState
        { fsCoeffs = fullFFT
        , fsOrigLen = origLen
        , fsPaddedLen = targetPaddedLen
        }

-- ============================================================================
-- DFT Vandermonde Matrix (for BLUE correction)
-- ============================================================================

-- Construct DFT Vandermonde matrix
-- F[n,k] = exp(+2*pi*i*k*n/paddedLen)  -- INVERSE DFT (reconstruction)
-- indices: which time points to evaluate at
-- paddedLen: padded FFT length (determines frequency spacing)
-- order: how many frequencies to include
-- Must match the frequency spacing in evaluateFourier!
-- Note: Does NOT include 1/N normalization (that's in evaluateFourier)
dftVandermonde :: [Int] -> Int -> Int -> [[Complex Double]]
dftVandermonde indices paddedLen order =
    [[ cis (2 * pi * fromIntegral k * fromIntegral idx / fromIntegral paddedLen)
     | k <- [0..order]]
    | idx <- indices]

-- ============================================================================
-- Fourier Evaluation (analogous to Taylor unroll)
-- ============================================================================

-- Evaluate Fourier series at a single time index
-- c: Fourier coefficients
-- n: time index  
-- origLen: not used (kept for compatibility)
-- paddedLen: padded length (for frequency spacing AND normalization)
-- Returns: (1/paddedLen) * sum of c[k] * exp(2*pi*i*k*n/paddedLen)
evaluateFourier :: FourierCoeffs -> Int -> Int -> Int -> Complex Double
evaluateFourier coeffs n _ paddedLen =
    let sumVal = sum $ zipWith (*) coeffs exponentials
        exponentials = [cis (2 * pi * fromIntegral k * fromIntegral n / fromIntegral paddedLen)
                       | k <- [0..length coeffs - 1]]
    in sumVal / (fromIntegral paddedLen :+ 0)

-- Evaluate Fourier series over a sequence of indices using FourierState
-- Returns real parts only (for real-valued signals)
-- Handles extrapolation beyond original data range
evaluateFourierState :: FourierState -> [Int] -> [Double]
evaluateFourierState state indices =
    let coeffs = fsCoeffs state
        origLen = fsOrigLen state
        paddedLen = fsPaddedLen state
        -- Scale indices to match the period of this FFT
        -- If evaluating beyond paddedLen, use periodicity
    in map realPart [evaluateFourier coeffs idx origLen paddedLen | idx <- indices]

-- Evaluate Fourier series over a sequence of indices
-- Returns real parts only (for real-valued signals)
evaluateFourierSequence :: FourierCoeffs -> [Int] -> Int -> Int -> [Double]
evaluateFourierSequence coeffs indices origLen paddedLen =
    map realPart [evaluateFourier coeffs idx origLen paddedLen | idx <- indices]
