{-# LANGUAGE ScopedTypeVariables #-}

module FourierCorrection
  ( FourierCVConfig(..)
  , FourierCorrectionMatrix
  , computeFourierCorrectionMatrix
  , applyFourierCorrection
  , applyFourierCorrectionIterative  -- NEW
  , testFourierCorrection
  , buildH
  ) where

import ComplexMatrix
import Fourier
import Data.Complex
import Numeric.FFT (fft, ifft)

-- =========================================================
-- Types and configuration
-- =========================================================

data FourierCVConfig = FourierCVConfig
  { fcvLength :: Int
  , fdeltaT  :: Double
  }

-- H pseudoinverse (128 x 256 for standard case)
type FourierCorrectionMatrix = [[Complex Double]]

-- =========================================================
-- Build H matrix: zero-padding transform (origLen -> targetLen)
-- =========================================================

buildH :: Int -> Int -> [[Complex Double]]
buildH origLen targetLen = 
    [ buildRow i | i <- [0..targetLen-1] ]
  where
    halfOrig = origLen `div` 2
    
    buildRow i
        | i < halfOrig =
            -- Positive frequencies: multiply by 2
            [ if j == i then (2.0 :+ 0) else (0 :+ 0) 
            | j <- [0..origLen-1] ]
        | i < targetLen - halfOrig =
            -- Middle zeros (new high frequencies)
            replicate origLen (0 :+ 0)
        | otherwise =
            -- Negative frequencies: multiply by 2
            let sourceIdx = i - (targetLen - origLen)
            in [ if j == sourceIdx then (2.0 :+ 0) else (0 :+ 0)
               | j <- [0..origLen-1] ]

-- =========================================================
-- Build H pseudoinverse (direct construction, no inversion)
-- Since H^H H = 4I, we have H^+ = (1/4) H^H = 0.5 H^H
-- =========================================================

buildHPseudoinverse :: Int -> Int -> [[Complex Double]]
buildHPseudoinverse origLen targetLen =
    [ buildRow i | i <- [0..origLen-1] ]
  where
    halfOrig = origLen `div` 2
    
    buildRow i
        | i < halfOrig =
            -- Extract from first halfOrig columns of C_targetLen
            [ if j == i then (0.5 :+ 0) else (0 :+ 0)
            | j <- [0..targetLen-1] ]
        | otherwise =
            -- Extract from last halfOrig columns of C_targetLen
            let targetIdx = i + (targetLen - origLen)
            in [ if j == targetIdx then (0.5 :+ 0) else (0 :+ 0)
               | j <- [0..targetLen-1] ]

-- =========================================================
-- Compute correction matrix (H pseudoinverse)
-- =========================================================

computeFourierCorrectionMatrix :: FourierState -> FourierState -> Maybe FourierCorrectionMatrix
computeFourierCorrectionMatrix stateReduced stateFull = 
    let reducedLen = fsOrigLen stateReduced
        fullLen = fsOrigLen stateFull
        -- Build H^+ directly (no matrix inversion needed)
        h_pinv = buildHPseudoinverse reducedLen fullLen
    in Just h_pinv

-- =========================================================
-- Apply BLUE correction in frequency space
-- =========================================================

applyFourierCorrection :: FourierState
                       -> FourierState
                       -> FourierCorrectionMatrix
                       -> FourierState
applyFourierCorrection stateReduced stateFull h_pinv =
    let c_reduced = fsCoeffs stateReduced
        c_full = fsCoeffs stateFull
        reducedLen = fsOrigLen stateReduced
        fullLen = fsOrigLen stateFull
        
        -- Build H matrix
        h = buildH reducedLen fullLen
        
        -- Compute residual in frequency space
        c_reduced_padded = cMatVec h c_reduced
        residual = zipWith (-) c_full c_reduced_padded
        
        -- BLUE correction: delta = H^+ * residual
        delta_c = cMatVec h_pinv residual
        
        -- Apply correction
        c_corrected = zipWith (+) c_reduced delta_c
        
    in stateReduced { fsCoeffs = c_corrected }

-- =========================================================
-- Test wrapper with cross-validation
-- =========================================================

testFourierCorrection :: FourierCVConfig -> [Double] -> Maybe ([Double], [Double], [Double])
testFourierCorrection config buffer = do
    let n = length buffer
        cvLen = fcvLength config
        reducedHistLen = n - cvLen
        
        -- Split data: training vs full
        trainingData = take reducedHistLen buffer
        
        -- Initialize Fourier states (using actual FFT, not padded length)
    let stateReduced = initializeFourierState trainingData
        stateFull = initializeFourierState buffer
        
    -- Compute correction matrix (H pseudoinverse)
    h_pinv <- computeFourierCorrectionMatrix stateReduced stateFull
    
    -- Apply correction
    let stateCorrected = applyFourierCorrection stateReduced stateFull h_pinv
        
        -- Evaluate in time domain (all indices)
        allIndices = [0..n-1]
        
        -- Uncorrected: zero-pad and IFFT
        h = buildH (fsOrigLen stateReduced) (fsOrigLen stateFull)
        c_uncorrected_padded = cMatVec h (fsCoeffs stateReduced)
        uncorrectedPreds = map realPart $ ifft c_uncorrected_padded
        
        -- Corrected: zero-pad corrected coefficients and IFFT
        c_corrected_padded = cMatVec h (fsCoeffs stateCorrected)
        correctedPreds = map realPart $ ifft c_corrected_padded
        
    return (buffer, uncorrectedPreds, correctedPreds)

-- =========================================================
-- Diagnostic: show errors for CV region only
-- =========================================================

testFourierCorrectionWithMetrics :: FourierCVConfig -> [Double] -> IO ()
testFourierCorrectionWithMetrics config buffer = do
    case testFourierCorrection config buffer of
        Nothing -> putStrLn "ERROR: Correction failed"
        Just (actual, uncorrected, corrected) -> do
            let n = length buffer
                cvLen = fcvLength config
                cvStart = n - cvLen
                
                -- CV region only
                cv_actual = drop cvStart actual
                cv_uncorrected = drop cvStart uncorrected
                cv_corrected = drop cvStart corrected
                
                error_uncorrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip cv_actual cv_uncorrected]
                error_corrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip cv_actual cv_corrected]
                improvement = (error_uncorrected - error_corrected) / error_uncorrected * 100
            
            putStrLn "Cross-Validation Region Errors:"
            putStrLn $ "  Uncorrected: " ++ show error_uncorrected
            putStrLn $ "  Corrected:   " ++ show error_corrected
            putStrLn $ "  Improvement: " ++ show improvement ++ "%"


-- =========================================================
-- Apply BLUE correction iteratively (for analysis)
-- Returns list of (iteration number, corrected state)
-- =========================================================

applyFourierCorrectionIterative :: Int 
                                -> FourierState
                                -> FourierState
                                -> FourierCorrectionMatrix
                                -> [(Int, FourierState)]
applyFourierCorrectionIterative numIterations stateReduced stateFull h_pinv =
    let c_full = fsCoeffs stateFull
        reducedLen = fsOrigLen stateReduced
        fullLen = fsOrigLen stateFull
        h = buildH reducedLen fullLen
        
        -- Single iteration step
        applyOneStep state iterNum =
            let c_current = fsCoeffs state
                c_padded = cMatVec h c_current
                residual = zipWith (-) c_full c_padded
                delta_c = cMatVec h_pinv residual
                c_corrected = zipWith (+) c_current delta_c
            in (iterNum, state { fsCoeffs = c_corrected })
        
        -- Generate iterations
        iterations = scanl (\(_, s) i -> applyOneStep s i) (0, stateReduced) [1..numIterations]
        
    in iterations