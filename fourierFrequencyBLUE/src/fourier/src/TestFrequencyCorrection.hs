{-# LANGUAGE ScopedTypeVariables #-}

module TestFrequencyCorrection where

import FourierFrequencyCorrection
import Data.Complex

-- ============================================================================
-- Run all tests
-- ============================================================================

go :: IO ()
go = do
    testFrequencyDrift
    testChirp
    testMultipleSinusoids

-- ============================================================================
-- Test Cases
-- ============================================================================

testFrequencyDrift :: IO ()
testFrequencyDrift = do
    putStrLn "=== Test: Frequency Drift ==="
    putStrLn "Signal: sine wave with slow frequency drift"
    putStrLn ""
    
    let n = 256
        -- Frequency drifts from 0.05 to 0.07 cycles/sample
        signal = [sin (2 * pi * (0.05 + 0.02 * fromIntegral t / fromIntegral n) * fromIntegral t)
                 | t <- [0..n-1]]
        cvLen = 128
        numPeaks = 1
    
    testFrequencyCorrection signal cvLen numPeaks
    putStrLn ""

testChirp :: IO ()
testChirp = do
    putStrLn "=== Test: Chirp Signal ==="
    putStrLn "Signal: quadratic chirp (accelerating frequency)"
    putStrLn ""
    
    let n = 256
        -- Quadratic chirp: f(t) = f0 + k*t^2
        signal = [let t_norm = fromIntegral t / fromIntegral n
                      freq = 0.03 + 0.1 * t_norm * t_norm
                  in sin (2 * pi * freq * fromIntegral t)
                 | t <- [0..n-1]]
        cvLen = 128
        numPeaks = 1
    
    testFrequencyCorrection signal cvLen numPeaks
    putStrLn ""

testMultipleSinusoids :: IO ()
testMultipleSinusoids = do
    putStrLn "=== Test: Multiple Sinusoids with Drift ==="
    putStrLn "Signal: two sine waves, both with frequency drift"
    putStrLn ""
    
    let n = 256
        -- Two sinusoids with different drift rates
        signal = [let t_norm = fromIntegral t / fromIntegral n
                      f1 = 0.04 + 0.01 * t_norm
                      f2 = 0.09 + 0.02 * t_norm
                      sine1 = sin (2 * pi * f1 * fromIntegral t)
                      sine2 = 0.6 * sin (2 * pi * f2 * fromIntegral t + 0.5)
                  in sine1 + sine2
                 | t <- [0..n-1]]
        cvLen = 128
        numPeaks = 2
    
    testFrequencyCorrection signal cvLen numPeaks
    putStrLn ""

-- ============================================================================
-- Detailed diagnostic test
-- ============================================================================

testDetailed :: IO ()
testDetailed = do
    putStrLn "=== Detailed Test: Frequency Drift with Diagnostics ==="
    putStrLn ""
    
    let n = 256
        cvLen = 128
        trainingLen = n - cvLen
        
        -- Frequency drifts from 0.05 to 0.08
        signal = [sin (2 * pi * (0.05 + 0.03 * fromIntegral t / fromIntegral n) * fromIntegral t)
                 | t <- [0..n-1]]
        
        trainingData = take trainingLen signal
        cvData = drop trainingLen signal
        cvIndices = [trainingLen .. n - 1]
        
        -- Extract initial peaks
        sinesInitial = extractSpectralPeaks 1 trainingData
    
    putStrLn "Training segment (samples 0-127):"
    putStrLn $ "  Actual frequency range: 0.05 -> 0.065 cycles/sample"
    putStrLn ""
    
    putStrLn "CV segment (samples 128-255):"
    putStrLn $ "  Actual frequency range: 0.065 -> 0.08 cycles/sample"
    putStrLn ""
    
    putStrLn "Extracted from training FFT:"
    mapM_ (\s -> putStrLn $ "  Frequency: " ++ show (sFreq s) ++ 
                            ", Amplitude: " ++ show (sAmp s) ++
                            ", Phase: " ++ show (sPhase s)) sinesInitial
    putStrLn ""
    
    -- Predictions
    let cvPredUncorrected = evaluateSinusoids sinesInitial cvIndices
        errorUncorrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip cvData cvPredUncorrected]
    
    putStrLn "Uncorrected prediction:"
    putStrLn $ "  Using constant frequency from training"
    putStrLn $ "  CV error: " ++ show errorUncorrected
    putStrLn ""
    
    -- Apply correction
    case applyFrequencyCorrection sinesInitial cvData cvIndices of
        Nothing -> putStrLn "ERROR: Correction failed"
        Just sinesCorrected -> do
            putStrLn "After frequency correction:"
            mapM_ (\s -> putStrLn $ "  Frequency: " ++ show (sFreq s) ++ 
                                    ", Amplitude: " ++ show (sAmp s) ++
                                    ", Phase: " ++ show (sPhase s)) sinesCorrected
            putStrLn ""
            
            let cvPredCorrected = evaluateSinusoids sinesCorrected cvIndices
                errorCorrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip cvData cvPredCorrected]
                improvement = (errorUncorrected - errorCorrected) / errorUncorrected * 100
            
            putStrLn $ "  CV error: " ++ show errorCorrected
            putStrLn $ "  Improvement: " ++ show improvement ++ "%"
            putStrLn ""
            
            putStrLn "Analysis:"
            let freqChange = sFreq (head sinesCorrected) - sFreq (head sinesInitial)
            putStrLn $ "  Frequency adjustment: " ++ show freqChange ++ " cycles/sample"
            putStrLn $ "  (Expected: positive, since frequency increases in CV)"