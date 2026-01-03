{-# LANGUAGE ScopedTypeVariables #-}

module TestHMatrix where

import ComplexMatrix
import Fourier
import Data.Complex
import Numeric.FFT (fft, ifft)


testIterativeBLUE :: IO ()
testIterativeBLUE = do
    putStrLn "=== Iterative BLUE Correction Test (Chirp) ==="
    putStrLn ""
    
    let n = 256
        testData = [ sin (2 * pi * (1 + fromIntegral t / fromIntegral n) * fromIntegral t / fromIntegral n)
                   | t <- [0..n-1] ]
        trainingData = take 128 testData
        c_128_original = fft $ map (:+ 0) trainingData
        c_256 = fft $ map (:+ 0) testData
        h = buildH 128 256
        h_pinv = buildHPseudoinverse 128 256
    
    putStrLn "Iterating BLUE correction..."
    putStrLn ""
    
    -- Iterate correction
    let iterations = 10
        results = iterate (applyOneIteration h h_pinv c_256) (c_128_original, 0)
        
    -- Show results for each iteration
    mapM_ (showIteration testData h) (take (iterations + 1) results)
    
    putStrLn "=== Test Complete ==="

-- Apply one iteration of BLUE correction
applyOneIteration :: [[Complex Double]] -> [[Complex Double]] -> [Complex Double] 
                  -> ([Complex Double], Int) -> ([Complex Double], Int)
applyOneIteration h h_pinv c_256 (c_128, iter) =
    let c_padded = cMatVec h c_128
        residual = zipWith (-) c_256 c_padded
        delta = cMatVec h_pinv residual
        c_new = zipWith (+) c_128 delta
    in (c_new, iter + 1)

-- Show metrics for one iteration
showIteration :: [Double] -> [[Complex Double]] -> ([Complex Double], Int) -> IO ()
showIteration testData h (c_128, iter) = do
    let c_padded = cMatVec h c_128
        time_domain = map realPart $ ifft c_padded
        second_half_actual = drop 128 testData
        second_half_predicted = drop 128 time_domain
        error_second = sqrt $ sum [(a - b)^2 | (a, b) <- zip second_half_actual second_half_predicted]
        
    putStrLn $ "Iteration " ++ show iter ++ ": Second half error = " ++ show error_second

-- Update the go function
go :: IO ()
go = do
    testBLUEWithDiscontinuity
    testBLUEWithChirp
    testBLUEWithStep
    testIterativeBLUE
-- ============================================================================
-- H Matrix Construction (128 -> 256 zero-padding)
-- ============================================================================

-- Forward transform: C_128 -> C_256 via zero-padding
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

-- ============================================================================
-- H Pseudoinverse Construction (direct, no matrix inversion needed)
-- ============================================================================

-- Since H^H H = 4I, we have H^+ = (1/4) H^H
buildHPseudoinverse :: Int -> Int -> [[Complex Double]]
buildHPseudoinverse origLen targetLen =
    [ buildRow i | i <- [0..origLen-1] ]
  where
    halfOrig = origLen `div` 2
    
    buildRow i
        | i < halfOrig =
            -- First 64 rows: extract from first 64 cols of C_256
            [ if j == i then (0.5 :+ 0) else (0 :+ 0)
            | j <- [0..targetLen-1] ]
            
        | otherwise =
            -- Last 64 rows: extract from last 64 cols of C_256
            let targetIdx = i + (targetLen - origLen)
            in [ if j == targetIdx then (0.5 :+ 0) else (0 :+ 0)
               | j <- [0..targetLen-1] ]

-- ============================================================================
-- Test with signal that has discontinuity at boundary
-- ============================================================================

testBLUEWithDiscontinuity :: IO ()
testBLUEWithDiscontinuity = do
    putStrLn "=== BLUE Correction Test (Discontinuous Signal) ==="
    putStrLn ""
    
    -- Test signal: linear ramp (always increasing, discontinuity at wrap)
    let n = 256
        testData = [ fromIntegral t / fromIntegral n | t <- [0..n-1] ]
    
    putStrLn "Generated 256-point linear ramp (0 to 1)"
    putStrLn $ "  Value at index 127: " ++ show (testData !! 127)
    putStrLn $ "  Value at index 128: " ++ show (testData !! 128)
    putStrLn $ "  (Clear discontinuity when periodic)"
    putStrLn ""
    
    runBLUETest testData

testBLUEWithChirp :: IO ()
testBLUEWithChirp = do
    putStrLn "=== BLUE Correction Test (Chirp Signal) ==="
    putStrLn ""
    
    -- Chirp signal (frequency increases linearly)
    let n = 256
        testData = [ sin (2 * pi * (1 + fromIntegral t / fromIntegral n) * fromIntegral t / fromIntegral n)
                   | t <- [0..n-1] ]
    
    putStrLn "Generated 256-point chirp (frequency 1 -> 2)"
    putStrLn ""
    
    runBLUETest testData

testBLUEWithStep :: IO ()
testBLUEWithStep = do
    putStrLn "=== BLUE Correction Test (Step Function) ==="
    putStrLn ""
    
    -- Step function at boundary
    let n = 256
        testData = [ if t < 128 then 0.0 else 1.0 | t <- [0..n-1] ]
    
    putStrLn "Generated 256-point step function (0->1 at index 128)"
    putStrLn ""
    
    runBLUETest testData

-- ============================================================================
-- Core BLUE test logic (using direct pseudoinverse)
-- ============================================================================

runBLUETest :: [Double] -> IO ()
runBLUETest testData = do
    let n = length testData
        trainingData = take 128 testData
    
    -- FFT both
    let c_128 = fft $ map (:+ 0) trainingData
        c_256 = fft $ map (:+ 0) testData
    
    putStrLn "Computed FFTs"
    putStrLn ""
    
    -- Build H and H^+
    let h = buildH 128 256
        h_pinv = buildHPseudoinverse 128 256
    
    putStrLn "Built H and H^+ (direct construction, no inversion)"
    putStrLn ""
    
    -- Baseline: H * C_128
    let c_128_padded = cMatVec h c_128
        residual = zipWith (-) c_256 c_128_padded
    
    -- BLUE correction: delta_C_128 = H^+ * residual
    let delta_c_128 = cMatVec h_pinv residual
        c_corrected_128 = zipWith (+) c_128 delta_c_128
        c_corrected_256 = cMatVec h c_corrected_128
    
    -- Errors (frequency domain)
    let error_freq_uncorrected = sqrt $ sum [magnitude r ^ 2 | r <- residual]
        residual_corrected = zipWith (-) c_256 c_corrected_256
        error_freq_corrected = sqrt $ sum [magnitude r ^ 2 | r <- residual_corrected]
    
    putStrLn "Frequency domain errors:"
    putStrLn $ "  Uncorrected: " ++ show error_freq_uncorrected
    putStrLn $ "  Corrected:   " ++ show error_freq_corrected
    putStrLn $ "  Improvement: " ++ show ((error_freq_uncorrected - error_freq_corrected) / error_freq_uncorrected * 100) ++ "%"
    putStrLn ""
    
    -- Errors (time domain)
    let time_uncorrected = map realPart $ ifft c_128_padded
        time_corrected = map realPart $ ifft c_corrected_256
        error_time_uncorrected = sqrt $ sum [(a - b)^2 | (a, b) <- zip testData time_uncorrected]
        error_time_corrected = sqrt $ sum [(a - b)^2 | (a, b) <- zip testData time_corrected]
    
    putStrLn "Time domain errors:"
    putStrLn $ "  Uncorrected: " ++ show error_time_uncorrected
    putStrLn $ "  Corrected:   " ++ show error_time_corrected
    putStrLn $ "  Improvement: " ++ show ((error_time_uncorrected - error_time_corrected) / error_time_uncorrected * 100) ++ "%"
    putStrLn ""
    
    -- Second half only
    let second_half_actual = drop 128 testData
        second_half_uncorrected = drop 128 time_uncorrected
        second_half_corrected = drop 128 time_corrected
        error_second_uncorrected = sqrt $ sum [(a - b)^2 | (a, b) <- zip second_half_actual second_half_uncorrected]
        error_second_corrected = sqrt $ sum [(a - b)^2 | (a, b) <- zip second_half_actual second_half_corrected]
    
    putStrLn "Second half (samples 128-255) errors:"
    putStrLn $ "  Uncorrected: " ++ show error_second_uncorrected
    putStrLn $ "  Corrected:   " ++ show error_second_corrected
    putStrLn $ "  Improvement: " ++ show ((error_second_uncorrected - error_second_corrected) / error_second_uncorrected * 100) ++ "%"
    putStrLn ""
    
    putStrLn "=== Test Complete ==="
    putStrLn ""