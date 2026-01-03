{-# LANGUAGE ScopedTypeVariables #-}

module TestFourierCorrection where

import FourierCorrection
import Fourier
import Data.Complex
import Numeric.FFT (fft, ifft) 
import ComplexMatrix

-- ============================================================================
-- Run all tests
-- ============================================================================


testComplexChirp :: IO ()
testComplexChirp = do
    putStrLn "=== Test: Complex Chirp (Quadratic + Amplitude Modulation) ==="
    let n = 256
        -- Quadratic chirp: frequency goes from 1 to 8 (faster acceleration)
        -- Plus amplitude decay and multiple harmonics
        complexChirp = [ let t_norm = fromIntegral t / fromIntegral n
                             freq = 1 + 7 * t_norm * t_norm  -- quadratic frequency sweep
                             amp = exp (negate (t_norm * 2))  -- exponential decay
                             fundamental = sin (2 * pi * freq * fromIntegral t / fromIntegral n)
                             harmonic2 = 0.3 * sin (2 * pi * 2 * freq * fromIntegral t / fromIntegral n)
                             harmonic3 = 0.1 * sin (2 * pi * 3 * freq * fromIntegral t / fromIntegral n)
                         in amp * (fundamental + harmonic2 + harmonic3)
                       | t <- [0..n-1] ]
        config = FourierCVConfig { fcvLength = 128, fdeltaT = 1.0 }
    testFourierCorrectionWithMetrics config complexChirp
    putStrLn ""

testIterativeComplexChirp :: IO ()
testIterativeComplexChirp = do
    putStrLn "=== Test: Iterative Correction (Complex Chirp) ==="
    putStrLn ""
    
    let n = 256
        complexChirp = [ let t_norm = fromIntegral t / fromIntegral n
                             freq = 1 + 7 * t_norm * t_norm
                             amp = exp (negate (t_norm * 2))
                             fundamental = sin (2 * pi * freq * fromIntegral t / fromIntegral n)
                             harmonic2 = 0.3 * sin (2 * pi * 2 * freq * fromIntegral t / fromIntegral n)
                             harmonic3 = 0.1 * sin (2 * pi * 3 * freq * fromIntegral t / fromIntegral n)
                         in amp * (fundamental + harmonic2 + harmonic3)
                       | t <- [0..n-1] ]
        trainingData = take 128 complexChirp
        
        stateReduced = initializeFourierState trainingData
        stateFull = initializeFourierState complexChirp
        
    case computeFourierCorrectionMatrix stateReduced stateFull of
        Nothing -> putStrLn "ERROR: Could not compute correction matrix"
        Just h_pinv -> do
            let iterations = applyFourierCorrectionIterative 10 stateReduced stateFull h_pinv
                reducedLen = fsOrigLen stateReduced
                fullLen = fsOrigLen stateFull
                h = buildH reducedLen fullLen
                
            putStrLn "Iteration | CV Error"
            putStrLn "----------|----------"
            
            mapM_ (showIterationError complexChirp h) iterations
            
            putStrLn ""
testAbsurdlyComplexChirp :: IO ()
testAbsurdlyComplexChirp = do
    putStrLn "=== Test: Absurdly Complex Chirp ==="
    putStrLn ""
    
    let n = 256
    let absurdChirp = [ 
                let t_norm = fromIntegral t / fromIntegral n
                    ti = fromIntegral t
                    -- Polynomial chirps of various orders
                    chirp1 = sin (2 * pi * (0.5 + 15 * t_norm^4 - 8 * t_norm^3) * ti / fromIntegral n)
                    chirp2 = 0.6 * cos (2 * pi * (12 * t_norm^5 - 20 * t_norm^2 + 3) * ti / fromIntegral n)
                    -- Chaotic frequency modulation (logistic map behavior)
                    chaos = let r = 3.9
                                x = 0.5 + 0.3 * sin (2 * pi * t_norm)
                            in sin (2 * pi * (r * x * (1 - x) * 10) * ti / fromIntegral n)
                    -- Self-similar fractal component
                    fractal = sum [ (1 / fromIntegral k) * sin (2 * pi * (2^k) * t_norm * ti / fromIntegral n)
                                  | k <- [1..5] ]
                    -- Discontinuous amplitude
                    amp_disc = if (floor (t_norm * 7) :: Int) `mod` 2 == 0 then 1.0 else 0.3
                    -- Multi-timescale amplitude modulation
                    amp_fast = 1 + 0.4 * sin (20 * pi * t_norm)
                    amp_slow = exp (negate (2 * t_norm)) * (1 + 0.5 * cos (pi * t_norm))
                    amp_chaos = 0.8 + 0.3 * sin (50 * pi * t_norm * (1 + 0.2 * sin (8 * pi * t_norm)))
                    -- FM synthesis with varying modulation index
                    mod_index = 5 * (1 + sin (3 * pi * t_norm))
                    carrier_freq = 6 + 4 * t_norm^2
                    modulator = sin (2 * pi * 1.3 * ti / fromIntegral n)
                    fm = sin (2 * pi * carrier_freq * ti / fromIntegral n + mod_index * modulator)
                    -- Waveshaping/distortion
                    distort x = tanh (2 * x)
                    -- Beat frequencies
                    beat1 = sin (2 * pi * 7.1 * ti / fromIntegral n)
                    beat2 = sin (2 * pi * 7.3 * ti / fromIntegral n)
                    -- Combine with nonlinear mixing
                    combined = chirp1 * amp_slow + chirp2 * amp_fast * amp_disc + 
                               0.5 * chaos + 0.3 * fractal + 
                               0.4 * fm * amp_chaos +
                               0.2 * distort (beat1 * beat2)
                in combined
              | t <- [0..n-1] ]
        trainingData = take 128 absurdChirp
        
        stateReduced = initializeFourierState trainingData
        stateFull = initializeFourierState absurdChirp
        
    case computeFourierCorrectionMatrix stateReduced stateFull of
        Nothing -> putStrLn "ERROR: Could not compute correction matrix"
        Just h_pinv -> do
            let iterations = applyFourierCorrectionIterative 10 stateReduced stateFull h_pinv
                reducedLen = fsOrigLen stateReduced
                fullLen = fsOrigLen stateFull
                h = buildH reducedLen fullLen
                
            putStrLn "Signal includes:"
            putStrLn "  - 3 interfering chirps (cubic, quadratic, sinusoidal frequency)"
            putStrLn "  - Multiple amplitude modulations (exponential decay, sinusoidal)"
            putStrLn "  - Phase modulation"
            putStrLn "  - FM synthesis"
            putStrLn ""
            putStrLn "Iteration | CV Error"
            putStrLn "----------|----------"
            
            mapM_ (showIterationError absurdChirp h) iterations
            
            putStrLn ""

-- Add to go
go :: IO ()
go = do
    testChirp
    testLinearRamp
    testStep
    testBeating
    testComplexChirp
    testIterative
    testIterativeComplexChirp
    testAbsurdlyComplexChirp  -- NEW

-- ============================================================================
-- Test helper function
-- ============================================================================

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

-- ============================================================================
-- Test Cases
-- ============================================================================

testChirp :: IO ()
testChirp = do
    putStrLn "=== Test: Chirp Signal ==="
    let n = 256
        chirp = [ sin (2 * pi * (1 + fromIntegral t / fromIntegral n) * fromIntegral t / fromIntegral n)
                | t <- [0..n-1] ]
        config = FourierCVConfig { fcvLength = 128, fdeltaT = 1.0 }
    testFourierCorrectionWithMetrics config chirp
    putStrLn ""

testLinearRamp :: IO ()
testLinearRamp = do
    putStrLn "=== Test: Linear Ramp ==="
    let n = 256
        ramp = [ fromIntegral t / fromIntegral n | t <- [0..n-1] ]
        config = FourierCVConfig { fcvLength = 128, fdeltaT = 1.0 }
    testFourierCorrectionWithMetrics config ramp
    putStrLn ""

testStep :: IO ()
testStep = do
    putStrLn "=== Test: Step Function ==="
    let n = 256
        step = [ if t < 128 then 0.0 else 1.0 | t <- [0..n-1] ]
        config = FourierCVConfig { fcvLength = 128, fdeltaT = 1.0 }
    testFourierCorrectionWithMetrics config step
    putStrLn ""

testBeating :: IO ()
testBeating = do
    putStrLn "=== Test: Beating Signal ==="
    let n = 256
        freq1 = 4.0
        freq2 = 4.1
        beating = [ 0.5 * sin (2 * pi * freq1 * fromIntegral t / fromIntegral n) 
                  + 0.5 * sin (2 * pi * freq2 * fromIntegral t / fromIntegral n)
                  | t <- [0..n-1] ]
        config = FourierCVConfig { fcvLength = 128, fdeltaT = 1.0 }
    testFourierCorrectionWithMetrics config beating
    putStrLn ""

-- ============================================================================
-- Detailed analysis test
-- ============================================================================

testDetailed :: IO ()
testDetailed = do
    putStrLn "=== Detailed Test: Chirp Signal ==="
    putStrLn ""
    
    let n = 256
        chirp = [ sin (2 * pi * (1 + fromIntegral t / fromIntegral n) * fromIntegral t / fromIntegral n)
                | t <- [0..n-1] ]
        config = FourierCVConfig { fcvLength = 128, fdeltaT = 1.0 }
    
    case testFourierCorrection config chirp of
        Nothing -> putStrLn "ERROR: Test failed"
        Just (actual, uncorrected, corrected) -> do
            let cvStart = 128
                
                -- Full signal errors
                error_full_uncorrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip actual uncorrected]
                error_full_corrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip actual corrected]
                
                -- First half (training region)
                first_actual = take cvStart actual
                first_uncorrected = take cvStart uncorrected
                first_corrected = take cvStart corrected
                error_first_uncorrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip first_actual first_uncorrected]
                error_first_corrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip first_actual first_corrected]
                
                -- Second half (CV region)
                second_actual = drop cvStart actual
                second_uncorrected = drop cvStart uncorrected
                second_corrected = drop cvStart corrected
                error_second_uncorrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip second_actual second_uncorrected]
                error_second_corrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip second_actual second_corrected]
            
            putStrLn "Full Signal:"
            putStrLn $ "  Uncorrected error: " ++ show error_full_uncorrected
            putStrLn $ "  Corrected error:   " ++ show error_full_corrected
            putStrLn $ "  Improvement: " ++ show ((error_full_uncorrected - error_full_corrected) / error_full_uncorrected * 100) ++ "%"
            putStrLn ""
            
            putStrLn "Training Region (samples 0-127):"
            putStrLn $ "  Uncorrected error: " ++ show error_first_uncorrected
            putStrLn $ "  Corrected error:   " ++ show error_first_corrected
            putStrLn $ "  Change: " ++ show ((error_first_uncorrected - error_first_corrected) / error_first_uncorrected * 100) ++ "%"
            putStrLn ""
            
            putStrLn "CV Region (samples 128-255):"
            putStrLn $ "  Uncorrected error: " ++ show error_second_uncorrected
            putStrLn $ "  Corrected error:   " ++ show error_second_corrected
            putStrLn $ "  Improvement: " ++ show ((error_second_uncorrected - error_second_corrected) / error_second_uncorrected * 100) ++ "%"
            putStrLn ""
            
            putStrLn "Sample values at boundary (indices 125-130):"
            mapM_ (\i -> putStrLn $ "  [" ++ show i ++ "] actual: " ++ show (actual !! i) ++
                                    ", uncorrected: " ++ show (uncorrected !! i) ++
                                    ", corrected: " ++ show (corrected !! i)) [125..130]


testIterative :: IO ()
testIterative = do
    putStrLn "=== Test: Iterative Correction (Chirp) ==="
    putStrLn ""
    
    let n = 256
        chirp = [ sin (2 * pi * (1 + fromIntegral t / fromIntegral n) * fromIntegral t / fromIntegral n)
                | t <- [0..n-1] ]
        config = FourierCVConfig { fcvLength = 128, fdeltaT = 1.0 }
        trainingData = take 128 chirp
        
        -- Build states
        stateReduced = initializeFourierState trainingData
        stateFull = initializeFourierState chirp
        
    case computeFourierCorrectionMatrix stateReduced stateFull of
        Nothing -> putStrLn "ERROR: Could not compute correction matrix"
        Just h_pinv -> do
            let iterations = applyFourierCorrectionIterative 10 stateReduced stateFull h_pinv
                reducedLen = fsOrigLen stateReduced
                fullLen = fsOrigLen stateFull
                h = buildH reducedLen fullLen
                
            putStrLn "Iteration | CV Error"
            putStrLn "----------|----------"
            
            mapM_ (showIterationError chirp h) iterations
            
            putStrLn ""

showIterationError :: [Double] -> [[Complex Double]] -> (Int, FourierState) -> IO ()
showIterationError testData h (iterNum, state) = do
    let c_padded = cMatVec h (fsCoeffs state)
        time_predicted = map realPart $ ifft c_padded
        cvStart = 128
        cv_actual = drop cvStart testData
        cv_predicted = drop cvStart time_predicted
        error_cv = sqrt $ sum [(a - p)^2 | (a, p) <- zip cv_actual cv_predicted]
        
    putStrLn $ "    " ++ show iterNum ++ "     | " ++ show error_cv
