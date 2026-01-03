{-# LANGUAGE ScopedTypeVariables #-}

module FourierFrequencyCorrection
  ( Sinusoid(..)
  , extractSpectralPeaks
  , evaluateSinusoids        -- ADD THIS
  , buildFrequencyJacobian
  , applyFrequencyCorrection
  , testFrequencyCorrection
  ) where

import ComplexMatrix
import Fourier
import Data.Complex
import Data.List (sortBy)
import Data.Ord (comparing)
import Numeric.FFT (fft)

-- =========================================================
-- Sinusoid representation
-- =========================================================

data Sinusoid = Sinusoid
  { sFreq :: Double      -- Frequency (cycles per sample)
  , sAmp :: Double       -- Amplitude
  , sPhase :: Double     -- Phase (radians)
  } deriving (Show, Eq)

-- =========================================================
-- Extract spectral peaks from FFT
-- =========================================================

-- Find top N spectral peaks using quadratic interpolation for sub-bin accuracy
extractSpectralPeaks :: Int -> [Double] -> [Sinusoid]
extractSpectralPeaks numPeaks signal = 
    let n = length signal
        coeffs = fft $ map (:+ 0) signal
        -- Only look at positive frequencies (up to Nyquist)
        halfN = n `div` 2
        
        -- Compute magnitude spectrum
        mags = [magnitude (coeffs !! k) | k <- [0..halfN]]
        
        -- Find peaks (local maxima)
        peaks = [(k, mags !! k) | k <- [1..halfN-1],
                                   mags !! k > mags !! (k-1),
                                   mags !! k > mags !! (k+1)]
        
        -- Sort by magnitude, take top N
        topPeaks = take numPeaks $ sortBy (flip $ comparing snd) peaks
        
        -- Refine each peak using quadratic interpolation
        refinePeak (k, _) =
            let alpha = magnitude (coeffs !! (k-1))
                beta  = magnitude (coeffs !! k)
                gamma = magnitude (coeffs !! (k+1))
                -- Quadratic interpolation for sub-bin accuracy
                delta = 0.5 * (alpha - gamma) / (alpha - 2*beta + gamma)
                freqBin = fromIntegral k + delta
                freq = freqBin / fromIntegral n
                
                -- Get amplitude and phase at refined frequency
                c = coeffs !! k
                amp = 2 * magnitude c / fromIntegral n  -- Factor of 2 for positive freqs only
                ph = phase c  -- Changed from 'phase' to 'ph'
            in Sinusoid freq amp ph  -- Changed from 'phase' to 'ph'
            
    in map refinePeak topPeaks

-- =========================================================
-- Evaluate sinusoidal model
-- =========================================================

evaluateSinusoids :: [Sinusoid] -> [Int] -> [Double]
evaluateSinusoids sines indices =
    [sum [sAmp s * sin (2 * pi * sFreq s * fromIntegral i + sPhase s) 
         | s <- sines]
    | i <- indices]

-- =========================================================
-- Build Frequency Jacobian Matrix
-- =========================================================

-- H[i,k] = ∂y(t_i)/∂f_k = 2πt_i · A_k · cos(2πf_k·t_i + φ_k)
buildFrequencyJacobian :: [Sinusoid] -> [Int] -> [[Double]]
buildFrequencyJacobian sines cvIndices =
    [[2 * pi * fromIntegral t * sAmp s * 
      cos (2 * pi * sFreq s * fromIntegral t + sPhase s)
     | s <- sines]
    | t <- cvIndices]

-- =========================================================
-- BLUE Frequency Correction
-- =========================================================

applyFrequencyCorrection :: [Sinusoid]   -- Initial sinusoids from reduced data
                        -> [Double]      -- Actual CV data
                        -> [Int]         -- CV indices
                        -> Maybe [Sinusoid]  -- Corrected sinusoids
applyFrequencyCorrection sinesInitial cvActual cvIndices = do
    -- Predict using initial frequencies
    let cvPredicted = evaluateSinusoids sinesInitial cvIndices
        
        -- Residual
        residual = zipWith (-) cvActual cvPredicted
        
        -- Build Jacobian
        h = buildFrequencyJacobian sinesInitial cvIndices
        ht = matTranspose h
        hth = matMul ht h
        
    -- Invert (may fail if ill-conditioned)
    hthInv <- matInvert hth
    
    let htr = matVec ht residual
        deltaFreq = matVec hthInv htr
        
        -- Apply corrections
        sinesCorrected = zipWith (\s df -> s { sFreq = sFreq s + df })
                                  sinesInitial
                                  deltaFreq
    
    return sinesCorrected

-- =========================================================
-- Helper: matrix-vector multiply
-- =========================================================

matVec :: [[Double]] -> [Double] -> [Double]
matVec m v = [sum $ zipWith (*) row v | row <- m]

matTranspose :: [[Double]] -> [[Double]]
matTranspose = transpose
  where transpose ([]:_) = []
        transpose m = map head m : transpose (map tail m)

matMul :: [[Double]] -> [[Double]] -> [[Double]]
matMul a b = [[sum $ zipWith (*) row col | col <- matTranspose b] | row <- a]

-- Simple matrix inversion (should use better algorithm for production)
matInvert :: [[Double]] -> Maybe [[Double]]
matInvert m = 
    let n = length m
        augmented = zipWith (++) m (identity n)
        identity n = [[if i == j then 1 else 0 | j <- [0..n-1]] | i <- [0..n-1]]
    in fmap (map (drop n)) (gaussJordan augmented)

gaussJordan :: [[Double]] -> Maybe [[Double]]
gaussJordan m = foldM reduce m [0..length m - 1]
  where
    reduce mat col = do
        -- Find pivot
        let pivot = maximum [abs (mat !! r !! col) | r <- [col..length mat - 1]]
        if pivot < 1e-10 then Nothing else do
            let pivotRow = head [r | r <- [col..length mat - 1], 
                                     abs (mat !! r !! col) == pivot]
            -- Swap rows
            let mat' = swapRows mat col pivotRow
            -- Scale pivot row
            let pivotVal = mat' !! col !! col
                mat'' = updateRow mat' col (map (/ pivotVal) (mat' !! col))
            -- Eliminate column
            return $ foldl (\m r -> if r == col then m 
                                     else let factor = m !! r !! col
                                          in updateRow m r (zipWith (\a b -> a - factor * b)
                                                                    (m !! r) (m !! col)))
                           mat'' [0..length mat - 1]
    
    swapRows m i j = [if k == i then m !! j 
                      else if k == j then m !! i 
                      else m !! k | k <- [0..length m - 1]]
    
    updateRow m i row = take i m ++ [row] ++ drop (i+1) m
    
    foldM f z [] = Just z
    foldM f z (x:xs) = case f z x of
                         Nothing -> Nothing
                         Just z' -> foldM f z' xs

-- =========================================================
-- Test wrapper
-- =========================================================

testFrequencyCorrection :: [Double] -> Int -> Int -> IO ()
testFrequencyCorrection signal cvLen numPeaks = do
    let n = length signal
        trainingData = take (n - cvLen) signal
        cvData = drop (n - cvLen) signal
        cvIndices = [n - cvLen .. n - 1]
        
    putStrLn "=== Frequency Correction Test ==="
    putStrLn ""
    
    -- Extract peaks from training data
    let sinesInitial = extractSpectralPeaks numPeaks trainingData
    
    putStrLn "Initial frequencies (from training data):"
    mapM_ (\(i, s) -> putStrLn $ "  Peak " ++ show i ++ ": f=" ++ 
                                  show (sFreq s) ++ ", A=" ++ show (sAmp s)) 
          (zip [1..] sinesInitial)
    putStrLn ""
    
    -- Uncorrected prediction
    let cvPredictedUncorrected = evaluateSinusoids sinesInitial cvIndices
        errorUncorrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip cvData cvPredictedUncorrected]
    
    putStrLn $ "Uncorrected CV error: " ++ show errorUncorrected
    putStrLn ""
    
    -- Apply frequency correction
    case applyFrequencyCorrection sinesInitial cvData cvIndices of
        Nothing -> putStrLn "ERROR: Could not compute frequency correction (singular matrix)"
        Just sinesCorrected -> do
            putStrLn "Corrected frequencies:"
            mapM_ (\(i, s) -> putStrLn $ "  Peak " ++ show i ++ ": f=" ++ 
                                          show (sFreq s) ++ ", A=" ++ show (sAmp s))
                  (zip [1..] sinesCorrected)
            putStrLn ""
            
            let cvPredictedCorrected = evaluateSinusoids sinesCorrected cvIndices
                errorCorrected = sqrt $ sum [(a - p)^2 | (a, p) <- zip cvData cvPredictedCorrected]
                improvement = (errorUncorrected - errorCorrected) / errorUncorrected * 100
            
            putStrLn $ "Corrected CV error: " ++ show errorCorrected
            putStrLn $ "Improvement: " ++ show improvement ++ "%"