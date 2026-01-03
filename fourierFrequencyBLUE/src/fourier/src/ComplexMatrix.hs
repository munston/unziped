{-# LANGUAGE ScopedTypeVariables #-}

module ComplexMatrix
    ( CMatrix
    , hermitian
    , cMatMul
    , cMatVec
    , cInvertMatrix
    , cPrintMatrix
    , complexNorm
    ) where

import Data.Complex

type CMatrix = [[Complex Double]]

-- ============================================================================
-- Complex Matrix Operations
-- ============================================================================

-- Matrix transpose (helper for Hermitian)
matTranspose :: [[a]] -> [[a]]
matTranspose = foldr (zipWith (:)) (repeat [])

-- Hermitian transpose (conjugate transpose)
hermitian :: CMatrix -> CMatrix
hermitian = matTranspose . map (map conjugate)

-- Complex matrix multiplication
cMatMul :: CMatrix -> CMatrix -> CMatrix
cMatMul a b =
    [[sum $ zipWith (*) row col | col <- matTranspose b] | row <- a]

-- Complex matrix-vector multiplication
cMatVec :: CMatrix -> [Complex Double] -> [Complex Double]
cMatVec m v = [sum $ zipWith (*) row v | row <- m]

-- Complex norm of a complex number (for column normalization)
complexNorm :: Complex Double -> Double
complexNorm c = sqrt $ realPart (conjugate c * c)

-- ============================================================================
-- Complex Matrix Inversion (Gauss-Jordan)
-- ============================================================================

-- Same algorithm as real version, but with Complex Double
cInvertMatrix :: CMatrix -> Maybe CMatrix
cInvertMatrix m
    | rows /= cols = Nothing
    | otherwise = cGaussJordan augmented
  where
    rows = length m
    cols = length (head m)
    identity = [[ if i == j then 1 :+ 0 else 0 :+ 0 
                | j <- [0..cols-1]] | i <- [0..rows-1]]
    augmented = zipWith (++) m identity

cGaussJordan :: CMatrix -> Maybe CMatrix
cGaussJordan m = do
    reduced <- cForwardElimination m 0
    result <- cBackSubstitution reduced (length reduced - 1)
    return $ map (drop (length result)) result

cForwardElimination :: CMatrix -> Int -> Maybe CMatrix
cForwardElimination m row
    | row >= length m = Just m
    | otherwise = do
        m' <- cPivotRow m row
        let m'' = cEliminateColumn m' row
        cForwardElimination m'' (row + 1)

cPivotRow :: CMatrix -> Int -> Maybe CMatrix
cPivotRow m row =
    case cFindPivot m row row of
        Nothing -> Nothing
        Just pivotIdx -> Just $ cSwapRows m row pivotIdx

cFindPivot :: CMatrix -> Int -> Int -> Maybe Int
cFindPivot m col row
    | row >= length m = Nothing
    | magnitude (m !! row !! col) > 1e-10 = Just row
    | otherwise = cFindPivot m col (row + 1)

cSwapRows :: CMatrix -> Int -> Int -> CMatrix
cSwapRows m i j =
    [ if k == i then m !! j
      else if k == j then m !! i
      else m !! k
    | k <- [0..length m - 1]]

cEliminateColumn :: CMatrix -> Int -> CMatrix
cEliminateColumn m row =
    [ if i <= row then m !! i
      else let pivot = m !! row !! row
               factor = (m !! i !! row) / pivot
           in zipWith (\a b -> a - factor * b) (m !! i) (m !! row)
    | i <- [0..length m - 1]]

cBackSubstitution :: CMatrix -> Int -> Maybe CMatrix
cBackSubstitution m row
    | row < 0 = Just m
    | magnitude (m !! row !! row) < 1e-10 = Nothing
    | otherwise = do
        let m' = cScaleRow m row
            m'' = cEliminateAbove m' row
        cBackSubstitution m'' (row - 1)

cScaleRow :: CMatrix -> Int -> CMatrix
cScaleRow m row =
    [ if i /= row then m !! i
      else let pivot = m !! row !! row
           in map (/ pivot) (m !! i)
    | i <- [0..length m - 1]]

cEliminateAbove :: CMatrix -> Int -> CMatrix
cEliminateAbove m row =
    [ if i >= row then m !! i
      else let factor = m !! i !! row
           in zipWith (\a b -> a - factor * b) (m !! i) (m !! row)
    | i <- [0..length m - 1]]

-- ============================================================================
-- Pretty Printing
-- ============================================================================

-- Pretty print complex matrix
cPrintMatrix :: CMatrix -> IO ()
cPrintMatrix m = mapM_ cPrintRow m
  where
    cPrintRow row = putStrLn $ unwords $ map formatComplex row
    formatComplex (r :+ i) = 
        let rStr = show r
            iStr = show i
        in "(" ++ take 8 (rStr ++ repeat ' ') ++ " + " ++ 
           take 8 (iStr ++ repeat ' ') ++ "i)"
