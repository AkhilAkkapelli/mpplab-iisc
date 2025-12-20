# Kartikeyan's Branch

## Results Till Now

(Kartikeyan)
1. Normalization Routine working properly (passes random test cases!) ✅
2. Creating new MPI from coefficient arrays working properly. ✅
3. New MPI from integer and MPI to integer work properly! ✅✅
4. Scalar multiplication tests all run successfully ✅

## Caveats and Changes

(Kartikeyan)
1. In testing scalar multiplication, scalar kept in integer(4) range, otherwise overflow can happen if scalar is very large!
2. Scalar Multiplication routine changed so that `ishft` does not encounter negative values. Sign resolution done after absolute value multiplication!


Comment (Akhil) : Resolved
