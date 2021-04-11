# Purpose
Determine how well the CMC method performs when we synthetically rotate the exact same scan or a very similar scan by a random angle. 
If the CMC method works well, then it should estimate the random rotation closely (e.g., within \pm 3 degrees)

# Experiments/Results
1. "Self-Comparison, Random Rotation" Compare a downsampled scan to a synthetically rotated (so that ground-truth is known) exact copy. 
  - Results: rotations are well-estimated.
2. "Copy Comparison, No Rotation" Compare a downsampled scan to a shifted copy of the same scan (changing the offset argument of x3p_sample). 
  - Restults: rotation (0 degrees here) is well-estimated.
3. "Copy Comparison, Random Rotation" Compare a downsampled scan to a shifted copy of the same scan that has been synthetically rotated. 
  - Results: rotations are well-estimated.

# Conclusions
The results of these experiments indicate that the CMC method performs well for self or near-self comparisons, even with synthetic rotations.
However, CMC method still doesn't seem to work as desired when matching cartridge case pairs are compared to each other (observed in other experiments).
There are realistically 2 courses of actions:
1. Change the preprocessing procedures to highlight similarities between matching (unclear how to exactly accomplish this)
2. Extract additional/alternative similarity features during the comparison procedure.

There appears to be a slight negative association between the rotation estimation error and the CMC count, indicating that cartridge case pairs for which the rotation estimation is off results in a lower CMC count. This makes intuitive sense. This negative relationship is more notable for the original method than the High CMC method. Might be indicative of the High CMC method being more robust to errors in the rotation estimation (more specifically, in a lack of consensus amongst the similarity features.)