Purpose:
Determine whether rotating a scan by an estimated "correct" angle would increase the CMC count value for known matches while leaving non-matches unaffected.

Experiments/Results:
1. Run CMC method without pre-rotating scans.
2. Run CMC method after estimating the optimal rotation using a polar domain transformation + CCF.

The CMC counts assigned to the matching scans didn't categorically increase. In particular, the CMC method assigned to some pairs decreased while for others it increased. The polar domain transformation does not seem to estimate the correct alignment for all known match pairs, which negatively impacts the CMC results. In particular, the CCF is thrown off if there is a large feature in one scan that isn't in the other (e.g., a chip in one breech face surface). The CCF seems to "attract" to such large features, even if they aren't shared between the two scans.

Conclusions/Insights:
The fact that the polar domain rotation estimation method doesn't work for all known matches indicates that the preprocessing procedure needs to change (e.g., to remove the large features unique to a particular cartridge case) or other features should be extracted (other than just CCF, translation, and rotation). In many cases, the similarity features simply do not reach a consensus as the CMC method (implicitly) assumes they should.

For some cartridge case scans, the decrease in CMC count after pre-rotating can be attributed to excessively large features in the scan that aren't shared with its known matches - take scan D3 from TopMatch, for instance. This causes the polar domain rotation estimation to be thrown off. I need to look more into other scans to determine whether there are any pairs such that, even after pre-rotating to approximately correct alignmnent, the CMC method still performed worse.