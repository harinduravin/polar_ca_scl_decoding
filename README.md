### Algorithm Details

Successive Cancellation List (SCL), Successive Cancellation (SC), and space-efficient SC decoding algorithms are based on [1]. One key difference in this implementation compared to the above paper is the incorporation of LLR values in the algorithm instead of a non-LLR implementation. Therefore, the LLR-based implementation is supported by an additional LLR path metric [2]. The polar code construction algorithm in this repository is based on the Gaussian Approximation (GA) method [3]. The use of CRC in polar SCL decoding is explained in [4]. Furthermore, [This GitHub repository](https://github.com/tavildar/Polar) implements CRC-aided SCL decoding using a CRC matrix instead of a CRC polynomial. Apart from this difference, its SCL decoding implementation aligns closely with the approach used in this repository.
### References

1. Tal, Ido, and Alexander Vardy. "List decoding of polar codes." *IEEE Transactions on Information Theory*, vol. 61, no. 5, 2015, pp. 2213-2226.
2. B. Yuan and K. K. Parhi, "Successive cancellation list polar decoder using log-likelihood ratios," *2014 48th Asilomar Conference on Signals, Systems and Computers*, Pacific Grove, CA, USA, 2014, pp. 548-552.
3. Wu, Daolong, Ying Li, and Yue Sun. "Construction and block error rate analysis of polar codes over AWGN channel based on Gaussian approximation." *IEEE Communications Letters*, vol. 18, no. 7, 2014, pp. 1099-1102.
4. K. Niu and K. Chen, "CRC-Aided Decoding of Polar Codes," *IEEE Communications Letters*, vol. 16, no. 10, pp. 1668-1671, October 2012.
