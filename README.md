# isorank
The IsoRank family of Network Alignment Algorithms

The IsoRank and IsoRank-N algorithms are designed for global alignment of multiple graphs against each other, especially when an independent measure of node similarity across the various graphs is available. They were designed for global alignment of protein-protein interaction (PPI) networks across species, with the sequence similarity of proteins being the node similarity measure. However, these are general algorithms that have been used for network alignment in other domains too (e.g., alignment of electron microscopy images). 

The executable is available at: 

- The original IsoRank algorithm was presented in Singh, Xu and Berger (RECOMB 2007 & PNAS 2008) and is made available here as IsoRank-Orig-2008

- Enhancements to it, as the IsoRank-N algorithm, were presented by Liao et al. (Bioinformatics 2009). Prof. Liao's lab has made improvements to the algorithm, with the aid of Kanghao Lu, Cheng-Yu Ma and Guan-Chung Chen. This source is made available here as IsoRankN-2018. 

- The executable available at http://isorank.csail.mit.edu corresponds to the IsoRank-N algorithm.

- A database of functional orthologs, as predicted by IsoRank, are available at http://isobase.csail.mit.edu

## References:
- R Singh, J Xu, B Berger. *Global alignment of multiple protein interaction networks with application to functional orthology detection*. Proceedings of the National Academy of Sciences 105 (35), 12763-12768

- CS Liao, K Lu, M Baym, R Singh, B Berger. *IsoRankN: spectral methods for global alignment of multiple protein networks*. Bioinformatics 25 (12), i253-i258

- R Singh, J Xu, B Berger. *Pairwise global alignment of protein interaction networks by matching neighborhood topology*. Annual International Conference on Research in Computational Molecular Biology (RECOMB), 2007

- D Park, R Singh, M Baym, CS Liao, B Berger. *IsoBase: a database of functionally related proteins across PPI networks*. Nucleic acids research 39 (suppl_1), D295-D300

- R Singh, J Xu, B Berger. *Global alignment of multiple protein interaction networks*. Pac Symp of Biocomputing 2008, 303-314
