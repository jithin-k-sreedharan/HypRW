# HypRW
Bayesian inference on network properties with partial crawl of the network using random walks.  
Software of the research paper "Inference in OSNs via Lightweight Partial Crawls" by Konstantin Avrachenkov, Bruno Ribeiro and Jithin K. Sreedharan to appear in the proceedings of ACM SIGMETRICS/PERFORMANCE 2016.

## Files
"friendster_community1_trimmed.edgelist" is a largest connected component extracted from one subgraph of Friendster. Data is collected from https://snap.stanford.edu/

The file "HypRW_friendster.py" is the script which does the random walk crawling and plots the approximate posterior distribution and the empirical distribution.

The Python script requires various modules like NetworkX, Scipy and NumPy.
