superpig-public
===============

Repository for publicly available code from the Superpig collaboration to search for ultra-short-period planets (USPs) in K2 data. Results have been published at: http://arxiv.org/abs/1603.06488

To run a search on a quarter of K2 data, follow these steps:

0. Download the data and make a file listing all data files (example: all_c00)

1. Run the EEBLS search for periodic signals using:
    python k2_superpig_search.py --fl allfiles > allfiles_results.txt

2. Use the iPython notebook Step2_K2_Identify_Candidates.ipynb to sift through for plausibly planet-like signals.

3. [Forthcoming] Information about the candidate stars

4. [Forthcoming] Transit light curve fitting 

5. [Forthcoming] Plots etc. used in Adams+ 2016