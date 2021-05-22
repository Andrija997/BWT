This is a GitHub repo containting novice implementation of BWT and FM index lookup. 
Project is tested on string lookup using Big Data genoms abailable at: https://ftp.ncbi.nlm.nih.gov/genomes/
Youtube video (in Serbian) explaining BWT and results of project: https://www.youtube.com/watch?v=OH0kSHD_uno&ab_channel=ivamil

To start this code you'll need to have access to special files containting sorted BW matrix or create your own one (using rotation_with_indices function -> This is the case when rotation matrix can't fit memory).

Command line exec assumes that you have all files needed in 'Data/';

Arguments needed:
	1. File name which contains sequnce eg. "13443_ref_Cara_1.0_chr1c.fa"
	2. Sequence for lookup eg. "CCTG"
	3. Suffix array offset eg. 1, 4, 16, 64 or 254
	4. Tally matrix offset eg. 1, 8, 32, 128 or  512


