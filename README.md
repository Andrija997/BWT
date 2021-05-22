This is a GitHub repo containting novice implementation of BWT and FM index lookup. <br/>
Project is tested on string lookup using Big Data genoms abailable at: https://ftp.ncbi.nlm.nih.gov/genomes/<br/>
Youtube video (in Serbian) explaining BWT and results of project: https://www.youtube.com/watch?v=OH0kSHD_uno&ab_channel=ivamil<br/>
<br/>
To start this code you'll need to have access to special files containting sorted BW matrix or create your own one (using rotation_with_indices function -> This is the case when rotation matrix can't fit memory).<br/>
<br/>
Command line exec assumes that you have all files needed in 'Data/';<br/>
<br/>
Arguments needed:<br/>
	1. File name which contains sequnce eg. "13443_ref_Cara_1.0_chr1c.fa"<br/>
	2. Sequence for lookup eg. "CCTG"<br/>
	3. Suffix array offset eg. 1, 4, 16, 64 or 254<br/>
	4. Tally matrix offset eg. 1, 8, 32, 128 or  512<br/>


