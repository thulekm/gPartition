# gPartition
<br>A fast and accurate method for partitioning genome alignments

<b>I.	About gPartition</b>
<br><b>gPartition</b> is a fast model-based program that automatically partition an alignment into disjoint subsets such that sites in the same subset follow the same evolutionary process. The gPartition program is orders of magnitudes faster than the existing partitioning programs and applicable to analyze genome alignments with million sites. The partition schemes generated from gPartition help construct better maximum likelihood trees. 
<br><br>
<b>II.	Versions</b>
<br>The <b>gPartition</b> program is available for both Linux and Windows operating system.
<br><br>
<b>III.	Setup</b>
<br>1.	The program requires several software
<br>  - python 2 or python 3 (with numpy library)
<br> - IQ-TREE software (http://www.iqtree.org/)
<br>2.	Run “setup.py”
<br>Note that the setup.py program will create a new config.py file. You can directly edit the path of IQ-TREE in the config.py file.

<b>IV.	Commands</b>
<br>python gPartition.py -f <i>alignment </i> -o <i>output</i>
<br><br>The resulting partition scheme will be written to the <i>alignment</i>.FINAL.nex in the <i>output</i> directory.
<br>  <i>-f</i>   A nucleotide alignment in the Phylip format.
<br>  <i>-o</i>   A folder that contains the result files.
<br><bR><b>Example:</b> <i>python gPartition.py -f sample.phy -o output</i>

