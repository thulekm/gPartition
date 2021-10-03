# gPartition
<br>A fast and accurate method for partitioning large alignments

<b>I.	About gPartition</b>
<br><b>gPartition</b> is a fast model-based program that automatically partition an alignment into subsets such that sites in a subset follow the same evolutionary process. The method aims at large alignments (especially genome alignments) which are hard to be partitioned by other methods. Beside consider the rate at each site, <b>gPartition</b> uses substitution model to properly group similar sites into subsets. Experiments on extremely large real datasets showed that <b>gPartition</b> was better than other partitioning methods tested. In addition, <b>gPartition</b> overcame the pitfall of site rates-based partitioning method that groups all invariant sites into one subset leading to incorrect trees.
<br>Using the <b>gPartition</b> method will make it possible to partitioning large or whole genome datasets. The result would help to enhances the accuracy of maximum likelihood tree inference process.

<b>II.	Versions</b>
<br>The <b>gPartition</b> program is available for Linux and Windows operating system.
<br><br>
<b>III.	Setup</b>
<br>1.	The program requires several software
<br>  - python 2 or python 3 (with numpy library)
<br>  - IQ-TREE (http://www.iqtree.org/)
<br>2.	Run “setup.py”
<br>Note that the setup.py program will create a new config.py file. You can directly edit the path of IQ-TREE in the config.py file.

<b>IV.	Commands</b>
<br><i>python gPartition.py -f [alignment] -o [output_folder] -k [minimum value of correlation between models]</i>
<br><br>The partition scheme will be written to the file [output_directory]/[alignment].FINAL.nex
<br>  <i>-f</i>   The path to a DNA alignment in Phylip format.
<br>  <i>-o</i>   The path to folder that contains the result files.
<br>  <i>-k</i>   The minimum value of correlation between models – optional. The default value is 0.9995.
<br><bR><b>Example:</b> <i>python gPartition.py -f Example.phy -o output</i>

