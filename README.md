# TAD-similarity

This repo contains the modified Go source code to run the different metrics based on "Quantifying the similarity of topological domains across normal and cancer human cell types" (Sauerwald and Kingsford), as well as a Python script to reproduce all statistics and figures used in the paper.
The source was taken from https://github.com/Kingsford-Group/localtadsim 

### EVALUATING TAD SET SIMILARITY

Build localdiff package: `go build go/src/localdiff`

Executing the binary : `./localdiff -h`

Usage of ./localdiff:
  * -gamma string
    	 if optimizing gamma, use 'opt,n' where n is the median TAD size to optimize for
  * -o string
     output filename
  * -res int
       resolution of Hi-C data (default 1)
  * -tad string
       comma-separated list of two TAD filenames or file patterns if optimizing gamma"""


We will go through an example, using the files provided in the go/example folder. These example TAD sets are random sample created to compare the different metrics. In order to compare these two TAD sets with different metrics, run the following command:

`TAD-similarity$ ./localdiff -tad=go/example/1.txt,go/example/2.txt -res=1 -o=out.csv`

Wrote output values to out.csv

### OUTPUT FILE DESCRIPTION
The output file is a csv which contains the following columns for each start and end of both TAD sets:
  * start : specifies the start boundary of the tad
  * end : specifies the end boundary of the tad
  * VI : Similarity score based on https://github.com/Kingsford-Group/localtadsim 
  * Jaccord : Jaccard similarity score
  * Overlap : Overlap similarity between Tad points
  * Dice : Dice similarity score
  * VI_boundaryless : VI similarity score excluding the boundary points


### ANALYSIS
