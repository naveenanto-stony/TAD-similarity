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
       comma-separated list of two TAD filenames


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


### ANALYSIS:

A python script can be found in the scripts folder to plot the different metrics.
<br>
<pre>usage: analysis.py [-h] -i -metrics

optional arguments:
  -h, --help        show this help message and exit
  -i                input file path (out.csv from go)
  -metrics          Comma separated metrics(Column names) to compare |
                    Accepted column names: 
                    'Jaccord','jacc_0.95_0.05',
                    'jacc_0.9_0.1', 'jacc_0.85_0.15','jacc_0.15_0.85',
                    'jacc_0.5_0.5','Overlap','Dice'
</pre>

The column name 'jacc_0.9_0.1' corresponds to 90% weightage to jaccard index (boundary) and 10% weightage to VI metric (interval)

The output of the analysis consist of a graph plotting the list of metrics provided. It also outputs the correlation with VI for each metric.

#### Sample Analysis:

Code: `python3 analysis.py -i ../out.csv -metrics=VI,jacc_0.15_0.85,Jaccord`

Command line output:
<pre>
Correlation with VI
VI                1.000000
jacc_0.15_0.85    0.606020
Jaccord           0.537286
</pre>
Image output: <br>
![Similarity Plot](https://user-images.githubusercontent.com/45582545/49527838-b91d2f00-f880-11e8-93db-66702f8d2bab.png)
