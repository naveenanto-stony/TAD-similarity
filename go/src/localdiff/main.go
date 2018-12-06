package main

import (
	"../hicutil"
	"bufio"
	"flag"
	"fmt"
	"math"
	"math/rand"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

type bdyvi struct {
	start                  int
	end                    int
	vi                     float64
	vi_without_boundary    float64
	jacc                   float64
	hamming                float64
	overlap                float64
	centroidSquareDistance float64
	centroidAbsDistance    float64
	dice                   float64
	pval                   float64
}

type fn func(a float64, b float64) float64

func main() {

	runtime.GOMAXPROCS(1)
	//numCPU := runtime.NumCPU()
	//pstart := time.Now()

	// inputs should be 2 TAD files, a resolution parameter, and the output file name
	tadin := flag.String("tad", "", "comma-separated list of two TAD filenames or file patterns if optimizing gamma")
	gammain := flag.String("gamma", "", "if optimizing gamma, use 'opt,n' where n is the median TAD size to optimize for")
	res := flag.Int("res", 1, "resolution of Hi-C data")
	outfile := flag.String("o", "", "output filename")
	//pcount := flag.Int("proc", numCPU, "number of processes to allow")

	flag.Parse()

	tadfilelist := strings.Split(*tadin, ",")

	var gammaopt bool
	medtadlen := 0.0
	var err error
	if len(*gammain) > 0 {
		gammaopt = true
		gammadata := strings.Split(*gammain, ",")
		medtadlen, err = strconv.ParseFloat(gammadata[1], 64)
		if err != nil {
			fmt.Println("Error: couldn't convert median TAD length value to float, make sure to input i.e. '-gamma=opt,100' ")
			os.Exit(1)
		}
	} else {
		gammaopt = false
	}
	/*readinputt := time.Now()
	readinputtime := readinputt.Sub(pstart)
	fmt.Println("time to read inputs:",readinputtime)*/

	// read TAD files and process TAD lists to fill in non-TAD regions
	tadlists := processTADLists(tadfilelist, res, gammaopt, medtadlen)
	fmt.Println("done processing TAD lists and choosing optimal gamma")
	/*processtads := time.Now()
	processtime := processtads.Sub(readinputt)
	fmt.Println("time to process tads:",processtime)*/

	// calculate VI values at boundaries (using DP)
	//bdyvis := calcVIatBdys(tadlists)
	bdyvis := calcNaiveMetrics(tadlists)
	//bdyvis := calcJaccatBdysNaive(tadlists)
	/*calcvis := time.Now()
	calcvitime := calcvis.Sub(processtads)
	fmt.Println("time to calculate VIs:",calcvitime)*/

	// calculate all p-values, select significant points
	//convergencecondition := 1e-5
	//numCPU = *pcount
	//runtime.GOMAXPROCS(numCPU)
	//sigpts := calcAllPvals(tadlists, bdyvis, numCPU, convergencecondition)
	//fmt.Println("done calculating all p-values")
	/*calcpvals := time.Now()
	calcpvaltime := calcpvals.Sub(calcvis)
	fmt.Println("time to calculate p-values:", calcpvaltime)*/

	// identify dominating points from significant ones
	//runtime.GOMAXPROCS(1)
	//dompts := findDomPts(sigpts)
	//fmt.Println("done finding dominating points")
	/*calcdompts := time.Now()
	calcdomptstime := calcdompts.Sub(calcpvals)
	fmt.Println("time to calculate dominating points:", calcdomptstime)*/

	// save results to a file
	writeOutputToFile(bdyvis, outfile)
	/*end := time.Now()
	totaltime := end.Sub(pstart)
	fmt.Println("total time:", totaltime)*/
}

func processTADLists(tadfilelist []string, res *int, gammaopt bool, medtadlen float64) [][][]int {

	tadlists := make([][][]int, 2)
	chrlength := 0
	var gamma float64
	for i := 0; i < 2; i++ {
		if gammaopt == true {
			tadlists[i], gamma = hicutil.ChooseGamma(medtadlen, tadfilelist[i], *res)
			_ = gamma
		} else {
			tadlists[i] = hicutil.ReadTADFile(tadfilelist[i], *res)
		}
		n := tadlists[i][len(tadlists[i])-1][1]
		if chrlength < n+1 {
			chrlength = n + 1
		}
	}
	for i := 0; i < 2; i++ {
		tadlists[i] = hicutil.FillinTADList(tadlists[i], chrlength)
	}

	return tadlists

}

type condentropy struct {
	vi     float64
	condh1 float64
	condh2 float64
}

func calcNaiveMetrics(tadlists [][][]int) []bdyvi {
	var bdyvilist []bdyvi
	var newbdyvi bdyvi

	var max_hamming float64 = 0

	for _, tadlist := range tadlists {
		for i, tadstart := range tadlist {
			for _, tadend := range tadlist[i:] {

				//if hicutil.ContainsHangingTAD(tadlists[0], tadstart[0], tadend[1]) || hicutil.ContainsHangingTAD(tadlists[1], tadstart[0], tadend[1]) {
				//	continue
				//}

				newbdyvi.start = tadstart[0]
				newbdyvi.end = tadend[1]
				//n1 := extendn(tadlists[0],tadend[1])
				//n2 := extendn(tadlists[1],tadend[1])

				var setACount, setBCount, unionCount, intersectionCount float64 = calcCardinality(tadlists, tadstart[0], tadend[1])

				newbdyvi.jacc = calcJaccordSimilarity(unionCount, intersectionCount)
				newbdyvi.hamming = calcHammingDistance(setACount, setBCount, unionCount, intersectionCount)
				newbdyvi.dice = calcTversky(setACount, setBCount, unionCount, intersectionCount)
				newbdyvi.overlap = calcOverlap(setACount, setBCount, intersectionCount)

				if newbdyvi.hamming > max_hamming {
					max_hamming = newbdyvi.hamming
				}
				//fmt.Println(newbdyvi.jacc,newbdyvi.hamming)
				newbdyvi.centroidSquareDistance = calcCentroidDistance(tadlists, tadstart[0], tadend[1], squareDistance)
				newbdyvi.centroidAbsDistance = calcCentroidDistance(tadlists, tadstart[0], tadend[1], absDistance)
				newbdyvi.vi = calcViNaive(tadlists, tadstart[0], tadend[1], false)
				newbdyvi.vi_without_boundary = calcViNaive(tadlists, tadstart[0], tadend[1], true)
				bdyvilist = append(bdyvilist, newbdyvi)
			}
		}
	}
	//normalizing hamming distance
	/*for i := range bdyvilist {
		//bdyvilist[i].hamming = (max_hamming - bdyvilist[i].hamming)
		fmt.Println(bdyvilist[i].hamming)
	}*/
	/*for _, vals := range bdyvilist {
		fmt.Println(vals.jacc,vals.hamming)
	}*/
	return bdyvilist

}

func squareDistance(a float64, b float64) float64 {
	return float64((a - b) * (a - b))
}

func absDistance(a float64, b float64) float64 {
	return math.Abs(float64(a - b))
}

/**
Finds the median of the intervals in both tad lists and then
tries to assign the nearest centroids of interval1 (interval from
tadset1) to interval2 (interval from tadset1).

The dissimilarity is measure of how much distance does the centroid
of tadset2 has to move in order to align with tadset1
*/
func calcCentroidDistance(tadlists [][][]int, tadstart int, tadend int, distance fn) float64 {
	intvl1 := hicutil.ProcessIntervals(tadlists[0], tadstart, tadend)
	intvl2 := hicutil.ProcessIntervals(tadlists[1], tadstart, tadend)

	// Calculate median of interval1 (tadset1)
	var medians []float64
	for i := range intvl1 {
		median := float64(intvl1[i][1]-intvl1[i][0]) / 2.0
		medians = append(medians, median)
	}
	// Find the nearest centroid for each interval in tadset2
	totalError := 0.0
	for i := range intvl2 {
		median := float64(intvl2[i][1]-intvl2[i][0]) / 2.0
		minDistance := math.MaxFloat64
		for j := range medians {
			// Distance to move centroid of tadset2 to align with tadset1
			distance := distance(medians[j], median)
			if distance < minDistance {
				minDistance = distance
			}
		}
		totalError += minDistance
	}

	// Calculate Root mean squared error
	rootMeanSquaredError := math.Sqrt(totalError / float64(len(intvl2)))
	// Normalize to get value between 0 to 1
	normalizedMeanSquaredError := rootMeanSquaredError / float64(tadend-tadstart)

	return normalizedMeanSquaredError
}

func calcViNaive(tadlists [][][]int, tadstart int, tadend int, exclude_boundary bool) float64 {
	var newbdyvi bdyvi
	var intvl1 [][]int
	var intvl2 [][]int
	var overlaps [][]int

	if !exclude_boundary {
		intvl1 = hicutil.ProcessIntervals(tadlists[0], tadstart, tadend)
		intvl2 = hicutil.ProcessIntervals(tadlists[1], tadstart, tadend)
		overlaps = hicutil.CalcOverlaps(intvl1, intvl2)
	} else {
		intvl1 = hicutil.ProcessIntervalsWithoutBoundary(tadlists[0], tadstart, tadend)
		intvl2 = hicutil.ProcessIntervalsWithoutBoundary(tadlists[1], tadstart, tadend)
		overlaps = hicutil.CalcOverlaps(intvl1, intvl2)
	}

	n := tadend - tadstart + 1

	clus1sizes := make([]int, len(intvl1))
	for c, clus := range intvl1 {
		clus1sizes[c] = clus[1] - clus[0] + 1
	}
	//newbdyvi.jacc = hicutil.CalcJaccSimilarity(overlaps, clus1sizes, n)

	clus2sizes := make([]int, len(intvl2))
	for c, clus := range intvl2 {
		clus2sizes[c] = clus[1] - clus[0] + 1
	}
	condh1 := hicutil.CalcCondEntropy(transpose(overlaps), clus1sizes, n)
	condh2 := hicutil.CalcCondEntropy(overlaps, clus2sizes, n)
	newbdyvi.vi = (condh1 + condh2) / math.Log(float64(n))

	return newbdyvi.vi

}
func calcCardinality(tadLists [][][]int, start int, end int) (float64, float64, float64, float64) {
	tadList1 := tadLists[0]
	tadList2 := tadLists[1]

	boundarySet := make(map[int]bool)
	boundaryIntersectionSet := make(map[int]bool)

	var setACount, setBCount float64 = 0, 0

	for _, tad := range tadList1 {
		if tad[0] >= start && tad[0] <= end && !boundarySet[tad[0]] {
			boundarySet[tad[0]] = true
			setACount += 1
		}
		if tad[1] >= start && tad[1] <= end && !boundarySet[tad[1]] {
			boundarySet[tad[1]] = true
			setACount += 1
		}

	}

	for _, tad := range tadList2 {
		if boundarySet[tad[0]] == true {
			boundaryIntersectionSet[tad[0]] = true
			setBCount += 1
		} else if tad[0] >= start && tad[0] <= end && !boundarySet[tad[0]] {
			boundarySet[tad[0]] = true
			setBCount += 1
		}

		if tad[0] != tad[1] && boundarySet[tad[1]] == true {
			boundaryIntersectionSet[tad[1]] = true
			setBCount += 1
		} else if tad[1] >= start && tad[1] <= end && !boundarySet[tad[1]] {
			boundarySet[tad[1]] = true
			setBCount += 1
		}
	}

	var unionCount float64 = float64(len(boundarySet))
	var intersectionCount float64 = float64(len(boundaryIntersectionSet))
	//fmt.Println(setACount,setBCount,unionCount,intersectionCount)
	return setACount, setBCount, unionCount, intersectionCount
}
func calcJaccordSimilarity(unionCount float64, intersectionCount float64) float64 {
	//var hammingDistance float64 = math.Abs((setACount - intersectionCount) + (setBCount-intersectionCount)) / unionCount
	//fmt.Println(hammingDistance,1.0 - float64(len(boundaryIntersectionSet)) / float64(len(boundarySet)) )
	return 1.0 - (intersectionCount / unionCount)
}

func calcOverlap(setACount float64, setBCount float64, intersectionCount float64) float64 {
	return 1.0 - intersectionCount/math.Min(setACount, setBCount)
}

func calcHammingDistance(setACount float64, setBCount float64, unionCount float64, intersectionCount float64) float64 {
	//fmt.Println(setACount,setBCount,unionCount,intersectionCount)
	var hammingDistance float64 = math.Abs((setACount-intersectionCount)+(setBCount-intersectionCount)) / unionCount
	//fmt.Println(hammingDistance,1.0 - float64(len(boundaryIntersectionSet)) / float64(len(boundarySet)) )
	//fmt.Println(hammingDistance)
	return hammingDistance
}

func calcTversky(setACount float64, setBCount float64, unionCount float64, intersectionCount float64) float64 {
	var setADiffSetB float64 = setACount - intersectionCount
	var setBDiffsetA float64 = setBCount - intersectionCount
	var alpha float64 = 0.5
	var beta float64 = 0.5
	/*var a,b float64 = 0,0
	if setADiffSetB > setBDiffsetA {
		b = setADiffSetB
		a = setBDiffsetA
	} else {
		a = setADiffSetB
		b = setBDiffsetA
	}
	var temp float64 = beta*( (a*alpha) + ((1-alpha)*b))

	var index float64 = intersectionCount / (intersectionCount + temp)*/

	var index2 float64 = intersectionCount / (intersectionCount + alpha*setADiffSetB + beta*setBDiffsetA)
	return 1.0 - index2
}

/*func scaleLastOverlap(overlaps [][]int, maxn int, lasttad []int) [][]int {

	scaleby := float64(maxn - lasttad[1] + 1)/float64(maxn - lasttad[0] + 1)
	overlaps[len(overlaps)-1][len(overlaps[0])-1] = int(float64(overlaps[len(overlaps)-1][len(overlaps[0])-1]) * scaleby)
	return overlaps
}*/

func calcVIatBdys(tadlists [][][]int) []bdyvi {

	var bdyvilist []bdyvi
	dpmap := make(map[int]map[int]condentropy) // keys should be [start][end]
	var currhvals condentropy
	// first initialize VI values of all single-TAD intervals
	for _, tadlist := range tadlists {
		for _, tads := range tadlist {
			var newbdyvi bdyvi
			newbdyvi.start = tads[0]
			newbdyvi.end = tads[1]
			//if tads[0] == tads[1] { continue }
			n := tads[1] - tads[0] + 1
			//newend := extendn(tadlists[int(math.Abs(float64(i-1)))], tads[1])
			//n := newend - tads[0] + 1
			intvl1 := hicutil.ProcessIntervals(tadlists[0], tads[0], tads[1])
			intvl2 := hicutil.ProcessIntervals(tadlists[1], tads[0], tads[1])
			// need cluster sizes and overlaps for conditional entropy calculations
			overlaps := hicutil.CalcOverlaps(intvl1, intvl2)
			clus1sizes := make([]int, len(intvl1))
			for c, clus := range intvl1 {
				clus1sizes[c] = clus[1] - clus[0] + 1
			}
			clus2sizes := make([]int, len(intvl2))
			for c, clus := range intvl2 {
				clus2sizes[c] = clus[1] - clus[0] + 1
			}
			currhvals.condh1 = hicutil.CalcCondEntropy(transpose(overlaps), clus1sizes, n)
			currhvals.condh2 = hicutil.CalcCondEntropy(overlaps, clus2sizes, n)
			currhvals.vi = currhvals.condh1 + currhvals.condh2
			if _, ok := dpmap[tads[0]]; !ok {
				dpmap[tads[0]] = make(map[int]condentropy)
			}
			dpmap[tads[0]][tads[1]] = currhvals
			if math.IsNaN(currhvals.vi) {
				fmt.Println("VI value of NaN in initialization")
				fmt.Println(currhvals)
				fmt.Println(tads)
				fmt.Println(intvl1)
				fmt.Println(tadlists[1])
				os.Exit(1)
			}
			newbdyvi.vi = currhvals.vi / math.Log(float64(tads[1]-tads[0]+1)) // divide by log(n) to normalize
			if tads[0] == tads[1] {
				continue
			}
			bdyvilist = append(bdyvilist, newbdyvi)
		}
	}
	for i, tadlist := range tadlists {
		var noti int
		if i == 0 {
			noti = 1
		} else {
			noti = 0
		}
		for numtads := 1; numtads < len(tadlist); numtads++ {
			for tadidx, starttad := range tadlist[:len(tadlist)-numtads] {
				endtad := tadlist[tadidx+numtads]
				n := endtad[1] - starttad[0] + 1
				prevn := endtad[0] - starttad[0]
				endtadlen := endtad[1] - endtad[0] + 1
				scale1 := float64(prevn) / float64(n)
				scale2 := float64(endtadlen) / float64(n)
				var newhvals condentropy
				var newbdyvis bdyvi
				newbdyvis.start = starttad[0]
				newbdyvis.end = endtad[1]
				if contains(tadlists[noti], endtad[0]-1) {
					// there are no overlapping TADs between the previous group and the new TAD so can just rescale and add the VIs/conditional entropies
					newhvals.vi = dpmap[starttad[0]][endtad[0]-1].vi*scale1 + dpmap[endtad[0]][endtad[1]].vi*scale2
					newhvals.condh1 = dpmap[starttad[0]][endtad[0]-1].condh1*scale1 + dpmap[endtad[0]][endtad[1]].condh1*scale2
					newhvals.condh2 = dpmap[starttad[0]][endtad[0]-1].condh2*scale1 + dpmap[endtad[0]][endtad[1]].condh2*scale2
				} else { // the other TAD set contains a TAD that crosses the boundary between endtad-1 and endtad, so we have to correct one of the conditional entropy terms
					if i == 0 {
						newhvals.condh1 = scale1*dpmap[starttad[0]][endtad[0]-1].condh1 + scale2*dpmap[endtad[0]][endtad[1]].condh1
						newhvals.condh2 = CorrectTADoverlap(tadlists, i, starttad[0], endtad, dpmap)
					} else {
						newhvals.condh2 = scale1*dpmap[starttad[0]][endtad[0]-1].condh2 + scale2*dpmap[endtad[0]][endtad[1]].condh2
						newhvals.condh1 = CorrectTADoverlap(tadlists, i, starttad[0], endtad, dpmap)
					}
					newhvals.vi = newhvals.condh1 + newhvals.condh2
				}
				if math.IsNaN(newhvals.vi) {
					fmt.Println("VI value of NaN in DP")
					fmt.Println(newhvals)
					fmt.Println(starttad, endtad)
					os.Exit(1)
				}
				dpmap[starttad[0]][endtad[1]] = newhvals
				newbdyvis.vi = newhvals.vi / math.Log(float64(endtad[1]-starttad[0]+1))
				// only add to list of intervals to consider if there is no hanging TAD
				if hicutil.ContainsHangingTAD(tadlists[0], starttad[0], endtad[1]) || hicutil.ContainsHangingTAD(tadlists[1], starttad[0], endtad[1]) {
					continue
				}
				bdyvilist = append(bdyvilist, newbdyvis)
			}
		}
	}
	return bdyvilist
}

func contains(a [][]int, x int) bool {
	for _, row := range a {
		if row[1] == x {
			return true
		}
	}
	return false
}

func transpose(a [][]int) [][]int {
	n := len(a)
	m := len(a[0])
	b := make([][]int, m)
	for j := 0; j < m; j++ {
		b[j] = make([]int, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			b[j][i] = a[i][j]
		}
	}
	return b
}

var wg sync.WaitGroup

func worker(tadlists [][][]int, job []bdyvi, result *[]bdyvi, convcond float64) {
	defer wg.Done()
	concurrentRand := rand.New(rand.NewSource(time.Now().UnixNano()))
	for i, querypt := range job {
		(*result)[i] = appendPval(tadlists, querypt, convcond, concurrentRand)
	}
}

func calcAllPvals(tadlists [][][]int, bdyvis []bdyvi, numCPU int, convcond float64) []bdyvi {

	var sigpts []bdyvi
	if len(bdyvis) < numCPU {
		numCPU = len(bdyvis)
	}
	bdyvis_pval := make([]bdyvi, len(bdyvis))
	allpvals := make([]float64, len(bdyvis))
	jobs := make([][]bdyvi, numCPU)
	results := make([][]bdyvi, numCPU)

	k := 0
	for i := 0; i < len(bdyvis); i++ {
		jobs[k] = append(jobs[k], bdyvis[i])
		k = (k + 1) % numCPU
	}

	for w, job := range jobs {
		wg.Add(1)
		results[w] = make([]bdyvi, len(job))
		go worker(tadlists, job, &results[w], convcond)
	}
	wg.Wait()

	i := 0
	for w := 0; w < numCPU; w++ {
		for _, res := range results[w] {
			bdyvis_pval[i] = res
			allpvals[i] = bdyvis_pval[i].pval
			i++
		}
	}

	bhidx := hicutil.MultHypTestBH(allpvals)
	if bhidx > -1 {
		print("Output Changed\n")
		sort.Slice(bdyvis_pval, func(i, j int) bool { return bdyvis_pval[i].pval < bdyvis_pval[j].pval })
		sigpts = bdyvis_pval[:bhidx+1]
	}
	return sigpts
}

/*func calcAllPvals(tadlists [][][]int, bdyvis []bdyvi, convcond float64) []bdyvi {

var sigpts []bdyvi
bdyvis_pval := make([]bdyvi, len(bdyvis))
allpvals := make([]float64,len(bdyvis_pval))
for i,querypt := range bdyvis {
	bdyvis_pval[i] = appendPval(tadlists, querypt, convcond)
	allpvals[i] = bdyvis_pval[i].pval
}
/*for _,query := range bdyvis_pval {
	if query.pval < 0.05 {
		sigpts = append(sigpts, query)
	}
}*/
/*bhidx := hicutil.MultHypTestBH(allpvals)
sort.Slice(bdyvis_pval, func(i,j int) bool {return bdyvis_pval[i].pval < bdyvis_pval[j].pval})
sigpts = bdyvis_pval[:bhidx+1]
return sigpts
}*/

func appendPval(tadlists [][][]int, querypt bdyvi, convcond float64, r *rand.Rand) bdyvi {

	intvl1 := hicutil.ProcessIntervals(tadlists[0], querypt.start, querypt.end)
	intvl2 := hicutil.ProcessIntervals(tadlists[1], querypt.start, querypt.end)
	n := querypt.end - querypt.start + 1
	p := hicutil.CalcPval(intvl1, intvl2, n, querypt.vi, convcond, r)
	if p < 0.05 && (len(intvl1) == 1 || len(intvl2) == 1) {
		fmt.Println(intvl1)
		fmt.Println(intvl2)
		fmt.Println(p)
		fmt.Println(n, querypt.vi)
		fmt.Println(querypt.start, querypt.end)
		//os.Exit(1)
	}
	querypt.pval = p
	return querypt
}

func findDomPts(sigpts []bdyvi) []bdyvi {

	//fmt.Println(sigpts)
	var dompts []bdyvi
	for i, querypt := range sigpts {
		isdom := true
		for j, comppt := range sigpts {
			if i == j {
				continue
			}
			/*if comppt.start <= querypt.start && comppt.end >= querypt.end {
				//fmt.Println("removing an interval:")
				//fmt.Println(querypt)
				//fmt.Println(comppt)
				isdom = false
				break
			}*/
			if comppt.start == querypt.start && comppt.end > querypt.end && comppt.vi < querypt.vi { // another interval starts at the same place but is strictly larger and has lower VI
				isdom = false
				break
			}
			if comppt.start < querypt.start && comppt.end == querypt.end && comppt.vi < querypt.vi { // another interval ends at the same place but is strictly larger and has lower VI
				isdom = false
				break
			}
			if comppt.start <= querypt.start {
				continue
			}
			if comppt.end >= querypt.end {
				continue
			}
			if comppt.vi < querypt.vi {
				isdom = false
				break
			}
		}
		if isdom == true {
			dompts = append(dompts, querypt)
		}
	}
	//if there are multiple dominating, significant points that start or end at the same place, remove the smaller interval
	var toremove []int
	for i, dompt1 := range dompts {
		for j, dompt2 := range dompts {
			if i == j {
				continue
			}
			if dompt1.start == dompt2.start && dompt1.end < dompt2.end {
				toremove = append(toremove, i)
				break
			}
			if dompt1.start > dompt2.start && dompt1.end == dompt2.end {
				toremove = append(toremove, i)
				break
			}
			if dompt1.start == dompt2.start && dompt1.end == dompt2.end && i < j {
				toremove = append(toremove, i)
				break
			}
		}
	}
	sort.Ints(toremove)
	for i, j := range toremove {
		dompts = append(dompts[:j-i], dompts[j-i+1:]...)
	}
	sort.Slice(dompts, func(i, j int) bool { return dompts[i].start < dompts[j].start })
	//fmt.Println(dompts)
	return dompts
}

func writeOutputToFile(domsigpts []bdyvi, outfile *string) {

	//write values to file
	f, err := os.Create(*outfile)
	if err != nil {
		panic(err)
	}
	//defer f.Close()

	w := bufio.NewWriter(f)
	labelline := []string{"start", "end", "VI", "Jaccord", "Overlap",
		"Dice", "CentroidSquareDistance", "CentroidAbsoluteDistance", "VI_boundaryless", "p-value"}
	//fmt.Println(strings.Join(labelline, "\t"))
	line1 := strings.Join(labelline, ",")
	//fmt.Println(line1)
	fmt.Fprintf(w, line1+"\n")

	for _, vals := range domsigpts {
		strvals := make([]string, 10)
		strvals[0] = strconv.Itoa(vals.start)
		strvals[1] = strconv.Itoa(vals.end)
		strvals[2] = strconv.FormatFloat(vals.vi, 'g', -1, 64)
		strvals[3] = strconv.FormatFloat(vals.jacc, 'g', -1, 64)
		strvals[4] = strconv.FormatFloat(vals.overlap, 'g', -1, 64)
		strvals[5] = strconv.FormatFloat(vals.dice, 'g', -1, 64)
		strvals[6] = strconv.FormatFloat(vals.centroidSquareDistance, 'g', -1, 64)
		strvals[7] = strconv.FormatFloat(vals.centroidAbsDistance, 'g', -1, 64)
		strvals[8] = strconv.FormatFloat(vals.vi_without_boundary, 'g', -1, 64)
		strvals[9] = strconv.FormatFloat(vals.pval, 'g', -1, 64)
		newline := strings.Join(strvals, ",")
		fmt.Fprintf(w, newline+"\n")
	}
	w.Flush()
	f.Close()
	fmt.Println("Wrote output values to", *outfile)
}
