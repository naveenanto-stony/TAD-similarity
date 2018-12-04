package hicutil

func CalcJaccSimilarity(overlaps [][]int, clusterSizes []int, n int) float64 {
	similarity := 0.0

	for i, _ := range clusterSizes {
		sum := 0.0
		for _, overlap := range overlaps[i] {
			sum +=  float64(overlap)/float64(n)
		}
		similarity += sum/float64(len(overlaps[i]))
	}

	return similarity/float64(len(clusterSizes))
}
