package hicutil

func OverlapInterval(tadlist [][]int, start int, end int)  []int {
	res := make([] int,2)
	res[0] = tadlist[findFloor(tadlist, start)][0]
	res[1] = tadlist[findCeil(tadlist, end)][1]
	return res
}
func findCeil(tadList [][]int, key int) int {
	left, right := 0, len(tadList)-1

	for ; left <= right ;  {
		mid := left + (right - left)/2
		if (mid == 0 || tadList[mid-1][1] < key) && tadList[mid][1] >= key  {
			return mid
		}
		if tadList[mid][0] < key {
			left = mid+1
		} else {
			right = mid-1
		}
	}

	return left
}
func findFloor(tadList [][]int, key int) int {
	left, right := 0, len(tadList)-1

	for ; left <= right ;  {
		mid := left + (right - left)/2
		if (mid == (len(tadList)-1) || tadList[mid+1][0] > key) && tadList[mid][0] <= key {
			return mid
		}
		if tadList[mid][0] < key {
			left = mid+1
		} else {
			right = mid-1
		}
	}

	return left
}

