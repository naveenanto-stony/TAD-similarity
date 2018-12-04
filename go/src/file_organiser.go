package src

import (
	"bufio"
	"os"
	)

func sortFile(filename string)  {
	file, err := os.Open(filename)
	defer file.Close()

	if err == nil {

	}
}
