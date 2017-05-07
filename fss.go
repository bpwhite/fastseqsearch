package main

import (
	"fmt"
)

func main() {
	string1 := "ATGCCGGGTTGGTGCGATTCAGGTGGACACGACGTTATAC"
	string2 := "ATGCCGGGAAGGCGCGGTTCATGTGGACACGATGTAATAC"

	n := len(string1)
	//b := make([]byte, n)
	matches := 0
	for i := 0; i < n; i++ {
		//b[i] =

		if int(string1[i]^string2[i]) == 0 {
			matches++
		}
		//fmt.Printf("%v\n", int(b[i]))

	}
	fmt.Printf("%d\n", matches)
}
