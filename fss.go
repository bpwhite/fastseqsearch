package main

import (
	"crypto/rand"
	"fmt"
	//"math"
	"math/big"
	//"time"
	"os"
	//"reflect"
	"sort"

	//"github.com/gonum/stat"
)

func main() {
	// CYTC005-12|Homo sapiens|COI-5P|HM771214
	string1 := "ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATATTGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATTGGCTTCCTGGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCATATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACTCCACGGAAGCAATAT"
	string2 := "ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATATTGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATTGGCTTCCTGGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCATATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACTCCACGGAAGCAATAT"

	// CYTC1000-12|Pan troglodytes schweinfurthii|COI-5P|JF727192
	//string2 := "CACAAAGATATTGGAACUUUACTATACCTACTATTCGGCGCATGGGCTGGAGTCCTGGGCACAGCCCTAAGTCTCCTTATTCGGGCTGAACTAGGCCAACCAGGCAACCTTCTAGGTAATGACCACATCTACAATGTCATCGTCACAGCCCATGCATTCGTAATAATCTTCTTCATAGTAATGCCTATCATAATCGGAGGCTTTGGCAACTGGCTAGTCCCCTTGATAATTGGTGCCCCCGACATGGCATTCCCCCGCATAAACAACATAAGCTTCTGACTCCTACCCCCTTCTCTCCTACTTCTACTTGCATCTGCCATAGTAGAAGCCGGCGCCGGAACAGGTTGAACGGTCTACCCTCCCTTAGCGGGAAACTACTCGCATCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTGCATCTGGCAGGCGTCTCCTCTATCCTAGGAGCCATTAACTTCATCACAACAATTATTAATATAAAACCTCCTGCCATAACCCAATACCAAACACCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCTTACTTCTCCTATCCCTCCCAGTCCTAGCTGCTGGCATCACCATACTATTGACAGATCGTAACCTCAACACTACCTTCTTCGATCCAGCCGGGGGAGGAGACCCTATTCTATATCAGCACTTATTCTGATTTTTTGGCCACCCCGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATTTCCCACATTGTAACTTATTACTCCGGAAAAAAAGAACCATTTGGATATATAGGCATGGTTTGAGCTATAATATCAATTGGTTTCCTAGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCCTATTTCACCTCCGCTACCATAATCATTGCTATTCCTACCGGCGTCAAAGTATTCAGCTGACTCGCTACACTTCACGGAAGC"

	//RemoveDuplicates(string1)

	a1_size := count_letters(string1)
	a1_letters := detect_letters(string1)
	fmt.Println("Alphabet 1 base size: ", a1_size, " Letters: ", a1_letters)
	a2_size := count_letters(string2)
	a2_letters := detect_letters(string2)
	fmt.Println("Alphabet 2 base size: ", a2_size, " Letters: ", a2_letters)

	generations := 20

	//mut_strength := 5
	//mut_rate := 15

	// search replicates
	//B := 1
	// number of sample hits to check
	sub_sample_size := 50
	// sequence slice size
	subn := 100
	//vary_dnastat := 1
	print_dists := 0

	// output
	//t := time.Now()
	//fmt.Println(t.Format("20060102150405"))
	//time_stamp := t.Format("20060102150405")
	//output_string := fmt.Sprint("output/", "seq_search_",B,"reps_", time_stamp, ".csv")
	//outp, _ := os.Create(output_string)
	// initial parameters

	n1 := len(string1)
	n2 := len(string2)

	for gen := 0; gen < generations; gen++ {

		var pdists = make([]float64, sub_sample_size)

		// sample data
		for sample := 0; sample < sub_sample_size; sample++ {

			// random substring size
			r1 := gen_cryp_num(int64(n1))
			r2 := gen_cryp_num(int64(n2))

			// end of substring
			end1 := r1 + int64(subn)
			end2 := r2 + int64(subn)

			// make sure substring is in range of string
			if end1 >= int64(n1) {
				sample--
				continue
			}
			if end2 >= int64(n2) {
				sample--
				continue
			}

			substr1 := string1[r1:end1]
			substr2 := string2[r2:end2]

			fmt.Println(substr1)
			rev_substr1 := reverseComplement(substr1)
			fmt.Println(rev_substr1)

			os.Exit(0)
			matches := 0
			for i := 0; i < subn; i++ {
				//b[i] =

				xor_compare := int(substr1[i] ^ substr2[i])

				if xor_compare != 0 {
					matches++
				}

			}
			pdist := 1.0 - float64(matches)/float64(subn)

			if print_dists == 1 {
				fmt.Printf("%d => %.2f\n", matches, pdists)
			}

			pdists[sample] = pdist

		}

		fmt.Println(pdists)
		//outp.WriteString(outp_string)
	}
	//outp.Sync()

	//fmt.Println(population)
}

func gen_cryp_num(input int64) (n int64) {
	nBig, err := rand.Int(rand.Reader, big.NewInt(input))
	if err != nil {
		panic(err)
	}
	n = nBig.Int64()
	//fmt.Printf("Here is a random %T in [0,27) : %d\n", n, n)
	return
}

func mutate(input float64, strength int64, rate int64) (mutated float64) {
	multiplier := 1.0
	if gen_cryp_num(100) <= rate {
		multiplier = float64(gen_cryp_num(strength))
		sign := 1.0
		if gen_cryp_num(100) < 50 {
			sign = -1.0
		}
		mutated = input + float64(multiplier)*sign
	} else {
		mutated = input
	}

	return mutated
}

func count_letters(input string) (num_letters int) {
	num_letters = 0
	letters := make(map[string]bool)
	letters = RemoveDuplicates(input)
	num_letters = len(letters)
	return num_letters
}

func detect_letters(input string) (letters []string) {
	letter_map := RemoveDuplicates(input)

	for key, _ := range letter_map {
		letters = append(letters, key)
	}
	sort.Strings(letters)
	return letters
}

func RemoveDuplicates(xs string) (found map[string]bool) {
	found = make(map[string]bool)

	for _, x := range xs {

		if !found[string(x)] {
			found[string(x)] = true
		}
	}
	return found
}

func reverseComplement(xs string) (rs string) {
	nuc_comp := make(map[string]string)
	nuc_comp["A"] = "T"
	nuc_comp["T"] = "A"
	nuc_comp["G"] = "C"
	nuc_comp["C"] = "G"

	for _, x := range xs {
		fmt.Println(x)
		os.Exit(0)
		//rs := fmt.Sprint(rs, nuc_comp[x])
	}

	//	rs = reverseString(rs)
	return rs
}

/*
func reverseString(s string) string {
	o := make([]int, utf8.RuneCountInString(s))
	i := len(o)
	for _, c := range s {
		i--
		o[i] = c
	}
	return string(o)
}
*/
