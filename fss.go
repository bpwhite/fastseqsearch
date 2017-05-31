package main

import (
	"crypto/rand"
	"fmt"
	//"math"
	"math/big"
	//"time"
	"os"
	"time"
	//"reflect"
	"sort"

	"github.com/gonum/stat"
)

func main() {
	// 16S rRNA homo sapiens
	//rRNA16S := "GCCAAACCTAGCCCCAAACCCACTCCACCTTACTACCAAACAACCTTAACCAAACCATTTACCCAAATAAAGTATAGGCGATAGAAATTGTAACCTGGCGCAATAGATATAGTACCGCAAGGGAAAGATGAAAAATTATAACCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAATTAACTAGAAATAACTTTGCAAGGAGAACCAAAGCTAAGACCCCCGAAACCAGACGAGCTACCTAAGAACAGCTAAAAGAGCACACCCGTCTATGTAGCAAAATAGTGGGAAGATTTATAGGTAGAGGCGACAAACCTACCGAGCCTGGTGATAGCTGGTTGTCCAAGATAGAATCTTAGTTCAACTTTAAATTTACCCACAGAACCCTCTAAATCCCCTTGTAAATTTAACTGTTAGTCCAAAGAGGAACAGCTCTTTGGACACTAGGAAAAAACCTTGTAGAGAGAGTAAAAAATTTAACACCCATAGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTTCAAGCTCAACACCCACTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATCACCTTATAGAAGAACTAATGTTAGTATAAGTAACATGAAAACATTCTCCTCCGCATAAGCCTGCGTCAGATTAAAACACTGAACTGACAATTAACAGCCCAATATCTACAATCAACCAACAGGCCATTATTACCCTCACTGTCAACCCAACACAGGCATGCTCATAAGGAAAGGTTAAAAAAAGTAAAAGGAACTCGGCAAATCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATTACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCATAATCACTTGTTCCTTAAATAGGGACCTGTATGAATGGCTCCACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCATAACACAGCAAGACGAGAAGACCCTATGGAGCTTTAATTTATTAATGCAAACAATACCTAACAAACCCACAGGTCCTAAACTACCAAACCTGCATTAAAAATTTCGGTTGGGGCGACCTCGGAGCACAACCCAACCTCCGAGCAGTACATGCTAAGACTTCACCAGTCAAAGCGAACTACCATACTCAATTGATCCAATAACTTGACCAACGGAACAAGTTACCCTAGGGATAACAGCGCAATCCTATTCCAGAGTCCATATCAACAATAGGGTTTACGACCTCGATGTTGGATCAGGACATCCCGATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAAAGTCCTACGTGATCTGAGTTCAGACCGGAGTAATCCAGGTCGGTTTCTATCTACTTCAAATTCCTCCCTGTACGAAAGGACAAGAGAAATAAGGCCTACTTCACAAAGCGCCTTCCCCCGTAAATGATATCATCTCAACTTAGTATTATACCCACACCCACCCAAGAACAGGGTTT"

	// CYTC005-12|Homo sapiens|COI-5P|HM771214
	COI5p := "ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATATTGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATTGGCTTCCTGGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCATATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACTCCACGGAAGCAATAT"

	// HBA1 homo sapiens
	//HBA1 := "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAGGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAA"

	string1 := COI5p
	string2 := COI5p

	// CYTC1000-12|Pan troglodytes schweinfurthii|COI-5P|JF727192
	//string2 := "CACAAAGATATTGGAACUUUACTATACCTACTATTCGGCGCATGGGCTGGAGTCCTGGGCACAGCCCTAAGTCTCCTTATTCGGGCTGAACTAGGCCAACCAGGCAACCTTCTAGGTAATGACCACATCTACAATGTCATCGTCACAGCCCATGCATTCGTAATAATCTTCTTCATAGTAATGCCTATCATAATCGGAGGCTTTGGCAACTGGCTAGTCCCCTTGATAATTGGTGCCCCCGACATGGCATTCCCCCGCATAAACAACATAAGCTTCTGACTCCTACCCCCTTCTCTCCTACTTCTACTTGCATCTGCCATAGTAGAAGCCGGCGCCGGAACAGGTTGAACGGTCTACCCTCCCTTAGCGGGAAACTACTCGCATCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTGCATCTGGCAGGCGTCTCCTCTATCCTAGGAGCCATTAACTTCATCACAACAATTATTAATATAAAACCTCCTGCCATAACCCAATACCAAACACCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCTTACTTCTCCTATCCCTCCCAGTCCTAGCTGCTGGCATCACCATACTATTGACAGATCGTAACCTCAACACTACCTTCTTCGATCCAGCCGGGGGAGGAGACCCTATTCTATATCAGCACTTATTCTGATTTTTTGGCCACCCCGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATTTCCCACATTGTAACTTATTACTCCGGAAAAAAAGAACCATTTGGATATATAGGCATGGTTTGAGCTATAATATCAATTGGTTTCCTAGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCCTATTTCACCTCCGCTACCATAATCATTGCTATTCCTACCGGCGTCAAAGTATTCAGCTGACTCGCTACACTTCACGGAAGC"

	//RemoveDuplicates(string1)

	a1_size := count_letters(string1)
	a1_letters := detect_letters(string1)
	fmt.Println("Alphabet 1 base size: ", a1_size, " Letters: ", a1_letters)
	a2_size := count_letters(string2)
	a2_letters := detect_letters(string2)
	fmt.Println("Alphabet 2 base size: ", a2_size, " Letters: ", a2_letters)

	//mut_strength := 5
	//mut_rate := 15

	// search replicates
	B := 50
	// number of sample hits to check
	sub_sample_size := 500
	// sequence slice size
	subn := 50
	//vary_dnastat := 1
	print_dists := 0

	// output
	t := time.Now()
	//fmt.Println(t.Format("20060102150405"))
	time_stamp := t.Format("20060102150405")
	output_file := fmt.Sprint("output/", "seq_search_", sub_sample_size, "sub_", B, "reps_", time_stamp, ".csv")
	outp, _ := os.Create(output_file)
	// initial parameters

	n1 := len(string1)
	n2 := len(string2)

	var medians = make([]float64, B)
	var skews = make([]float64, B)
	var b_weights = make([]float64, B)

	for rep := 0; rep < B; rep++ {

		var pdists = make([]float64, sub_sample_size)
		var rev_pdists = make([]float64, sub_sample_size)
		var h_comps = make([]float64, sub_sample_size)
		var weights = make([]float64, sub_sample_size)
		b_weights[rep] = 1.0

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

			//fmt.Println("Forward: ", substr1)
			//rev_substr1 := reverseComplement(substr1)
			//fmt.Println("Reverse: ", rev_substr1)
			rev_substr2 := reverseString(substr2)

			matches := countMatch(substr1, substr2)
			rev_matches := countMatch(substr1, rev_substr2)
			num_comp := countComplementary(substr1, substr2)
			rev_num_comp := countComplementary(substr1, rev_substr2)

			pdist := 1.0 - float64(matches)/float64(subn)
			rev_pdist := 1.0 - float64(rev_matches)/float64(subn)

			// select lower of 2 distances
			l_dist := pdist
			if pdist > rev_pdist {
				l_dist = rev_pdist
			}

			// select higher of 2 complementary
			h_comp := num_comp
			if rev_num_comp > num_comp {
				h_comp = rev_num_comp
			}

			if print_dists == 1 {
				fmt.Printf("%d => %.2f \t %d => %.2f | %.2f | %d\n", matches, pdist, rev_matches, rev_pdist, l_dist, h_comp)
			}
			output_dists := fmt.Sprintf("%d,%.2f,%d,%.2f,%.2f,%d\n", matches, pdist, rev_matches, rev_pdist, l_dist, h_comp)
			outp.WriteString(output_dists)
			pdists[sample] = pdist
			rev_pdists[sample] = rev_pdist

			// Complementary stats
			h_comps[sample] = h_comp
			weights[sample] = 1.0
		}
		sort.Float64s(h_comps)
		median := stat.Quantile(0.5, 1, h_comps, weights)
		skew := stat.Skew(h_comps, weights)
		medians[rep] = median
		skews[rep] = skew
	}
	outp.Sync()
	//fmt.Println(h_comps)

	//fmt.Println(medians)
	sort.Float64s(skews)
	//fmt.Println(skews)
	mean_skew := stat.Mean(skews, b_weights)
	fmt.Println("Mean skew: ", mean_skew)
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
		rs = fmt.Sprint(rs, nuc_comp[string(x)])
	}

	rs = reverseString(rs)

	return rs
}

func reverseString(s string) string {
	chars := []rune(s)
	for i, j := 0, len(chars)-1; i < j; i, j = i+1, j-1 {
		chars[i], chars[j] = chars[j], chars[i]
	}
	return string(chars)
}

func countMatch(substr1 string, substr2 string) (matches int) {
	matches = 0
	len1 := len(substr1)
	len2 := len(substr2)
	if len1 != len2 {
		fmt.Println("Substring length mismatch.")
	}
	for i := 0; i < len1; i++ {

		xor_compare := int(substr1[i] ^ substr2[i])

		if xor_compare != 0 {
			matches++
		}

	}
	return matches
}

func countComplementary(substr1 string, substr2 string) (complements float64) {
	complements = 0.0

	nuc_comp := make(map[string]string)
	nuc_comp["A"] = "T"
	nuc_comp["T"] = "A"
	nuc_comp["G"] = "C"
	nuc_comp["C"] = "G"

	len1 := len(substr1)
	len2 := len(substr2)
	if len1 != len2 {
		fmt.Println("Substring length mismatch.")
	}
	for i := 0; i < len1; i++ {

		if nuc_comp[string(substr1[i])] == string(substr2[i]) {
			complements++
		}

	}
	return complements
}
