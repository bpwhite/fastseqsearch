package main

import (
	"crypto/rand"
	"fmt"
	"math"
	"math/big"
	//"time"
	"os"
	//"reflect"
	"sort"

	"github.com/gonum/stat"
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

	init_params1 := make(map[string]float64)
	init_params2 := make(map[string]float64)

	init_dims1 := make(map[string]int)
	init_dims2 := make(map[string]int)

	for i := 0; i < len(a1_letters); i++ {
		init_params1[a1_letters[i]] = 1.0
		init_dims1[a1_letters[i]] = 0
	}
	for i := 0; i < len(a2_letters); i++ {
		init_params2[a2_letters[i]] = 1.0
		init_dims2[a2_letters[i]] = 0
	}

	fmt.Println(init_params1, init_dims1)
	fmt.Println(init_params2, init_dims2)

	generations := 20
	population_size := 10

	//mut_strength := 5
	//mut_rate := 15

	// search replicates
	//B := 1
	// number of sample hits to check
	sub_sample_size := 1
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

	var population = make([]string, population_size)

	for gen := 0; gen < generations; gen++ {

		dims1 := make(map[string]int)
		dims2 := make(map[string]int)

		params1 := make(map[string]float64)
		params2 := make(map[string]float64)

		if gen == 0 {
			dims1 = init_dims1
			dims2 = init_dims2
			params1 = init_params1
			params2 = init_params2
		}

		var correlations = make([]float64, generations)
		var pdists = make([]float64, sub_sample_size)
		var eucdists = make([]float64, sub_sample_size)
		var weights = make([]float64, sub_sample_size)

		fmt.Println("Dims1: ", dims1, " Dims2: ", " Params1: ", params1, " Params2: ", params2)

		// sample data
		for sample := 0; sample < sub_sample_size; sample++ {

			// hold current position in c-space
			var s1_pt = make([]float64, len(dims1))
			var s2_pt = make([]float64, len(dims2))

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

			matches := 0
			for i := 0; i < subn; i++ {
				//b[i] =

				xor_compare := int(substr1[i] ^ substr2[i])

				if xor_compare != 0 {
					matches++
				}

				//fmt.Println(string(substr1[i]))
				s1_pt[dims1[string(substr1[i])]] += params1[string(substr1[i])]
				s2_pt[dims2[string(substr2[i])]] += params2[string(substr2[i])]

				fmt.Println(s1_pt, " => ", s2_pt)

			}
			os.Exit(0)
			pdist := 1.0 - float64(matches)/float64(subn)

			//fmt.Printf("s1: %v, s2: %v\n", s1_pt, s2_pt)

			//radius1 := math.Sqrt(s1_pt[0]*s1_pt[0] + s1_pt[1]*s1_pt[1] + s1_pt[2]*s1_pt[2])
			//radius2 := math.Sqrt(s2_pt[0]*s2_pt[0] + s2_pt[1]*s2_pt[1] + s2_pt[2]*s2_pt[2])
			//radist := math.Abs(radius1 - radius2)

			eucdist := math.Sqrt(math.Pow(s1_pt[0]-s2_pt[0], 2) + math.Pow(s1_pt[1]-s2_pt[1], 2) + math.Pow(s1_pt[2]-s2_pt[2], 2))
			if print_dists == 1 {
				fmt.Printf("%d => %.2f, => %v => %v\n", matches, pdist, eucdist)
			}

			pdists[sample] = pdist
			eucdists[sample] = eucdist
			weights[sample] = 1.0

		}

		c := stat.Correlation(pdists, eucdists, weights)
		//fmt.Printf("Correlation is %.5f\n", c)
		correlations[gen] = c
		//final_cor = c

		// rep, a, t, g, c, dna_6, dna_23, dna_4, dna_19, dna_2, dna_21, correlation
		//outp_string := fmt.Sprintf("%d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.5f\n", gen, dna_a, dna_t, dna_g, dna_c, dna_6, dna_23, dna_4, dna_19, dna_2, dna_21, c)
		//fmt.Print(outp_string)

		//outp.WriteString(outp_string)
	}
	//population[pop] = fmt.Sprintf("%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.5f\n", dna_a, dna_t, dna_g, dna_c, dna_6, dna_23, dna_4, dna_19, dna_2, dna_21, final_cor)
	//outp.Sync()

	fmt.Println(population)
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
	//fmt.Println(found)
	return found
}
