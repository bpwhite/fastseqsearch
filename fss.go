package main

import (
	"crypto/rand"
	"fmt"
	"math"
	"math/big"
	//"time"
	"github.com/gonum/stat"
	//"os"
)


func main() {
	// CYTC005-12|Homo sapiens|COI-5P|HM771214
	string1 := "ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATATTGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATTGGCTTCCTGGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCATATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACTCCACGGAAGCAATAT"
	// CYTC1000-12|Pan troglodytes schweinfurthii|COI-5P|JF727192
	string2 := "CACAAAGATATTGGAACACTATACCTACTATTCGGCGCATGGGCTGGAGTCCTGGGCACAGCCCTAAGTCTCCTTATTCGGGCTGAACTAGGCCAACCAGGCAACCTTCTAGGTAATGACCACATCTACAATGTCATCGTCACAGCCCATGCATTCGTAATAATCTTCTTCATAGTAATGCCTATCATAATCGGAGGCTTTGGCAACTGGCTAGTCCCCTTGATAATTGGTGCCCCCGACATGGCATTCCCCCGCATAAACAACATAAGCTTCTGACTCCTACCCCCTTCTCTCCTACTTCTACTTGCATCTGCCATAGTAGAAGCCGGCGCCGGAACAGGTTGAACGGTCTACCCTCCCTTAGCGGGAAACTACTCGCATCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTGCATCTGGCAGGCGTCTCCTCTATCCTAGGAGCCATTAACTTCATCACAACAATTATTAATATAAAACCTCCTGCCATAACCCAATACCAAACACCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCTTACTTCTCCTATCCCTCCCAGTCCTAGCTGCTGGCATCACCATACTATTGACAGATCGTAACCTCAACACTACCTTCTTCGATCCAGCCGGGGGAGGAGACCCTATTCTATATCAGCACTTATTCTGATTTTTTGGCCACCCCGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATTTCCCACATTGTAACTTATTACTCCGGAAAAAAAGAACCATTTGGATATATAGGCATGGTTTGAGCTATAATATCAATTGGTTTCCTAGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCCTATTTCACCTCCGCTACCATAATCATTGCTATTCCTACCGGCGTCAAAGTATTCAGCTGACTCGCTACACTTCACGGAAGC"

	// G -> C = 4 transversion
	// C -> G = 4 transversion
	// G -> T = 19 transversion
	// T -> G = 19 transversion
	// C -> T = 23 transition
	// T -> C = 23 transition
	// A -> C = 2 transversion
	// C -> A = 2 transversion
	// A -> T = 21 transversion
	// T -> A = 21 transversion
	// A -> G = 6 transition
	// G -> A = 6 transition
	
	generations := 20
	population_size := 10
	
	// search replicates
	//B := 1
	// number of sample hits to check
	sub_sample_size := 500
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
	dna_a := 1.0
	dna_t := 1.0
	dna_g := 1.0
	dna_c := 1.0

	dna_6 	:= 1.0 
	dna_23 	:= 1.0
	dna_4 	:= 1.0
	dna_19 	:= 1.0
	dna_2 	:= 1.0
	dna_21 	:= 1.0
		
	n1 := len(string1)
	n2 := len(string2)
	
	var population = make([]string, population_size)
	
	last_dna_a := dna_a
	last_dna_t := dna_t
	last_dna_g := dna_g
	last_dna_c := dna_c
	
	last_dna_6 := dna_6
	last_dna_23 := dna_23
	last_dna_4 := dna_4
	last_dna_19	:= dna_19
	last_dna_2 := dna_2
	last_dna_21 := dna_21
		
	for pop := 0; pop < population_size; pop++ {
		
		final_cor := 0.0
		

		for gen := 0; gen < generations; gen++ {			

			var correlations = make([]float64, generations)
			var pdists = make([]float64, sub_sample_size)
			var eucdists = make([]float64, sub_sample_size)
			var weights = make([]float64, sub_sample_size)
			
			if gen == 0 {
				dna_a = last_dna_a
				dna_t = last_dna_t
				dna_g = last_dna_g
				dna_c = last_dna_c
				
				dna_6 = last_dna_6
				dna_23 = last_dna_23
				dna_4 = last_dna_4
				dna_19	= last_dna_19
				dna_2 = last_dna_2
				dna_21 = last_dna_21
			} else {
				dna_a 	= mutate(dna_a, 5, 15)
				dna_t 	= mutate(dna_t, 5, 15)
				dna_g 	= mutate(dna_g, 5, 15)
				dna_c 	= mutate(dna_c, 5, 15)

				dna_6	= mutate(dna_6, 5, 15)
				dna_23 	= mutate(dna_23, 5, 15)
				dna_4 	= mutate(dna_4, 5, 15)
				dna_19 	= mutate(dna_19, 5, 15)
				dna_2 	= mutate(dna_2, 5, 15)
				dna_21 	= mutate(dna_21, 5, 15)	
			}

			// sample data
			for sample := 0; sample < sub_sample_size; sample++ {

				r1 := gen_cryp_num(int64(n1))
				r2 := gen_cryp_num(int64(n2))

				end1 := r1 + int64(subn)
				end2 := r2 + int64(subn)

				if end1 >= int64(n1) {
					sample--
					continue
				}
				if end2 >= int64(n2) {
					sample--
					continue
				}

				//		os.Exit(1)

				substr1 := string1[r1:end1]
				substr2 := string2[r2:end2]

				s1_pt := [3]float64{0, 0, 0}
				s2_pt := [3]float64{0, 0, 0}

				matches := 0
				for i := 0; i < subn; i++ {
					//b[i] =

					xor_compare := int(substr1[i] ^ substr2[i])

					if xor_compare != 0 {
						matches++
					}

					if i > 0 {
						xor_prev1 := int(substr1[i] ^ substr1[i-1])
						xor_prev2 := int(substr2[i] ^ substr2[i-1])

						switch xor_prev1 {
						case 6:
							s1_pt[2] += dna_6
						case 23:
							s1_pt[2] += dna_23
						case 4:
							s1_pt[2] += dna_4
						case 19:
							s1_pt[2] += dna_19
						case 2:
							s1_pt[2] += dna_2
						case 21:
							s1_pt[2] += dna_21
						}
						switch xor_prev2 {
						case 6:
							s2_pt[2] += dna_6
						case 23:
							s2_pt[2] += dna_23
						case 4:
							s2_pt[2] += dna_4
						case 19:
							s2_pt[2] += dna_19
						case 2:
							s2_pt[2] += dna_2
						case 21:
							s2_pt[2] += dna_21
						}
					}
					//continue

					switch string(substr1[i]) {
					case "A":
						s1_pt[0] += dna_a
						s1_pt[1] += dna_a
					case "T":
						s1_pt[0] += dna_t
						s1_pt[1] += dna_t
					case "G":
						s1_pt[0] += dna_g
						s1_pt[1] += dna_g
					case "C":
						s1_pt[0] += dna_c
						s1_pt[1] += dna_c
					}

					switch string(substr2[i]) {
					case "A":
						s2_pt[0] += dna_a
						s2_pt[1] += dna_a
					case "T":
						s2_pt[0] += dna_t
						s2_pt[1] += dna_t
					case "G":
						s2_pt[0] += dna_g
						s2_pt[1] += dna_g
					case "C":
						s2_pt[0] += dna_c
						s2_pt[1] += dna_c
					}

				}
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
			final_cor = c
			
			// rep, a, t, g, c, dna_6, dna_23, dna_4, dna_19, dna_2, dna_21, correlation	
			//outp_string := fmt.Sprintf("%d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.5f\n", gen, dna_a, dna_t, dna_g, dna_c, dna_6, dna_23, dna_4, dna_19, dna_2, dna_21, c)
			//fmt.Print(outp_string)
			

			//outp.WriteString(outp_string)
		}
		population[pop] = fmt.Sprintf("%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.5f\n", dna_a, dna_t, dna_g, dna_c, dna_6, dna_23, dna_4, dna_19, dna_2, dna_21, final_cor)
		//outp.Sync()
	}
	
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
	if (gen_cryp_num(100) <= rate) {
		multiplier = float64(gen_cryp_num(strength))
		sign := 1.0
		if (gen_cryp_num(100) < 50) {
			sign = -1.0
		}
		mutated = input + float64(multiplier)*sign
	} else {
		mutated = input
	}
	
	return mutated
}
