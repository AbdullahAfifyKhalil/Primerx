/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package primerDesigning;

import java.util.Scanner;
import java.util.Vector;

/**
 *
 * @author 7asa
 */
public class Primer_DesignerV1 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Scanner input = new Scanner(System.in);

        System.out.print("Enter Your Sequence: ");
        
        String sequence;
        boolean DNA_flag = true;
        sequence = input.next();
        DeoxyRibose_NA DNA;
        for (int i = 0; i < sequence.length(); i++) {
            if (sequence.charAt(i) == 'U') {
                DNA_flag = false;
            }
        }
        if (DNA_flag == true) {
            DNA = new DeoxyRibose_NA(sequence);
        } else {
            String newSequence = null;
            for (int i = 0; i < sequence.length(); i++) {
                if (sequence.charAt(i) == 'U') {
                    newSequence = sequence.substring(0, i) + 'T' + sequence.substring(i + 1);
                }
            }
            DNA = new DeoxyRibose_NA(newSequence);
        }
        
        DNA.Sequencing();
        Primer p = new Primer();
        p.get_all_Forward_Primers(131); //index al bracket 
        System.out.println("=====================================");
        p.get_all_Reverse_Primers(190);


    }
    //GACGAGGGAAGTCCCAATCACTTATGCCCATCATAATTCCCCCCGCGCCTTCCTTTTATAGGTGCCACGACGAGCCCCATTCAGATTATGCGTTACTGGTCTCCTCCTCACCCACTTCCTCCAGATCGTCGCTCATATTGTACGTCAGGGCCACGACGACTAGGCCGGGCCTTGAACTTCGAGCTTTCCCTCGTTTCCGTTGCACCTTGCTCAACATCTGCTCCAACAATGTTGAGATCACCCTTACGGACCATTCAAGAACCGACAATCCCGCCATCGCCCTTTCTCCCCCCCCCATAATAGGAAATCAGCCGGATAGTGTCATGGGTATATCTCTTGGACTCTTCTCTCTTGTTTCTTCCGGGGTTATAAAAAAGCTTACGGGGGTTGACTTTTTCTTACTTTTGTTAATTAGTTCTTAACATCAACACAAGTAAGCTAGGAGATCTAGAGATATAATAAATTTTTAAATGTCTTTTTTTTCATGTAGTTGTCTAGTTTT
//1310 primer from 0 to 131

}
