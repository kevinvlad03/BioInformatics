package org.example;

import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class Main {

    // Static map to hold the genetic code (RNA codon -> Amino Acid)
    private static final Map<String, String> geneticCode = new HashMap<>();

    // Static initializer block to populate the genetic code map
    static {
        // Phenylalanine (Phe)
        geneticCode.put("UUU", "Phe");
        geneticCode.put("UUC", "Phe");
        // Leucine (Leu)
        geneticCode.put("UUA", "Leu");
        geneticCode.put("UUG", "Leu");
        geneticCode.put("CUU", "Leu");
        geneticCode.put("CUC", "Leu");
        geneticCode.put("CUA", "Leu");
        geneticCode.put("CUG", "Leu");
        // Isoleucine (Ile)
        geneticCode.put("AUU", "Ile");
        geneticCode.put("AUC", "Ile");
        geneticCode.put("AUA", "Ile");
        // Methionine (Met) / START
        geneticCode.put("AUG", "Met");
        // Valine (Val)
        geneticCode.put("GUU", "Val");
        geneticCode.put("GUC", "Val");
        geneticCode.put("GUA", "Val");
        geneticCode.put("GUG", "Val");
        // Serine (Ser)
        geneticCode.put("UCU", "Ser");
        geneticCode.put("UCC", "Ser");
        geneticCode.put("UCA", "Ser");
        geneticCode.put("UCG", "Ser");
        geneticCode.put("AGU", "Ser");
        geneticCode.put("AGC", "Ser");
        // Proline (Pro)
        geneticCode.put("CCU", "Pro");
        geneticCode.put("CCC", "Pro");
        geneticCode.put("CCA", "Pro");
        geneticCode.put("CCG", "Pro");
        // Threonine (Thr)
        geneticCode.put("ACU", "Thr");
        geneticCode.put("ACC", "Thr");
        geneticCode.put("ACA", "Thr");
        geneticCode.put("ACG", "Thr");
        // Alanine (Ala)
        geneticCode.put("GCU", "Ala");
        geneticCode.put("GCC", "Ala");
        geneticCode.put("GCA", "Ala");
        geneticCode.put("GCG", "Ala");
        // Tyrosine (Tyr)
        geneticCode.put("UAU", "Tyr");
        geneticCode.put("UAC", "Tyr");
        // Histidine (His)
        geneticCode.put("CAU", "His");
        geneticCode.put("CAC", "His");
        // Glutamine (Gln)
        geneticCode.put("CAA", "Gln");
        geneticCode.put("CAG", "Gln");
        // Asparagine (Asn)
        geneticCode.put("AAU", "Asn");
        geneticCode.put("AAC", "Asn");
        // Lysine (Lys)
        geneticCode.put("AAA", "Lys");
        geneticCode.put("AAG", "Lys");
        // Aspartic Acid (Asp)
        geneticCode.put("GAU", "Asp");
        geneticCode.put("GAC", "Asp");
        // Glutamic Acid (Glu)
        geneticCode.put("GAA", "Glu");
        geneticCode.put("GAG", "Glu");
        // Cysteine (Cys)
        geneticCode.put("UGU", "Cys");
        geneticCode.put("UGC", "Cys");
        // Tryptophan (Trp)
        geneticCode.put("UGG", "Trp");
        // Arginine (Arg)
        geneticCode.put("CGU", "Arg");
        geneticCode.put("CGC", "Arg");
        geneticCode.put("CGA", "Arg");
        geneticCode.put("CGG", "Arg");
        geneticCode.put("AGA", "Arg");
        geneticCode.put("AGG", "Arg");
        // Glycine (Gly)
        geneticCode.put("GGU", "Gly");
        geneticCode.put("GGC", "Gly");
        geneticCode.put("GGA", "Gly");
        geneticCode.put("GGG", "Gly");
        // STOP codons
        geneticCode.put("UAA", "Stop");
        geneticCode.put("UAG", "Stop");
        geneticCode.put("UGA", "Stop");
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.println("Enter the DNA coding sequence:");
        String dnaSequence = scanner.nextLine();
        scanner.close();

        if (dnaSequence == null || dnaSequence.isEmpty()) {
            System.err.println("Error: Input DNA sequence cannot be empty.");
            return;
        }

        System.out.println("\nInput DNA: " + dnaSequence);

        String proteinSequence = translateDnaToProtein(dnaSequence);

        if (proteinSequence != null) {
            System.out.println("Translated Protein: " + proteinSequence);
        } else {
            System.out.println("Translation failed (e.g., no start codon found or invalid sequence).");
        }
    }


    public static String translateDnaToProtein(String dna) {
        // 1. Convert DNA to uppercase and then to RNA (T -> U)
        String rna = dna.toUpperCase().replace('T', 'U');
        System.out.println("RNA sequence:   " + rna);

        // 2. Find the first occurrence of the start codon "AUG"
        int startIndex = rna.indexOf("AUG");
        if (startIndex == -1) {
            System.err.println("Error: Start codon 'AUG' not found.");
            return null;
        }

        // 3. Translate codon by codon starting after the start codon index
        StringBuilder protein = new StringBuilder();
        boolean stopCodonFound = false;

        // Iterate through the RNA sequence in steps of 3 (codons)
        for (int i = startIndex; i < rna.length() - 2; i += 3) {
            // Extract the 3-letter codon
            String codon = rna.substring(i, i + 3);

            // Look up the amino acid in the genetic code map
            String aminoAcid = geneticCode.get(codon);

            if (aminoAcid == null) {
                System.err.println("Warning: Invalid codon '" + codon + "' encountered at position " + i + ". Stopping translation.");
                return protein.length() > 0 ? protein.toString() : null;
            }

            if (aminoAcid.equals("Stop")) {
                stopCodonFound = true;
                break;
            }

            // Append the amino acid abbreviation to the protein sequence
            protein.append(aminoAcid);
        }

        // Check if the loop ended because of reaching the end without a stop codon
        if (!stopCodonFound && (rna.length() - startIndex) % 3 != 0) {
            System.err.println("Warning: Sequence ended before a stop codon was found, and length is not a multiple of 3 after start codon.");
        } else if (!stopCodonFound) {
            System.err.println("Warning: Sequence ended before a stop codon was found.");
        }


        return protein.toString();
    }
}