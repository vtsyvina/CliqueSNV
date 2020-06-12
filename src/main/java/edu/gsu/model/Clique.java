package edu.gsu.model;

import edu.gsu.algorithm.AbstractSNV;
import edu.gsu.algorithm.SNVPacBioMethod;
import edu.gsu.util.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.Set;

public class Clique {
    public String minors;
    public List<Integer> snps;
    public Set<Integer> splittedSnps;
    // only for VC
    public List<Double> frequencies;

    public Clique(String haplo, String consensus){
        this.snps = new ArrayList<>();
        StringBuilder str = new StringBuilder();
        for (int i = 0; i < haplo.length(); i++) {
            if (haplo.charAt(i) != consensus.charAt(i)){
                snps.add(i);
                str.append(haplo.charAt(i));
            }
        }
        this.minors = str.toString();
    }

    public Clique(Set<Integer> splittedSnps, String consensus) {
        this.snps = new ArrayList<>();
        StringBuilder str = new StringBuilder();
        splittedSnps.stream().sorted().forEach(i -> {
            int pos = i / SNVPacBioMethod.minorCount;
            int allele1 = i % SNVPacBioMethod.minorCount >= Utils.getMajorAllele(consensus, AbstractSNV.al, pos) ? i % SNVPacBioMethod.minorCount + 1 : i % SNVPacBioMethod.minorCount;
            char m1 = SNVPacBioMethod.al.charAt(allele1);
            this.snps.add(pos);
            str.append(m1);
        });
        this.minors = str.toString();
        this.splittedSnps = splittedSnps;
    }

    public Clique(Set<Integer> splittedSnps, String consensus, List<Double> frequencies){
        this(splittedSnps, consensus);
        this.frequencies = frequencies;
    }

    @Override
    public String toString() {
        return minors + snps;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Clique clique = (Clique) o;

        if (!Objects.equals(minors, clique.minors)) return false;
        if (!Objects.equals(snps, clique.snps)) return false;
        return Objects.equals(splittedSnps, clique.splittedSnps);
    }

    @Override
    public int hashCode() {
        int result = minors != null ? minors.hashCode() : 0;
        result = 31 * result + (snps != null ? snps.hashCode() : 0);
        result = 31 * result + (splittedSnps != null ? splittedSnps.hashCode() : 0);
        return result;
    }
}