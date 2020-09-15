package edu.gsu.model;

import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Class contains haplotype obtained in SNV method and some additional
 * human-friendly fields to understand and analyse results
 */
public class SNVResultContainer {

    public Map<Integer, String> pacBioCluster;
    public List<PairEndRead> illuminaCluster;
    public Clique haploClique;
    public String haplotype;
    public String clusterString;
    public Clique sourceClique;
    public double frequency;


    public SNVResultContainer(String clusterString, Map<Integer, String> pacBioCluster, Clique haplotypeClique, String haplotype) {
        this.pacBioCluster = pacBioCluster;
        this.haploClique = haplotypeClique;
        this.haplotype = haplotype;
        this.clusterString = clusterString;
    }

    public SNVResultContainer(String clusterString, Map<Integer, String> pacBioCluster, Clique haplotypeClique, String haplotype, double frequency) {
        this.pacBioCluster = pacBioCluster;
        this.haploClique = haplotypeClique;
        this.haplotype = haplotype;
        this.clusterString = clusterString;
        this.frequency = frequency;
    }



    public SNVResultContainer(String clusterString,  Clique haplotypeClique, String haplotype, List<PairEndRead> illuminaCluster) {
        this.illuminaCluster = illuminaCluster;
        this.haploClique = haplotypeClique;
        this.haplotype = haplotype;
        this.clusterString = clusterString;
    }

    public SNVResultContainer(String clusterString,  Clique haplotypeClique, Clique sourceClique, String haplotype, List<PairEndRead> illuminaCluster,Map<Integer, String> pacBioCluster,  double frequency) {
        this.illuminaCluster = illuminaCluster;
        this.pacBioCluster = pacBioCluster;
        this.sourceClique = sourceClique;
        this.haploClique = haplotypeClique;
        this.haplotype = haplotype;
        this.clusterString = clusterString;
        this.frequency = frequency;
    }

    public SNVResultContainer copy(){
        return new SNVResultContainer(clusterString, haploClique, sourceClique, haplotype, illuminaCluster, pacBioCluster, frequency);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        SNVResultContainer that = (SNVResultContainer) o;
        return Objects.equals(haploClique, that.haploClique);
    }

    @Override
    public int hashCode() {
        return Objects.hash(haploClique);
    }

    @Override
    public String toString() {
        return "{\n" +
                " snps=" + haploClique +
                ",\n source clique='" + sourceClique +
                ",\n frequency= " + frequency +
                ",\n haplotype='" + haplotype + "\n\'" +
                "}";
    }
}
