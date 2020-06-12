package edu.gsu.algorithm.imputation;

import edu.gsu.algorithm.SNVPacBioMethod;
import edu.gsu.model.SNVResultContainer;
import edu.gsu.model.Sample;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This class implements algorithm for imputation problem based in CliqueSNV algorithm
 */
public class Imputation {

     private SNVPacBioMethod pacBio;
     private Sample sample;
     private boolean log;

    /**
     * @param sample       Given input reads with N in the beginning or end of read to make them all equal size
     * @param minThreshold threshold for O22 value to consider alleles as SNPs
     * @param minFreq      minimum frequency (relative to reads' coverage) for O22 value to consider alleles as SNP
     */
    public Imputation(Sample sample, int minThreshold, double minFreq, boolean log) {
        pacBio = new SNVPacBioMethod(sample, minThreshold, minFreq, log);
        pacBio.CLOSE_POSITIONS_THRESHOLD = 0;
        pacBio.MINIMUM_READS_NUMBER_IN_CLUSTER = 0;
        pacBio.READS_PORTION_TO_FILTER = 0;
        pacBio.setForImputation(true);
        this.sample = sample;
        this.log = log;
    }

    public Sample getImputedHaplotypes(){
        Sample result = new Sample();
        result.name = sample.name+"_imputed";
        List<SNVResultContainer> haplotypes = pacBio.getHaplotypes();
        if (log) System.out.println("Haplotypes "+haplotypes.size());
        Map<Integer, String> imputedMap = new HashMap<>();
        // haplotypes are sorted by frequency, so no problems to put first occurrence
        haplotypes.forEach( h -> h.pacBioCluster.forEach((key, value) -> {
            if (!imputedMap.containsKey(key)) {
                StringBuilder str = new StringBuilder();
                for (int i = 0; i < value.length(); i++) {
                    if (value.charAt(i) == 'N') {
                        str.append(h.haplotype.charAt(i));
                    } else {
                        str.append(value.charAt(i));
                    }
                }
                imputedMap.put(key, str.toString());
            }
        }));
        result.reads = new String[sample.reads.length];
        for (int i = 0; i < sample.reads.length; i++) {
            if (!imputedMap.containsKey(i)){
                System.out.println("Error! Imputed haplotype wasn't found for index "+i);
            } else {
                result.reads[i] = imputedMap.get(i);
            }

        }
        return result;
    }
}
