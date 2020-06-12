package edu.gsu.algorithm.em;

import edu.gsu.model.Sample;

import java.util.List;

public class PacBioEM extends AbstractEM {

    private static final double e = 0.01;

    public PacBioEM() {
        init();
    }

    public PacBioEM(int max_read_length) {
        this.MAX_READ_LENGTH = max_read_length;
        init();
    }

    public List<Double> frequencies(List<String> haplotypes, Sample sample) {
        double[][] h = new double[haplotypes.size()][sample.reads.length];
        int f = 0;
        for (String haplotype : haplotypes) {
            int s = 0;
            for (String read : sample.reads) {
                int misses = 0;
                for (int i = 0; i < read.length(); i++) {
                    if (haplotype.charAt(i) != read.charAt(i) && read.charAt(i) != 'N') {
                        misses++;
                    }
                }
                if (misses >= MAX_MISSES) {
                    misses = MAX_MISSES - 1;
                }
                h[f][s++] = pre[read.length()][misses];
            }
            f++;
        }
        return calculateFrequencies(haplotypes.size(), sample.reads.length, h, false);
    }

    @Override
    protected double getE() {
        return e;
    }
}
