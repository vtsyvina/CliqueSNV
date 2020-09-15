package edu.gsu.algorithm.em;

import edu.gsu.model.IlluminaSNVSample;
import edu.gsu.model.PairEndRead;
import edu.gsu.start.Start;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

public class IlluminaEM extends AbstractEM {

    private static final double e = 0.001;

    public IlluminaEM() {
        init();
    }

    @Override
    protected double getE() {
        return e;
    }

    public List<Double> frequencies(List<String> haplotypes, IlluminaSNVSample sample) {
        double[] frequencies = new double[haplotypes.size()];
        for (int i = 0; i < frequencies.length; i++) {
            frequencies[i] = 1 / (double) frequencies.length;
        }
        double[][] h = new double[haplotypes.size()][sample.reads.size()];
        List<Callable<Boolean>> tasks = new ArrayList<>();
        for (int f = 0; f < haplotypes.size(); f++) {
            int finalF = f;
            tasks.add(() -> {
                String haplotype = haplotypes.get(finalF);
                for (int s = 0; s < sample.reads.size(); s++) {
                    int misses = 0;
                    PairEndRead read = sample.reads.get(s);
                    for (int i = 0; i < read.l.length(); i++) {
                        if (haplotype.charAt(read.lOffset + i) != read.l.charAt(i)) {
                            misses++;
                        }
                    }
                    for (int i = 0; i < read.r.length(); i++) {
                        if (haplotype.charAt(read.rOffset + i) != read.r.charAt(i)) {
                            misses++;
                        }
                    }
                    if (misses >= MAX_MISSES){
                        misses = MAX_MISSES-1;
                    }
                    int length = Math.min(MAX_READ_LENGTH-1, read.l.length() + read.r.length());
                    h[finalF][s] = pre[length][misses];
                }
                return true;
            });
        }
        List<Future<Boolean>> futures = null;
        try {
            futures = Start.executor.invokeAll(tasks);
            for (Future<Boolean> f : futures) {
                f.get();
            }
        } catch (InterruptedException | ExecutionException ex) {
            ex.printStackTrace();
        }

        return calculateFrequencies(haplotypes.size(), sample.reads.size(), h, true);
    }
}
