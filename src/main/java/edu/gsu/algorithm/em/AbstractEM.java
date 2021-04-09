package edu.gsu.algorithm.em;

import edu.gsu.start.Start;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

public abstract class AbstractEM {
    protected int MAX_EM_ITER = 100;
    protected int MAX_READ_LENGTH = 3000;
    protected final int MAX_MISSES = 1000;
    private final double e = getE();
    protected static final double eps = 0.001;

    protected double[][] pre;

    protected double[][] pij;

    protected void init() {
        pre = new double[MAX_READ_LENGTH][MAX_MISSES];
        pre[1][0] = 1 - e;
        pre[1][1] = e / 3;
        for (int i = 2; i < MAX_READ_LENGTH; i++) {
            pre[i][0] = pre[i - 1][0] * (1 - e);
            for (int j = 1; j <= i && j < MAX_MISSES; j++) {
                pre[i][j] = pre[i - 1][j - 1] * (e / 3);
            }
        }
    }

    protected abstract double getE();


    protected double euclidDistance(double[] x, double[] y) {
        double d = 0;
        for (int i = 0; i < x.length; i++) {
            d += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(d);
    }

    protected List<Double> calculateFrequencies(int haplotypesCount, int readsCount, double[][] h, boolean log) {
        pij = new double[haplotypesCount][readsCount];
        double[] frequencies = new double[haplotypesCount];
        Arrays.fill(frequencies, 1 / (double) frequencies.length);
        double[] oldFrequencies;
        if(log) System.out.println("Start iterations EM");
        int iter = 0;
        do {
            oldFrequencies = frequencies;
            frequencies = new double[haplotypesCount];
            double[] denominators = new double[readsCount];
            for (int i = 0; i < readsCount; i++) {
                double d = 0;
                for (int j = 0; j < haplotypesCount; j++) {
                    d += oldFrequencies[j] * h[j][i];
                }
                denominators[i] = d;
            }
            double[] m = new double[haplotypesCount];
            double sum = 0;
            List<Callable<Double>> tasks = new ArrayList<>();

            for (int j = 0; j < haplotypesCount; j++) {
                int finalJ = j;
                double[] finalOldFrequencies = oldFrequencies;
                tasks.add(() -> {
                    double result = 0;
                    for (int i = 0; i < readsCount; i++) {
                        if (denominators[i] == 0.0) {
                            continue;
                        }
                        pij[finalJ][i] = finalOldFrequencies[finalJ] * h[finalJ][i] / denominators[i];
                        result += pij[finalJ][i];
                    }
                    return result;
                });

            }
            try {
                List<Future<Double>> futures = Start.executor.invokeAll(tasks);
                for (int j = 0; j < haplotypesCount; j++) {
                    m[j] = futures.get(j).get();
                    sum+=m[j];
                }
            } catch (InterruptedException | ExecutionException ex) {
                ex.printStackTrace();
            }

            for (int j = 0; j < haplotypesCount; j++) {
                frequencies[j] = m[j] / sum;
            }
            //if (log) System.out.print("\r"+iter++);
        } while (euclidDistance(oldFrequencies, frequencies) > eps && iter < MAX_EM_ITER);
        if(log) System.out.println();
        return Arrays.stream(frequencies)
                .boxed()
                .collect(Collectors.toList());
    }

    public double[][] getPij(){
        return pij;
    }

}
