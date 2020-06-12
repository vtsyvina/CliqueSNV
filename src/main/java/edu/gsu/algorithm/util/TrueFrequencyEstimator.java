package edu.gsu.algorithm.util;

import java.util.Arrays;
import java.util.Comparator;

public class TrueFrequencyEstimator {
    public static final double ILLUMINA_EPS = 0.002;
    public static final double PACBIO_EPS = 0.002;

    /**
     * Calculates true frequencies based on observed and given transposition probability
     * by solving linear equations system
     *
     * @param observed observed frequencies
     * @param allele
     * @return array of
     */
    public static double[] estimateFrequency(double[] observed, int alleles, int allele, double epsLimit) {

        C[] obsSorted = new C[observed.length];
        for (int i = 0; i < observed.length; i++) {
            obsSorted[i] = new C(i, observed[i]);
        }
        Arrays.sort(obsSorted, Comparator.comparingDouble(c -> ((C) c).value).reversed());
        /*
            in case when some allele that is not an SNP has greater frequency then given allele
            Example: allele frequency
                        A      50%
                        -       9%
                        C       2%
            Algorithm says that there are only two alleles - A, C. Then we need to think that there are 3 alleles to properly calculate frequencies
         */
        while (obsSorted[alleles-1].value > observed[allele]) {
            alleles++;
        }

        if (alleles > 4) {
            return observed;
        }
        double[] result = new double[observed.length];
        if (alleles == 2) {
            double eps = Math.min((obsSorted[2].value + obsSorted[3].value + obsSorted[4].value) / 3, epsLimit);
            result[obsSorted[0].index] = (obsSorted[0].value - eps) / (1 - 5 * eps);
            result[obsSorted[1].index] = (obsSorted[1].value - eps) / (1 - 5 * eps);
        }
        if (alleles == 3) {
            double eps = Math.min((obsSorted[3].value + obsSorted[4].value) / 2, epsLimit);
            result[obsSorted[0].index] = (obsSorted[0].value - eps) / (1 - 5 * eps);
            result[obsSorted[1].index] = (obsSorted[1].value - eps) / (1 - 5 * eps);
            result[obsSorted[2].index] = (obsSorted[2].value - eps) / (1 - 5 * eps);
        }
        if (alleles == 4) {
            double eps = Math.min(obsSorted[4].value, epsLimit);
            result[obsSorted[0].index] = (obsSorted[0].value - eps) / (1 - 5 * eps);
            result[obsSorted[1].index] = (obsSorted[1].value - eps) / (1 - 5 * eps);
            result[obsSorted[2].index] = (obsSorted[2].value - eps) / (1 - 5 * eps);
            result[obsSorted[3].index] = (obsSorted[3].value - eps) / (1 - 5 * eps);
        }
        return result;
    }

}

class C {
    int index;
    double value;

    C(int index, double value) {
        this.index = index;
        this.value = value;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("C{");
        sb.append("index=").append(index);
        sb.append(", value=").append(value);
        sb.append('}');
        return sb.toString();
    }
}
