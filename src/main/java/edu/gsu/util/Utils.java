package edu.gsu.util;

import edu.gsu.model.IlluminaSNVSample;
import edu.gsu.model.PairEndRead;
import edu.gsu.model.Sample;
import org.apache.commons.math3.distribution.BinomialDistribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.function.Predicate;

/**
 * Just class with some useful functions
 */
public class Utils {
    private static final List<String> impossibleCharacters = new ArrayList<>();
    private static final char impossibleChar = '%';
    public static final String DEFAULT_ALPHABET = "ACGT";

    public static List<String> numbers = new ArrayList<>();

    static {
        for (int i = 0; i < 100_000; i++) {
            numbers.add(String.valueOf(i));
        }
    }

    static {
        impossibleCharacters.add("");
        for (int i = 1; i < 3000; i++) {
            impossibleCharacters.add(impossibleCharacters.get(i - 1) + impossibleChar);
        }
    }

    public static int convertLetterToDigit(char c) {
        return convertLetterToDigit(c, DEFAULT_ALPHABET);
    }

    public static int convertLetterToDigit(char c, String alphabet) {
        int i = alphabet.indexOf(c);
        // 5 == N, if there is some trash letters for some reason
        return i == -1? 5 : i;
    }

    public static char convertIntToLetter(int i) {
        switch (i) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            case 4:
                return '-';
        }
        return 'N';
    }

    public static String consensus(String[] sequences) {
        return consensus(sequences, DEFAULT_ALPHABET);
    }

    public static String consensus(String[] sequences, String alphabet) {
        if (sequences.length == 0) {
            return "";
        }
        int l = sequences[0].length();
        int[][] count = new int[alphabet.length()][l];
        for (String s : sequences) {
            for (int i = 0; i < s.length(); i++) {
                count[convertLetterToDigit(s.charAt(i), alphabet)][i]++;
            }
        }
        StringBuilder str = new StringBuilder();
        int alphabetLength = alphabet.indexOf('N') == -1 ? alphabet.length() : alphabet.length() - 1;
        for (int i = 0; i < l; i++) {
            int max = 0;
            for (int j = 1; j < alphabetLength; j++) {
                if (count[j][i] > count[max][i]) {
                    max = j;
                }
            }
            if (count[max][i] == 0) {
                max = -1;
            }
            str.append(convertIntToLetter(max));
        }
        return str.toString();
    }

    public static String consensus(Sample sample) {
        return consensus(sample.reads);
    }

    /**
     * Return consensus for given profile. Put $ if all profile values equals 0
     *
     * @param profile  sample profile
     * @param alphabet alphabet for profile
     * @return consensus
     */
    public static String consensus(double[][] profile, String alphabet) {
        char[] majors = new char[profile[0].length];
        for (int j = 0; j < profile[0].length; j++) {
            int majorAllele = Utils.getMajorAllele(profile, j);
            majors[j] = majorAllele == -1 ? '-' : alphabet.charAt(majorAllele);
        }
        return new String(majors);
    }

    public static double[][] profile(Sample sample) {
        return profile(sample, DEFAULT_ALPHABET);
    }

    public static double[][] profile(Sample sample, String alphabet) {
        String[] sequences = sample.reads;
        if (sequences.length == 0) {
            return new double[0][alphabet.length()];
        }
        int l = sequences[0].length();
        int[][] count = new int[alphabet.length()][l];
        for (String s : sequences) {
            for (int i = 0; i < s.length(); i++) {
                count[convertLetterToDigit(s.charAt(i), alphabet)][i]++;
            }
        }
        double[][] result = new double[alphabet.length()][l];
        for (int i = 0; i < alphabet.length(); i++) {
            for (int j = 0; j < l; j++) {
                result[i][j] = count[i][j] / (double) sequences.length;
            }
        }
        return result;
    }

    public static int[][] countCoverage(IlluminaSNVSample sample, String alphabet){
        return countCoverage(sample.reads, alphabet, sample.referenceLength);
    }

    public static int[][] countCoverage(List<PairEndRead> reads, String alphabet, int referenceLength){
        if (reads.size() == 0) {
            return new int[0][alphabet.length()];
        }
        int[][] count = new int[alphabet.length()][referenceLength];
        reads.forEach(r -> {
            for (int i = 0; i < r.l.length(); i++) {
                count[convertLetterToDigit(r.l.charAt(i), alphabet)][i + r.lOffset]++;
            }
            for (int i = 0; i < r.r.length(); i++) {
                count[convertLetterToDigit(r.r.charAt(i), alphabet)][i + r.rOffset]++;
            }
        });
        return count;
    }

    public static int[][] countCoverage(Sample sample, String alphabet){
        if (sample.reads.length == 0) {
            return new int[0][alphabet.length()];
        }
        int[][] count = new int[alphabet.length()][sample.reads[0].length()];
        for (String read : sample.reads) {
            for (int i = 0; i < read.length(); i++) {
                int idx = convertLetterToDigit(read.charAt(i), alphabet);
                if (idx == 5){
                    continue;
                }
                count[idx][i]++;
            }

        }
        return count;
    }

    public static double[][] profile(IlluminaSNVSample sample, String alphabet) {
        List<PairEndRead> reads = sample.reads;
        if (reads.size() == 0) {
            return new double[0][alphabet.length()];
        }
        int[][] count = new int[alphabet.length()][sample.referenceLength];
        reads.forEach(r -> {
            for (int i = 0; i < r.l.length(); i++) {
                count[convertLetterToDigit(r.l.charAt(i), alphabet)][i + r.lOffset]++;
            }
            for (int i = 0; i < r.r.length(); i++) {
                count[convertLetterToDigit(r.r.charAt(i), alphabet)][i + r.rOffset]++;
            }
        });
        double[][] result = new double[alphabet.length()][sample.referenceLength];
        for (int i = 0; i < count[0].length; i++) {
            int sum = 0;
            for (int[] aCount : count) {
                sum += aCount[i];
            }
            if (sum == 0) {
                sum = 1;
            }
            for (int j = 0; j < count.length; j++) {
                result[j][i] = count[j][i] / (double) sum;
            }
        }
        return result;
    }

    /**
     * Append missing characters to string so they have the same size
     */
    public static String[] stringsForHamming(String[] sequences) {
        int max = Arrays.stream(sequences).mapToInt(String::length).max().getAsInt();
        String[] result = new String[sequences.length];
        for (int i = 0; i < sequences.length; i++) {
            result[i] = sequences[i].length() < max ?
                    sequences[i] + impossibleCharacters.get(max - sequences[i].length()) :
                    sequences[i];
        }
        return result;
    }

    public static String[] stringsForHamming(String[] sequences, char character) {
        List<String> chars = new ArrayList<>();
        chars.add("");
        for (int i = 1; i < 3000; i++) {
            chars.add(chars.get(i - 1) + character);
        }
        int max = Arrays.stream(sequences).mapToInt(String::length).max().getAsInt();
        String[] result = new String[sequences.length];
        for (int i = 0; i < sequences.length; i++) {
            result[i] = sequences[i].length() < max ?
                    sequences[i] + chars.get(max - sequences[i].length()) :
                    sequences[i];
        }
        return result;
    }

    public static void expandNumbers(int size) {
        if (size > numbers.size()) {
            for (int i = numbers.size(); i < size; i++) {
                numbers.add(String.valueOf(i));
            }
        }
    }

    public static int getMajorAllele(double[][] profile, int i) {
        int major = 0;
        //don't count N
        for (int j = 1; j < profile.length - 1; j++) {
            if (profile[j][i] > profile[major][i]) {
                major = j;
            }
        }
        //in case there is no reads covering this position
        return profile[major][i] < 0.001 ? -1 : major;
    }

    public static int getMajorAllele(String consensus, String alphabet, int i) {
        return alphabet.indexOf(consensus.charAt(i));
    }

    public static String byteArrayToString(byte[] arr) {
        StringBuilder str = new StringBuilder();
        for (byte b : arr) {
            str.append((char) b);
        }
        return str.toString();
    }

    public static double binomialPvalue(int s, double p, int n) {
        return 1 - new BinomialDistribution(n, p).cumulativeProbability(s);
    }

    public static double binomialLogPvalue(int s, double p, int n) {
        return -new BinomialDistribution(n, p).logProbability(s);
    }

    public static <T> Predicate<T> distinctByKey(Function<? super T, ?> keyExtractor) {
        Map<Object, Boolean> seen = new ConcurrentHashMap<>();
        return t -> seen.putIfAbsent(keyExtractor.apply(t), Boolean.TRUE) == null;
    }

    /**
     * Used to convert frequencies to String. Leaves at least 3 digits if x < 0.1, leaves 4 digits otherwise
     * 0.853434347 -> 0.853
     * 0.001112 -> 0.00111
     * 0.000034394 -> 0.0000344
     *
     * @param x
     * @return
     */
    public static String smartDoubleToString(double x) {
        int i = 0;
        double eps = 1;
        while (x < eps && i < 10) {
            i++;
            eps /= 10;
        }
        if (i == 1) {
            i = 3;
        } else if (i < 10) {
            i += 2;
        } else {
            i = 1;
        }
        return String.format("%." + i + "f", x);
    }

    public static int splittedPosition(int pos, char minor, int minorCount, String al, String consensus) {
        int r = pos * minorCount;
        int allele = al.indexOf(minor) >= Utils.getMajorAllele(consensus, al, pos) ? al.indexOf(minor) - 1 : al.indexOf(minor);
        return r + allele;
    }

    public static double[] normalize(double[] x) {
        double sum = 0;
        double[] res = new double[x.length];
        for (int i = 0; i < res.length; i++) {
            res[i] = x[i];
            sum += x[i];
        }
        if (sum == 0) {
            return res;
        }
        for (int i = 0; i < res.length; i++) {
            res[i] /= sum;
        }
        return res;
    }

    public static String humanReadableSI(long bytes) {
        String s = bytes < 0 ? "-" : "";
        long b = bytes == Long.MIN_VALUE ? Long.MAX_VALUE : Math.abs(bytes);
        return b < 1000L ? bytes +""
                : b < 999_950L ? String.format("%s%.1fk", s, b / 1e3)
                : (b /= 1000) < 999_950L ? String.format("%s%.1fM", s, b / 1e3)
                : (b /= 1000) < 999_950L ? String.format("%s%.1fG", s, b / 1e3)
                : (b /= 1000) < 999_950L ? String.format("%s%.1fT", s, b / 1e3)
                : (b /= 1000) < 999_950L ? String.format("%s%.1fP", s, b / 1e3)
                : String.format("%s%.1f EB", s, b / 1e6);
    }

    public static int[][] rotateMatrix(int[][] x){
        int[][] result = new int[x[0].length][x.length];
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                result[i][j] = x[j][i];
            }
        }
        return  result;
    }

    public static int hammingDistance(String x, String y){
        int c = 0;
        for (int i = 0; i < x.length(); i++) {
            if(x.charAt(i) != y.charAt(i)){
                c++;
            }
        }
        return c;
    }
}