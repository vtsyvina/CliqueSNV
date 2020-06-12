package edu.gsu.algorithm.util;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;

public class EpsilonEstimator {

    /**
     * Estimates probability of transaction based on observed and true frequencies based on quadratic optimization
     * @param observed array {O11, O12, O21, O22} of observed frequencies
     * @param trueFrequency array {T11, T12, T21, T22} of true frequencies
     * @return array {eA12, eA21, eB12, eB21} of transition probabilities that best describe observed frequencies
     */
    double[] estimateEpsilon(double[] observed, double[] trueFrequency){
        double[][] m = {

        };
        Array2DRowRealMatrix array2DRowRealMatrix = new Array2DRowRealMatrix();
        //LUDecomposition luDecomposition = new LUDecomposition();
        return null;
    }
}
