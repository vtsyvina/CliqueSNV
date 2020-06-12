package edu.gsu.algorithm.em;

import java.util.List;

public class SNPEM extends AbstractEM{

    public double[] frequencies(double eA12, double eA21, double eB12, double eB21){
        double[] f = {3./8, 1./4, 1./4, 1./8};
        double eA11 = 1 - eA21;
        double eA22 = 1 - eA12;
        double eB11 = 1 - eB21;
        double eB22 = 1 - eB12;

        //I really hope that I didn't make any typo
        double[][] h = {
                {eA11*eB11, eA11*eB21, eA21*eB11, eA21*eB21}, //O11
                {eA11*eB12, eA11*eB22, eA21*eB12, eA21*eB22}, //O12
                {eA12*eB11, eA12*eB21, eA22*eB11, eA22*eB21}, //O21
                {eA12*eB12, eA12*eB22, eA22*eB12, eA22*eB22} //O22
        };
        return calculateFrequencies(4, 4, h, false).stream().mapToDouble(Double::doubleValue).toArray();
    }

    @Override
    protected double getE() {
        //we won't really use it
        return 0;
    }
}
