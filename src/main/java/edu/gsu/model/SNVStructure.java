package edu.gsu.model;

public class SNVStructure {
    public double[][] profile;
    public int[][] rowMinors;//splitted
    public int[][] rowN;//not splitted
    public int[][] colMinors;//splitted
    public int[] majorsInRow;//not splitted
    //for Illumina only
    public int[][] readsAtPosition;//not splitted
}
