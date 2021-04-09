package edu.gsu.model;

import java.util.List;

public class IlluminaSNVSample {
    public String name;
    public List<PairEndRead> reads;
    public int referenceLength;
    public boolean singleRead; // tells if the sample consists mainly of single-end reads

    public IlluminaSNVSample(){}

    public IlluminaSNVSample(String name, List<PairEndRead> reads, int referenceLength, boolean singleRead) {
        this.name = name;
        this.reads = reads;
        this.referenceLength = referenceLength;
        this.singleRead = singleRead;
    }
}
