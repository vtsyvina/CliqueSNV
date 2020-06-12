package edu.gsu.model;

import java.util.List;

public class IlluminaSNVSample {
    public String name;
    public List<PairEndRead> reads;
    public int referenceLength;

    public IlluminaSNVSample(){}

    public IlluminaSNVSample(String name, List<PairEndRead> reads, int referenceLength) {
        this.name = name;
        this.reads = reads;
        this.referenceLength = referenceLength;
    }
}
