package edu.gsu.model;

/**
 * Data container to store sample data
 */
public class Sample {
    public String name;
    /**
     * Contains list of all reads
     */
    public String[] reads;

    public String[] readNames;

    public Sample() {
    }

    public Sample(String name, String[] sequences) {
        this.name = name;
        this.reads = sequences;
    }

    public Sample(String name, String[] sequences, String[] readNames) {
        this.name = name;
        this.reads = sequences;
        this.readNames = readNames;
    }
}

