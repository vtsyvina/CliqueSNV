package edu.gsu.model;

public class PairEndRead {
    public String l;
    public String r;
    public int lOffset;
    public int rOffset;
    public String name;

    public PairEndRead(){

    }

    public PairEndRead(String l, String r, int lOffset, int rOffset, String name) {
        this.l = l;
        this.r = r;
        this.lOffset = lOffset;
        this.rOffset = rOffset;
        this.name = name;
    }
}
