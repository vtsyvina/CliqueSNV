package edu.gsu.model;

import java.util.Objects;

public class Interval implements Comparable<Interval> {
    public int s;
    public int e;

    public Interval() {
    }

    public Interval(int s, int e) {
        this.s = s;
        this.e = e;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Interval interval = (Interval) o;
        return s == interval.s && e == interval.e;
    }

    @Override
    public int hashCode() {
        return Objects.hash(s, e);
    }

    @Override
    public String toString() {
        return "[" + s + ", " + e + "]";
    }

    @Override
    public int compareTo(Interval o) {
        return Integer.compare(s, o.s);
    }
}
