package edu.gsu.util;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.Interval;

import java.util.function.Function;
import java.util.function.Predicate;

/**
 BEGIN_DOC

 ## Motivation

 "Sequence logo ( http://weblogo.berkeley.edu/logo.cgi ) for different alleles or generated from SAM/BAM" http://www.biostars.org/p/73021

 ![ScreenShot](https://raw.github.com/lindenb/jvarkit/master/doc/sam2weblogo.png)


 ## Example

 ```bash
 $ java -jar dist/sam4weblogo.jar -r seq1:80-110  sorted.bam  2> /dev/null | head -n 50
 >B7_593:4:106:316:452/1
 TGTTG--------------------------
 >B7_593:4:106:316:452a/1
 TGTTG--------------------------
 >B7_593:4:106:316:452b/1
 TGTTG--------------------------
 >B7_589:8:113:968:19/2
 TGGGG--------------------------
 >B7_589:8:113:968:19a/2
 TGGGG--------------------------
 >B7_589:8:113:968:19b/2
 TGGGG--------------------------
 >EAS54_65:3:321:311:983/1
 TGTGGG-------------------------
 >EAS54_65:3:321:311:983a/1
 TGTGGG-------------------------
 >EAS54_65:3:321:311:983b/1
 TGTGGG-------------------------
 >B7_591:6:155:12:674/2
 TGTGGGGG-----------------------
 >B7_591:6:155:12:674a/2
 TGTGGGGG-----------------------
 >B7_591:6:155:12:674b/2
 TGTGGGGG-----------------------
 >EAS219_FC30151:7:51:1429:1043/2
 TGTGGGGGGCGCCG-----------------
 >EAS219_FC30151:7:51:1429:1043a/2
 TGTGGGGGGCGCCG-----------------
 >EAS219_FC30151:7:51:1429:1043b/2
 TGTGGGGGGCGCCG-----------------
 >B7_591:5:42:540:501/1
 TGTGGGGGCCGCAGTG---------------
 >EAS192_3:5:223:142:410/1
 TGGGGGGGGCGCAGT----------------
 >B7_591:5:42:540:501a/1
 TGTGGGGGCCGCAGTG---------------
 >EAS192_3:5:223:142:410a/1
 TGGGGGGGGCGCAGT----------------
 >B7_591:5:42:540:501b/1
 TGTGGGGGCCGCAGTG---------------
 >EAS192_3:5:223:142:410b/1
 TGGGGGGGGCGCAGT----------------
 ```

 ## See also

 * https://www.biostars.org/p/103052/
 * http://www.sciencedirect.com/science/article/pii/S1874778715300210

 END_DOC

 */

public class SAM4WebLogo {
    private boolean useClip = false;

    private final Function<SAMRecord, Integer> readStart = rec ->
            useClip ? rec.getUnclippedStart() : rec.getAlignmentStart();

    private final Function<SAMRecord, Integer> readEnd = rec ->
            useClip ? rec.getUnclippedEnd() : rec.getAlignmentEnd();
    private Interval interval;

    public SAM4WebLogo(SamReader samReader) {
        final IntervalParser intervalParser = new IntervalParser().setFixContigName(true);
        this.interval = intervalParser.
                setDictionary(samReader.getFileHeader().getSequenceDictionary()).
                parse(samReader.getFileHeader().getSequenceDictionary().getSequences().get(0).getSequenceName() + ":0-" + samReader.getFileHeader().getSequenceDictionary().getSequences().get(0).getSequenceLength());
    }

    public String printRead(
            final SAMRecord rec
    ) {

        final Cigar cigar = rec.getCigar();
        final Function<Integer, Character> read2base = IDX -> {
            byte bases[] = rec.getReadBases();
            if (SAMRecord.NULL_SEQUENCE.equals(bases)) return '?';
            if (IDX < 0 || IDX >= bases.length) return '?';
            return (char) bases[IDX];
        };

        final Predicate<Integer> inInterval = IDX -> IDX >= interval.getStart() && IDX <= interval.getEnd();

        final StringBuilder seq = new StringBuilder(interval.length());

        int refPos = Math.min(
                interval.getStart(),
                rec.getUnclippedStart()
        );

        while (refPos < rec.getUnclippedStart()) {
            if (inInterval.test(refPos)) {
                seq.append('-');
            }
            ++refPos;
        }

        int readPos = 0;
        for (int i = 0; i < cigar.numCigarElements(); ++i) {
            final CigarElement ce = cigar.getCigarElement(i);
            final CigarOperator op = ce.getOperator();
            switch (op) {
                case P:
                    break;
                case I: {
                    readPos += ce.getLength();
                    break;
                }
                case D:
                case N: {
                    for (int j = 0; j < ce.getLength() && refPos <= interval.getEnd(); ++j) {
                        if (inInterval.test(refPos)) {
                            seq.append('-');
                        }
                        refPos++;
                    }
                    break;
                }
                case H: {
                    for (int j = 0; j < ce.getLength() && refPos <= interval.getEnd(); ++j) {
                        if (inInterval.test(refPos)) {
                            seq.append(useClip ? 'n' : '-');
                        }
                        refPos++;
                    }
                    break;
                }
                case S: {
                    for (int j = 0; j < ce.getLength() && refPos <= interval.getEnd(); ++j) {
                        if (inInterval.test(refPos)) {
                            if (useClip) {
                                seq.append(Character.toLowerCase(read2base.apply(readPos)));
                            } else {
                                seq.append('-');
                            }
                        }
                        readPos++;
                        refPos++;
                    }
                    break;
                }
                case M:
                case X:
                case EQ: {
                    for (int j = 0; j < ce.getLength() && refPos <= interval.getEnd(); ++j) {
                        if (inInterval.test(refPos)) {
                            seq.append(read2base.apply(readPos));
                        }
                        readPos++;
                        refPos++;
                    }
                    break;
                }
                default:
                    throw new IllegalStateException("Not handled. op:" + op);
            }
        }

        while(refPos<= interval.getEnd())
        {
            seq.append('-');
            refPos++;
        }
        return seq.toString();
    }
}