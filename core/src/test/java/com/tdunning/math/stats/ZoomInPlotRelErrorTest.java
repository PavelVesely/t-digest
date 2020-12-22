/*
 * Licensed to Ted Dunning under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.tdunning.math.stats;

import org.apache.mahout.common.RandomUtils;
import org.junit.*;

import java.io.*;
import java.lang.Math;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;


/**
 *
 */
//@Ignore
public class ZoomInPlotRelErrorTest extends AbstractTest {

    private static final int[] CompressionsForTesting = {
        100}; //100, 200, 500, 1000, 5000, 10000}; // param delta of t-digest
    private static final int NumberOfPoints = 200; // number of points where we probe the rank estimates

    // params for generating the input
    private static final int N = 50000; // 
    private static final int K = 200; // N / K should be roughly at most 200 (otherwise, we don't have enough precision)
    private static final int PrefixSize = 0; // number of points below zero in each iteration
    private static final int NumberOfRepeats = 100; // N * NumberOfRepeats is (approx.) the number of points in the final instance
    private static final String InputStreamFileName = "t-digest-genInput";
    private static final String InputStreamFileDir = "../../../../data/inputs/"; // CHANGE AS APPROPRIATE
    private static final String DigestStatsFileName = "t-digest-results";
    private static final String DigestStatsDir = "../../../../data/results/"; // CHANGE AS APPROPRIATE
    private static final String FileSuffix = ".csv";

    @BeforeClass
    public static void freezeSeed() {
        RandomUtils.useTestSeed();
    }

    @Ignore
    public void testZoomIn() throws FileNotFoundException, IOException {

        List<Double> sortedData = new ArrayList<Double>();
        List<Double> sortedDataPart = new ArrayList<Double>();
        List<Double> data = new ArrayList<Double>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        PrintWriter w = new PrintWriter(InputStreamFileDir + InputStreamFileName + FileSuffix);
        int n = 0;
        for (int r = 0; r < NumberOfRepeats; r++) {  // same instance is essentially repeated
            double max = Double.MAX_VALUE / (100 * K + r); // initial interval max
            double min = 0; //Double.MIN_VALUE;
            int nn = 0;
            for (int i = -PrefixSize + 1; i <= K - PrefixSize;
                i++) {//for (int i = 0; i <= K; i++) {
                double item = min + (i / (double) (K + 1)) * (max - min);
                data.add(item);
                sortedData.add(item);
                sortedDataPart.add(item);
                w.println(String.valueOf(item));
                nn++;
                n++;
            }
            int maxREindex = -1;
            while (nn < N) {
                maxREindex = -Collections.binarySearch(sortedDataPart, Double.MIN_NORMAL) - 1;
                min = 0; //maxREindex > 0 ? sortedDataPart.get(maxREindex - 1) : 0;
                max = sortedDataPart.get(maxREindex);
                //System.out.printf(String.format("phase maxRErTrue %d;\tmin %s;\tmax %s;\t means %d\n",
                //		maxREindex, String.valueOf(min), String.valueOf(max), digest.centroidCount()));
                if ((max - min) / (double) (K + 1) < Double.MIN_VALUE) {
                    System.out.printf("TOO SMALL max - min\n");
                }
                for (int i = -PrefixSize; i <= K - PrefixSize; i++) {
                    double item = min + (i / (double) (K + 1)) * (max - min);
                    //if (item != 0) { // try to avoid duplicates
                    data.add(item);
                    sortedData.add(item);
                    sortedDataPart.add(item);
                    w.println(String.valueOf(item));
                    n++;
                    nn++;
                    //}
                }
                Collections.sort(sortedDataPart);
            }
            sortedDataPart.clear();
        }
        Collections.sort(sortedData);
        w.close();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
            for (double item : data) {
                digest.add(item);
            }
            digest.compress();

            writeResults(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_K=" + String.valueOf(K) + "_N=" + String
                    .valueOf(N) + "_PS=" + String.valueOf(PrefixSize) + "_repeats=" + String
                    .valueOf(NumberOfRepeats) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
    }


    @Test
    public void carefulNested() throws Exception {

        BigDecimal EPSILON = BigDecimal.valueOf(Double.MIN_VALUE);

        double delta = 500;
        //double delta = 1000;

        List<BigDecimal> data = new ArrayList<>();
        MergingDigest digest = new MergingDigest(delta);
        digest.setScaleFunction(ScaleFunction.K_0);
        digest.setUseAlternatingSort(false);

        int initializingHalfBatchSize = (int) Math.floor(delta * 10);

        // delta=500, denom=100000 seems to work okay..

        double denom = 100000d;
        //denom *= 100;
        BigDecimal infty = BigDecimal
            .valueOf((Double.MAX_VALUE) / denom); // so we can safely average

        BigDecimal increment = infty.divide(BigDecimal.valueOf(initializingHalfBatchSize));

        for (int i = 0; i < 2 * initializingHalfBatchSize; i++) {
            BigDecimal point = BigDecimal.valueOf(-1).multiply(infty);
            for (int j = 0; j < i; j++) {
                point = point.add(increment);
            }
            data.add(point);
            digest.add(point);
        }
        digest.compress();

        // add a glob of zeroes, aim for error there?!?!?!?
        for (int j = 0; j < 5 * initializingHalfBatchSize; j++) {
            data.add(BigDecimal.valueOf(0d));
            digest.add(BigDecimal.valueOf(0d));
        }
        digest.compress();

        Centroid centroidToAttack = belowValue(BigDecimal.ZERO, digest.centroids());
        System.out.println(centroidToAttack);
        BigDecimal centerOfAttack = centroidToAttack.bigMean();
        int weightToRight = 0;
        for (Centroid centroid : digest.centroids()) {
            if (centroid.bigMean().compareTo(centerOfAttack) > 0) {
                weightToRight += centroid.count();
            }
        }

        double weightToLeft;
        int its = 0;
        int weightGoal;
        int currentDeficit;
        int weightOfAttacked;

        Collections.sort(data);
        BigDecimal nextStreamValue = nextValue(centerOfAttack, data);

        BigDecimal previous_c;
        BigDecimal previous_y;

        while (true) {
            its++;

            previous_c = centerOfAttack;
            previous_y = nextStreamValue;

            centerOfAttack = centroidToAttack.bigMean();
            weightOfAttacked = centroidToAttack.count();

            weightToLeft = digest.size() - weightOfAttacked - weightToRight;

            Collections.sort(data);
            nextStreamValue = nextValue(centerOfAttack, data);
            if (!(centerOfAttack.compareTo(nextStreamValue) < 0)) {
                break;
            }

            // weight of the centroid we will fabricate
            weightGoal = (int) Math.ceil((weightToLeft + weightToRight) / ((delta / 2d) - 2d));
            currentDeficit = weightGoal - weightOfAttacked;
            assert currentDeficit >= 0;

            // fill up the old one
            for (int v = 0; v < currentDeficit; v++) {
                BigDecimal point = centerOfAttack.add(EPSILON);
                digest.add(point, 1);
                data.add(point);
            }

            // make the new one
            BigDecimal anotherPoint = centerOfAttack.add(EPSILON).add(EPSILON);
            digest.add(anotherPoint, 1);
            data.add(anotherPoint);

            BigDecimal leftEdge = nextStreamValue.subtract(EPSILON);
            for (int v = 0; v < weightGoal - 2; v++) {
                digest.add(leftEdge, 1);
                data.add(leftEdge);
            }

            digest.compress();
            centroidToAttack = aboveValue(anotherPoint, digest.centroids());

            BigDecimal bad_point = (centerOfAttack.add(nextStreamValue))
                .divide(BigDecimal.valueOf(2d));
            Collections.sort(data);

            System.out.println("finished iteration: " + its);
            System.out.println("td " + digest.cdf(bad_point));
            System.out
                .println("truth " + countBelow(bad_point, data) / (double) data.size() + "\n");

        }

        BigDecimal bad_point = (previous_c.add(previous_y)).divide(BigDecimal.valueOf(2d));
        Collections.sort(data);

        System.out.println("td " + digest.cdf(bad_point));
        System.out.println("truth " + countBelow(bad_point, data) / (double) data.size());
        System.out.println("iterations" + its);
    }

    // assume centroids in ascending order
    private Centroid aboveValue(BigDecimal c, Collection<Centroid> centroids) throws Exception {
        Iterator<Centroid> centroidIterator = centroids.iterator();
        Centroid centroid = centroidIterator.next();
        while (centroidIterator.hasNext()) {
            if (c.compareTo(centroid.bigMean()) < 0) {
                return centroid;
            } else {
                centroid = centroidIterator.next();
            }
        }
        throw new Exception("couldn't find a centroid there");
    }

    // assume centroids in ascending order
    private Centroid belowValue(BigDecimal c, Collection<Centroid> centroids) throws Exception {
        Iterator<Centroid> centroidIterator = centroids.iterator();
        Centroid centroid = centroidIterator.next();
        Centroid previous = null;
        while (centroidIterator.hasNext()) {
            if (centroid.bigMean().compareTo(c) >= 0) {
                return previous;
                // centroid = centroidIterator.next();
            } else {
                previous = centroid;
                centroid = centroidIterator.next();
            }
        }
        throw new Exception("couldn't find a centroid there");
    }

    // assume data is sorted
    private BigDecimal nextValue(BigDecimal c, List<BigDecimal> sortedData) {
        int index = 0;
        while (c.compareTo(sortedData.get(index)) >= 0) {
            index++;
        }
        assert sortedData.get(index).compareTo(c) > 0;
        assert sortedData.get(index - 1).compareTo(c) <= 0;
        return sortedData.get(index);
    }

    private int countBelow(BigDecimal c, List<BigDecimal> sortedData) {
        int index = 0;
        while (c.compareTo(sortedData.get(index)) > 0) {
            index++;
        }
        // assert sortedData.get(index) > c;
        return index;
    }


    private static void writeResults(int compr, int size, TDigest digest,
        List<Double> sortedData, String digestStatsDir, String outName) throws
        IOException {
        Files.createDirectories(Paths.get(DigestStatsDir));
        System.out.printf("stats file:" + outName + "\n");
        File fout = new File(outName);
        fout.createNewFile();
        FileWriter fwout = new FileWriter(fout);

        //System.out.printf("computing rel. errors\n");
        //System.out.flush();

        fwout.write("true quantile;true rank;est. rank;rel. error;abs. error;item\n");
        for (int t = 0; t < NumberOfPoints; t++) {
            //THE FOLLOWING IS EXTREMELY SLOW: Dist.cdf(item, sortedData);
            int rTrue = (int) Math.ceil(t / (float) NumberOfPoints * size) + 1;
            double item = sortedData.get(rTrue - 1);
            // handling duplicate values -- rank is then rather an interval
            int rTrueMin = rTrue;
            int rTrueMax = rTrue;
            while (rTrueMin >= 2 && item == sortedData.get(rTrueMin - 2)) {
                rTrueMin--;
            }
            while (rTrueMax < sortedData.size() && item == sortedData.get(rTrueMax)) {
                rTrueMax++;
            }
            double rEst = digest.cdf(item) * size + 0.5;
            double relErr = 0;
            double addErr = 0;
            if (rEst < rTrueMin) {
                relErr = Math.abs(rTrueMin - rEst) / rTrue;
                addErr = (rEst - rTrueMin) / size;
            }
            if (rEst > rTrueMax) {
                relErr = Math.abs(rTrueMax - rEst) / rTrue;
                addErr = (rEst - rTrueMax) / size;
            }
            fwout.write(String
                .format("%.6f;%d;%.6f;%.6f;%.6f;%s\n", rTrue / (float) size, (int) rTrue, rEst,
                    relErr, addErr, String.valueOf(item)));
        }
        fwout.write("\n");
        fwout.write(String.format("n=%d\n", size));
        fwout.write(String.format("scale func. = %s\n", digest.scale.toString()));
        fwout.write(String.format("delta = %d (compression param of t-digest)\n", compr));
        fwout.write(String.format("# of centroids = %d\n", digest.centroids().size()));
        fwout.write(String.format("size in bytes = %d\n", digest.byteSize()));
        fwout.write("\nCentroids:\n");
        for (Centroid centr : digest.centroids()) {
            fwout.write(centr.toString() + "\n");
        }
        fwout.close();
        System.out.flush();

    }
}
