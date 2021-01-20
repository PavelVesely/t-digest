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

    //mvn test -Dtests.seed=1639B3E7B594EF9 -DrunSlowTests=false -Dtests.class=com.tdunning.math.stats.ZoomInPlotRelErrorTest
    //      	at __randomizedtesting.SeedInfo.seed([1639B3E7B594EF9:39272979715260C0]:0)
    @Test
    public void carefulNested() throws Exception {

        BigDecimal EPSILON = BigDecimal.valueOf(0.1); //BigDecimal.valueOf(Double.MIN_VALUE);

        double delta = 500;
        //double delta = 1000;

        List<BigDecimal> data = new ArrayList<>();
        MergingDigest digest = new MergingDigest(delta);
        //AVLTreeDigest digest = new AVLTreeDigest(delta);
        digest.setScaleFunction(ScaleFunction.K_0);
        digest.setUseAlternatingSort(false);

        int initializingHalfBatchSize = (int) Math.floor(delta * 20);

        // delta=500, denom=100000 seems to work okay..

        double denom = 1000000000d;
        //denom *= 100;
        BigDecimal infty = BigDecimal
            .valueOf((Double.MAX_VALUE) / denom); // so we can safely average
        BigDecimal increment = infty.divide(BigDecimal.valueOf(initializingHalfBatchSize));

        for (int i = 0; i < initializingHalfBatchSize; i++) {
            //BigDecimal point = BigDecimal.valueOf(-1).multiply(infty);
            BigDecimal point = BigDecimal.TEN.multiply(BigDecimal.valueOf(-1d));
//            for (int j = 0; j < i; j++) {
//                point = point.add(increment);
//            }
            point = point.add(increment.multiply(BigDecimal.valueOf(i)));
            data.add(point);
            digest.add(point);
        }
        // add a glob of zeroes, aim for error there.
        for (int j = 0; j < 5 * initializingHalfBatchSize; j++) {
            data.add(BigDecimal.ZERO);
            digest.add(BigDecimal.ZERO);
        }
        digest.compress();

        Centroid centroidToAttack = belowValue(BigDecimal.TEN, digest.centroids());
        Centroid rightNeighbor = aboveValue(centroidToAttack.bigMean(), digest.centroids());

        BigDecimal centerOfAttack = centroidToAttack.bigMean();
        BigDecimal centerOfRightNeighbor = rightNeighbor.bigMean();

        double weightToRight;
        double weightToLeft;
        int its = 0;
        int weightGoal;
        int currentDeficit;
        int weightOfAttacked;

        int currentDeficitRight;
        int weightOfRightNeighbor;

        Collections.sort(data);
        BigDecimal nextStreamValue = nextValue(centerOfAttack, data);

        BigDecimal previous_c;
        BigDecimal previous_y;

        BigDecimal badPoint;

        double maximalError = 0;
        int indexError = 0;

        while (true) {
            its++;

            previous_c = centerOfAttack;
            previous_y = nextStreamValue;

            centerOfAttack = centroidToAttack.bigMean();
            weightOfAttacked = centroidToAttack.count();

            centerOfRightNeighbor = rightNeighbor.bigMean();
            weightOfRightNeighbor = rightNeighbor.count();

            weightToRight = 0;
            for (Centroid centroid : digest.centroids()) {
                if (centroid.bigMean().compareTo(centerOfRightNeighbor) > 0) {
                    weightToRight += centroid.count();
                }
            }

            weightToLeft = digest.size() - weightOfAttacked - weightOfRightNeighbor - weightToRight;

            Collections.sort(data);
            nextStreamValue = nextValue(centerOfAttack, data);
            countCentroids(centerOfAttack, nextStreamValue, digest.centroids(), "center to nextValue (before add)");
            if (!(centerOfAttack.compareTo(nextStreamValue) < 0)) {
                break;
            }

            if (its > 1) {
                if (previous_c.compareTo(centerOfAttack) > 0) {
                    System.out.println(String
                        .format("previous center: %f\n current: %f\n", previous_c, centerOfAttack));
                    throw new Exception("wrongly ordered, you probably ran out of precision");
                }
                if (nextStreamValue.compareTo(previous_y) > 0) {
                    System.out.println(String
                        .format("current next: %f\nprevious: %f\n", nextStreamValue, previous_y));
                    throw new Exception("wrongly ordered, you probably ran out of precision");
                }
            }

            // weight of the centroid we will fabricate
            // this is the formula for K_0
            // we also maintain centroid to the right of the attack
            weightGoal = (int) Math.floor((weightToLeft + weightToRight) / ((delta / 2d) - 3d));

            System.out.println(
                "centroids exceeding goal at start: " + centroidsExceedingCount(digest.centroids(),
                    weightGoal));

            currentDeficit = weightGoal - weightOfAttacked;
            //assert currentDeficit >= 0;

            double EPS = .001;
            double EPS_2 = .01;

            // fill up the old one
            BigDecimal point = convexCombination(nextStreamValue, centerOfAttack, EPSILON);
            for (int v = 0; v < currentDeficit; v++) {
                //double point = centerOfAttack + (nextStreamValue - centerOfAttack) / 10000; // * onePlus; //no Plus Epsilon
                BigDecimal alpha = BigDecimal.valueOf(((double) v) / currentDeficit);
                BigDecimal point_ = convexCombination(point, centerOfAttack, alpha);
                digest.add(point_, 1);
                data.add(point_);
            }
            digest.compress();

            // make the new one
            BigDecimal anotherPoint = convexCombination(nextStreamValue, centerOfAttack,
                EPSILON.multiply(BigDecimal.valueOf(2d)));
            digest.add(anotherPoint, 1);
            data.add(anotherPoint);

            BigDecimal leftEdge = convexCombination(centerOfAttack, nextStreamValue, EPSILON);
            for (int v = 0; v < weightGoal - 1; v++) {
                BigDecimal alpha = BigDecimal.valueOf(((double) v) / (weightGoal - 1));
                BigDecimal point_ = convexCombination(leftEdge, nextStreamValue, alpha);
                digest.add(point_, 1);
                data.add(point_);
            }

//            digest.compress();


            // fill up the centroid to the right
            currentDeficitRight = weightGoal - rightNeighbor.count();
           // assert currentDeficit >= 0;
            //Centroid rightCentroid = aboveValue(leftEdge, digest.centroids());
            //int deficit = weightGoal - rightCentroid.count();
            BigDecimal rightCentroidVal = rightNeighbor.bigMean();
            for (int pp = 0; pp < currentDeficitRight; pp++) {
                digest.add(rightCentroidVal);
                data.add(rightCentroidVal);
            }
//
//            System.out.println(String
//                .format("%s\n\n%s\n\n%f\n\n%f\n\n%f\n\n", centerOfAttack, point, anotherPoint, leftEdge,
//                    nextStreamValue));

            digest.compress();

            printComparison(centerOfAttack, point);
            printComparison(point, anotherPoint);
            printComparison(anotherPoint, leftEdge);
            printComparison(leftEdge, nextStreamValue);

//            for (Centroid c : digest.centroids()) {
//                System.out.println(c.count() + "; ");
//            }


            System.out.println(
                "centroids exceeding goal at end: " + centroidsExceedingCount(digest.centroids(),
                    weightGoal - 1));

            int a = countCentroids(centerOfAttack, point, digest.centroids(), "center to point");
            int U = countCentroids(centerOfAttack, leftEdge, digest.centroids(), "center to leftEdge");
            int b = countCentroids(anotherPoint, leftEdge, digest.centroids(), "another to left");
            int c = countCentroids(leftEdge, rightCentroidVal, digest.centroids(), "left to right");
            int t = countCentroids(centerOfAttack, rightCentroidVal, digest.centroids(), "center to right");
            int d = countCentroids(rightCentroidVal, rightCentroidVal.multiply(BigDecimal.TEN),
                digest.centroids(), "right to 10*right, whatevs");

            int maxCount = 0;
            for (Centroid centroid : digest.centroids()) {
                if (centroid.count() > maxCount) {
                    maxCount = centroid.count();
                }
            }
            System.out.println("weight goal was: " + weightGoal + "; max weight is: " + maxCount);

            System.out.println("oldCenter to point (expect 1): " + a);
            System.out.println("oldCenter to leftEdge (expect 2?): " + U);
            System.out.println("another to leftEdge (expect 1): " + b);
            System.out.println("leftTo to oldRight (expect 1):" + c);
            System.out.println("oldCenter to oldRight (total, 3?):" + t);
            System.out.println("to the right: " + d);

            badPoint = convexCombination(centerOfAttack, nextStreamValue, EPSILON);
            Collections.sort(data);

            System.out.println("finished iteration: " + its);
            System.out.println("td " + digest.cdf(badPoint));
            System.out
                .println("truth " + countBelow(badPoint, data) / (double) data.size());

            double error = Math
                .abs(digest.cdf(badPoint) - countBelow(badPoint, data) / (double) data.size());

            if (error > maximalError) {
                maximalError = error;
                indexError = its;
            }
            System.out.println("maximal error so far: " + maximalError + " on " + indexError);
            System.out.println("num centroids: " + digest.centroids().size() + "\n");

            // hmmm ....
            //centroidToAttack = aboveValue(convexCombination(anotherPoint, leftEdge, 0.99999),
            //    digest.centroids());
            //centroidToAttack = belowValue(nextStreamValue, digest.centroids());
            centroidToAttack = largestCentroidBetweenValues(centerOfAttack, nextStreamValue,
                digest.centroids());


            System.out.println("size of centroidTOAttack next: " + centroidToAttack.count());
            if (nextValue(centroidToAttack.bigMean(), data).compareTo(nextStreamValue) > 0) {
                System.out.println("odd next val orderings");
                centroidToAttack = belowValue(
                    centroidToAttack.bigMean().subtract(BigDecimal.valueOf(Double.MIN_VALUE)),
                    digest.centroids());
            }

            if (centroidToAttack.bigMean().compareTo(centerOfAttack) < 0) {
                System.out.println("wtf");
            }
            rightNeighbor = aboveValue(centroidToAttack.bigMean(), digest.centroids());

            detectDupes(digest.centroids());
            System.out.println("done with this; NEXT\n\n\n\n");
        }



        badPoint = convexCombination(previous_c, previous_y,
            0.5); //previous_c + previous_y) / 2d;
        Collections.sort(data);

        double error = Math
            .abs(digest.cdf(badPoint) - countBelow(badPoint, data) / (double) data.size());

        if (error > maximalError) {
            maximalError = error;
            indexError = its;
        }
        System.out.println("maximal error so far: " + maximalError + " on " + indexError);

        System.out.println("td " + digest.cdf(badPoint));
        System.out.println("truth " + countBelow(badPoint, data) / (double) data.size());
        System.out.println("iterations" + its);
    }

    private void detectDupes(Collection<Centroid> centroids) {
        Iterator<Centroid> centroidIterator = centroids.iterator();
        Centroid centroid = centroidIterator.next();
        Centroid previous = null;
        while (centroidIterator.hasNext()) {
            if (previous != null) {
                if (previous.bigMean() == centroid.bigMean()) {
                    System.out.println("dupe!");
                }
            }
            previous = centroid;
            centroid = centroidIterator.next();
        }
    }


    private Centroid largestCentroidBetweenValues(BigDecimal lowerBound, BigDecimal upperBound,
        Collection<Centroid> centroids) {
        Centroid centroid = null;
        BigDecimal largest = BigDecimal.valueOf(Double.MAX_VALUE * -1);
        for (Centroid c : centroids) {
            if (c.bigMean().compareTo(lowerBound) >= 0 && c.bigMean().compareTo(upperBound) <= 0) {
                if (c.bigMean().compareTo(largest) >= 0) {
                    centroid = c;
                    largest = c.bigMean();
                }
            }
        }
        return centroid;
    }

    private void printComparison(BigDecimal v1, BigDecimal v2) throws Exception {
        if (v1.compareTo(v2) > 0) {
            System.out.println(String
                .format("should be smaller: %f\n\n should be larger: %f\n\n", v1, v2));
            throw new Exception("wrongly ordered, you probably ran out of precision");
        }
    }

    private BigDecimal convexCombination(BigDecimal val1, BigDecimal val2, BigDecimal alpha) {
        return val1.multiply(alpha).add(val2.multiply(BigDecimal.ONE.subtract(alpha)));
    }

    private BigDecimal convexCombination(BigDecimal val1, BigDecimal val2, double alphaDouble) {
        BigDecimal alpha = BigDecimal.valueOf(alphaDouble);
        return val1.multiply(alpha).add(val2.multiply(BigDecimal.ONE.subtract(alpha)));
    }

    private int centroidsExceedingCount(Collection<Centroid> centroids, int countThreshold) {
        int j = 0;
        for (Centroid c : centroids) {
            if (c.count() >= countThreshold) {
                j++;
            }
        }
        return j;
    }

    private int countCentroids(BigDecimal lowerBound, BigDecimal upperBound,
        Collection<Centroid> centroids, String out) {
        int count = 0;
        for (Centroid c : centroids) {
            if (c.bigMean().compareTo(lowerBound) >= 0 && c.bigMean().compareTo(upperBound) <= 0) {
                count++;
                System.out.println(out + "; count: " + c.count());
            }
        }
        System.out.println("\n");
        return count;
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
