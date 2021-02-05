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

import java.lang.Math;
import java.util.*;
import java.io.*;

import org.junit.Ignore;

//import com.tdunning.math.stats.SpeedComparison.getProperty;

/**
 *
 */
public class CarefulAttack {

    private int weightGoal(ScaleFunction scaleFunction, double delta, double weightToLeft,
        double weightToRight, long size, double deltaMult) {
        if (scaleFunction == ScaleFunction.K_0) {
            return (int) Math.ceil((weightToLeft + weightToRight) / ((delta * deltaMult) - 3d));
        } else if (scaleFunction == ScaleFunction.K_3) {
            double norm = 4 * Math.log(size / delta) + 21;
            //double norm = 10d;
            return (int) ((Math.exp(norm / delta) - 1d) * Math
                .min(weightToRight, weightToLeft));
        } else if (scaleFunction == ScaleFunction.K_3_NO_NORM) {
            //double norm = 4 * Math.log(1d / delta);
            double norm = 1d;
            return (int) ((Math.exp(norm / delta) - 1d) * Math
                .min(weightToRight, weightToLeft));
        } else {
            return 0;
        }
    }

    //@ Test  - WIP
    public List<Double> carefulNestedAroundZeroK_3() throws Exception {
        return carefulNestedAroundZero(ScaleFunction.K_3_NO_NORM, 100, "tree",
            true, 44, true,
            false, new NestedInputParams(0.001, 0.1, 0.0001, true, 0.5, false, 0.2, 10),
            RandomUtils.getRandom().nextLong(), false, false, "reqsketch");
    }

    //@Test
    public List<Double> carefulNestedAroundZeroK_0() throws Exception {
        return carefulNestedAroundZero(ScaleFunction.K_0, 500, "merging",  // merging or tree
            false, 1000, true,
            false, new NestedInputParams(0.000000001, 0.2, 0.0000000001, false, 0.5, false, 0.2, 10),
            RandomUtils.getRandom().nextLong(), false, true, "reqsketch");
    }


    @Test // resulting data files used by error curves notebook
    public void writeCarefulNestedAroundZeroK_0() throws Exception {
        carefulNestedAroundZero(ScaleFunction.K_0, 500, "merging",
            false, 1000, true,
            true, new NestedInputParams(0.0000001, 0.21, 0.0000000001, false, 0.5, false, 0.2, 10),
            981198271346L, true, true, "kll");
        carefulNestedAroundZero(ScaleFunction.K_0, 500, "merging",
            true, 1000, true,
            true, new NestedInputParams(0.0000001, 0.21, 0.0000000001, false, 0.5, false, 0.2, 10),
            981198271346L, true, true, "kll");
        carefulNestedAroundZero(ScaleFunction.K_0, 500, "tree",
            false, 1300, true,
            true, new NestedInputParams(0.000000001, 0.26, 0.0000000001, false, 1/1.48d, true, 0.18, 8),
            981198271346L, true, false, "kll");
    }
    
    int reqK, kllK, NumberOfPoints;
    String DigestStatsDir, FileSuffix;
    
    public CarefulAttack(String configFile) throws Exception {
        System.out.println("processing config file: " + configFile);
        Properties prop = new Properties();
        FileInputStream instream = new FileInputStream(configFile);
        prop.load(instream);

        // load properties
        String DigestImpl = SpeedComparison.getProperty(prop, "DigestImpl");
        NumberOfPoints = Integer.parseInt(
            SpeedComparison.getProperty(prop, "NumberOfPoints")); // number of points where we probe the rank estimates // TODO not used now
        boolean WriteCentroidData = Boolean.parseBoolean(SpeedComparison.getProperty(prop, "WriteCentroidData"));
        boolean useAlternatingSort = Boolean.parseBoolean(SpeedComparison.getProperty(prop, "useAlternatingSort"));
        int Compression = Integer.parseInt(SpeedComparison.getProperty(prop, "Compression")); // delta for t-digest
        ScaleFunction scale = ScaleFunction.valueOf(SpeedComparison.getProperty(prop, "ScaleFunction")); // ScaleFunction for t-digest
        DigestStatsDir = SpeedComparison.getProperty(prop, "DigestStatsDir");
        FileSuffix = SpeedComparison.getProperty(prop, "FileSuffix");
        reqK = Integer.parseInt(SpeedComparison.getProperty(prop, "ReqK"));
        kllK = Integer.parseInt(SpeedComparison.getProperty(prop, "KllK"));
        int iterations = Integer.parseInt(SpeedComparison.getProperty(prop, "iterations"));
        long seed = Long.parseLong(SpeedComparison.getProperty(prop, "seed"));

        carefulNestedAroundZero(scale, Compression, DigestImpl,
            useAlternatingSort, iterations, true,
            WriteCentroidData, new NestedInputParams(Double.parseDouble(SpeedComparison.getProperty(prop, "newCentroidNextMinusCenterCoeff")), 
                                                     Double.parseDouble(SpeedComparison.getProperty(prop, "newCentroidNextMultiplier")),
                                                     Double.parseDouble(SpeedComparison.getProperty(prop, "rightCentroidNudge")),
                                                     Boolean.parseBoolean(SpeedComparison.getProperty(prop, "useConvexCombination")),
                                                     Double.parseDouble(SpeedComparison.getProperty(prop, "deltaMult")),
                                                     Boolean.parseBoolean(SpeedComparison.getProperty(prop, "compressBeforePositiveUpdates")),
                                                     Double.parseDouble(SpeedComparison.getProperty(prop, "fracNegativeUpdates")),
                                                     Double.parseDouble(SpeedComparison.getProperty(prop, "initBatchSizeMult"))),
            seed, Boolean.parseBoolean(SpeedComparison.getProperty(prop, "CompareToSorted")),
                  Boolean.parseBoolean(SpeedComparison.getProperty(prop, "maintainRightCentroid")),
                  SpeedComparison.getProperty(prop, "compareTo"));
    }

    private class NestedInputParams {

        private final double newCentroidNextMinusCenterCoeff;
        private final double newCentroidNextMultiplier;
        private final double rightCentroidNudge;
        private final boolean useConvexCombination;
        private final double deltaMult;
        private final boolean compressBeforePositiveUpdates;
        private final double fracNegativeUpdates;
        private final double initBatchSizeMult;

        private NestedInputParams(double p1, double p2, double p3, boolean p4, double p5, boolean p6, double p7, double p8) {
            newCentroidNextMinusCenterCoeff = p1;
            newCentroidNextMultiplier = p2;
            rightCentroidNudge = p3;
            useConvexCombination = p4;
            deltaMult = p5;
            compressBeforePositiveUpdates = p6;
            fracNegativeUpdates = p7;
            initBatchSizeMult = p8;
        }
    }

    /*
    iterations - will eventually run out of heap space, this caps it
     */
    public List<Double> carefulNestedAroundZero(ScaleFunction scaleFunction, double delta,
        String implementation, boolean useAlternatingSort, int iterations, boolean writeResults,
        boolean writeCentroidData, NestedInputParams params, long seed, boolean compareToSorted,
        boolean maintainRightCentroid, String compareTo) throws Exception {

        List<Double> errors = new ArrayList<>();

        List<Double> data = new ArrayList<>();
        List<Double> sortedData = new ArrayList<>();
        TDigest digest = CarefulAttackAuxiliary.digest(delta, implementation, scaleFunction, useAlternatingSort, seed);
        System.out.println(
                "delta:\t" + delta);

        int initializingHalfBatchSize = (int) Math.floor(delta * params.initBatchSizeMult);

        double denom = 100000000d;
        double infty = (Double.MAX_VALUE) / denom; // so we can safely average

        double increment = infty / initializingHalfBatchSize;
        System.out.println(
            "infty:\t" + infty + "\ninitializingHalfBatchSize: \t" + initializingHalfBatchSize);

        // try to place the attack in the right tail
        if (!(scaleFunction == ScaleFunction.K_0)) {
            int init = 0;
            while (init < 250 * initializingHalfBatchSize) {
                double p = -infty * (2d - ((double) init / 250 / initializingHalfBatchSize));
                digest.add(p);
                data.add(p);
                sortedData.add(p);
                init++;
            }
        }

        for (int i = 0; i < initializingHalfBatchSize; i++) {
            double p = -infty * (1d - ((double) i / 2 / initializingHalfBatchSize));
            sortedData.add(p);
            digest.add(p);
            data.add(p);

            double q = infty * (1d - ((double) i / 2 / initializingHalfBatchSize));
            sortedData.add(q);
            digest.add(q);
            data.add(q);
        }

        if (!(scaleFunction == ScaleFunction.K_0)) {
            long sz = digest.size();
            double val = digest.min - 1d;
            for (long i = 0; i < sz * 8; i++) {
                digest.add(val);
                data.add(val);
                sortedData.add(val);
            }
        }

        digest.compress();

        Centroid centroidToAttack = belowValue(0d, digest.centroids());
        Centroid rightNeighbor = aboveValue(centroidToAttack.mean(), digest.centroids());

        double centerOfAttack = centroidToAttack.mean();
        double centerOfRightNeighbor = rightNeighbor.mean();

        double weightToRight;
        double weightToLeft;
        int iteration = 0;
        int weightGoal;
        int currentDeficit;
        int weightOfAttacked;

        int currentDeficitRight;
        int weightOfRightNeighbor;

        Collections.sort(sortedData);

        //  1> maximal error so far: 0.3265282585631959 on 1244
        //detectDupes(data, "initial batch");
        double nextStreamValue = nextValue(centerOfAttack, sortedData);

        double previous_c;
        double previous_y;

        double maximalError = 0;
        int indexError = 0;

        while (iteration < iterations) {
            iteration++;

            previous_c = centerOfAttack;
            previous_y = nextStreamValue;

            centerOfAttack = centroidToAttack.mean();
            weightOfAttacked = centroidToAttack.count();

            centerOfRightNeighbor = rightNeighbor.mean();
            weightOfRightNeighbor = rightNeighbor.count();

            weightToRight = 0;
            for (Centroid centroid : digest.centroids()) {
                if (centroid.mean() > centerOfRightNeighbor) {
                    weightToRight += centroid.count();
                }
            }

            weightToLeft = digest.size() - weightOfAttacked - weightOfRightNeighbor - weightToRight;

            Collections.sort(sortedData);
            nextStreamValue = nextValue(0, sortedData);   //centerOfAttack
            //nextStreamValue = nextValue(centerOfAttack, data); //<-- doesn't work very well

            if (nextStreamValue < 100 * Double.MIN_VALUE) {
                System.out.println(String
                    .format("too small nextStreamValue:\t%s", nextStreamValue));
                break;
            }

            if (!(centerOfAttack < nextStreamValue)) {
                System.out.println(String
                    .format(
                        "centerOfAttack < nextStreamValue: nextStreamValue:\t%s\ncenterOfAttack:\t%s\n",
                        nextStreamValue, centerOfAttack));
                break;
            }

            if (iteration > 1) {
                if ((previous_c > centerOfAttack)) {
                    System.out.println(String
                        .format("previous_c:\t%s\ncenterOfAttack:\t%s\n", previous_c,
                            centerOfAttack));
//                    throw new Exception(
//                        "previous_c > centerOfAttack: wrongly ordered, you probably ran out of precision");
                }
                if (nextStreamValue > previous_y) {
                    System.out.println(String
                        .format("nextStreamValue:\t%s\nprevious_y:\t%s\n", nextStreamValue,
                            previous_y));
//                    throw new Exception(
//                        "nextStreamValue > previous_y: wrongly ordered, you probably ran out of precision");
                }
            }

            // weight of the centroid we will fabricate
            System.out.println("centroids count: " + digest.centroids().size());
            weightGoal = weightGoal(scaleFunction, delta, weightToLeft, weightToRight,
                digest.size(), params.deltaMult);
            //(int) Math.ceil((weightToLeft + weightToRight) / ((delta / 2d) - 3d));
            System.out.println("weightGoal: " + weightGoal);

            // this information is relevant mainly for K_0, where the behavior is uniform
            System.out.println(
                "centroids exceeding goal at start: " + centroidsExceedingCount(digest.centroids(),
                    weightGoal));
            currentDeficit = weightGoal - weightOfAttacked;
            // if (its > 2) {assert currentDeficit >= 0;}

            // fill up the old one
            //double point = centerOfAttack; // * (1d - EPS) + nextStreamValue * EPS;
            for (int v = 0; v < currentDeficit; v++) {
                //   / 1000000); // * onePlus; //no Plus Epsilon
                double point = centerOfAttack;
                digest.add(point, 1);
                data.add(point);
                sortedData.add(point);
            }

            // make the new one
            int v = 0;

            double anotherPoint = (centerOfAttack + params.newCentroidNextMinusCenterCoeff * (
                nextStreamValue - centerOfAttack));

            //double anotherPoint = centerOfAttack * 0.1;

            double newCentroidMainValue = nextStreamValue * params.newCentroidNextMultiplier;

            //centerOfAttack * (1d - EPS_2) + nextStreamValue * EPS_2; nextStreamValue * 0.1; //centerOfAttack * EPS +

            for (; v < weightGoal * params.fracNegativeUpdates; v++) {
                double toAdd = anotherPoint; //+ (newCentroidMainValue - anotherPoint) * v  / (double) weightGoal / 10000;
                if (params.useConvexCombination) {
                    toAdd += params.newCentroidNextMinusCenterCoeff * anotherPoint;
                }
                digest.add(toAdd);
                data.add(toAdd);
                sortedData.add(toAdd);
            }
            if (params.compressBeforePositiveUpdates)
                digest.compress();

            for (; v < weightGoal; v++) {
                double toAdd = newCentroidMainValue; //+ (anotherPoint - newCentroidMainValue) * v / (double) weightGoal / 10000;
                if (params.useConvexCombination) {
                    toAdd -= params.newCentroidNextMinusCenterCoeff * toAdd;
                }
                digest.add(toAdd, 1);
                data.add(toAdd);
                sortedData.add(toAdd);
            }

            if (maintainRightCentroid) {
                // we also maintain centroid to the right of the attack
                currentDeficitRight = weightGoal - rightNeighbor.count();
                //assert currentDeficit >= 0;

                //Centroid rightCentroid = aboveValue(newCentroidMainValue, digest.centroids());
                //int deficit = weightGoal - rightCentroid.count();
                double rightCentroidVal = rightNeighbor.mean();
                for (int pp = 0; pp < currentDeficitRight; pp++) {
                    rightCentroidVal = rightCentroidVal + params.rightCentroidNudge;
                    digest.add(rightCentroidVal);
                    data.add(rightCentroidVal);
                    sortedData.add(rightCentroidVal);
                }
            }
            System.out.println(String
                .format(
                    "centerOfAttack:\t%s\npoint:\t\t%s\nanotherPoint:\t%s\nnewCentroidMainValue:  \t%s\nnextStreamValue:\t%s",
                    centerOfAttack, centerOfAttack, anotherPoint, newCentroidMainValue,
                    nextStreamValue));

            digest.compress();

            System.out.println(
                "centroids exceeding goal at end: " + centroidsExceedingCount(digest.centroids(),
                    weightGoal - 1));

            double bad_point = 0;
            //double bad_point = newCentroidMainValue;
            Centroid belowZeroC = belowValue(0, digest.centroids());
            Centroid aboveZeroC = aboveValue(0, digest.centroids());
            //double bad_point = belowZeroC.mean();
            Collections.sort(sortedData);

            // detectDupes(data, "iteration " + its);
            // detectDupes(digest.centroids());

            System.out.println("finished iteration: " + iteration);
            System.out.println("belowZeroC mean: " + belowZeroC.mean());
            System.out.println("aboveZeroC mean: " + aboveZeroC.mean());
            //System.out.println("bad point: " + bad_point);

            if (iteration > 1) {
                System.out.println("td " + digest.cdf(bad_point));
                double truth = countBelow(bad_point, sortedData) / (double) data.size();
                System.out
                    .println("truth " + truth);
                double error = Math.abs(digest.cdf(bad_point) - truth);

                errors.add(error);

                if (error > maximalError) {
                    maximalError = error;
                    indexError = iteration;
                }
            }
            System.out.println("maximal error so far: " + maximalError + " on " + indexError);
            System.out.println("num centroids: " + digest.centroids().size() + "\n");

            //double yetAnotherPoint = centerOfAttack * 0.9 + nextStreamValue * ;
            centroidToAttack = belowZeroC; //aboveValue(anotherPoint, digest.centroids()); //centerOfAttack //belowZeroC; //
            // ^ works better for AVL?!

            /////centroidToAttack = belowValue(0, digest.centroids()); //
            //
            // centerOfAttack

            //centroidToAttack = belowValue(newCentroidMainValue, digest.centroids()); //centerOfAttack
            //            centroidToAttack = aboveValue(centerOfAttack, digest.centroids()); //centerOfAttack

            // the good one for Merging?
            //centroidToAttack = aboveValue(anotherPoint, digest.centroids()); //centerOfAttack

            if (centroidToAttack.mean() < centerOfAttack) {
                System.out.println("wtf: centroidToAttack.mean()=" + centroidToAttack.mean()
                    + " < centerOfAttack");
            }
            //if (centroidToAttack.mean() > 0) {
            //  System.out.println("above zero: centroidToAttack.mean()=" + centroidToAttack.mean());
            //  centroidToAttack = belowZeroC;
            //}
            rightNeighbor = aboveValue(centroidToAttack.mean(), digest.centroids());
        }

        double bad_point = 0; // (previous_c + previous_y) / 2d;
        Collections.sort(sortedData);

        //System.out.println("td " + digest.cdf(bad_point));
        System.out.println("truth " + countBelow(bad_point, sortedData) / (double) data.size());
        System.out.println("iterations" + iteration);

        String altSort =
            (implementation.equals("merging")) ? "_alt_" + String.valueOf(useAlternatingSort) : "";
        if (writeResults) {
            CarefulAttackAuxiliary.writeResults((int) delta, data.size(), digest, data, sortedData, compareTo,
                DigestStatsDir,
                DigestStatsDir + String.format(
                    "careful_iterations=%d_samples=%d_scalefunc=%s_delta=%d_centroids=%d_sizeBytes=%d_impl=%s",
                    iteration,
                    data.size(), digest.scale.toString(), (int) delta, digest.centroidCount(),
                    digest.byteSize(), implementation) + altSort + FileSuffix, false);
        }
        if (writeCentroidData) {
            CarefulAttackAuxiliary.writeCentroidData(digest,
                DigestStatsDir + String.format(
                    "centroids_careful_iterations=%d_samples=%d_scalefunc=%s_delta=%d_centroids=%d_sizeBytes=%d_impl=%s",
                    iteration,
                    data.size(), digest.scale.toString(), (int) delta, digest.centroidCount(),
                    digest.byteSize(), implementation) + altSort + FileSuffix);
        }
        if (compareToSorted) {
            TDigest tDigest = CarefulAttackAuxiliary.rerunOnSortedInput(digest, sortedData);
            if (writeResults) {
                CarefulAttackAuxiliary.writeResults((int) delta, data.size(), tDigest, data, sortedData, compareTo,
                    DigestStatsDir,
                    DigestStatsDir + String.format(
                        "careful_iterations=%d_samples=%d_scalefunc=%s_delta=%d_centroids=%d_sizeBytes=%d_impl=%s",
                        iteration,
                        data.size(), tDigest.scale.toString(), (int) delta, tDigest.centroidCount(),
                        tDigest.byteSize(), implementation) + altSort + "_sorted" + FileSuffix,
                    false);
            }
            if (writeCentroidData) {
                CarefulAttackAuxiliary.writeCentroidData(tDigest,
                    DigestStatsDir + String.format(
                        "centroids_careful_iterations=%d_samples=%d_scalefunc=%s_delta=%d_centroids=%d_sizeBytes=%d_impl=%s",
                        iteration,
                        data.size(), tDigest.scale.toString(), (int) delta, tDigest.centroidCount(),
                        tDigest.byteSize(), implementation) + altSort + "_sorted" +  FileSuffix);
            }
        }
        return errors;
    }

    private void detectDupes(List<Double> inputs, String context) {
        for (int index = 0; index < inputs.size() - 1; index++) {
            if (inputs.get(index).equals(inputs.get(index + 1))) {
                System.out.println("input dupe detected: " + context);
            }
        }
    }

    // assume sorted...
    private void detectDupes(Collection<Centroid> centroids) {
        Iterator<Centroid> centroidIterator = centroids.iterator();
        Centroid centroid = centroidIterator.next();
        Centroid previous = null;
        while (centroidIterator.hasNext()) {
            if (previous != null) {
                if (centroid.mean() == previous.mean()) {
                    System.out.println("centroid dupe detected");
                }
            }
            previous = centroid;
            centroid = centroidIterator.next();
        }
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

    /*
    Find the smallest centroid above the given threshold
     */
    private Centroid aboveValue(double c, Collection<Centroid> centroids) throws Exception {
        Iterator<Centroid> centroidIterator = centroids.iterator();
        Centroid centroid = centroidIterator.next();
        Centroid smallestCentroid = null;
        double smallestMean = Double.POSITIVE_INFINITY;
        while (centroidIterator.hasNext()) {
            if (c < centroid.mean() && centroid.mean() < smallestMean) {
                smallestCentroid = centroid;
                smallestMean = centroid.mean();
            }
            centroid = centroidIterator.next();
        }
        if (smallestCentroid != null) {
            return smallestCentroid;
        } else {
            throw new Exception("couldn't find a centroid above threshold");
        }
    }

    /*
    Find the largest centroid below the given threshold
     */
    private Centroid belowValue(double c, Collection<Centroid> centroids) throws Exception {
        Iterator<Centroid> centroidIterator = centroids.iterator();
        Centroid centroid = centroidIterator.next();
        Centroid largestCentroid = null;
        double largestMean = Double.NEGATIVE_INFINITY;
        while (centroidIterator.hasNext()) {
            if (c > centroid.mean() && centroid.mean() > largestMean) {
                largestCentroid = centroid;
                largestMean = centroid.mean();
            }
            centroid = centroidIterator.next();
        }
        if (largestCentroid != null) {
            return largestCentroid;
        } else {
            throw new Exception("couldn't find a centroid below threshold");
        }
    }

//
//    // assume centroids in ascending order
//    private Centroid belowValue(double c, Collection<Centroid> centroids) throws Exception {
//        Iterator<Centroid> centroidIterator = centroids.iterator();
//        Centroid centroid = centroidIterator.next();
//        Centroid previous = null;
//        while (centroidIterator.hasNext()) {
//            if (centroid.mean() >= c) {
//                return previous;
//                // centroid = centroidIterator.next();
//            } else {
//                previous = centroid;
//                centroid = centroidIterator.next();
//            }
//        }
//        throw new Exception("couldn't find a centroid there");
//    }

    // assume data is sorted
    private double nextValue(double c, List<Double> sortedData) {
        int index = 0;
        while (c >= sortedData.get(index)) {
            index++;
        }
        assert sortedData.get(index) > c;
        assert index == 0 || sortedData.get(index - 1) <= c;
        return sortedData.get(index);
    }

    private int countBelow(double c, List<Double> sortedData) {
        int index = 0;
        while (c > sortedData.get(index)) {
            index++;
        }
        // assert sortedData.get(index) > c;
        return index;
    }

    @SuppressWarnings("unused")
    public static void main(final String[] args) throws Exception {
        for (int j = 0; j < args.length; j++) {
            new CarefulAttack(args[j]);
        }
    }

}


