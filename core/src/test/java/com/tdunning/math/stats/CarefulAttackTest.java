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

import com.google.common.collect.Lists;

import org.apache.mahout.common.RandomUtils;
import org.apache.mahout.math.jet.random.AbstractContinousDistribution;
import org.apache.mahout.math.jet.random.Gamma;
import org.apache.mahout.math.jet.random.Normal;
import org.apache.mahout.math.jet.random.Uniform;
import org.junit.*;

import java.io.*;
import java.lang.Math;
import java.nio.ByteBuffer;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import org.junit.Ignore;


/**
 *
 */
//@Ignore
public class CarefulAttackTest extends AbstractTest {
    
    @BeforeClass
    public static void freezeSeed() {
        RandomUtils.useTestSeed();
    }


    @Ignore
    public void carefulNested() throws Exception {

        double EPSILON = Double.MIN_VALUE;

        double delta = 500;
        //double delta = 1000;

        List<Double> data = new ArrayList<>();
        MergingDigest digest = new MergingDigest(delta);
        digest.setScaleFunction(ScaleFunction.K_0);
        digest.setUseAlternatingSort(false);

        int initializingHalfBatchSize = (int) Math.floor(delta * 10);

        // delta=500, denom=100000 seems to work okay..

        double denom = 10000000d;
        //denom *= 100;
        double infty = (Double.MAX_VALUE) / denom; // so we can safely average

        double increment = infty / initializingHalfBatchSize;

        // add a glob of zeroes, aim for error there.
        for (int j = 0; j < 5 * initializingHalfBatchSize; j++) {
            data.add(0d);
            digest.add(0d);
        }

        for (int i = 0; i < 2 * initializingHalfBatchSize; i++) {
            double point = -infty;
            for (int j = 0; j < i; j++) {
                point += increment;
            }
            data.add(point);
            digest.add(point);
        }
        digest.compress();

        Centroid centroidToAttack = belowValue(0d, digest.centroids());
        Centroid rightNeighbor = aboveValue(centroidToAttack.mean(), digest.centroids());

        double centerOfAttack = centroidToAttack.mean();
        double centerOfRightNeighbor = rightNeighbor.mean();

        double weightToRight;
        double weightToLeft;
        int its = 0;
        int weightGoal;
        int currentDeficit;
        int weightOfAttacked;

        int currentDeficitRight;
        int weightOfRightNeighbor;

        Collections.sort(data);
        double nextStreamValue = nextValue(centerOfAttack, data);

        double previous_c;
        double previous_y;

        double maximalError = 0;
        int indexError = 0;

        while (true) {
            its++;

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

            Collections.sort(data);
            nextStreamValue = nextValue(centerOfAttack, data);
            if (!(centerOfAttack < nextStreamValue)) {
                break;
            }

            if (its > 1) {
                if ((previous_c > centerOfAttack)) {
                    System.out.println(String
                        .format("%f\n%f\n", previous_c, centerOfAttack));
                    throw new Exception("wrongly ordered, you probably ran out of precision");
                }
                if (nextStreamValue > previous_y) {
                    System.out.println(String
                        .format("%f\n%f\n", nextStreamValue, previous_y));
                    throw new Exception("wrongly ordered, you probably ran out of precision");
                }
            }

            // weight of the centroid we will fabricate
            // this is the formula for K_0
            // we also maintain centroid to the right of the attack
            weightGoal = (int) Math.ceil((weightToLeft + weightToRight) / ((delta / 2d) - 3d));

            System.out.println(
                "centroids exceeding goal at start: " + centroidsExceedingCount(digest.centroids(),
                    weightGoal));
            currentDeficit = weightGoal - weightOfAttacked;
            assert currentDeficit >= 0;

            double EPS = .001;
            double EPS_2 = .01;

            // fill up the old one
            double point = centerOfAttack * (1d - EPS) + nextStreamValue * EPS;
            for (int v = 0; v < currentDeficit; v++) {
                //double point = centerOfAttack + (nextStreamValue - centerOfAttack) / 10000; // * onePlus; //no Plus Epsilon
                digest.add(point, 1);
                data.add(point);
            }

            // make the new one
            double anotherPoint = centerOfAttack * (1d - EPS_2) + nextStreamValue * EPS_2;
            digest.add(anotherPoint, 3);
            data.add(anotherPoint);

            double leftEdge = centerOfAttack * EPS + nextStreamValue * (1 - EPS);
            for (int v = 0; v < weightGoal - 3; v++) {
                digest.add(leftEdge, 1);
                data.add(leftEdge);
            }

            // fill up the centroid to the right
            currentDeficitRight = weightGoal - rightNeighbor.count();
            assert currentDeficit >= 0;
            //Centroid rightCentroid = aboveValue(leftEdge, digest.centroids());
            //int deficit = weightGoal - rightCentroid.count();
            double rightCentroidVal = rightNeighbor.mean();
            for (int pp = 0; pp < currentDeficitRight; pp++) {
                digest.add(rightCentroidVal);
                data.add(rightCentroidVal);
            }

            System.out.println(String
                .format("%s\n%s\n%f\n%f\n%f\n", centerOfAttack, point, anotherPoint, leftEdge,
                    nextStreamValue));

            digest.compress();

            System.out.println(
                "centroids exceeding goal at end: " + centroidsExceedingCount(digest.centroids(),
                    weightGoal - 1));

            double bad_point = centerOfAttack * .00000000001 + nextStreamValue * (1 - .00000000001);
            Collections.sort(data);

            System.out.println("finished iteration: " + its);
            System.out.println("td " + digest.cdf(bad_point));
            System.out
                .println("truth " + countBelow(bad_point, data) / (double) data.size());

            double error = Math
                .abs(digest.cdf(bad_point) - countBelow(bad_point, data) / (double) data.size());

            if (error > maximalError) {
                maximalError = error;
                indexError = its;
            }
            System.out.println("maximal error so far: " + maximalError + " on " + indexError);
            System.out.println("num centroids: " + digest.centroids().size() + "\n");

            centroidToAttack = aboveValue(centerOfAttack, digest.centroids());
            if (centroidToAttack.mean() < centerOfAttack) {
                System.out.println("wtf");
            }
            rightNeighbor = aboveValue(centroidToAttack.mean(), digest.centroids());
        }

        double bad_point = (previous_c + previous_y) / 2d;
        Collections.sort(data);

        System.out.println("td " + digest.cdf(bad_point));
        System.out.println("truth " + countBelow(bad_point, data) / (double) data.size());
        System.out.println("iterations" + its);
    }
    

    @Test
    public void carefulNestedAroundZero() throws Exception {

        double EPSILON = Double.MIN_VALUE;

        double delta = 500;
        //double delta = 1000;

        List<Double> data = new ArrayList<>();
        MergingDigest digest = new MergingDigest(delta);
        digest.setScaleFunction(ScaleFunction.K_0);
        digest.setUseAlternatingSort(false);

        int initializingHalfBatchSize = (int) Math.floor(delta * 10);

        // delta=500, denom=100000 seems to work okay..

        double denom = 10000000000d;
        //denom *= 100;
        double infty = (Double.MAX_VALUE) / denom; // so we can safely average

        double increment = infty / initializingHalfBatchSize;
        System.out.println("infty:\t" + infty + "\ninitializingHalfBatchSize: \t" + initializingHalfBatchSize);
        
        // add a glob of zeroes, aim for error there.
//        for (int j = 0; j < 5 * initializingHalfBatchSize; j++) {
//            data.add(0d);
//            digest.add(0d);
//        }
//
//        for (int i = 0; i < 2 * initializingHalfBatchSize; i++) {
//            double point = -infty;
//            for (int j = 0; j < i; j++) {
//                point += increment;
//            }
//            data.add(point);
//            digest.add(point);
//        }
        for (int i = 0; i < 2 * initializingHalfBatchSize; i++) {
          data.add(-infty);
          digest.add(-infty);
          data.add(infty);
          digest.add(infty);
        }
        digest.compress();

        Centroid centroidToAttack = belowValue(0d, digest.centroids());
        Centroid rightNeighbor = aboveValue(centroidToAttack.mean(), digest.centroids());

        double centerOfAttack = centroidToAttack.mean();
        double centerOfRightNeighbor = rightNeighbor.mean();

        double weightToRight;
        double weightToLeft;
        int its = 0;
        int weightGoal;
        int currentDeficit;
        int weightOfAttacked;

        int currentDeficitRight;
        int weightOfRightNeighbor;

        Collections.sort(data);
        double nextStreamValue = nextValue(centerOfAttack, data);

        double previous_c;
        double previous_y;

        double maximalError = 0;
        int indexError = 0;

        while (true) {
            its++;

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

            Collections.sort(data);
            nextStreamValue = nextValue(0, data); //centerOfAttack
            
            if (nextStreamValue < 100*Double.MIN_VALUE) {
                System.out.println(String
                  .format("too small nextStreamValue:\t%s", nextStreamValue));
                break;
            }
            
            if (!(centerOfAttack < nextStreamValue)) {
                System.out.println(String
                  .format("centerOfAttack >= nextStreamValue:\nnextStreamValue:\t%s\ncenterOfAttack:\t%s\n", nextStreamValue, centerOfAttack));
                break;
            }

            if (its > 1) {
                if ((previous_c > centerOfAttack)) {
                    System.out.println(String
                        .format("previous_c:\t%s\ncenterOfAttack:\t%s\n", previous_c, centerOfAttack));
                    throw new Exception("previous_c > centerOfAttac: wrongly ordered, you probably ran out of precision");
                }
                if (nextStreamValue > previous_y) {
                    System.out.println(String
                        .format("nextStreamValue:\t%s\nprevious_y:\t%s\n", nextStreamValue, previous_y));
                    throw new Exception("nextStreamValue > previous_y: wrongly ordered, you probably ran out of precision");
                }
            }

            // weight of the centroid we will fabricate
            // this is the formula for K_0
            // we also maintain centroid to the right of the attack
            System.out.println("centroids count: " + digest.centroids().size());
            weightGoal = (int) Math.ceil((weightToLeft + weightToRight) / ((delta / 2d) - 3d));
            System.out.println("weightGoal: " + weightGoal);

            System.out.println(
                "centroids exceeding goal at start: " + centroidsExceedingCount(digest.centroids(),
                    weightGoal));
            currentDeficit = weightGoal - weightOfAttacked;
            assert currentDeficit >= 0;

            double EPS = .001;
            double EPS_2 = .01;

            // fill up the old one
            double point = centerOfAttack; // * (1d - EPS) + nextStreamValue * EPS;
            for (int v = 0; v < currentDeficit; v++) {
                //double point = centerOfAttack + (nextStreamValue - centerOfAttack) / 10000; // * onePlus; //no Plus Epsilon
                digest.add(point, 1);
                data.add(point);
            }


            // make the new one
            int v = 0;
            double anotherPoint = centerOfAttack * (1d - EPS_2) + nextStreamValue * EPS_2;
            
            double leftEdge = nextStreamValue * 0.1; //centerOfAttack * EPS + 
            
            for (; v < weightGoal / 5; v++) {
                digest.add(anotherPoint);
                data.add(anotherPoint);
            }            
            for (; v < weightGoal; v++) {
                digest.add(leftEdge, 1);
                data.add(leftEdge);
            }
            

            // fill up the centroid to the right
            currentDeficitRight = weightGoal - rightNeighbor.count();
            assert currentDeficit >= 0;
            //Centroid rightCentroid = aboveValue(leftEdge, digest.centroids());
            //int deficit = weightGoal - rightCentroid.count();
            double rightCentroidVal = rightNeighbor.mean();
            for (int pp = 0; pp < currentDeficitRight; pp++) {
                digest.add(rightCentroidVal);
                data.add(rightCentroidVal);
            }

            System.out.println(String
                .format("centerOfAttack:\t%s\npoint:\t\t%s\nanotherPoint:\t%s\nleftEdge:  \t%s\nnextStreamValue:\t%s", centerOfAttack, point, anotherPoint, leftEdge,
                    nextStreamValue));

            digest.compress();

            System.out.println(
                "centroids exceeding goal at end: " + centroidsExceedingCount(digest.centroids(),
                    weightGoal - 1));

            double bad_point = 0; //centerOfAttack * .00000000001 + nextStreamValue * (1 - .00000000001);
            Centroid belowZeroC = belowValue(0, digest.centroids());
            Centroid aboveZeroC = aboveValue(0, digest.centroids());
            //double bad_point = belowZeroC.mean();
            Collections.sort(data);

            System.out.println("finished iteration: " + its);
            System.out.println("belowZeroC mean: " + belowZeroC.mean());
            System.out.println("aboveZeroC mean: " + aboveZeroC.mean());
            //System.out.println("bad point: " + bad_point);
            System.out.println("td " + digest.cdf(bad_point));
            System.out
                .println("truth " + countBelow(bad_point, data) / (double) data.size());

            double error = Math
                .abs(digest.cdf(bad_point) - countBelow(bad_point, data) / (double) data.size());

            if (error > maximalError) {
                maximalError = error;
                indexError = its;
            }
            System.out.println("maximal error so far: " + maximalError + " on " + indexError);
            System.out.println("num centroids: " + digest.centroids().size() + "\n");

            //double yetAnotherPoint = centerOfAttack * 0.9 + nextStreamValue * ;
            centroidToAttack = aboveValue(anotherPoint, digest.centroids()); //centerOfAttack //belowZeroC; //
            if (centroidToAttack.mean() < centerOfAttack) {
                System.out.println("wtf: centroidToAttack.mean()=" + centroidToAttack.mean() + " < centerOfAttack");
            }
            //if (centroidToAttack.mean() > 0) {
            //  System.out.println("above zero: centroidToAttack.mean()=" + centroidToAttack.mean());
            //  centroidToAttack = belowZeroC;
            //}
            rightNeighbor = aboveValue(centroidToAttack.mean(), digest.centroids());
        }

        double bad_point = (previous_c + previous_y) / 2d;
        Collections.sort(data);

        System.out.println("td " + digest.cdf(bad_point));
        System.out.println("truth " + countBelow(bad_point, data) / (double) data.size());
        System.out.println("iterations" + its);
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

    // assume centroids in ascending order
    private Centroid aboveValue(double c, Collection<Centroid> centroids) throws Exception {
        Iterator<Centroid> centroidIterator = centroids.iterator();
        Centroid centroid = centroidIterator.next();
        while (centroidIterator.hasNext()) {
            if (c < (centroid.mean())) {
                return centroid;
            } else {
                centroid = centroidIterator.next();
            }
        }
        throw new Exception("couldn't find a centroid above threshold");
    }

    // assume centroids in ascending order
    private Centroid belowValue(double c, Collection<Centroid> centroids) throws Exception {
        Iterator<Centroid> centroidIterator = centroids.iterator();
        Centroid centroid = centroidIterator.next();
        Centroid previous = null;
        while (centroidIterator.hasNext()) {
            if (centroid.mean() >= c) {
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

}


