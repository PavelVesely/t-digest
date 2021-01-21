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
public class ZoomInPlotRelErrorTest extends AbstractTest {

    private static final int[] CompressionsForTesting = {
        100}; //100, 200, 500, 1000, 5000, 10000}; // param delta of t-digest
    private static final int NumberOfPoints = 200; // number of points where we probe the rank estimates

    // params for generating the input
    private static final int N = 10000000; // stream length 
    private static final int K = 2; // N / K should be roughly at most 200 (otherwise, we don't have enough precision)
    private static final int PrefixSize = 0; // number of points below zero in each iteration
    private static final int NumberOfRepeats = 400; // for zoom in generator; N * NumberOfRepeats is (approx.) the number of items in the final instance
    private static final Boolean NegativeNumbers = true; // if false, zoom in generates positive numbers only
    private static final Boolean WriteCentroidData = false; 
    private static final double lambda = 0.000001; // for exponential distribution
    private static final String InputStreamFileName = "t-digest-genInput";
    private static final String InputStreamFileDir = "../../../../data/inputs/";
    //private static final String InputStreamFileDir = "/aux/vesely/TD-inputs/"; // // CHANGE AS APPROPRIATE
    private static final String DigestStatsFileName = "t-digest-results";
    private static final String DigestStatsDir = "../../../../data/results/";
    //private static final String DigestStatsDir = "../../../../TD-stats/"; // CHANGE AS APPROPRIATE
    private static final String FileSuffix = ".csv";
    //private static final String StatsFileDir = "/home/vesely/research/biasedQuantiles/TD-stats/"; // CHANGE AS APPROPRIATE
    
    
    
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
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_K=" + String.valueOf(K) + "_N=" + String
                .valueOf(N) + "_PS=" + String.valueOf(PrefixSize) + "_repeats=" + String
                .valueOf(NumberOfRepeats) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
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
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
            for (double item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            writeResults(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_K=" + String.valueOf(K) + "_N=" + String
                    .valueOf(N) + "_PS=" + String.valueOf(PrefixSize) + "_repeats=" + String
                    .valueOf(NumberOfRepeats) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
    }

    @Ignore
    public void testZoomInIIDgenerator() throws FileNotFoundException, IOException {
        List<Double> sortedData = new ArrayList<Double>();
        List<Double> data = new ArrayList<Double>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_zoominIIDitems" + (NegativeNumbers
                ? "_wNegativeNumbers" : "_PositiveOnly") + "_N=" + String.valueOf(N) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        Random rand = new Random();
        int n;
        final int maxExp = (int) (Math.log(Double.MAX_VALUE / 100) / Math.log(2));
        for (n = 0; n < N; n++) {
            double item = Math.pow(2, (rand.nextDouble() - 0.5) * 2 * maxExp);
            if (NegativeNumbers && rand.nextDouble() < 0.5) {
                item = -item;
            }
            data.add(item);
            sortedData.add(item);
            w.println(String.valueOf(item));
        }
        Collections.sort(sortedData);
        w.close();
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
            digest.recordAllData(); // tracks centroids
            for (double item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            if (WriteCentroidData)
                writeCentroidData(digest,
                      DigestStatsDir + DigestStatsFileName + "_zoominIIDitems" + "_N=" + String.valueOf(N)
                        + "-centroids-PP_" + String.valueOf(NumberOfPoints)
                        + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);

            writeResults(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_zoominIIDitems" + (
                    NegativeNumbers ? "_wNegativeNumbers" : "_PositiveOnly") + "_N=" + String
                    .valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);

        }
    }

    @Ignore
    public void testUniformIIDgenerator() throws FileNotFoundException, IOException {
        List<Double> sortedData = new ArrayList<Double>();
        List<Double> data = new ArrayList<Double>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_uniformIIDitems" + "_N=" + String
                .valueOf(N) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        Random rand = new Random();
        int n;
        for (n = 0; n < N; n++) {
            double item = rand.nextDouble() * Double.MAX_VALUE / N;
            data.add(item);
            sortedData.add(item);
            w.println(String.valueOf(item));
        }
        Collections.sort(sortedData);
        w.close();
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
            for (double item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            writeResults(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_uniformIIDitems" + "_N=" + String
                    .valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
    }

    @Ignore
    public void test2DistribUniformIIDgenerator() throws FileNotFoundException, IOException {
        List<Double> sortedData = new ArrayList<Double>();
        List<Double> data = new ArrayList<Double>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_2DistribUniformIIDitems" + "_N=" + String
                .valueOf(N) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        Random rand = new Random();
        int n;
        for (n = 0; n < N / 2; n++) {
            double item = rand.nextDouble() * Double.MAX_VALUE / N;
            data.add(item);
            sortedData.add(item);
            w.println(String.valueOf(item));
        }
        for (n = N / 2; n < N; n++) {
            double item = rand.nextDouble() * N * Double.MIN_VALUE;
            data.add(item);
            sortedData.add(item);
            w.println(String.valueOf(item));
        }
        Collections.sort(sortedData);
        w.close();
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
            for (double item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            writeResults(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_2DistribUniformIIDitems" + "_N=" + String
                    .valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
    }


    @Ignore
    public void test2valuesIIDgenerator() throws FileNotFoundException, IOException {
        List<Double> sortedData = new ArrayList<Double>();
        List<Double> data = new ArrayList<Double>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_2valuesIIDitems" + "_N=" + String
                .valueOf(N) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        Random rand = new Random();
        int n;
        for (n = 0; n < N; n++) {
            double item = (0.99 + rand.nextDouble() / 100);
            if (rand.nextDouble() < 0.5) {
                item = -item;
            } else {
                data.add(item);
            }
            sortedData.add(item);
            w.println(String.valueOf(item));
        }
        Collections.sort(sortedData);
        w.close();
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
            for (double item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            writeResults(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_2valuesIIDitems" + "_N=" + String
                    .valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
    }

    @Ignore
    public void test3valuesIIDgenerator() throws FileNotFoundException, IOException {
        List<Double> sortedData = new ArrayList<Double>();
        List<Double> data = new ArrayList<Double>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_3valuesIIDitems" + "_N=" + String
                .valueOf(N) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        Random rand = new Random();
        int n;
        final int maxExp = (int) Math.log10(Double.MAX_VALUE / 100);
        for (n = 0; n < N; n++) {
            double item = (0.99 + rand.nextDouble() / 100);
            double r = rand.nextDouble();
            if (r < 0.333) {
                item = -item * Double.MAX_VALUE / N;
            } else if (r > 0.666) {
                item = item * Double.MAX_VALUE / N;
            } else {
                item = item * Math.pow(10, -maxExp);
            }
            data.add(item);
            sortedData.add(item);
            w.println(String.valueOf(item));
        }
        Collections.sort(sortedData);
        w.close();
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
            for (double item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            writeResults(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_3valuesIIDitems" + "_N=" + String
                    .valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
    }

    @Ignore
    public void testClusteredIIDgenerator() throws FileNotFoundException, IOException {
        List<Double> sortedData = new ArrayList<Double>();
        List<Double> data = new ArrayList<Double>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_clusteredIIDitems" + "_K=" + String
                .valueOf(K) + "_N=" + String.valueOf(N) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        Random rand = new Random();
        int n;
        final int maxExp = (int) Math.log10(Double.MAX_VALUE / 100);
        for (n = 0; n < N; n++) {
            double item =
                (0.9999 + rand.nextDouble() / 100000000) * (rand.nextInt(K) + 1) * (Double.MAX_VALUE / (N
                    * K));
            if (rand.nextDouble() < 0.5) {
                item = -item;
            }
            data.add(item);
            sortedData.add(item);
            w.println(String.valueOf(item));
        }
        Collections.sort(sortedData);
        w.close();
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
            for (double item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            writeResults(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_clusteredIIDitems" + "_K=" + String
                    .valueOf(K) + "_N=" + String.valueOf(N) + "-stats-PP_" + String
                    .valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
    }

    @Ignore
    public void testExpDistribIIDgenerator() throws FileNotFoundException, IOException {
        List<Double> sortedData = new ArrayList<Double>();
        List<Double> data = new ArrayList<Double>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "ExpDistrib" + "_lambda=" + String
                .valueOf(lambda) + "_N=" + String.valueOf(N) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        Random rand = new Random();
        int n;
        for (n = 0; n < N; n++) {
            double item = Math.log(1 - rand.nextDouble()) / (-lambda);
            if (rand.nextDouble() < 0.5) {
                item = -item;
            }
            data.add(item);
            sortedData.add(item);
            w.println(String.valueOf(item));
        }
        Collections.sort(sortedData);
        w.close();
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
            for (double item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            writeResults(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "ExpDistrib" + "_lambda=" + String
                    .valueOf(lambda) + "_N=" + String.valueOf(N) + "-stats-PP_" + String
                    .valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
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
        for (int t = 0; t <= NumberOfPoints; t++) {
            //THE FOLLOWING IS EXTREMELY SLOW: Dist.cdf(item, sortedData);
            int rTrue = (int) Math.ceil(t / (float) NumberOfPoints * size) + 1;
            if (rTrue > size) rTrue--;
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

    @Ignore
    public void testZoomInFloat() throws FileNotFoundException, IOException {

        List<Float> sortedData = new ArrayList<Float>();
        List<Float> sortedDataPart = new ArrayList<Float>();
        List<Float> data = new ArrayList<Float>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_Float" + "_K=" + String.valueOf(K) + "_N="
                + String
                .valueOf(N) + "_PS=" + String.valueOf(PrefixSize) + "_repeats=" + String
                .valueOf(NumberOfRepeats) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        int n = 0;
        for (int r = 0; r < NumberOfRepeats; r++) {  // same instance is essentially repeated
            float max = Float.MAX_VALUE / (10 * K + r); // initial interval max
            float min = 0; //Float.MIN_VALUE;
            int nn = 0;
            for (int i = -PrefixSize + 1; i <= K - PrefixSize;
                i++) {//for (int i = 0; i <= K; i++) {
                float item = min + (i / (float) (K + 1)) * (max - min);
                data.add(item);
                sortedData.add(item);
                sortedDataPart.add(item);
                w.println(String.valueOf(item));
                nn++;
                n++;
            }
            int maxREindex = -1;
            while (nn < N) {
                maxREindex = -Collections.binarySearch(sortedDataPart, Float.MIN_NORMAL) - 1;
                min = 0; //maxREindex > 0 ? sortedDataPart.get(maxREindex - 1) : 0;
                max = sortedDataPart.get(maxREindex);
                //System.out.printf(String.format("phase maxRErTrue %d;\tmin %s;\tmax %s;\t means %d\n",
                //		maxREindex, String.valueOf(min), String.valueOf(max), digest.centroidCount()));
                if ((max - min) / (float) (K + 1) < Float.MIN_VALUE) {
                    System.out.printf("TOO SMALL max - min\n");
                }
                for (int i = -PrefixSize; i <= K - PrefixSize; i++) {
                    float item = min + (i / (float) (K + 1)) * (max - min);
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
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
            for (float item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            writeResultsFloat(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_Float" + "_K=" + String.valueOf(K) + "_N="
                    + String
                    .valueOf(N) + "_PS=" + String.valueOf(PrefixSize) + "_repeats=" + String
                    .valueOf(NumberOfRepeats) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
    }


    @Ignore
    public void testZoomInExp2IIDgeneratorFloat() throws FileNotFoundException, IOException {
        List<Float> sortedData = new ArrayList<Float>();
        List<Float> data = new ArrayList<Float>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_ZoomInExp2IIDitems_minusInfty_Float"
                + "_N=" + String.valueOf(N) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        Random rand = new Random();
        int n;
        final int maxExp = (int) Math.log(Float.MAX_VALUE / 10);
        for (n = 0; n < N; n++) {
            float item = (0.5f + (float) rand.nextDouble() / 2) * (float) Math
                .pow(2, rand.nextInt(2 * maxExp) - maxExp);
            if (rand.nextDouble() < 0.5) {
                item = -item;
            }
            data.add(item);
            sortedData.add(item);
            w.println(String.valueOf(item));
        }
        Collections.sort(sortedData);
        w.close();
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();

        for (int compr : CompressionsForTesting) {
            TDigest digest = new MergingDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_3); //_GLUED;
            for (float item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            writeResultsFloat(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_ZoomInExp2IIDitems_minusInfty_Float"
                    + "_N=" + String.valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
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
                  .format("centerOfAttack < nextStreamValue: nextStreamValue:\t%s\ncenterOfAttack:\t%s\n", nextStreamValue, centerOfAttack));
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
            centroidToAttack = aboveValue(anotherPoint, digest.centroids()); //centerOfAttack
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

    private static void writeResultsFloat(int compr, int size, TDigest digest,
        List<Float> sortedData, String digestStatsDir, String outName) throws
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
            float item = sortedData.get(rTrue - 1);
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

    private static void writeCentroidData(TDigest digest, String outName) throws IOException {
        Files.createDirectories(Paths.get(DigestStatsDir));

        System.out.printf("stats file:" + outName + "\n");
        File fout = new File(outName);
        fout.createNewFile();
        FileWriter fwout = new FileWriter(fout);

        int index = 0;
        for (Centroid c : digest.centroids()) {

            fwout.write(index + ",");
            for (Double d : c.data()) {
                fwout.write(d + ",");
            }
            index++;
            fwout.write("\n");
        }
        fwout.close();
        System.out.flush();
    }

}


