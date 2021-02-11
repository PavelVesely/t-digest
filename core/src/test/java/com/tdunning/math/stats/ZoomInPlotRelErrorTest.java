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

import org.junit.*;

import java.io.*;
import java.lang.Math;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

import org.junit.Ignore;

/**
 *
 */
@Ignore
public class ZoomInPlotRelErrorTest extends AdversarialAttackTest {

    @Ignore
    public void testZoomIn() throws FileNotFoundException, IOException, Exception {

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

            CarefulAttackAuxiliary.writeResults(compr, n, digest, data, sortedData, "DO NOT COMPARE TO ANYTHING", DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_K=" + String.valueOf(K) + "_N=" + String
                    .valueOf(N) + "_PS=" + String.valueOf(PrefixSize) + "_repeats=" + String
                    .valueOf(NumberOfRepeats) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix,
                true);
        }
    }

    @Test
    public void testZoomInIIDgenerator() throws FileNotFoundException, IOException, Exception {
        List<Double> sortedData = new ArrayList<Double>();
        List<Double> data = new ArrayList<Double>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_AVL_zoominIIDitems" + (NegativeNumbers
                ? "_wNegativeNumbers" : "_PositiveOnly") + "_N=" + String.valueOf(N) + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        Random rand = new Random();
        int n;
        final int maxExp = (int) (Math.log(Double.MAX_VALUE / 10000) / Math.log(2));
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
            TDigest digest = new AVLTreeDigest(compr);
            digest.setScaleFunction(ScaleFunction.K_2_GLUED); //_GLUED
            digest.recordAllData(); // tracks centroids
            for (double item : data) {
                digest.add(item);
            }
            digest.compress();
            System.out
                .println("processing by t-digest done for compression =" + String.valueOf(compr));
            System.out.flush();

            if (WriteCentroidData) {
                CarefulAttackAuxiliary.writeCentroidData(digest,
                    DigestStatsDir + DigestStatsFileName + "_zoominIIDitems" + "_N=" + String
                        .valueOf(N)
                        + "-centroids-PP_" + String.valueOf(NumberOfPoints)
                        + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
            }

            CarefulAttackAuxiliary.writeResults(compr, n, digest, data, sortedData, "DO NOT COMPARE TO ANYTHING", DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_AVL_zoominIIDitems" + (
                    NegativeNumbers ? "_wNegativeNumbers" : "_PositiveOnly") + "_N=" + String
                    .valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix,
                true);

        }
    }

    @Ignore
    public void testUniformIIDgenerator() throws FileNotFoundException, IOException, Exception {
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

            CarefulAttackAuxiliary.writeResults(compr, n, digest, data, sortedData, "DO NOT COMPARE TO ANYTHING", DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_uniformIIDitems" + "_N=" + String
                    .valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix,
                true);
        }
    }

    @Ignore
    public void test2DistribUniformIIDgenerator() throws FileNotFoundException, IOException, Exception {
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

            CarefulAttackAuxiliary.writeResults(compr, n, digest, data, sortedData, "DO NOT COMPARE TO ANYTHING", DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_2DistribUniformIIDitems" + "_N=" + String
                    .valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix,
                true);
        }
    }

    @Ignore
    public void test2valuesIIDgenerator() throws FileNotFoundException, IOException, Exception {
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

            CarefulAttackAuxiliary.writeResults(compr, n, digest, data, sortedData, "DO NOT COMPARE TO ANYTHING", DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_2valuesIIDitems" + "_N=" + String
                    .valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix,
                true);
        }
    }

    @Ignore
    public void test3valuesIIDgenerator() throws FileNotFoundException, IOException, Exception {
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

            CarefulAttackAuxiliary.writeResults(compr, n, digest, data, sortedData, "DO NOT COMPARE TO ANYTHING", DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_3valuesIIDitems" + "_N=" + String
                    .valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix,
                true);
        }
    }

    @Ignore
    public void testClusteredIIDgenerator() throws FileNotFoundException, IOException, Exception {
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
                (0.9999 + rand.nextDouble() / 100000000) * (rand.nextInt(K) + 1) * (Double.MAX_VALUE
                    / (N
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

            CarefulAttackAuxiliary.writeResults(compr, n, digest, data, sortedData, "DO NOT COMPARE TO ANYTHING", DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_clusteredIIDitems" + "_K=" + String
                    .valueOf(K) + "_N=" + String.valueOf(N) + "-stats-PP_" + String
                    .valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix,
                true);
        }
    }

    @Ignore
    public void testExpDistribIIDgenerator() throws FileNotFoundException, IOException, Exception {
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

            CarefulAttackAuxiliary.writeResults(compr, n, digest, data, sortedData, "DO NOT COMPARE TO ANYTHING", DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "ExpDistrib" + "_lambda=" + String
                    .valueOf(lambda) + "_N=" + String.valueOf(N) + "-stats-PP_" + String
                    .valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix,
                true);
        }
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

            CarefulAttackAuxiliary.writeResultsFloat(compr, n, digest, sortedData, DigestStatsDir,
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

            CarefulAttackAuxiliary.writeResultsFloat(compr, n, digest, sortedData, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + "_ZoomInExp2IIDitems_minusInfty_Float"
                    + "_N=" + String.valueOf(N) + "-stats-PP_" + String.valueOf(NumberOfPoints)
                    + "_compr_" + String.valueOf(compr) + digest.scale.toString() + FileSuffix);
        }
    }


}


