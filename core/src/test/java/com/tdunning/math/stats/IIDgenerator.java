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

import java.time.format.DateTimeFormatter;
import java.time.LocalDateTime;

import java.io.*;
import java.lang.Math;
import java.nio.ByteBuffer;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.time.*;

import org.junit.Ignore;


/**
 *
 */

public class IIDgenerator {

    // confidence intervals of normal distrib. 
    public static final double M4SD = 0.0000316712418331; //minus 4 StdDev
    public static final double M3SD = 0.0013498980316301; //minus 3 StdDev
    public static final double M2SD = 0.0227501319481792; //minus 2 StdDev
    public static final double M1SD = 0.1586552539314570; //minus 1 StdDev
    public static final double P1SD = 0.8413447460685430; //plus  1 StdDev
    public static final double P2SD = 0.9772498680518210; //plus  2 StdDev
    public static final double P3SD = 0.9986501019683700; //plus  3 StdDev
    public static final double P4SD = 0.9999683287581670; //plus  4 StdDev

    // properties from config file
    final String Distribution;
    final int N; // the number of i.i.d. items in the generated input
    final int T; // the number of trials
    final int LgN; // log_2 of the number of i.i.d. items in the generated input
    final int LgT; // log_2 of the number of trials
    final int NumClusters; // for clustered distribution
    final int NumberOfPoints; // number of points where we probe the rank estimates
    final Boolean NegativeNumbers; // if false, generate positive numbers only
    final Boolean WriteCentroidData;
    final double lambda; // parameter for exponential distribution
    final String DigestImpl; // "Merging" or "AVLTree"
    final int Compression; // delta for t-digest
    final ScaleFunction scale; // ScaleFunction for t-digest
    //final String InputStreamFileName;
    //final String InputStreamFileDir;
    final String DigestStatsFileName;
    final String DigestStatsDir;
    final String FileSuffix;

    // vars for experiments
    Random rand;
    Properties prop;
    int maxExpBase2;
    int n;


    public IIDgenerator(final String configFile) throws Exception {
        Instant startTime = Instant.now();

        System.out.println("processing config file: " + configFile);
        prop = new Properties();
        FileInputStream instream = new FileInputStream(configFile);
        prop.load(instream);

        // load properties
        Distribution = getProperty("Distribution");
        DigestImpl = getProperty("DigestImpl");
        LgN = Integer.parseInt(getProperty("LgN")); // base 2
        N = 1 << LgN;
        LgT = Integer.parseInt(getProperty("LgT")); // base 2
        T = 1 << LgT;
        NumClusters = Integer.parseInt(getProperty("NumClusters"));
        NumberOfPoints = Integer.parseInt(
            getProperty("NumberOfPoints")); // number of points where we probe the rank estimates
        NegativeNumbers = Boolean.parseBoolean(
            getProperty("NegativeNumbers")); // if false, generate positive numbers only
        WriteCentroidData = Boolean.parseBoolean(getProperty("WriteCentroidData"));
        lambda = Double.parseDouble(getProperty("Lambda")); // for exponential distribution
        Compression = Integer.parseInt(getProperty("Compression")); // delta for t-digest
        scale = ScaleFunction.valueOf(getProperty("ScaleFunction")); // ScaleFunction for t-digest
//        InputStreamFileName = getProperty("InputStreamFileName");
//        InputStreamFileDir = getProperty("InputStreamFileDir");
        DigestStatsFileName = getProperty("DigestStatsFileName");
        DigestStatsDir = getProperty("DigestStatsDir");
        FileSuffix = getProperty("FileSuffix");

        List<Double> data = new ArrayList<Double>();
        //Files.createDirectories(Paths.get(InputStreamFileDir));
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss");
        LocalDateTime now = LocalDateTime.now();
        String fileNamePart =
            "_" + dtf.format(now) + "_" + Distribution + (NegativeNumbers ? "" : "_PositiveOnly")
                + "_lgN=" + String.valueOf(LgN);
//        String inputFilePath =
//            InputStreamFileDir + InputStreamFileName + "_" + fileNamePart + FileSuffix;
//        PrintWriter w = new PrintWriter(inputFilePath);
        rand = new Random();

        // intialize t-digests to collect errors of t-digest:
        TDigest[] errorDigests = new TDigest[NumberOfPoints
            + 1]; // TODO use t-digest or something else?
        for (int t = 0; t <= NumberOfPoints; t++) {
            errorDigests[t] = new MergingDigest(200);
            errorDigests[t].setScaleFunction(ScaleFunction.K_2); // we do not need extreme quantiles
        }

        maxExpBase2 = (int) (Math.log(Double.MAX_VALUE / N) / Math.log(2));

//        for (n = 0; n < N; n++) {
//            double item = generateItem();
//            data.add(item);
//            sortedData.add(item);
//            //w.println(String.valueOf(item));
//        }
//        Collections.sort(sortedData);
//        w.close();
//        System.out.println("file with generated input: " + inputFilePath);
//        System.out.flush();
        TDigest digest = null;
        List<Double> sortedData = null;
        for (int trial = 0; trial < T; trial++) { // trials
            if (DigestImpl.equals("Merging")) {
                digest = new MergingDigest(Compression);
            } else if (DigestImpl.equals("AVLTree")) {
                digest = new AVLTreeDigest(Compression);
            } else {
                throw new Exception("unknown digest implementation: '" + DigestImpl + "'");
            }
            assert digest.size() == 0;
            digest.setScaleFunction(scale);
            if (WriteCentroidData && trial == T - 1) {
                digest.recordAllData(); // tracks centroids during the last trial
            }
            sortedData = new ArrayList<Double>();
            for (n = 0; n < N; n++) {
                double item = generateItem();
                sortedData.add(item);
                //if (trial == T-1) data.add(item);
                digest.add(item);
            }
            Collections.sort(sortedData);
            digest.compress();
            // extract error
            for (int t = 0; t <= NumberOfPoints; t++) {
                //THE FOLLOWING IS EXTREMELY SLOW: Dist.cdf(item, sortedData);
                int rTrue = (int) Math.ceil(t / (float) NumberOfPoints * N) + 1;
                if (rTrue > N) {
                    rTrue--;
                }
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
                double rEst = digest.cdf(item) * N + 0.5; // TODO: should the +0.5 be here?
                //double relErr = 0;
                double addErr = 0;
                if (rEst < rTrueMin) {
                    //relErr = Math.abs(rTrueMin - rEst) / (N - rTrue + 1);
                    addErr = (rEst - rTrueMin) / N;
                }
                if (rEst > rTrueMax) {
                    //relErr = Math.abs(rTrueMax - rEst) / (N - rTrue + 1);
                    addErr = (rEst - rTrueMax) / N;
                }
                errorDigests[t].add(addErr);
            }
        }
        //Collections.sort(data);

        System.out
            .println("processing by t-digest done for compression =" + String.valueOf(Compression));
        System.out.flush();

        if (WriteCentroidData) // from the last trial
        {
            writeCentroidData(digest, DigestStatsDir,
                DigestStatsDir + DigestStatsFileName + fileNamePart + "_" + DigestImpl
                    + "_compr=" + String.valueOf(Compression) + digest.scale.toString()
                    + FileSuffix);
        }

        writeResults(Compression, n, NumberOfPoints, prop, digest, errorDigests, sortedData,
            startTime, DigestStatsDir,
            DigestStatsDir + DigestStatsFileName + fileNamePart + "_" + DigestImpl + "-stats-PP="
                + String.valueOf(NumberOfPoints)
                + "_compr=" + String.valueOf(Compression) + "_" + digest.scale.toString()
                + FileSuffix);
        // TODO write all properties into results
    }

    // the following does an additional trim and removes possible comment
    String getProperty(String propName) {
        String value = prop.getProperty(propName);
        int inx = value.indexOf('#');
        if (inx >= 0) {
            value = value.substring(0, inx - 1);
        }
        return value.trim();
    }

    double generateItem() throws Exception {
        double item = 0;
        switch (Distribution) {
            case "loguniform":
                item = Math.pow(2, (rand.nextDouble() - 0.5) * 2 * maxExpBase2);
                break;
            case "exponential":
                item = Math.log(1 - rand.nextDouble()) / (-lambda);
                break;
            case "clustered":
                item =
                    (0.9999 + rand.nextDouble() / 100000000) * (rand.nextInt(NumClusters / 2) + 1)
                        * (Double.MAX_VALUE / (N * NumClusters));
                break;
            case "uniform":
                item = rand.nextDouble();
                break;
            case "normal":
                item = rand.nextGaussian();
                break;
            default:
                throw new Exception("Distribution '" + Distribution + "' undefined");
        }
        if (NegativeNumbers && rand.nextDouble() < 0.5) {
            item = -item;
        }
        return item;
    }


    public static void writeResults(int compr, int size, int numPoints, Properties prop,
        TDigest digest, TDigest[] errorDigests,
        List<Double> sortedData, Instant startTime, String digestStatsDir, String outName) throws
        IOException {
        Files.createDirectories(Paths.get(digestStatsDir));
        System.out.printf("stats file:" + outName + "\n");
        File fout = new File(outName);
        fout.createNewFile();
        FileWriter fwout = new FileWriter(fout);

        //System.out.printf("computing rel. errors\n");
        //System.out.flush();

        fwout.write(
            "true quantile;-2 std. dev. of error; median error;+2 std. dev. of error;item\n");
        for (int t = 0; t <= numPoints; t++) {
            int rTrue = (int) Math.ceil(t / (float) numPoints * size) + 1;
            if (rTrue > size) {
                rTrue--;
            }
            double item = sortedData.get(rTrue - 1); // in the last trial
            double addErrM2SD = errorDigests[t].quantile(M2SD);
            double addErrMed = errorDigests[t].quantile(0.5);
            double addErrP2SD = errorDigests[t].quantile(P2SD);

            //relErr = Math.abs(rTrueMax - rEst) / (size - rTrue + 1);
            fwout.write(String
                .format("%.6f;%.6f;%.6f;%.6f;%s\n", rTrue / (float) size, addErrM2SD, addErrMed,
                    addErrP2SD, String.valueOf(item)));
        }
        fwout.write("\n");
        fwout.write(String.format("n=%d\n", size));
        fwout.write(String.format("scale func. = %s\n", digest.scale.toString()));
        fwout.write(String.format("delta = %d (compression param of t-digest)\n", compr));
        fwout.write(String.format("# of centroids = %d\n", digest.centroids().size()));
        fwout.write(String.format("size in bytes = %d\n", digest.byteSize()));
        Duration diff = Duration.between(startTime, Instant.now());
        String hms = String.format("%d:%02d:%02d",
            diff.toHours());//,
        //diff.toMinutesPart(),
        //diff.toSecondsPart());
        fwout.write(String.format("time taken = %s\n", hms));

        fwout.write("\nProperties:\n");
        for (Object key : prop.keySet()) {
            fwout.write(key + " = " + prop.getProperty(key.toString()) + "\n");
        }

        fwout.write("\nCentroids:\n");
        for (Centroid centr : digest.centroids()) {
            fwout.write(centr.toString() + "\n");
        }
        fwout.close();
        System.out.flush();

    }

    public static void writeCentroidData(TDigest digest, String digestStatsDir, String outName)
        throws IOException {
        Files.createDirectories(Paths.get(digestStatsDir));

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


    @SuppressWarnings("unused")
    public static void main(final String[] args) throws Exception {
        for (int j = 0; j < args.length; j++) {
            new IIDgenerator(args[j]);
        }
    }

}


