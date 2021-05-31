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

import com.tdunning.math.stats.datasketches.req.ReqSketch;
import com.tdunning.math.stats.datasketches.kll.KllDoublesSketch;

/**
 *
 */

public class LoguniformWithVaryingMaxExpGenerator {


    // properties from config file
    final String Distribution;
    final int N; // the number of i.i.d. items in the generated input
    final int T; // the number of trials
    final int LgN; // log_2 of the number of i.i.d. items in the generated input
    final int LgT; // log_2 of the number of trials
    //final int NumClusters; // for clustered distribution
    final int NumberOfPoints; // number of points where we probe the rank estimates
    final Boolean NegativeNumbers; // if false, generate positive numbers only
    //final Boolean WriteCentroidData;
    //final double lambda; // parameter for exponential distribution
    final String DigestImpl; // "Merging" or "AVLTree"
    final int Compression; // delta for t-digest
    final ScaleFunction scale; // ScaleFunction for t-digest
    //final String InputStreamFileName;
    //final String InputStreamFileDir;
    final String DigestStatsFileName;
    final String DigestStatsDir;
    final String FileSuffix;
    final int reqK;
    int maxExpStep;

    // vars for experiments
    Random rand;
    Properties prop;
    int maxExp, maxExpForDoubleBase10;
    int n;
    

    double[] maxRelErrorClustering;
    double[] maxRelErrorMerging;
    double[] maxRelErrorRS_2SD;
    double[] maxRelErrorRS_median;
    double[] avgRelErrorClustering;
    double[] avgRelErrorMerging;
    double[] avgRelErrorRS_2SD;
    double[] avgRelErrorRS_median;

    public LoguniformWithVaryingMaxExpGenerator(final String configFile) throws Exception {
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
        NumberOfPoints = Integer.parseInt(
            getProperty("NumberOfPoints")); // number of points where we probe the rank estimates
        NegativeNumbers = Boolean.parseBoolean(
            getProperty("NegativeNumbers")); // if false, generate positive numbers only
        //WriteCentroidData = Boolean.parseBoolean(getProperty("WriteCentroidData")); // disable in this experiment for simplicity (run IIDgenerator instead)
        Compression = Integer.parseInt(getProperty("Compression")); // delta for t-digest
        scale = ScaleFunction.valueOf(getProperty("ScaleFunction")); // ScaleFunction for t-digest
//        InputStreamFileName = getProperty("InputStreamFileName");
//        InputStreamFileDir = getProperty("InputStreamFileDir");
        DigestStatsFileName = getProperty("DigestStatsFileName");
        DigestStatsDir = getProperty("DigestStatsDir");
        FileSuffix = getProperty("FileSuffix");
        reqK = Integer.parseInt(getProperty("ReqK"));
        maxExpStep = Integer.parseInt(getProperty("MaxExpStep"));

        if (DigestImpl.trim().length() > 0 && !DigestImpl.equals("AVLTree") && !DigestImpl
            .equals("Merging")) {
            throw new Exception("unknown digest implementation: '" + DigestImpl + "'");
        }

        List<Double> data = new ArrayList<Double>();
        //Files.createDirectories(Paths.get(InputStreamFileDir));
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss");
        LocalDateTime now = LocalDateTime.now();
        String fileNamePart = "_" + Distribution + (NegativeNumbers ? "" : "_PositiveOnly")
            + "_lgN=" + String.valueOf(LgN) + "_lgT=" + String.valueOf(LgT) + "_maxExpStep=" + String.valueOf(maxExpStep);
//        String inputFilePath =
//            InputStreamFileDir + InputStreamFileName + "_" + fileNamePart + FileSuffix;
//        PrintWriter w = new PrintWriter(inputFilePath);
        rand = new Random();

        maxExpForDoubleBase10 = (int) (Math.log(Double.MAX_VALUE / N) / Math.log(10));        
        int numMaxExps = maxExpForDoubleBase10 / maxExpStep;
        maxRelErrorClustering = new double[numMaxExps];
        maxRelErrorMerging = new double[numMaxExps];
        maxRelErrorRS_2SD = new double[numMaxExps];
        maxRelErrorRS_median = new double[numMaxExps];
        avgRelErrorClustering = new double[numMaxExps];
        avgRelErrorMerging = new double[numMaxExps];
        avgRelErrorRS_2SD = new double[numMaxExps];
        avgRelErrorRS_median = new double[numMaxExps];


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
        //TDigest digest = null;
        MergingDigest merging = null;
        AVLTreeDigest clustering = null;
        ReqSketch reqsk = null;
        //ReqSketchBuilder reqskBuilder = new ReqSketchBuilder();
        //reqskBuilder.setK(reqK);
        //reqskBuilder.setLessThanOrEqual(true);
        List<Double> sortedData = null;
        int maxExpIter = 0;
        for (maxExp = maxExpStep; maxExp <= maxExpForDoubleBase10; maxExp+=maxExpStep) { // varying maxExp// intialize klls to collect errors of t-digest and of ReqSketch:
            KllDoublesSketch[] errorKllsClustering = new KllDoublesSketch[NumberOfPoints + 1];
            KllDoublesSketch[] errorKllsMerging = new KllDoublesSketch[NumberOfPoints + 1];
            KllDoublesSketch[] errorKllsRS = new KllDoublesSketch[NumberOfPoints + 1];
            for (int t = 0; t <= NumberOfPoints; t++) {
                errorKllsClustering[t] = new KllDoublesSketch(200);  // we do not need extreme quantiles
                errorKllsMerging[t] = new KllDoublesSketch(200);
                errorKllsRS[t] = new KllDoublesSketch(200);
            }
            // the following is similar to the IID generator
            for (int trial = 0; trial < T; trial++) { // trials
                if (DigestImpl.equals("Merging") || DigestImpl.trim().length() == 0) {
                    merging = new MergingDigest(Compression);
                    merging.setScaleFunction(scale);
                    if (scale == ScaleFunction.K_2_GLUED || scale == ScaleFunction.K_3_GLUED) {
                        merging.setUseAlternatingSort(
                            false); // fixing an issue with asymmetric scale functions (too few centroid in the end)
                    }
                }
                if (DigestImpl.equals("AVLTree") || DigestImpl.trim().length() == 0) {
                    clustering = new AVLTreeDigest(Compression);
                    clustering.setScaleFunction(scale);
                }
                reqsk = new ReqSketch(reqK, true, null);//reqskBuilder.build(); //
                reqsk.setLessThanOrEqual(true);
                /*if (WriteCentroidData && trial == T - 1) {
                    clustering.recordAllData(); // tracks centroids during the last trial
                }*/
                sortedData = new ArrayList<Double>();
                for (n = 0; n < N; n++) {
                    double item = generateItem();
                    sortedData.add(item);
                    //if (trial == T-1) data.add(item);
                    if (merging != null) {
                        merging.add(item);
                    }
                    if (clustering != null) {
                        clustering.add(item);
                    }
                    reqsk.update(item);
                }
                Collections.sort(sortedData);
                if (merging != null) {
                    merging.compress();
                }
                if (clustering != null) {
                    clustering.compress();
                }
                reqsk.compress();
                // extract error from t-digest
                for (int t = 0; t <= NumberOfPoints; t++) {
                    // t-digest error
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
                    // merging error
                    if (merging != null) {
                        double rEstM = merging.cdf(item) * N + 0.5;
                        //double relErr = 0;
                        double addErrM = 0;
                        if (rEstM < rTrueMin) {
                            //relErr = Math.abs(rTrueMin - rEst) / (N - rTrue + 1);
                            addErrM = (rEstM - rTrueMin) / N;
                        }
                        if (rEstM > rTrueMax) {
                            //relErr = Math.abs(rTrueMax - rEst) / (N - rTrue + 1);
                            addErrM = (rEstM - rTrueMax) / N;
                        }
                        errorKllsMerging[t].update(addErrM);
                    }
    
                    // clustering error
                    if (clustering != null) {
                        double rEstM = clustering.cdf(item) * N + 0.5;
                        //double relErr = 0;
                        double addErrM = 0;
                        if (rEstM < rTrueMin) {
                            //relErr = Math.abs(rTrueMin - rEst) / (N - rTrue + 1);
                            addErrM = (rEstM - rTrueMin) / N;
                        }
                        if (rEstM > rTrueMax) {
                            //relErr = Math.abs(rTrueMax - rEst) / (N - rTrue + 1);
                            addErrM = (rEstM - rTrueMax) / N;
                        }
                        errorKllsClustering[t].update(addErrM);
                    }
    
                    // ReqSketch error
                    double rEstRS = reqsk.getRank(item) * N;
                    double addErrRS = 0;
                    if (rEstRS < rTrueMin) {
                        addErrRS = (rEstRS - rTrueMin) / N;
                    }
                    if (rEstRS > rTrueMax) {
                        addErrRS = (rEstRS - rTrueMax) / N;
                    }
                    errorKllsRS[t].update(addErrRS);
                }
            }
            // aggregate err
            maxRelErrorRS_2SD[maxExpIter] = 0;
            maxRelErrorRS_median[maxExpIter] = 0;
            maxRelErrorMerging[maxExpIter] = 0;
            maxRelErrorClustering[maxExpIter] = 0;
            avgRelErrorRS_2SD[maxExpIter] = 0;
            avgRelErrorRS_median[maxExpIter] = 0;
            avgRelErrorMerging[maxExpIter] = 0;
            avgRelErrorClustering[maxExpIter] = 0;
            int numAvg = 0;
            for (int t = 0; t < NumberOfPoints; t++) { // no need to do it for t == NumberOfPoints
                double relErr =  Math.abs(errorKllsRS[t].getQuantile(0.5)) * NumberOfPoints / (NumberOfPoints - t);
                if (relErr > maxRelErrorRS_median[maxExpIter]) maxRelErrorRS_median[maxExpIter] = relErr;
                //if (t > NumberOfPoints * 0.6)
                avgRelErrorRS_median[maxExpIter] += relErr;
                relErr =  Math.abs(errorKllsRS[t].getQuantile(IIDgenerator.P2SD)) * NumberOfPoints / (NumberOfPoints - t);
                if (relErr > maxRelErrorRS_2SD[maxExpIter]) maxRelErrorRS_2SD[maxExpIter] = relErr;
                //if (t > NumberOfPoints * 0.6)
                avgRelErrorRS_2SD[maxExpIter] += relErr;
                if (merging != null) {
                    relErr =  Math.abs(errorKllsMerging[t].getQuantile(0.5)) * NumberOfPoints / (NumberOfPoints - t);
                    if (relErr > maxRelErrorMerging[maxExpIter]) maxRelErrorMerging[maxExpIter] = relErr;
                    avgRelErrorMerging[maxExpIter] += relErr;
                }
                if (clustering != null) {
                    relErr =  Math.abs(errorKllsClustering[t].getQuantile(0.5)) * NumberOfPoints / (NumberOfPoints - t);
                    if (relErr > maxRelErrorClustering[maxExpIter]) maxRelErrorClustering[maxExpIter] = relErr;
                    avgRelErrorClustering[maxExpIter] += relErr;
                }
            }
            avgRelErrorRS_2SD[maxExpIter] /= NumberOfPoints;
            avgRelErrorRS_median[maxExpIter] /= NumberOfPoints;
            avgRelErrorClustering[maxExpIter] /= NumberOfPoints;
            avgRelErrorMerging[maxExpIter] /= NumberOfPoints;
            System.out
                .println("finished experiment with maxExp = " + String.valueOf(maxExp));
            System.out.flush();
            
            maxExpIter++;
        }
        //Collections.sort(data);

        writeAggregatedResults(startTime, merging, clustering, reqsk, DigestStatsDir,
            DigestStatsDir + DigestStatsFileName + fileNamePart + "_" + DigestImpl
                + "_compr=" + String.valueOf(Compression) + "_" + scale.toString()
                + FileSuffix);
        // TODO write all properties into results
    }

    // the following does an additional trim and removes possible comment
    String getProperty(String propName) {
        String value = prop.getProperty(propName);
        int inx = value.indexOf('#');
        if (inx >= 0) {
            value = value.substring(0, inx);
        }
        return value.trim();
    }

    double generateItem() throws Exception { // only those with maxExp
        double item = 0;
        switch (Distribution) {
            case "loguniform":
                item = Math.pow(10, (rand.nextDouble() - 0.5) * 2 * maxExp);
                break;
            case "loguniform2":
                item = Math.pow(10, (Math.pow(rand.nextDouble(), 2) - 0.5) * 2 * maxExp);
                break;
            default:
                throw new Exception("Distribution '" + Distribution + "' undefined");
        }
        if (NegativeNumbers && rand.nextDouble() < 0.5) {
            item = -item;
        }
        return item;
    }


    public void writeAggregatedResults(Instant startTime, MergingDigest merging, AVLTreeDigest clustering, ReqSketch reqsk,
        String digestStatsDir, String outName) throws IOException {
        Files.createDirectories(Paths.get(digestStatsDir));
        System.out.printf("stats file:" + outName + "\n");

        File fout = new File(outName);
        fout.createNewFile();
        FileWriter fwout = new FileWriter(fout);

        //System.out.printf("computing rel. errors\n");
        //System.out.flush();

        fwout.write(
            "MaxExp;" + (merging != null ? "Merging MaxRE;" : "") + (clustering != null
                ? "Clustering MaxRE;" : "") + "ReqSketch MaxRE (median);ReqSketch MaxRE (+2SD);" + (merging != null ? "Merging AvgRE;" : "") + (clustering != null
                ? "Clustering AvgRE;" : "") + "ReqSketch AvgRE (median);ReqSketch AvgRE (+2SD)\n");
        int maxExpIter = 0;
        for (maxExp = maxExpStep; maxExp <= maxExpForDoubleBase10; maxExp+=maxExpStep) {
            if (merging != null && clustering != null) {
                fwout.write(String
                    .format("%d;%.6f;%.6f;%.6f;%.6f;%.6f;%.6f;%.6f;%.6f\n", maxExp,
                            maxRelErrorMerging[maxExpIter], maxRelErrorClustering[maxExpIter], maxRelErrorRS_median[maxExpIter], maxRelErrorRS_2SD[maxExpIter],
                            avgRelErrorMerging[maxExpIter], avgRelErrorClustering[maxExpIter], avgRelErrorRS_median[maxExpIter], avgRelErrorRS_2SD[maxExpIter]));
            }
            if (merging != null && clustering == null) {
                fwout.write(String
                        .format("%d;%.6f;%.6f;%.6f;%.6f;%.6f;%.6f\n", maxExp,
                                maxRelErrorMerging[maxExpIter], maxRelErrorRS_median[maxExpIter], maxRelErrorRS_2SD[maxExpIter],
                                avgRelErrorMerging[maxExpIter], avgRelErrorRS_median[maxExpIter], avgRelErrorRS_2SD[maxExpIter]));
            }
            if (merging == null && clustering != null) {
                fwout.write(String
                        .format("%d;%.6f;%.6f;%.6f;%.6f;%.6f;%.6f\n", maxExp,
                                maxRelErrorClustering[maxExpIter], maxRelErrorRS_median[maxExpIter], maxRelErrorRS_2SD[maxExpIter],
                                avgRelErrorClustering[maxExpIter], avgRelErrorRS_median[maxExpIter], avgRelErrorRS_2SD[maxExpIter]));
            }

            maxExpIter++;
        }
        fwout.close();

        if (merging != null) {
            fout = new File(outName + "_merging_centroids");
            System.out.printf("centroids file:" + outName + "_merging_centroids\n");
            fout.createNewFile();
            fwout = new FileWriter(fout);
            fwout.write(String.format("n=%d\n", N));
            fwout.write(String.format("scale func. = %s\n", merging.scale.toString()));
            fwout.write(String.format("delta = %d (compression param of t-digest)\n", Compression));
            fwout.write(String.format("# of centroids = %d\n", merging.centroids().size()));
            fwout.write(String.format("t-digest size in bytes = %d\n", merging.byteSize()));
            fwout.write(String.format("ReqSketch w/ k=%d size in bytes = %d\n", reqsk.getK(),
                reqsk.getSerializationBytes()));
            //Duration diff = Duration.between(startTime, Instant.now());
            /*String hms = String.format("%d:%02d:%02d", diff.toHours(), // doesn't work in JDK 1.8
                (int) (diff.toMinutes() % 60),
                (int) (diff.toSeconds() % 60));
            fwout.write(String.format("time taken = %s\n", hms));*/

            fwout.write("\nProperties:\n");
            for (Object key : prop.keySet()) {
                fwout.write(key + " = " + prop.getProperty(key.toString()) + "\n");
            }

            fwout.write("\nCentroids:\n");
            for (Centroid centr : merging.centroids()) {
                fwout.write(centr.toString() + "\n");
            }
            fwout.close();
        }
        if (clustering != null) {
            fout = new File(outName + "_clustering_centroids");
            System.out.printf("centroids file:" + outName + "_clustering_centroids\n");
            fout.createNewFile();
            fwout = new FileWriter(fout);
            fwout.write(String.format("n=%d\n", N));
            fwout.write(String.format("scale func. = %s\n", clustering.scale.toString()));
            fwout.write(String.format("delta = %d (compression param of t-digest)\n", Compression));
            fwout.write(String.format("# of centroids = %d\n", clustering.centroids().size()));
            fwout.write(String.format("t-digest size in bytes = %d\n", clustering.byteSize()));
            fwout.write(String.format("ReqSketch w/ k=%d size in bytes = %d\n", reqsk.getK(),
                reqsk.getSerializationBytes()));
//            Duration diff = Duration.between(startTime, Instant.now());
//            String hms = String.format("%d:%02d:%02d", diff.toHours(), % doesn't work in JDK 1.8
//                (int) (diff.toMinutes() % 60),
//                (int) (diff.toSeconds() % 60));
//            fwout.write(String.format("time taken = %s\n", hms));

            fwout.write("\nProperties:\n");
            for (Object key : prop.keySet()) {
                fwout.write(key + " = " + prop.getProperty(key.toString()) + "\n");
            }

            fwout.write("\nCentroids:\n");
            for (Centroid centr : clustering.centroids()) {
                fwout.write(centr.toString() + "\n");
            }
            fwout.close();
        }
    }


    @SuppressWarnings("unused")
    public static void main(final String[] args) throws Exception {
        for (int j = 0; j < args.length; j++) {
            new LoguniformWithVaryingMaxExpGenerator(args[j]);
        }
    }

}


