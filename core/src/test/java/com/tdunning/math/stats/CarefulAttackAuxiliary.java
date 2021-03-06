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

import java.io.*;
import java.lang.Math;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

import com.tdunning.math.stats.datasketches.req.ReqSketch;
import com.tdunning.math.stats.datasketches.kll.KllDoublesSketch;

/**
 *
 */
public class CarefulAttackAuxiliary {

    protected static final int NumberOfPoints = 200; // number of points where we probe the rank estimates
    protected static final int NumberOfTrials = 1000; // number of executions of ReqSketch/KLL on data
    protected static final int reqK = 4; // accuracy param. of ReqSketch
    protected static final int kllK = 100; // accuracy param. of KLL
    protected static final String DigestStatsDir = "../../../data/results/";
    protected static final String FileSuffix = ".csv";

    protected static TDigest digest(final double delta, String implementation,
        ScaleFunction scaleFunction,
        boolean useAlternatingSort, boolean recordAllData, long seed) {
        TDigest digest;
        switch (implementation.toLowerCase()) {
            case "merging":
                digest = new MergingDigest(delta);
                ((MergingDigest) digest).setUseAlternatingSort(useAlternatingSort);
                break;
            case "tree":
                digest = new AVLTreeDigest(delta);
                ((AVLTreeDigest) digest).gen.setSeed(seed);
                break;
            default:
                digest = new AVLTreeDigest(-1d);
        }
        try {
            digest.setScaleFunction(scaleFunction);
        } catch (IllegalArgumentException e) {
            digest.setUnnormalizedScaleFunction(scaleFunction);
        }
        if (recordAllData) {
            digest.recordAllData();
        }
        return digest;
    }

    // compareTo should be ReqSketch or KLL
    protected static void writeResults(int compr, int size, TDigest digest,
        List<Double> data, List<Double> sortedData, String compareTo, String digestStatsDir,
        String outName, boolean writeCentroids)
        throws
        IOException, Exception {
        Files.createDirectories(Paths.get(digestStatsDir));
        System.out.printf("stats file:" + outName + "\n");
        File fout = new File(outName);
        fout.createNewFile();
        FileWriter fwout = new FileWriter(fout);
        KllDoublesSketch[] errorKlls = null;
        if (compareTo.equalsIgnoreCase("reqsketch")) {
            System.out.printf("running ReqSketch on data\n");
            errorKlls = runReqSketchOnData(data, sortedData);
        }
        if (compareTo.equalsIgnoreCase("kll")) {
            System.out.printf("running KLL on data\n");
            errorKlls = runKllOnData(data, sortedData);
        }
        //else System.out.printf("nothing to compare to: '" + compareTo + "'\n"); 

        System.out.printf("\n");
        System.out.printf(String.format("n=%d\n", size));
        System.out.printf(String.format("scale func. = %s\n", digest.scale.toString()));
        System.out.printf(String.format("delta = %d\n", compr));
        System.out.printf(String.format("# of centroids = %d\n", digest.centroids().size()));
        System.out.printf(String.format("size in bytes = %d\n", digest.byteSize()));

        System.out.printf("computing errors\n");
        //System.out.flush();
        if (errorKlls == null) {
            fwout.write("true quantile;true rank;est. rank;rel. error;abs. error;item\n");
        } else {
            fwout.write(String
                .format("true quantile;TD abs. error;%s -2SD error; %s +2SD error;item\n",
                    compareTo, compareTo));
        }
        for (int t = 0; t <= NumberOfPoints; t++) {
            //THE FOLLOWING IS EXTREMELY SLOW: Dist.cdf(item, sortedData);
            int rTrue = (int) Math.ceil(t / (float) NumberOfPoints * size) + 1;
            if (rTrue > size) {
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
            if (errorKlls == null) {
                fwout.write(String
                    .format("%.6f;%d;%.6f;%.6f;%.6f;%s\n", rTrue / (float) size, (int) rTrue, rEst,
                        relErr, addErr, String.valueOf(item)));
            } else {
                double addErrRSM2SD = errorKlls[t].getQuantile(IIDgenerator.M2SD);
                double addErrRSP2SD = errorKlls[t].getQuantile(IIDgenerator.P2SD);
                fwout.write(String
                    .format("%.6f;%.6f;%.6f;%.6f;%s\n", rTrue / (float) size, addErr, addErrRSM2SD,
                        addErrRSP2SD, String.valueOf(item)));
            }
        }

        if (writeCentroids) {
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
        }
        fwout.close();
        System.out.flush();

    }

    protected static KllDoublesSketch[] runReqSketchOnData(List<Double> data,
        List<Double> sortedData) throws Exception {
        KllDoublesSketch[] errorKlls = new KllDoublesSketch[NumberOfPoints
            + 1]; // TODO use t-digest or something else?
        for (int t = 0; t <= NumberOfPoints; t++) {
            errorKlls[t] = new KllDoublesSketch(200); // we do not need extreme quantiles
        }
        int N = data.size();
        ReqSketch reqsk = null;
        for (int trial = 0; trial < NumberOfTrials; trial++) { // trials
            reqsk = new ReqSketch(reqK, true, null);//reqskBuilder.build(); //
            reqsk.setLessThanOrEqual(true);
            for (double item : data) {
                reqsk.update(item);
            }
            reqsk.compress();
            for (int t = 0; t <= NumberOfPoints; t++) {
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
                // ReqSketch error
                double rEstRS = reqsk.getRank(item) * N;
                double addErrRS = 0;
                if (rEstRS < rTrueMin) {
                    addErrRS = (rEstRS - rTrueMin) / N;
                }
                if (rEstRS > rTrueMax) {
                    addErrRS = (rEstRS - rTrueMax) / N;
                }
                errorKlls[t].update(addErrRS);
            }
        }
        System.out.printf(String.format("ReqSketch w/ k=%d size in bytes = %d\n", reqsk.getK(),
            reqsk.getSerializationBytes()));
        return errorKlls;
    }

    protected static KllDoublesSketch[] runKllOnData(List<Double> data, List<Double> sortedData)
        throws Exception {
        KllDoublesSketch[] errorKlls = new KllDoublesSketch[NumberOfPoints
            + 1]; // TODO use t-digest or something else?
        for (int t = 0; t <= NumberOfPoints; t++) {
            errorKlls[t] = new KllDoublesSketch(200);
        }
        int N = data.size();
        KllDoublesSketch kll = null;
        for (int trial = 0; trial < NumberOfTrials; trial++) { // trials
            kll = new KllDoublesSketch(kllK);
            for (double item : data) {
                kll.update(item);
            }
            //kll.compress(); // no such method available
            for (int t = 0; t <= NumberOfPoints; t++) {
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
                // KLL error
                double rEstRS = kll.getRank(item) * N;
                double addErrRS = 0;
                if (rEstRS < rTrueMin) {
                    addErrRS = (rEstRS - rTrueMin) / N;
                }
                if (rEstRS > rTrueMax) {
                    addErrRS = (rEstRS - rTrueMax) / N;
                }
                errorKlls[t].update(addErrRS);
            }
        }
        System.out.printf(String
            .format("KLL w/ k=%d size in bytes = %d\n", kll.getK(), kll.getSerializedSizeBytes()));
        return errorKlls;
    }

    protected static void writeResultsFloat(int compr, int size, TDigest digest,
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

    protected static void writeCentroidData(TDigest digest, String outName) throws IOException {
        Files.createDirectories(Paths.get(DigestStatsDir));

        System.out.printf("stats file:" + outName + "\n");
        File fout = new File(outName);
        fout.createNewFile();
        FileWriter fwout = new FileWriter(fout);

        int index = 0;
        for (Centroid c : digest.centroids()) {
            if (c.data() != null) {
                fwout.write(index + ",");
                for (Double d : c.data()) {
                    fwout.write(d + ",");
                }
                index++;
                fwout.write("\n");
            } else {
                fwout.write(c.mean() + ";" + c.count() + "\n");
            }
        }
        fwout.close();
        System.out.flush();
    }

    protected static TDigest rerunOnSortedInput(TDigest digest, List<Double> sortedData)
        throws Exception {

        TDigest freshDigest;

        double delta = digest.compression();
        if (digest instanceof MergingDigest) {
            freshDigest = new MergingDigest(delta);
            ((MergingDigest) freshDigest)
                .setUseAlternatingSort(((MergingDigest) digest).useAlternatingSort);
        } else if (digest instanceof AVLTreeDigest) {
            freshDigest = new AVLTreeDigest(delta);
        } else {
            throw new Exception("digest has no implementation");
        }

        freshDigest.setScaleFunction(digest.scale);

        if (digest.isRecording()) {
            freshDigest.recordAllData();
        }

        for (double item : sortedData) {
            freshDigest.add(item);
        }
        freshDigest.compress();

        return freshDigest;

    }

}


