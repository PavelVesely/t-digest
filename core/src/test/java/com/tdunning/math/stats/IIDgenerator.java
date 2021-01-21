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

import org.junit.Ignore;


/**
 *
 */

public class IIDgenerator  {

    // properties from config file
    final String Distribution;
    final int N;
    final int LgN; // base 2
    final int NumClusters; // for clustered distribution
    final int NumberOfPoints; // number of points where we probe the rank estimates
    final Boolean NegativeNumbers; // if false, generate positive numbers only
    final Boolean WriteCentroidData; 
    final double lambda; // for exponential distribution
    final int Compression; // delta for t-digest
    final ScaleFunction scale; // ScaleFunction for t-digest
    final String InputStreamFileName;
    final String InputStreamFileDir;
    final String DigestStatsFileName;
    final String DigestStatsDir;
    final String FileSuffix;
    
    // vars for experiments
    Random rand;
    Properties prop;
    int maxExpBase2;
    int n;

    
    public IIDgenerator(final String configFile) throws Exception {
        System.out.println("processing config file: " + configFile);
        prop = new Properties();
        FileInputStream instream = new FileInputStream(configFile);
        prop.load(instream);
        
        // load properties
        Distribution = getProperty("Distribution");
        LgN = Integer.parseInt(getProperty("LgN")); // base 2
        N = 1 << LgN;
        NumClusters = Integer.parseInt(getProperty("NumClusters"));
        NumberOfPoints = Integer.parseInt(getProperty("NumberOfPoints")); // number of points where we probe the rank estimates
        NegativeNumbers = Boolean.parseBoolean(getProperty("NegativeNumbers")); // if false, generate positive numbers only
        WriteCentroidData = Boolean.parseBoolean(getProperty("WriteCentroidData")); 
        lambda = Double.parseDouble(getProperty("Lambda")); // for exponential distribution
        Compression = Integer.parseInt(getProperty("Compression")); // delta for t-digest
        scale = ScaleFunction.valueOf(getProperty("ScaleFunction")); // ScaleFunction for t-digest
        InputStreamFileName = getProperty("InputStreamFileName");
        InputStreamFileDir = getProperty("InputStreamFileDir");
        DigestStatsFileName = getProperty("DigestStatsFileName");
        DigestStatsDir = getProperty("DigestStatsDir");
        FileSuffix = getProperty("FileSuffix");

        List<Double> sortedData = new ArrayList<Double>();
        List<Double> data = new ArrayList<Double>();
        Files.createDirectories(Paths.get(InputStreamFileDir));
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss");  
        LocalDateTime now = LocalDateTime.now();
        String fileNamePart = "_" + dtf.format(now) + "_" + Distribution + (NegativeNumbers ? "_wNegativeNumbers" : "_PositiveOnly") + "_lgN=" + String.valueOf(LgN);
        String inputFilePath =
            InputStreamFileDir + InputStreamFileName + "_" + fileNamePart + FileSuffix;
        PrintWriter w = new PrintWriter(inputFilePath);
        rand = new Random();
        
        maxExpBase2 = (int) (Math.log(Double.MAX_VALUE / 10000) / Math.log(2));
        for (n = 0; n < N; n++) {
            double item = generateItem();
            data.add(item);
            sortedData.add(item);
            w.println(String.valueOf(item));
        }
        Collections.sort(sortedData);
        w.close();
        System.out.println("file with generated input: " + inputFilePath);
        System.out.flush();
      
        TDigest digest = new AVLTreeDigest(Compression);
        digest.setScaleFunction(scale);
        if (WriteCentroidData) 
            digest.recordAllData(); // tracks centroids
        for (double item : data) {
            digest.add(item);
        }
        digest.compress();
        System.out
            .println("processing by t-digest done for compression =" + String.valueOf(Compression));
        System.out.flush();

        if (WriteCentroidData)
            writeCentroidData(digest, DigestStatsDir,
                  DigestStatsDir + DigestStatsFileName + fileNamePart
                    + "_compr=" + String.valueOf(Compression) + digest.scale.toString() + FileSuffix);

        writeResults(Compression, n, NumberOfPoints, prop, digest, sortedData, DigestStatsDir,
            DigestStatsDir + DigestStatsFileName + fileNamePart + "-stats-PP=" + String.valueOf(NumberOfPoints)
                + "_compr=" + String.valueOf(Compression) + "_" + digest.scale.toString() + FileSuffix);
        // TODO write all properties into results
    }
    
    // the following does an additional trim and removes possible comment
    String getProperty(String propName) {
        String value = prop.getProperty(propName);
        int inx = value.indexOf('#');
        if (inx >= 0)
             value = value.substring(0, inx-1);
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
                item = (0.9999 + rand.nextDouble() / 100000000) * (rand.nextInt(NumClusters/2) + 1) * (Double.MAX_VALUE / (N*NumClusters));
                break;
            case "uniform":
                item = rand.nextDouble();
                break;
            default:
                throw new Exception("Distribution '" + Distribution + "' undefined");
        }
        if (NegativeNumbers && rand.nextDouble() < 0.5) {
            item = -item;
        }
        return item;
    }
    

            
    public static void writeResults(int compr, int size, int numPoints, Properties prop, TDigest digest,
        List<Double> sortedData, String digestStatsDir, String outName) throws
        IOException {
        Files.createDirectories(Paths.get(digestStatsDir));
        System.out.printf("stats file:" + outName + "\n");
        File fout = new File(outName);
        fout.createNewFile();
        FileWriter fwout = new FileWriter(fout);

        //System.out.printf("computing rel. errors\n");
        //System.out.flush();

        fwout.write("true quantile;true rank;est. rank;rel. error;abs. error;item\n");
        for (int t = 0; t <= numPoints; t++) {
            //THE FOLLOWING IS EXTREMELY SLOW: Dist.cdf(item, sortedData);
            int rTrue = (int) Math.ceil(t / (float) numPoints * size) + 1;
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
                relErr = Math.abs(rTrueMin - rEst) / (size - rTrue + 1);
                addErr = (rEst - rTrueMin) / size;
            }
            if (rEst > rTrueMax) {
                relErr = Math.abs(rTrueMax - rEst) / (size - rTrue + 1);
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
        fwout.write("\nProperties:\n");
        for (Object key: prop.keySet()) {
            fwout.write(key + " = " + prop.getProperty(key.toString()) + "\n");
        }
        
        fwout.write("\nCentroids:\n");
        for (Centroid centr : digest.centroids()) {
            fwout.write(centr.toString() + "\n");
        }
        fwout.close();
        System.out.flush();

    }

    public static void writeCentroidData(TDigest digest, String digestStatsDir, String outName) throws IOException {
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

