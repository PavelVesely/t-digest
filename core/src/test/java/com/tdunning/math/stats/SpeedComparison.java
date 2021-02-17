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

import com.tdunning.math.stats.datasketches.req.ReqSketch;
import com.tdunning.math.stats.datasketches.kll.KllDoublesSketch;

/**
 *
 */

public class SpeedComparison {

    public static void runSpeedComparison(String configFile) throws Exception {
        System.out.println("processing config file: " + configFile);
        Properties prop = new Properties();
        FileInputStream instream = new FileInputStream(configFile);
        prop.load(instream);

        // load properties
        int LgNmin = Integer.parseInt(getProperty(prop, "LgNmin")); // base 2
        int LgNmax = Integer.parseInt(getProperty(prop, "LgNmax")); // base 2
        int Compression = Integer.parseInt(getProperty(prop, "Compression")); // delta for t-digest
        ScaleFunction scale = ScaleFunction
            .valueOf(getProperty(prop, "ScaleFunction")); // ScaleFunction for t-digest
        String DigestStatsDir = getProperty(prop, "DigestStatsDir");
        String FileSuffix = getProperty(prop, "FileSuffix");
        int reqK = Integer.parseInt(getProperty(prop, "ReqK"));
        int kllK = Integer.parseInt(getProperty(prop, "KllK"));
        int reqKmax = Integer.parseInt(getProperty(prop,
            "ReqKmax")); // if reqKmax > reqK, try all even values of k in [reqK, reqKmax] on input of size 2^{lgNmax}

        Files.createDirectories(Paths.get(DigestStatsDir));
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss");
        //LocalDateTime now = LocalDateTime.now();
        String outName = DigestStatsDir + "speedComparison_LgNmin=" + LgNmin + "_LgNmax=" + LgNmax
            + "_Compr=" + Compression + "_Scale=" + scale + "_ReqK=" + reqK + FileSuffix;
        System.out.printf("stats file:" + outName + "\n");

        File fout = new File(outName);
        fout.createNewFile();
        FileWriter fwout = new FileWriter(fout);
        fwout.write("lgN;merging;clustering;ReqSketch;KLL\n");

        Random rand = new Random();
        for (int lgN = LgNmin - 1; lgN <= LgNmax;
            lgN++) { // starting from lgN = LgNmin - 1 as the first iteration is somewhat slower for some reason
            long N = 1l << lgN;
            MergingDigest merging = new MergingDigest(Compression);
            merging.setScaleFunction(scale);
            Instant startTime = Instant.now();
            for (long i = 0; i < N; i++) {
                merging.add(rand.nextDouble());
            }
            double mergingNs = Duration.between(startTime, Instant.now()).toNanos() / (double) N;

            AVLTreeDigest tree = new AVLTreeDigest(Compression);
            tree.setScaleFunction(scale);
            startTime = Instant.now();
            for (long i = 0; i < N; i++) {
                tree.add(rand.nextDouble());
            }
            double treeNs = Duration.between(startTime, Instant.now()).toNanos() / (double) N;

            ReqSketch reqsk = new ReqSketch(reqK, true, null);
            startTime = Instant.now();
            for (long i = 0; i < N; i++) {
                reqsk.update(rand.nextDouble());
            }
            double reqskNs = Duration.between(startTime, Instant.now()).toNanos() / (double) N;

            ReqSketch reqskLazy = new ReqSketch(reqK, true, null);
            reqskLazy.setFullLaziness(true);
            startTime = Instant.now();
            for (long i = 0; i < N; i++) {
                reqskLazy.update(rand.nextDouble());
            }
            double reqskLazyNs = Duration.between(startTime, Instant.now()).toNanos() / (double) N;

            KllDoublesSketch kll = new KllDoublesSketch(kllK);
            startTime = Instant.now();
            for (long i = 0; i < N; i++) {
                kll.update(rand.nextDouble());
            }
            double kllNs = Duration.between(startTime, Instant.now()).toNanos() / (double) N;

            if (lgN >= LgNmin) {
                fwout.write(String
                    .format("%d;%.2f;%.2f;%.2f;%.2f\n", lgN, mergingNs, treeNs, reqskNs, kllNs));
            }
        }
        fwout.close();
        if (reqKmax > reqK) {
            outName = DigestStatsDir + "speedComparison-laziness_LgNmax=" + LgNmax
                + "_ReqK=" + reqK + "_ReqKmax=" + reqKmax + FileSuffix;
            fout = new File(outName);
            fout.createNewFile();
            fwout = new FileWriter(fout);
            //fwout.write("\nReqSketch:\n");
            fwout.write("k;Partial laziness;Full laziness:\n");
            for (int k = reqK; k <= reqKmax; k += 2) {
                long N = 1l << LgNmax;
                ReqSketch reqsk = new ReqSketch(k, true, null);
                Instant startTime = Instant.now();
                for (long i = 0; i < N; i++) {
                    reqsk.update(rand.nextDouble());
                }
                double reqskNs = Duration.between(startTime, Instant.now()).toNanos() / (double) N;

                ReqSketch reqskLazy = new ReqSketch(k, true, null);
                reqskLazy.setFullLaziness(true);
                startTime = Instant.now();
                for (long i = 0; i < N; i++) {
                    reqskLazy.update(rand.nextDouble());
                }
                double reqskLazyNs =
                    Duration.between(startTime, Instant.now()).toNanos() / (double) N;

                fwout.write(String.format("%d;%.2f;%.2f\n", k, reqskNs, reqskLazyNs));
            }
            fwout.close();
        }

        //fwout.write("\nProperties:\n");
        //for (Object key : prop.keySet()) {
        //    fwout.write(key + " = " + prop.getProperty(key.toString()) + "\n");
        //}
        //fwout.close();
    }

    // the following does an additional trim and removes possible comment
    public static String getProperty(Properties prop, String propName) {
        String value = prop.getProperty(propName);
        int inx = value.indexOf('#');
        if (inx >= 0) {
            value = value.substring(0, inx - 1);
        }
        return value.trim();
    }

    @SuppressWarnings("unused")
    public static void main(final String[] args) throws Exception {
        for (int j = 0; j < args.length; j++) {
            runSpeedComparison(args[j]);
        }
    }


}


