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

/**
 *
 */

public class SpeedComparisonTest extends AbstractTest {
    static final int LgNmin = 20;
    static final int LgNmax = 30;
    
    @BeforeClass
    public static void freezeSeed() {
        RandomUtils.useTestSeed();
    }

    @Test
    public void speedComparison() throws Exception {
        
        System.out.println("lgN;merging;tree;reqSketch");
        Random rand = new Random();
        for (int lgN = LgNmin; lgN <= LgNmax; lgN++) {
            long N = 1l << lgN;
            MergingDigest merging = new MergingDigest(500);
            Instant startTime = Instant.now();
            for (long i = 0; i < N; i++) {
                merging.add(rand.nextDouble());
            }
            double mergingNs = Duration.between(startTime, Instant.now()).toNanos() / (double)N;

            AVLTreeDigest tree = new AVLTreeDigest(500);
            startTime = Instant.now();
            for (long i = 0; i < N; i++) {
                tree.add(rand.nextDouble());
            }
            double treeNs = Duration.between(startTime, Instant.now()).toNanos() / (double)N;
            
            ReqSketch reqsk = new ReqSketch(4, true, null);
            startTime = Instant.now();
            for (long i = 0; i < N; i++) {
                reqsk.update(rand.nextDouble());
            }
            double reqskNs = Duration.between(startTime, Instant.now()).toNanos() / (double)N;
            
            System.out.println(String.format("%d;%.2f;%.2f;%.2f", lgN, mergingNs, treeNs, reqskNs));
            System.out.flush();
        }
    }


}


