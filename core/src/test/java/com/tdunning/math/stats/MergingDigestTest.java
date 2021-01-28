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

import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.util.Pair;
import org.apache.mahout.common.RandomUtils;
import org.apache.mahout.math.jet.random.AbstractContinousDistribution;
import org.apache.mahout.math.jet.random.Exponential;
import org.apache.mahout.math.jet.random.Uniform;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.Ignore;

import com.carrotsearch.randomizedtesting.annotations.Seed;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

//to freeze the tests with a particular seed, put the seed on the next line
//@Seed("84527677CF03B566:A6FF596BDDB2D59D")
@Seed("1CD6F48E8CA53BD1:379C5BDEB3A02ACB")
@Ignore
public class MergingDigestTest extends TDigestTest {
    @BeforeClass
    public static void setup() throws IOException {
        TDigestTest.setup("merge");
    }

    protected DigestFactory factory(final double compression) {
        return new DigestFactory() {
            @Override
            public TDigest create() {
                return new MergingDigest(compression);
            }
        };
    }

    @Before
    public void testSetUp() {
        RandomUtils.useTestSeed();
    }

    @Override
    protected TDigest fromBytes(ByteBuffer bytes) {
        return MergingDigest.fromBytes(bytes);
    }

    @Test
    public void writeUniformAsymmetricScaleFunctionResults() {
        try {
            writeAsymmetricScaleFunctionResults(Distribution.UNIFORM);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void writeExponentialAsymmetricScaleFunctionResults() {
        try {
            writeAsymmetricScaleFunctionResults(Distribution.EXPONENTIAL);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void writeAsymmetricScaleFunctionResults(Distribution distribution) throws Exception {

        List<ScaleFunction> scaleFcns = Arrays.asList(ScaleFunction.K_0, ScaleFunction.K_1,
                ScaleFunction.K_2, ScaleFunction.K_3, ScaleFunction.K_1_GLUED,
                ScaleFunction.K_2_GLUED, ScaleFunction.K_3_GLUED, ScaleFunction.K_QUADRATIC);

        int numTrials = 100;

        Map<String, Pair<ScaleFunction, Boolean>> digestParams = new HashMap<>();

        for (ScaleFunction fcn : scaleFcns) {
            if (fcn.toString().endsWith("GLUED") || fcn.toString().endsWith("QUADRATIC")) {
                digestParams.put(fcn.toString(), new Pair<>(fcn, false));
            } else {
                digestParams.put(fcn.toString() + "_USUAL", new Pair<>(fcn, false));
            }
        }
        writeSeveralDigestUniformResults(digestParams, numTrials, distribution,
            "../docs/asymmetric/data/merging/" + distribution.name() + "/");
    }

    private void writeSeveralDigestUniformResults(Map<String, Pair<ScaleFunction, Boolean>> digestParams,
        int numTrials, Distribution distribution, String writeLocation) throws Exception {

        int trialSize = 1_000_000;
        double compression = 100;
        double[] quants = new double[]{0.00001, 0.0001, 0.001, 0.01, 0.1,
                0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999};

        Map<String, List<Integer>> centroidCounts= new HashMap<>();

        Map<String, List<List<Integer>>> centroidSequences= new HashMap<>();


        for (Map.Entry<String, Pair<ScaleFunction, Boolean>> entry : digestParams.entrySet()) {
            centroidCounts.put(entry.getKey(), new ArrayList<Integer>());
            centroidSequences.put(entry.getKey(), new ArrayList<List<Integer>>());
            try {
                Map<Double, List<String>> records = new HashMap<>();
                for (double q : quants) {
                    records.put(q, new ArrayList<String>());
                }
                for (int j = 0; j < numTrials; j++) {
                        MergingDigest digest = (MergingDigest) factory(compression).create();
                        digest.setScaleFunction(entry.getValue().getFirst());
                        digest.setUseAlternatingSort(entry.getValue().getSecond());
                        Random rand = new Random();
                        AbstractContinousDistribution gen;
                        if (distribution.equals(Distribution.UNIFORM)) {
                            gen = new Uniform(50, 51, rand);
                        } else if (distribution.equals(Distribution.EXPONENTIAL)) {
                            gen = new Exponential(5, rand);
                        } else throw new Exception("distribution not specified");
                        double[] data = new double[trialSize];
                        for (int i = 0; i < trialSize; i++) {
                            data[i] = gen.nextDouble();
                            digest.add(data[i]);
                        }
                        Arrays.sort(data);
                        digest.compress();
                        for (double q : quants) {
                            double x1 = Dist.quantile(q, data);
                            double q1 = Dist.cdf(x1, data);
                            double q2 = digest.cdf(x1);
                            records.get(q).add(String.valueOf(Math.abs(q1 - q2)) + "," +
                                    String.valueOf(Math.abs(q1 - q2) / Math.min(q, 1 - q)) + "\n");
                        }
                        centroidCounts.get(entry.getKey()).add(digest.centroids().size());

                        List<Integer> seq = new ArrayList<>();
                        for (Centroid c : digest.centroids()) {
                            seq.add(c.count());
                        }
                        centroidSequences.get(entry.getKey()).add(seq);
                }
                for (double q : quants) {
                    FileWriter csvWriter = new FileWriter(writeLocation + entry.getKey() + "_" + String.valueOf(q) + ".csv");
                    csvWriter.append("error_q,norm_error_q\n");
                    for (String obs : records.get(q)) {
                        csvWriter.append(obs);
                    }
                    csvWriter.flush();
                    csvWriter.close();
                }

                FileWriter csvWriter = new FileWriter(writeLocation + entry.getKey()  + "_centroid_counts.csv");
                csvWriter.append("centroid_count\n");
                for (Integer ct : centroidCounts.get(entry.getKey())) {
                    csvWriter.append(ct.toString()).append("\n");
                }
                csvWriter.flush();
                csvWriter.close();


                FileWriter csvWriter2 = new FileWriter(writeLocation + entry.getKey()  + "_centroid_sizes.csv");
                for (List<Integer> ct : centroidSequences.get(entry.getKey())) {
                    for (Integer c : ct) {
                        csvWriter2.append(c.toString()).append(",");
                    }
                    csvWriter2.append("\n");
                }
                csvWriter2.flush();
                csvWriter2.close();

            } catch (IOException e) {
                System.out.println(e.toString());
                return;
            }
        }
    }
    // This test came from PR#145 by github user pulver
    @Test
    public void testNanDueToBadInitialization() {
        int compression = 30;
        int factor = 5;
        MergingDigest md = new MergingDigest(compression, (factor + 1) * compression, compression);

        final int M = 10;
        List<MergingDigest> mds = new ArrayList<>();
        for (int i = 0; i < M; ++i) {
            mds.add(new MergingDigest(compression, (factor + 1) * compression, compression));
        }

        // Fill all digests with values (0,10,20,...,80).
        List<Double> raw = new ArrayList<>();
        for (int i = 0; i < 9; ++i) {
            double x = 10 * i;
            md.add(x);
            raw.add(x);
            for (int j = 0; j < M; ++j) {
                mds.get(j).add(x);
                raw.add(x);
            }
        }
        Collections.sort(raw);

        // Merge all mds one at a time into md.
        for (int i = 0; i < M; ++i) {
            List<MergingDigest> singleton = new ArrayList<>();
            singleton.add(mds.get(i));
            md.add(singleton);
        }
//        md.add(mds);

//        Assert.assertFalse(Double.isNaN(md.quantile(0.01)));
        // Output
        System.out.printf("%4s\t%10s\t%10s\t%10s\t%10s\n", "q", "estimated", "actual", "error_cdf", "error_q");
        String dashes = "==========";

        System.out.printf("%4s\t%10s\t%10s\t%10s\t%10s\n", dashes.substring(0, 4), dashes, dashes, dashes, dashes);
        for (double q : new double[]{0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99}) {
            double est = md.quantile(q);
            double actual = Dist.quantile(q, raw);
            double qx = md.cdf(actual);
            Assert.assertEquals(q, qx, 0.08);
            Assert.assertEquals(est, actual, 3.5);
            System.out.printf("%4.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\n", q, est, actual, Math.abs(est - actual), Math.abs(qx - q));
        }
    }


    /**
     * Verifies interpolation between a singleton and a larger centroid.
     */
    @Test
    public void singleMultiRange() {
        TDigest digest = factory(50).create();
        digest.setScaleFunction(ScaleFunction.K_0);
        for (int i = 0; i < 100; i++) {
            digest.add(1);
            digest.add(2);
            digest.add(3);
        }
        // this check is, of course true, but it also forces merging before we change scale
        assertTrue(digest.centroidCount() < 300);
        digest.add(0);
        // we now have a digest with a singleton first, then a heavier centroid next
        Iterator<Centroid> ix = digest.centroids().iterator();
        Centroid first = ix.next();
        Centroid second = ix.next();
        assertEquals(1, first.count());
        assertEquals(0, first.mean(), 0);
//        assertTrue(second.count() > 1);
        assertEquals(1.0, second.mean(), 0);

        assertEquals(0.5 / digest.size(), digest.cdf(0), 0);
        assertEquals(1.0 / digest.size(), digest.cdf(1e-10), 1e-10);
        assertEquals(1.0 / digest.size(), digest.cdf(0.25), 1e-10);
    }

    /**
     * Make sure that the first and last centroids have unit weight
     */
    @Test
    public void testSingletonsAtEnds() {
        TDigest d = new MergingDigest(50);
        d.recordAllData();
        Random gen = new Random(1);
        double[] data = new double[100];
        for (int i = 0; i < data.length; i++) {
            data[i] = Math.floor(gen.nextGaussian() * 3);
        }
        for (int i = 0; i < 100; i++) {
            for (double x : data) {
                d.add(x);
            }
        }
        int last = 0;
        for (Centroid centroid : d.centroids()) {
            if (last == 0) {
                assertEquals(1, centroid.count());
            }
            last = centroid.count();
        }
        assertEquals(1, last);
    }
}
