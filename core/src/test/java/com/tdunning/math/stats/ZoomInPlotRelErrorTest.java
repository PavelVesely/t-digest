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
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import org.junit.Ignore;


/**
 * 
 */
//@Ignore
public class ZoomInPlotRelErrorTest extends AbstractTest { 
    private static final int[] CompressionsForTesting = {100}; //100, 200, 500, 1000, 5000, 10000}; // param delta of t-digest
    private static final int NumberOfPoints = 200; // number of points where we probe the rank estimates
    
    // params for generating the input
    private static final int N = 50000; // 
    private static final int K = 200; // N / K should be roughly at most 200 (otherwise, we don't have enough precision)
    private static final int PrefixSize = 0; // number of points below zero in each iteration
    private static final int NumberOfRepeats = 100; // N * NumberOfRepeats is (approx.) the number of points in the final instance
    private static final String OutputFileName = "t-digest-genInput";
    private static final String OutputFileDir = "/tmp/"; // CHANGE AS APPROPRIATE 
    private static final String StatsFileDir = "/tmp/"; // CHANGE AS APPROPRIATE
    
    @BeforeClass
    public static void freezeSeed() {
        RandomUtils.useTestSeed();
    }

    @Test
    public void testZoomIn() throws FileNotFoundException, IOException {
        
      List<Double> sortedData = new ArrayList<Double>();
      List<Double> sortedDataPart = new ArrayList<Double>();
      List<Double> data = new ArrayList<Double>();
      PrintWriter w = new PrintWriter(OutputFileDir + OutputFileName);
      int n = 0;
      for (int r = 0; r < NumberOfRepeats; r++) {  // same instance is essentially repeated
        double max = Double.MAX_VALUE / (100*K + r); // initial interval max
        double min = 0; //Double.MIN_VALUE;
        int nn = 0;
        for (int i = -PrefixSize+1; i <= K - PrefixSize; i++) {//for (int i = 0; i <= K; i++) {
    		    double item = min + (i/(double)(K+1))*(max - min);
    		    data.add(item);
            sortedData.add(item);
            sortedDataPart.add(item);
            w.println(String.valueOf(item));
            nn++; n++;
    	  }
        int maxREindex = -1;
        while (nn < N) {
        	maxREindex = -Collections.binarySearch(sortedDataPart, Double.MIN_NORMAL) - 1;
        	min = 0; //maxREindex > 0 ? sortedDataPart.get(maxREindex - 1) : 0;
        	max = sortedDataPart.get(maxREindex);
        	//System.out.printf(String.format("phase maxRErTrue %d;\tmin %s;\tmax %s;\t means %d\n",
        	//		maxREindex, String.valueOf(min), String.valueOf(max), digest.centroidCount()));
        	if ((max - min) / (double)(K+1) < Double.MIN_VALUE) {
        		System.out.printf("TOO SMALL max - min\n");
        	}
        	for (int i = -PrefixSize; i <= K - PrefixSize; i++) {
        		double item = min + (i/(double)(K+1))*(max - min);
        		//if (item != 0) { // try to avoid duplicates
		    		data.add(item);
		        sortedData.add(item);
            sortedDataPart.add(item);
		        w.println(String.valueOf(item));
		        n++; nn++;
        		//}
        	}
          Collections.sort(sortedDataPart);
        }
        sortedDataPart.clear();
      }
      Collections.sort(sortedData);
      w.close();

            
            for (int compr : CompressionsForTesting) {  
                TDigest digest = new MergingDigest(compr);
                digest.setScaleFunction(ScaleFunction.K_2); //_GLUED;
                String outName = StatsFileDir + OutputFileName + "_K=" + String.valueOf(K) + "_N=" + String.valueOf(N)  +  "_PS=" + String.valueOf(PrefixSize) + "_repeats=" + String.valueOf(NumberOfRepeats) + "-stats-PP_" + String.valueOf(NumberOfPoints) + "_compr_" +  String.valueOf(compr) + "_" + digest.scale.toString();
                System.out.printf("stats file:" + outName + "\n");
                File fout = new File(outName);
                fout.createNewFile();
                FileWriter fwout = new FileWriter(fout);
                for (double item : data) {
                    digest.add(item);
                }
    	          digest.compress();
                //System.out.printf("computing rel. errors\n");
                //System.out.flush();
                
                fwout.write("true quantile;true rank;est. rank;rel. error;abs. error;item\n");
                for (int t = 0; t < NumberOfPoints; t++) {
                    //THE FOLLOWING IS EXTREMELY SLOW: Dist.cdf(item, sortedData);
                    int rTrue = (int)Math.ceil(t / (float)NumberOfPoints * n) + 1;
                    double item = sortedData.get(rTrue - 1);
                    // handling duplicate values -- rank is then rather an interval
                    int rTrueMin = rTrue;
                    int rTrueMax = rTrue;
                    while (rTrueMin >= 2 && item == sortedData.get(rTrueMin-2)) rTrueMin--;
                    while (rTrueMax < sortedData.size() && item == sortedData.get(rTrueMax)) rTrueMax++; 
                    double rEst = digest.cdf(item) * n + 0.5;
                    double relErr = 0;
                    double addErr = 0;
                    if (rEst < rTrueMin) {
                      relErr = Math.abs(rTrueMin - rEst) / rTrue;
                      addErr = (rEst - rTrueMin) / n;
                    }
                    if (rEst > rTrueMax) {
                      relErr = Math.abs(rTrueMax - rEst) / rTrue;
                      addErr = (rEst - rTrueMax) / n;                    
                    }
                    fwout.write(String.format("%.6f;%d;%.6f;%.6f;%.6f;%s\n", rTrue / (float)n, (int)rTrue, rEst, relErr, addErr, String.valueOf(item)));
                }
                fwout.write("\n");
                fwout.write(String.format("n=%d\n", n)); 
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
    }

}
