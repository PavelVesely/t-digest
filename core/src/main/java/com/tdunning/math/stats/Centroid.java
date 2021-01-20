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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * A single centroid which represents a number of data points.
 */
public class Centroid implements Comparable<Centroid>, Serializable {

    private static final AtomicInteger uniqueCount = new AtomicInteger(1);

    //private double centroid = 0;
    private BigDecimal centroid = BigDecimal.valueOf(0d);
    private int count = 0;

    // The ID is transient because it must be unique within a given JVM. A new
    // ID should be generated from uniqueCount when a Centroid is deserialized.
    private transient int id;

    //private List<Double> actualData = null;
    private List<BigDecimal> actualData = null;

    private Centroid(boolean record) {
        id = uniqueCount.getAndIncrement();
        if (record) {
            actualData = new ArrayList<>();
        }
    }

    private static List<Double> convertToDouble(List<BigDecimal> d) {
        List<Double> _ = new ArrayList<>(d.size());
        for (BigDecimal __ : d) {
            _.add(__.doubleValue());
        }
        return _;
    }

    private static List<BigDecimal> convertToBig(List<Double> d) {
        List<BigDecimal> _ = new ArrayList<>(d.size());
        for (Double __ : d) {
            _.add(BigDecimal.valueOf(__));
        }
        return _;
    }


    // conuld convert these to use Big constructor
    public Centroid(double x) {
        this(false);
        start(x, 1, uniqueCount.getAndIncrement());
    }

    public Centroid(double x, int w) {
        this(false);
        start(x, w, uniqueCount.getAndIncrement());
    }

    public Centroid(double x, int w, int id) {
        this(false);
        start(x, w, id);
    }

    public Centroid(double x, int id, boolean record) {
        this(record);
        start(x, 1, id);
    }

    Centroid(double x, int w, List<Double> data) {
        this(x, w);
        actualData = convertToBig(data);
    }


    public Centroid(BigDecimal x) {
        this(false);
        start(x, 1, uniqueCount.getAndIncrement());
    }

    public Centroid(BigDecimal x, int w) {
        this(false);
        start(x, w, uniqueCount.getAndIncrement());
    }

    public Centroid(BigDecimal x, int w, int id) {
        this(false);
        start(x, w, id);
    }

    public Centroid(BigDecimal x, int id, boolean record) {
        this(record);
        start(x, 1, id);
    }


    Centroid(BigDecimal x, int w, List<BigDecimal> data) {
        this(x, w);
        actualData = data;
    }


    private void start(double x, int w, int id) {
        this.id = id;
        add(x, w);
    }


    private void start(BigDecimal x, int w, int id) {
        this.id = id;
        add(x, w);
    }

    public void add(double x, int w) {
        add(BigDecimal.valueOf(x), w);
    }


    public void add(BigDecimal x, int w) {
        if (actualData != null) {
            actualData.add(x);
        }
        count += w;
        //centroid += w * (x - centroid) / count;
        //System.out.println("\n x: " + x.toString());
        //System.out.println(toString());
        x = x.subtract(centroid);
        //System.out.println("x after subtraction: " + x.toString());
        //System.out.println(toString());
        BigDecimal incr = BigDecimal.valueOf(w).multiply(x).divide(BigDecimal.valueOf(count));
        //System.out.println("incr: " + incr.toString());
        centroid = centroid.add(incr);
        //System.out.println("add: " + add.toString());
        //System.out.println(toString() + "\n");

    }

    public double mean() {
        return centroid.doubleValue();
    }

    public BigDecimal bigMean() {
        return centroid;
    }

    public int count() {
        return count;
    }

    public int id() {
        return id;
    }

    @Override
    public String toString() {
        return "Centroid{" +
            "centroid=" + centroid +
            ", count=" + count +
            '}';
    }

    @Override
    public int hashCode() {
        return id;
    }

    @Override
    public int compareTo(@SuppressWarnings("NullableProblems") Centroid o) {
        //int r = Double.compare(centroid, o.centroid);
//
//        int r = centroid.compareTo(o.centroid);
//        if (r == 0) {
//            r = id - o.id;
//        }
//        return r;
        int r = centroid.compareTo(o.centroid);
        if ((centroid.subtract(o.centroid).abs().compareTo(BigDecimal.valueOf(Double.MIN_VALUE)) < 0)) {
            System.out.println("randomness intensifies");
            r = id - o.id;
        }
        return r;
    }


    public List<Double> data() {
        return convertToDouble(actualData);
    }


    public List<BigDecimal> bigData() {
        return actualData;
    }

    @SuppressWarnings("WeakerAccess")
    public void insertData(BigDecimal x) {
        if (actualData == null) {
            actualData = new ArrayList<>();
        }
        actualData.add(x);
    }

    @SuppressWarnings("WeakerAccess")
    public void insertData(double x) {
        insertData(BigDecimal.valueOf(x));
    }


    public static Centroid createWeighted(double x, int w, Iterable<? extends Double> data) {
        Centroid r = new Centroid(data != null);
        r.add(x, w, data);
        return r;
    }

    public void add(double x, int w, Iterable<? extends Double> data) { ///aiyyy?A??A?A
//        Iterable<? extends BigDecimal> bigData = new Iterable<BigDecimal>()
//        {
//
//            @Override
//            public Iterator<BigDecimal> iterator() {
//                return data_;
//            }
//        };
        add(BigDecimal.valueOf(x), w, data);
    }

    public void add(BigDecimal x, int w, Iterable<? extends Double> data) {
        if (actualData != null) {
            if (data != null) {
                for (Double old : data) {
                    actualData.add(BigDecimal.valueOf(old));
                }
            } else {
                actualData.add(x);
            }
        }
        centroid = AbstractTDigest.weightedAverage(centroid, count, x, w);
        count += w;
    }


    private void readObject(ObjectInputStream in) throws ClassNotFoundException, IOException {
        in.defaultReadObject();
        id = uniqueCount.getAndIncrement();
    }
}
