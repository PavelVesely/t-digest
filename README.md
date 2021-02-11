# Experiments with t-digest and ReqSketch

Both [t-digest](https://arxiv.org/abs/1902.04023) and [ReqSketch](https://arxiv.org/abs/2004.01668) are recently proposed algorithms for estimating quantiles of a dataset, or equivalently, for approximating the cumulative distribution function. While there many more quantile estimation algorithms (such the Greenwald-Khanna algorithm or the KLL sketch), the advantage of these two is that they aim for more accurate results for extreme quantile queries, such as the 99.5th percentile, thus giving a better understanding of the tails of the distribution, compared to many previous algorithms. These two algorithms are nevertheless very different: ReqSketch comes with formal mathematical analysis, ensuring that such "relative-error" guarantees hold even on worst-case inputs at the cost of the sketch size depending logarithmically on the input stream length. On the other hand, the space used by t-digest is constant (which is set a priori), but there are no accuracy guanrantees for it. Furthermore, t-digest works with numerical data only, while ReqSketch can be applied to any data with a total ordering, such as lexicographically ordered strings (as ReqSketch applies comparisons to items only, while t-digest relies on linear interpolation to produce accurate rank estimates).

The t-digest is already known to perform very well on data drawn from uniform or normal distribution (because linear interpolation works well on such data) -- see [the original paper](https://arxiv.org/abs/1902.04023) or [a comparison with the KLL sketch by the DataSketches library](https://datasketches.apache.org/docs/Quantiles/KllSketchVsTDigest.html). The purpose of this repository is to provide a broader comparison. Namely, we provide the following three experiments with the two algorithms:
1. A careful construction of a hard input (aka attack), which shows that t-digest may suffer a very large error. This implies that indeed there can be no formal accuracy worst-case guanrantees for t-digest.
2. A generator of i.i.d. samples from a specified distribution, e.g.: uniform, normal, log-uniform, and log-uniform^2
3. A comparison of the runtime of t-digest (for both the merging and clustering variants), ReqSketch and KLL sketch.

These experiments can be used to reproduce results in the paper 'Theory meets Practice: worst case behavior of quantile algorithms' (to be available at arXiv soon).

This repository is a clone of the [Ted Dunning's repository for t-digest](https://github.com/tdunning/t-digest), also incorporating [asymmetric scale functions](https://arxiv.org/abs/2005.09599) from [signal-fx/t-digest repository](https://github.com/signalfx/t-digest/tree/asymmetric). We included the implementations of ReqSketch and the KLL sketch from the [Apache DataSketches library](https://datasketches.apache.org/). Finally, we developed new generators of instances for testing accuracy of these two algorithms and one test for comparing the speed, which are described below.

## Running experiments

The first step is to compile the whole repository using [Apache Maven](https://maven.apache.org/) in the `core/` directory via, for example (possibly, one can skip running the unit tests using `-DskipTests=true`):

    $ mvn clean install

Before running any of the aforementioned three experiments, we need to modify the CLASSPATH variable for Java:

    $ export CLASSPATH="$CLASSPATH:./target/classes:./target/test-classes:./target/classes/org/apache/datasketches/req/ReqSketch:../../datasketches/target/classes:../../../.m2/repository/org/apache/datasketches/datasketches-java/1.3.0-incubating/datasketches-java-1.3.0-incubating.jar"
Next, we run an experiment implemented as class `[Class]` (see below) with parameters specified in a configuration file using:

    $ java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.[Class] [configuration files]...

We provide prepared configuration files in the `core/resources/` directory (options for each parameter are listed in the comment). The classes for the three tests are:

1. `com.tdunning.math.stats.CarefulAttack` for the careful construction of a hard input for t-digest. There are three configuration files available:
    - `core/resources/CarefulAttack_k_0_merging.conf` -- for scale function k_0 and the merging variant of t-digest,
    - `core/resources/CarefulAttack_k_0_clustering.conf` -- for scale function k_0 and the clustering variant of t-digest, and
    - `core/resources/CarefulAttack_k_3_clustering.conf` -- for scale function k_3 and the clustering variant of t-digest (this one was not used in the paper).
2. `com.tdunning.math.stats.IIDgenerator` for the generator of i.i.d. samples from a specified distribution; the prepared configuration file is `core/resources/IIDgenerator.conf`.
3. `com.tdunning.math.stats.SpeedComparison` for the comparison of the runtime of t-digest (using either the merging or clustering variant), ReqSketch and KLL sketch.

Each of these tests outputs a CSV file with generated results to a specified directory (which is `data/results/` by default).

To generate the plots and tables in the paper, the following experiments should be run:
- scenario 1. above with both `CarefulAttack_k_0_merging.conf` and `CarefulAttack_k_0_clustering.conf`
- scenario 2. with:
  - `IIDgenerator.conf` as given
  - `IIDgenerator.conf` modified to use `Distribution=loguniform`
  - `IIDgenerator.conf` modified to use `MaxExp=10` (keeping `Distribution=loguniform2`)
- scenario 3. with `resources/SpeedComparison.conf` as given to reproduce Table 1 and then with `ReqKmax=50` to reproduce Figure 5
        
Note this can be sped up considerably by reducing the number of trials (parameter `LgT`).

Then from `docs/python/adversarial_plots`, run `make install` and then `make notebook`. Then run the entire [error_plots](docs/python/adversarial_plots/notebooks/error_plots.ipynb) notebook. This will both render the plots in the notebook, and save image files to `docs/python/adversarial_plots/images/`.

Run also the [speed_comparison](docs/python/adversarial_plots/notebooks/speed_comparison.ipynb) notebook, which plots the data resulting from scenario 3. above.

## Additional plots

As the weak ordering of centroids plays a role in the t-digest error, we provide some plots to understand the nature of the weak ordering in the various scenarios.

Run scenario 1. (the careful attacks) with `WriteCentroidData=true, CompareToSorted=true`. (This will overwrite outputs with identical files if you have already run the experiment with the given configuration.)
Run scenario 2. (the IID case) with `writeCentroidData=true, DigestImpl=`.

These generate outputs which are used by both the [overlap_computation](docs/python/adversarial_plots/notebooks/overlap_computation.ipynb) notebook and a cell in the [error_plots](docs/python/adversarial_plots/notebooks/error_plots.ipynb) notebook demonstrating the (theoretically clear) result that the carefully constructed input is not difficult for the t-digest when presented in sorted order.

The "local overlap" plots do not appear in the paper.

## Remarks

- Merged all commits till Jan 28, 2021 from the main branch of [Ted Dunning's repository for t-digest](https://github.com/tdunning/t-digest).
- Code taken from the DataSketches library is as of Jan 29, 2021.
- ReqSketch and KLL are implemented by the DataSketches library using the `float` data type, however, t-digest uses `double`. As we needed to work with the `double` type in our experiments, we switched the implementations of these two algorithms to `double` (essentially just using find/replace).
