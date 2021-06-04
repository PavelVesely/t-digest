# Experiments with t-digest and ReqSketch

Both [t-digest](https://arxiv.org/abs/1902.04023) and [ReqSketch](https://arxiv.org/abs/2004.01668) are recently proposed algorithms for estimating quantiles of a dataset, or equivalently, for approximating the cumulative distribution function. While there are many more quantile estimation algorithms (such as the Greenwald-Khanna algorithm or the KLL sketch), the advantage of these two is that they aim for more accurate results for extreme quantile queries, such as the 99.5th percentile, thus giving a better understanding of the tails of the distribution, compared to many previous algorithms. These two algorithms are nevertheless very different: ReqSketch comes with formal mathematical analysis, ensuring that such "relative-error" guarantees hold even on worst-case inputs at the cost of the sketch size depending logarithmically on the input stream length. On the other hand, the space used by t-digest is constant (which is set a priori), but there are no accuracy guarantees for it. Furthermore, t-digest works with numerical data only, while ReqSketch can be applied to any data with a total ordering, such as lexicographically ordered strings (as ReqSketch applies comparisons to items only, while t-digest relies on averaging of centroids as data arrives, and on linear interpolation to produce accurate rank estimates).

The t-digest is already known to perform very well on data drawn from uniform or normal distribution (because linear interpolation works well on such data) -- see [the original paper](https://arxiv.org/abs/1902.04023) or [a comparison with the KLL sketch by the DataSketches library](https://datasketches.apache.org/docs/Quantiles/KllSketchVsTDigest.html). The purpose of this repository is to provide a broader comparison. Namely, we provide the following four experiments with the two algorithms:
1. A careful construction of a hard input (aka attack), which shows that t-digest may suffer a very large error. This implies that indeed there can be no formal accuracy worst-case guarantees for t-digest.
2. A generator of i.i.d. samples from a specified distribution, e.g.: uniform, normal, log-uniform, and log-uniform^2.
3. A study of the relative error of the two approaches for the log-uniform and log-uniform^2 distributions, as a key parameter (`maxExp`) varies.
4. A comparison of the runtime of t-digest (for both the merging and clustering variants), ReqSketch and KLL sketch.

These experiments can be used to reproduce results in the paper 'Theory meets Practice at the Median: a worst case comparison of relative error quantile algorithms' to appear in KDD 2021 (Applied Data Science Track); [arXiv preprint](https://arxiv.org/abs/2102.09299).

This repository is a clone of the [Ted Dunning's repository for t-digest](https://github.com/tdunning/t-digest), also incorporating [asymmetric scale functions](https://arxiv.org/abs/2005.09599) from [signal-fx/t-digest repository](https://github.com/signalfx/t-digest/tree/asymmetric). We included the implementations of ReqSketch and the KLL sketch from the [Apache DataSketches library](https://datasketches.apache.org/). Finally, we developed new generators of instances for testing accuracy of these two algorithms and one test for comparing the speed, which are described below.

## Running experiments

The first step is to compile the whole repository using [Apache Maven](https://maven.apache.org/) in the `core/` directory via, for example:

    $ mvn clean install

To generate all required data to reproduce all plots and tables in the paper, run the following script in the `core/` directory (see below for generating the plots):

    $ ./RunAllExperiments.sh

Note that some experiments take up to several hours to finish, and only after that one can [generate the plots](#creating-plots).
Experiment can be configured to require less time by reducing the number of trials or the input size, as described below.

### Adjusting experiments

Alternatively, one can run all experiments one by one, with a possibility to adjust parameters.
Before running any of the aforementioned experiments, we need to modify the CLASSPATH variable for Java:

    $ export CLASSPATH="$CLASSPATH:./target/classes:./target/test-classes"
Next, we run an experiment implemented as class `[Class]` (see below) with parameters specified in a configuration file using:

    $ java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.[Class] [configuration files]...

We provide prepared configuration files in the `core/resources/` directory (description of each parameter is inside the configuration file). The classes for the four tests are:

1. `com.tdunning.math.stats.CarefulAttack` for the careful construction of a hard input for t-digest. There are three configuration files available:
    - `core/resources/CarefulAttack_k_0_merging.conf` -- for scale function k_0 and the merging variant of t-digest,
    - `core/resources/CarefulAttack_k_0_clustering.conf` -- for scale function k_0 and the clustering variant of t-digest, and
    - `core/resources/CarefulAttack_k_3_clustering.conf` -- for scale function k_3 and the clustering variant of t-digest (this one was not used in the paper).
2. `com.tdunning.math.stats.IIDgenerator` for the generator of i.i.d. samples from a specified distribution; the prepared configuration file is `core/resources/IIDgenerator.conf`.
3. `com.tdunning.math.stats.LoguniformWithVaryingMaxExpGenerator` to study the behavior of the relative error as a key parameter (`maxExp`) of the log-uniform or log-uniform^2 distributions varies; the prepared configuration file is ` core/resources/LoguniformWithVaryingMaxExpGenerator.conf`.
4. `com.tdunning.math.stats.SpeedComparison` for the comparison of the runtime of t-digest (using either the merging or clustering variant), ReqSketch and KLL sketch.

Each of these tests outputs a CSV file with generated results to a specified directory (which is `data/results/` by default).

To generate the plots and tables in the paper, the following experiments should be run:
- scenario 1. above with both `CarefulAttack_k_0_merging.conf` and `CarefulAttack_k_0_clustering.conf` to reproduce Figure 2
- scenario 2. with:
  - `IIDgenerator.conf` as given to reproduce Figure 3
  - `IIDgenerator.conf` modified to use `Distribution=loguniform` to reproduce Figure 4
- scenario 3. with `LoguniformWithVaryingMaxExpGenerator.conf` as given to reproduce Figure 5
- scenario 4. with `resources/SpeedComparison.conf` as given to reproduce Table 1 and Figure 6
        
Note this can be sped up considerably by reducing the number of trials (parameter `LgT`) in the first two scenarios and reducing the input size (parameter `LgNmax`) in scenario 4,
which may, however, cause the experimental results to be less accurate.

## Creating plots

After generating results of the experiments as outlined above, from `docs/python/adversarial_plots`, run `make install` and then `make notebook`. Then run the entire [error_plots](docs/python/adversarial_plots/notebooks/error_plots.ipynb) notebook. This will both render the plots in the notebook, and save image files to `docs/python/adversarial_plots/images/`. Note that `RunAllExperiments.sh` will take some time to complete, so the data required by the notebooks will not be immediately available.

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
- ReqSketch and KLL are implemented by the DataSketches library using the `float` data type, however, t-digest uses `double`. As we needed to work with the `double` type in our experiments, we switched the implementations of these two algorithms to `double` (essentially just using find/replace). Also, we recommend using JDK 8 since the DataSketches repo requires it (and t-digest is compatible with it).
