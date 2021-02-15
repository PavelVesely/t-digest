#!/bin/sh

# runs all experiments mentioned in README (scenarios 1., 2. and 3.) to reproduce the results in the paper 

export CLASSPATH="$CLASSPATH:./target/classes:./target/test-classes:.target/classes/org/apache/datasketches/req/ReqSketch:../../datasketches/target/classes:../../../.m2/repository/org/apache/datasketches/datasketches-java/1.3.0-incubating/datasketches-java-1.3.0-incubating.jar"

# scenario 1.
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.CarefulAttack resources/CarefulAttack_k_0_merging.conf &
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.CarefulAttack resources/CarefulAttack_k_0_clustering.conf &

# scenario 2.
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.IIDgenerator resources/IIDgenerator.conf &
sed 's/Distribution=loguniform2/Distribution=loguniform/' resources/IIDgenerator.conf >/tmp/IIDgenerator-loguniform.conf
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.IIDgenerator /tmp/IIDgenerator-loguniform.conf &
sed 's/MaxExp=0/MaxExp=10/' resources/IIDgenerator.conf >/tmp/IIDgenerator-MaxExp.conf
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.IIDgenerator /tmp/IIDgenerator-MaxExp.conf &

# scenario 3.
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.SpeedComparison resources/SpeedComparison.conf &