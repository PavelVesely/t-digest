#!/bin/sh

# runs all experiments mentioned in README (scenarios 1., 2., 3., and 4.) to reproduce the results in the paper 

export CLASSPATH="$CLASSPATH:./target/classes:./target/test-classes"

# scenario 1.
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.CarefulAttack resources/CarefulAttack_k_0_merging.conf &
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.CarefulAttack resources/CarefulAttack_k_0_clustering.conf &

# scenario 2.
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.IIDgenerator resources/IIDgenerator.conf &
sed 's/Distribution=loguniform2/Distribution=loguniform/' resources/IIDgenerator.conf >/tmp/IIDgenerator-loguniform.conf
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.IIDgenerator /tmp/IIDgenerator-loguniform.conf &

# scenario 3. (varying maxExp scenario)
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.LoguniformWithVaryingMaxExpGenerator resources/LoguniformWithVaryingMaxExpGenerator.conf &

# scenario 4.
java -ea -Dfile.encoding=UTF-8 com.tdunning.math.stats.SpeedComparison resources/SpeedComparison.conf &
