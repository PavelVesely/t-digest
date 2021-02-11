#!/usr/bin/python3
# Licensed to the Apache Software Foundation (ASF) under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The ASF licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Authors: Joseph Ross and Pavel Vesely

import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser(description='simple plot generator for update time')
parser.add_argument('-f', help='csv file to process')
args = parser.parse_args()
filepath = args.f
with open(filepath) as f:
    df = pd.read_csv(f, header=0, sep=';')
x_axis = df.columns[0]
#for col in df.columns[1:-1]:
#    plt.plot(df[x_axis], df[col], label=col)
plt.plot(df[x_axis], df[df.columns[1]], linestyle='solid', label=df.columns[1])
plt.plot(df[x_axis], df[df.columns[2]], linestyle='dashed', label=df.columns[2])
plt.xlabel('k')
plt.ylabel('update time [ns]')
plt.legend()
plt.savefig('plot' + os.path.basename(filepath) + '.png', dpi=300)
