#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
import matplotlib.pyplot as plt


# get_ipython().run_line_magic('matplotlib', 'inline')


def get_readable_name(cname):
    if '_' in cname:
        values = cname.split('_')[1:]
        return 'Jaccard : ' + values[0] + ' VI : ' + values[1]
    return cname


plt.rcParams['figure.figsize'] = [12, 7]

parser = argparse.ArgumentParser()

parser.add_argument('-i', help='input file path (out.csv from go)', required=True)
parser.add_argument('-metrics', help="Comma separated metrics(Column names) to compare  |"
                                     " Accepted column names: 'Jaccord','jacc_0.95_0.05', 'jacc_0.9_0.1', 'jacc_0.85_0.15',"
                                     "'jacc_0.15_0.85', 'jacc_0.5_0.5','Overlap','Dice','CentroidDistance'", required=True)
# parser.add_argument('-o', help='output file that contains the prediction of test file', required=True)
# parser.add_argument('-a', help='alpha value used for calculating the probability', required=False, default='1.0')

args = vars(parser.parse_args())

df = pd.read_csv(args["i"])
metrics = args["metrics"].split(",")

# In[33]:


boundary_score = [0.95, 0.9, 0.85, 0.15, 0.5]
interval_score = [0.05, 0.1, 0.15, 0.85, 0.5]

# In[34]:


df.fillna(0, inplace=True)

# In[42]:


for i in range(len(boundary_score)):
    df['jacc_' + str(boundary_score[i]) + '_' + str(interval_score[i])] = boundary_score[i] * df['Jaccord'] + \
                                                                          interval_score[i] * df['VI']

print("Correlation with VI")
print(df[metrics].corrwith(df['VI']))

# In[63]:


ax = df.plot(y=metrics)
ax.legend(list(map(lambda f: get_readable_name(f), metrics)))
ax.set_ylabel('Similarity')
plt.show()

# In[ ]:
