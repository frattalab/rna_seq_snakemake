import re, string, pandas as pd, numpy as np, yaml


from librarySize import *

project_folder = "/SAN/vyplab/alb_projects/data/4su_full_ward_tdp_kd_ipsc/"
sampleCSVpath = project_folder + "ward_time_course_samples.csv"
star_output_folder = project_folder + "STAR_aligned/"
feature_counts_output_folder = project_folder + "feature_counts/"

mappingStats = project_folder + "mappingStats.csv"

SAMPLES = pd.read_csv(sampleCSVpath, sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

allStarLogs = [star_output_folder + name + ".Log.final.out" \
               for name in SAMPLE_NAMES]

allFeatureLogs = [feature_counts_output_folder + name + "_featureCounts_results.txt.summary" \
                  for name in SAMPLE_NAMES]


libs = list(map(getLibSize, allStarLogs))
starMapped = list(map(getStarMapped, allStarLogs))
featureMapped = list(map(getFeatureCountsMapped, allFeatureLogs))


stats = {'sample': SAMPLE_NAMES, \
'librarySize': libs,\
'percentMapped':starMapped,\
'featureCountsMappedToGene': featureMapped}
df = pd.DataFrame(stats)
print(df)

df.to_csv(mappingStats, sep=',')
