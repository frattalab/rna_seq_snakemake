import re, string, pandas as pd, numpy as np
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"
include: "librarySize.py"
SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

allStarLogs = expand(config['star_output_folder'] + "{name}.Log.final.out", name = SAMPLE_NAMES)
allFeatureLogs = expand(config['feature_counts_output_folder'] + "{name}_featureCounts_results.txt.summary", name = SAMPLE_NAMES)



starMapped = list(map(getStarMapped, allStarLogs))
featureMapped = list(map(getFeatureCountsMapped, allFeatureLogs))
stats = {'sample': SAMPLE_NAMES, \
'librarySize': libs,\
'percentMapped':lstarMapped,\
'featureCountsMappedToGene': featureMapped}
df = pd.DataFrame(stats)
print(df)
