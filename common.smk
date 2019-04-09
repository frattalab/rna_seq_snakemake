#this rule is for running fastqc across multiple fastqc, it simply collects all of the fastqc in the sample sheet

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])
def get.all.fastqfiles(DATA):
	#