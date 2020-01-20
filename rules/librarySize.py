import re, string

def getLibSize(STARLogFile):
    # read in the logFile from star and split the lines
    lines = open(STARLogFile).read().splitlines()
    # filter for the unique mappers
    uni_map = list(filter(lambda x:'Uniquely mapped reads number' in x, lines))[0]
    # filter for the multi mappers
    multi_map = list(filter(lambda x:'Number of reads mapped to multiple loci' in x, lines))[0]
    # parse them down and turn them into numbers
    num_unique = int(uni_map.split('\t')[1])
    num_multi = int(multi_map.split('\t')[1])
    libSize = num_unique + num_multi
    return(libSize)

def getStarMapped(STARLogFile):
    # read in the logFile from star and split the lines
    lines = open(STARLogFile).read().splitlines()
    # filter for the unique mappers
    map_percent = list(filter(lambda x:'Uniquely mapped reads %' in x, lines))[0]
    map_percent_string = float(re.sub('%', '', map_percent.split('\t')[1]))
    return(map_percent_string)

def getFeatureCountsMapped(featureCountsLogFile):
    lines = open(featureCountsLogFile).read().splitlines()
    temp = [x.split('\t') for x in lines]
    info = { k[0]: k[1] for k in temp }
    # bet there's a better way of doing this but I'm not playing code golf today
    # Loop through all the keys, add everything that isn't 'Status' or 'Assigned'
    # to the denominator
    for key, value in info.items():
        unassigned = 0
        if key == 'Assigned':
            assigned = int(value)
        elif key == 'Status':
            continue
        else:
            unassigned += int(value)
    # return the fraction of reads assigned to genes
    return(assigned / unassigned)
