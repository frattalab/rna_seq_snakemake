
import json
import re
import getopt, sys, os.path
import os
import pandas as pd


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["help", "input=","output="])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    infile = None
    outfile = None

    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            if os.path.isfile(arg):
                infile = arg
        elif opt in ("-o", "--output"):
            outfile = arg
        else:
            print(opt)
            assert False, "Unhandled option"

    if infile is not None and outfile is not None:
        parsed_df = parse_fastp(infile)
        parsed_df.to_csv(outfile)
    else:
        usage()


def usage():
    print("\nUsage: python parse_fastp_json [options] <mandatory>")
    print("Used for turning fastp's json output into a table of sample insert metrics")
    print("Options:")
    print("\t-h, --help:\n\t\t show this help message and exit")
    print("Mandatory:")
    print("\t-i, --input:\n\t\t File with the regions in bed format")
    print("\t-o, --output:\n\t\t Name of the gtf file output file. Directory where the file will be created should exist!")

def parse_fastp(fastp_file):
    sample_name = re.sub('_fastp.json', '',fastp_file)
    sample_name = os.path.basename(sample_name).split('_', 1)[1]
    #load in the data as a json
    with open(fastp_file) as f:
        data = json.load(f)
    #take the summary stats
    temp = pd.DataFrame.from_dict(data['summary'])
    #remove the index
    temp.reset_index(level=0, inplace=True)
    #get the peak of the insert size
    insert_peak = data['insert_size']['peak']
    #build new rows, and then insert them
    newrow = {"index" : "insert_size_peak","before_filtering": insert_peak,"after_filtering" : insert_peak}
    newrow2 = {"index" : "sample_name","before_filtering": sample_name,"after_filtering" : sample_name}
    temp = temp.append(newrow,ignore_index=True)
    temp = temp.append(newrow2,ignore_index=True)
    #transpose the df and return it
    temp2 = temp.transpose()
    temp2.columns = temp2.iloc[0]
    temp2 = temp2.drop('index')

    return(temp2)




if __name__ == "__main__":
    main()

