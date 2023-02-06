
import os
from glob import glob
from slugify import slugify
import numpy as np
from scipy.stats import ttest_ind

GROUP_MIN = 8

def get_genome_list(subgroup):
    output = open(f"../data/genome-list/{subgroup}.tsv").read().split("\n")
    output = [a.split("\t") for a in output if len(a) > 0]
    output = [f"{os.path.basename(a[4])}-{slugify(a[0])}" for a in output]
    return output

sex_chromosomes="XYZW"
if __name__ == '__main__':
    microsatellite_list=open("../data/ssr-stats/exhaustive-stats/exhaustive-microsatellite-list.txt").read().split("\n")
    subgroups=glob("../data/genome-list/*")
    subgroups=[os.path.basename(a).split(".")[0] for a in subgroups]
    report = [["Sex Chromosome", "Microsatellite", "Subgroup", "In-Group genome analyzed", "In-Group fold diff. Mean", "In-Group fold diff. STD", "Other group fold diff. Mean", "Other group fold diff. STD", "T-Score", "P-value"]]
    for sex_chromosome in sex_chromosomes:

        genome_chromosome_list = glob("../data/ssr-stats/exhaustive-stats/highlight-sex-chromosomes/*.tsv")
        genome_chromosome_list = [os.path.basename(a).split(".")[0] for a in genome_chromosome_list]
        genome_chromosome_list = [(a[:-2], a[-1]) for a in genome_chromosome_list]
        genome_chromosome_list = [a[0] for a in genome_chromosome_list if a[1] == sex_chromosome]

        fold_diff_exhaustive_list = {a: {b: None for b in genome_chromosome_list} for a in microsatellite_list}

        for genome_chromosome in genome_chromosome_list:
            ms_fold_diff = open(f"../data/ssr-stats/exhaustive-stats/highlight-sex-chromosomes/{genome_chromosome}-{sex_chromosome}.tsv").read().split("\n")[1:]
            ms_fold_diff = [a.split("\t") for a in ms_fold_diff]
            ms_fold_diff = [a for a in ms_fold_diff]
            for row in ms_fold_diff:
                fold_diff_exhaustive_list[row[2]][row[0]] = float(row[5])

        # Check if ONE microsatellite has a tendency have fold diff. in one subgroup
        for subgroup in subgroups:
            h0_subgroups=[a for a in subgroups if a != subgroup]
            ha_subgroups=[a for a in subgroups if a in subgroup]
            h0_genome_list=sum([get_genome_list(a) for a in h0_subgroups], [])
            ha_genome_list=sum([get_genome_list(a) for a in ha_subgroups], [])

            for microsatellite in microsatellite_list:
                print(f"Analyzing ({microsatellite})n")
                h0_fold_diff_list={a: fold_diff_exhaustive_list[microsatellite][a] for a in h0_genome_list if a in fold_diff_exhaustive_list[microsatellite].keys()}
                ha_fold_diff_list={a: fold_diff_exhaustive_list[microsatellite][a] for a in ha_genome_list if a in fold_diff_exhaustive_list[microsatellite].keys()}
                h0_fold_diff_list=np.array([h0_fold_diff_list[a] for a in h0_fold_diff_list.keys()])
                ha_fold_diff_list=np.array([ha_fold_diff_list[a] for a in ha_fold_diff_list.keys()])

                stat, p=ttest_ind(h0_fold_diff_list, ha_fold_diff_list)
                # High T-Stats: more difference between groups
                # Low T-Stats: more similarity between groups

                row = [
                    sex_chromosome,
                    microsatellite,
                    subgroup,
                    len(ha_fold_diff_list),
                    float(ha_fold_diff_list.mean()),
                    float(ha_fold_diff_list.std()),
                    float(h0_fold_diff_list.mean()),
                    float(h0_fold_diff_list.std()),
                    float(stat),
                    float(p)
                ]
                # p-value listing should be done, in group length should be analyzed also
                if len(ha_fold_diff_list) >= GROUP_MIN:
                    if(row[4] > row[6]) and (row[4] > 1): # A high fold diff is constated in the studied group
                        report.append(row)

    head=report[0]
    report=report[1:]
    report.sort(key=lambda x:x[-2], reverse=False)
    report=[head]+report
    with open("../data/ssr-stats/exhaustive-stats/p-value.tsv", "w") as f:
        f.write("\n".join(["\t".join([str(b) for b in a]) for a in report]))
    print("Done")