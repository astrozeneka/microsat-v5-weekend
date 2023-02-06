
from PyPDF2 import PdfMerger
import pandas as pd
import os
from glob import glob
from slugify import slugify
import numpy as np
import matplotlib.pyplot as plt
pastel_colors = ['#FD8A8A', '#A8D1D1', '#9EA1D4', '#ADA2FF', '#FFCAC8', "#DBA39A"] # One colors for each microsatellites


sex_chromosomes="XYZW"
def get_genome_list(subgroup):
    output = open(f"../data/genome-list/{subgroup}.tsv").read().split("\n")
    output = [a.split("\t") for a in output if len(a) > 0]
    output = [f"{os.path.basename(a[4])}-{slugify(a[0])}" for a in output]
    return output

if __name__ == '__main__':

    filename_index=1
    microsatellite_list=open("../data/ssr-stats/exhaustive-stats/exhaustive-microsatellite-list.txt").read().split("\n")

    for sex_chromosome in sex_chromosomes:
        genome_list = glob(f"../data/ssr-stats/exhaustive-stats/highlight-sex-chromosomes/*-{sex_chromosome}.tsv")
        genome_list = [os.path.basename(a).split('.')[0] for a in genome_list]
        genome_list  = [a[:-2] for a in genome_list]

        subgroups = glob("../data/genome-list/*")
        subgroups = [os.path.basename(a).split(".")[0] for a in subgroups]
        subgroup_genomes = {}
        for subgroup in subgroups:
            _list = open(f"../data/genome-list/{subgroup}.tsv").read().split("\n")
            _list = [a.split("\t") for a in _list if len(a.strip()) != 0]
            _list = [f"{os.path.basename(a[4])}-{slugify(a[0])}" for a in _list]
            subgroup_genomes[subgroup] = list(set(_list) & set(genome_list))

        fold_diff_exhaustive_list = {a: {b: None for b in genome_list} for a in microsatellite_list}
        for genome in genome_list:
            ms_fold_diff = open(f"../data/ssr-stats/exhaustive-stats/highlight-sex-chromosomes/{genome}-{sex_chromosome}.tsv").read().split("\n")[1:]
            ms_fold_diff = [a.split("\t") for a in ms_fold_diff]
            ms_fold_diff = [a for a in ms_fold_diff]
            for row in ms_fold_diff:
                fold_diff_exhaustive_list[row[2]][row[0]] = float(row[5])

        data = open("../data/ssr-stats/exhaustive-stats/p-value.tsv").read().split("\n")
        data = [a.split("\t") for a in data if a.strip() != ""][1:]
        data = [a for a in data if float(a[-1]) < 0.05]  # P-value < alpha
        pv_table = [[a[0], a[1], a[2], a[-1]] for a in data if a[0] == sex_chromosome]

        for entry in pv_table:
            microsatellite=entry[1]
            subgroup=entry[2]
            h0_subgroups=[a for a in subgroups if a != subgroup]
            ha_subgroups=[a for a in subgroups if a in subgroup]
            h0_genome_list=sum([get_genome_list(a) for a in h0_subgroups], [])
            ha_genome_list=sum([get_genome_list(a) for a in ha_subgroups], [])
            h0_genome_list=[a for a in h0_genome_list if a in genome_list]
            ha_genome_list=[a for a in ha_genome_list if a in genome_list]

            h0_fold_diff_list={a: fold_diff_exhaustive_list[microsatellite][a] for a in h0_genome_list}
            ha_fold_diff_list={a: fold_diff_exhaustive_list[microsatellite][a] for a in ha_genome_list}
            h0_fold_diff_list = np.array([h0_fold_diff_list[a] for a in h0_fold_diff_list.keys()])
            ha_fold_diff_list = np.array([ha_fold_diff_list[a] for a in ha_fold_diff_list.keys()])

            plt.figure()
            x = np.array( + ha_fold_diff_list.shape[0] * [1])
            y = np.array(list(h0_fold_diff_list) + list(ha_fold_diff_list))
            plt.xticks([-1, 0, 1, 2], ["", "Other groups", subgroup, ""])
            plt.scatter(h0_fold_diff_list.shape[0] * [0]
                        , h0_fold_diff_list, s=3, color=pastel_colors[0])
            plt.scatter(ha_fold_diff_list.shape[0] * [1]
                        , ha_fold_diff_list, s=3, color=pastel_colors[1])
            plt.scatter([-1], [1], s=0)
            plt.scatter([2], [1], s=0)
            plt.title(f"({microsatellite})N density fold diff. on {sex_chromosome} chromosome of {subgroup} \n p={round(float(entry[-1]), 4)}")
            plt.savefig(f"../docs/plots-highlight-sex-chromosome/{str(filename_index).zfill(3)}-{sex_chromosome}-{subgroup}-{microsatellite}.pdf")
            print(f"{str(filename_index).zfill(3)}-{sex_chromosome}-{subgroup}-{microsatellite}.pdf")
            filename_index+=1
            print()