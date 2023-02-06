
import os
from glob import glob
from slugify import slugify


_GENOME_LENGTH={}
def get_genome_length(genome):
    if len(_GENOME_LENGTH.keys()) != 0:
        return _GENOME_LENGTH[genome]
    subgroups=glob("../data/genome-list/*")
    for subgroup in subgroups:
        genome_list=open(subgroup).read().split("\n")
        genome_list=[a.split("\t") for a in genome_list if len(a.strip()) != 0]
        for a in genome_list:
            _GENOME_LENGTH[f"{os.path.basename(a[4])}-{slugify(a[0])}"] = int(float(a[1].replace(',', '.'))*(10**6))
    return _GENOME_LENGTH[genome]
sex_chromosomes="XYZW"
if __name__ == '__main__':
    microsatellite_list=open("../data/ssr-stats/exhaustive-stats/exhaustive-microsatellite-list.txt").read().split("\n")

    genome_list=glob("../data/ssr-stats/ssr-quantity/chromosome-wide/*")
    genome_list=[os.path.basename(a) for a in genome_list]
    for genome  in genome_list:
        chromosome_list=glob(f"../data/ssr-stats/ssr-quantity/chromosome-wide/{genome}/*")
        chromosome_list=[os.path.basename(a).split(".")[0] for a in chromosome_list]

        genome_exhaustive_loci_length_list = {a: 0 for a in microsatellite_list}
        chromosome_exhaustive_loci_length_list = {a: {} for a in chromosome_list}
        genome_length=get_genome_length(genome)
        chromosome_length={a:0 for a in chromosome_list}

        for chromosome in chromosome_list:
            chromosome_exhaustive_loci_length_list[chromosome]={a:0 for a in microsatellite_list}
            ssr_loci_list=open(f"../data/ssr-stats/ssr-quantity/chromosome-wide/{genome}/{chromosome}.tsv").read().split("\n")
            ssr_loci_list=[a.split("\t") for a in ssr_loci_list]
            ssr_loci_list=[a for a in ssr_loci_list if a[0] in microsatellite_list]
            for a in ssr_loci_list:
                chromosome_exhaustive_loci_length_list[chromosome][a[0]]=int(a[2])
                genome_exhaustive_loci_length_list[a[0]]+=int(a[2])

        # IMPORTANT, Chromosome length has been faked to to lack of data
        _chromosome_ssr_total_length={a:sum([chromosome_exhaustive_loci_length_list[a][b] for b in list(chromosome_exhaustive_loci_length_list[a])]) for a in chromosome_list}
        _genome_ssr_total_length=sum([_chromosome_ssr_total_length[a] for a in _chromosome_ssr_total_length.keys()])
        _chromosome_estimated_length={a: int(_chromosome_ssr_total_length[a] * genome_length / _genome_ssr_total_length) for a in chromosome_list}

        for sex_chromosome in set(chromosome_list) & set(sex_chromosomes):
            # For each microsatellites, calculate :
            # 1) chromosome-wide density
            # 2) genome-wide density
            # 3) fold difference
            report=[["Genome", "Chromosome", "Microsatellite", "Chromosome-wide density", "Genome-wide density", "Fold diff."]]
            for microsatellite in microsatellite_list:
                chromosome_wide_density=chromosome_exhaustive_loci_length_list[sex_chromosome][microsatellite]/_chromosome_estimated_length[sex_chromosome]
                genome_wide_density=genome_exhaustive_loci_length_list[microsatellite]/genome_length
                fold_diff=1 if genome_wide_density == 0 else chromosome_wide_density/genome_wide_density
                report.append([genome, sex_chromosome, microsatellite, chromosome_wide_density, genome_wide_density, fold_diff])
            with open(f"../data/ssr-stats/exhaustive-stats/highlight-sex-chromosomes/{genome}-{sex_chromosome}.tsv", "w") as f:
                f.write("\n".join(["\t".join([str(b) for b in a]) for a in report]))
            print(f"Write {genome}-{sex_chromosome}.tsv")