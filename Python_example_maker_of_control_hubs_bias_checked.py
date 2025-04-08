import random

    """
    This script is used to create a Control or "dummy" file that's based on an experimental file that contains key positional information 
    for "Adaptive hubs" which is a name we're giving to chromatin loops that contain genes that are key to adaptation (Ie: contain signatures of 
    strong selection, as calculated by CLR and Tajima's D measures, contain differentially expressed genes, between hindwing tissue with red color patterns and hindwing tissue with no red color pattern, and which are bound by Optix TFs, known
    to control red color patterns that are well evidenced to be under strong selection from field studies)
    In Particular this generates control hubs that are based on an input file (of adaptive hubs) but whose positions will be randomised. 
    It also takes a chromosome sizes file so that it knows the real genomic bounderies within which to provide the random positions for control hubs.
    In downstream analysess this control file will be used to analyse differences in interchromosomal Linkage Disequilibrium between adaptive hubs and between control hubs.
    """

def read_chromosome_sizes(file_path):
    """
    This'll read a file containing chromosome sizes and returns a dictionary.
    The file format must be be: chromosome_name\tsize 
    i.e: Herato01 19290100
    """
    chromosome_sizes = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            chromosome_sizes[parts[0]] = int(parts[1])
    return chromosome_sizes

def read_experimental_data(file_path):
    """
    Reads the experimental data file and returns detailed range information.
    The file format should be: chromosome\tstart\tend\tother_columns.. 
    i.e: Herato01 10000   23000   bla bla bla
    """
    detailed_ranges = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            detailed_ranges.append((chrom, start, end))
    return detailed_ranges

def generate_control_data(chromosome_sizes, experimental_data):
    """
    Generates control data and "Bias Checks", ensuring each randomly selected chromosome 
    has the same number of representations as a corresponding chromosome in the experimental data.
    A previous version of this script did not ensure that if a given chromosome in the experimental file 
    was represented a given number of times, that the output control file gave a random chromosome with 
    the same number of representations. This could drive a bias in LD values being higher in the file which
    had more repetitions of the same chromosome for obvious LD reasons. 
    Therefore, if experimental file shows that chromosome 8 has 3 instances, then a random chromosome (ie Chr21) 
    will be given 3 instances as well. 
    """
    control_data = []
    chromosomes = list(chromosome_sizes.keys())

    # Count the number of ranges for each chromosome in the experimental data
    ranges_per_chromosome = {}
    for chrom, start, end in experimental_data:
        if chrom not in ranges_per_chromosome:
            ranges_per_chromosome[chrom] = []
        ranges_per_chromosome[chrom].append(end - start)

    # Assign a random chromosome to each chromosome in the experimental data
    random_chromosomes = random.sample(chromosomes, len(ranges_per_chromosome))

    # Generate control data
    for (exp_chrom, range_sizes), rand_chrom in zip(ranges_per_chromosome.items(), random_chromosomes):
        for size in range_sizes:
            max_start_pos = chromosome_sizes[rand_chrom] - size
            if max_start_pos > 0:
                random_start = random.randint(0, max_start_pos)
                random_end = random_start + size
                control_data.append((rand_chrom, random_start, random_end, "control_dataset"))

    return control_data

def write_control_data(file_path, control_data):
    """
    Writes the control data to a file, sorted by chromosome and start position.
    """
    control_data.sort(key=lambda x: (x[0], x[1], x[2]))
    with open(file_path, 'w') as file:
        for entry in control_data:
            file.write('\t'.join(map(str, entry)) + '\n')

# Example usage
chromosome_sizes_file = "/home/rpapa/jogilvie/work/amalfreda/LD_SVB/making_control_loci_for_LD/chrom.sizes"
experimental_data_file = "/home/rpapa/jogilvie/work/easley/LD_easley/making_control_loci_for_LD/adaptive_hubs_all.amalfreda_edited.bed"
control_data_file = "/home/rpapa/jogilvie/work/easley/LD_easley/making_control_loci_for_LD/amalfreda_control_hubs.bed"

chromosome_sizes = read_chromosome_sizes(chromosome_sizes_file)
range_sizes = read_experimental_data(experimental_data_file)
control_data = generate_control_data(chromosome_sizes, range_sizes)
write_control_data(control_data_file, control_data)

