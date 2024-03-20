#TO DO add efficient way to store header names, so last rule is not a bottleneck
from snakemake.io import directory, expand, temp
from pathlib import Path
import logging

# logging.basicConfig(level=logging.ERROR)
ROOT = Path(config['output'])
LOGS = ROOT / 'logs'
BENCH = ROOT / 'bench'
PYTHON_SCRIPTS = config['python_scripts']
CLASSIFIERS_LOCATIONS = config['locations']
TAXONOMY_DB = CLASSIFIERS_LOCATIONS['taxonomy']['db_path']

CLASS_TYPE = ['reads']
CLASS_LEVEL = ['genus', 'species']
# reverse bool
filter_reads = not config['no_filter']
USE_CLASSIFIERS = config['classifiers']
truth_file = config['truth_file']

ENV = Path('/scratch/alvanuffelen/env/')


for classifier in USE_CLASSIFIERS:
    if classifier in ['ccmetagen', 'bracken']:
        # It is already checked in run_map_classifiers.py if kma and kraken2 are present for ccmetagen and bracken, respectively.
        continue
    include: f'{classifier}.smk'

rule all:
    input:
        # expand(ROOT / 'output' / '{class_type}' / '{class_level}' / 'output_{classifier}.tsv', class_type=CLASS_TYPE, class_level=CLASS_LEVEL, classifier=use_classifiers),
        # expand(ROOT / 'output' / '{class_type}' / '{class_level}' / 'output_kaiju.tsv', class_type=CLASS_TYPE, class_level=CLASS_LEVEL),
        # ROOT/ 'logs' / 'merge_check_output_tax.log',
        expand(Path(truth_file).parent / '{ground_truth_type}_updated.yml', ground_truth_type=['ground_truth', 'ground_truth_tax'])

rule filter_quality_length:
    """
    Filters reads based on quality and length
    """
    input:
        config['input']['fastq']
    output:
        ROOT / 'input' / (Path(config['input']['fastq']).stem + '_HQ' + Path(config['input']['fastq']).suffix)
    params:
        length = 1000,
        quality = 7
    shell:
        """
        ml seqkit/2.3.1; \
        seqkit seq \
        -m {params.length} \
        -Q {params.quality} \
        {input} > {output}
        """

rule stats_filtered_reads:
    """
    Get stats on the filtered reads
    """
    input:
        rules.filter_quality_length.output  if filter_reads else config['input']['fastq']
    output:
        TSV = ROOT / 'input' / 'HQ.stats'
    shell:
        """
        ml seqkit/2.3.1; \
        seqkit stats -aT {input} > {output}
        """


rule canu_assembly:
    """
    Assembles reads with Canu
    """
    input:
        rules.filter_quality_length.output if filter_reads else config['input']['fastq']
    output:
        ROOT / 'canu' / 'hackaton.contigs.fasta'
    params:
        dir_out = directory(ROOT / 'canu'),
        prefix = 'hackaton'
    threads: 24
    shell:
        """
        ml load canu/2.2;
        canu -p {params.prefix} -d {params.dir_out} genomeSize=1m stopOnLowCoverage=0.000001 -nanopore-raw {input} -maxThreads={threads}
        """

FILENAME= dict(assembly = rules.canu_assembly.output, reads = rules.canu_assembly.input)

rule convert_seq_to_tax_abundace:
    input:
        truth_file
    output:
        Path(truth_file).parent / 'ground_truth_tax.yml'
    log:
        ROOT / "logs" / "ground_truth_tax.log"
    script:
        "../scripts/convert_seq_to_tax_abun.py"


rule update_ground_truth:
    """
    The taxonomy of the output from the classifiers will be determined by means of the provided taxdump.
    However, due to name changes, this can differ from the taxonomy of the ground truth.
    This rule checks if there are any older names in the ground truth.
    It renames the taxonomy of given ground truth file and creates an updated version, if necessary
    """
    input:
        lambda wc: {'ground_truth': truth_file, 'ground_truth_tax': rules.convert_seq_to_tax_abundace.output}[wc.ground_truth_type]
    output:
        Path(truth_file).parent / '{ground_truth_type}_updated.yml'
    params:
        TAXONOMY = TAXONOMY_DB
    log:
        ROOT / "logs" / "update_{ground_truth_type}.log"
    shell:
        """
        ml taxonkit/0.15.0
        export TAXONKIT_DB={params.TAXONOMY}
        # Create table with old_name, found taxid and new_name
        NEW_NAMES=$(cut -f 1 -d ':' {input} | taxonkit reformat --format "{{s}}" -I 1)
        # Print is a name will be updated
        echo "taxid\tnew_name" > {log}
        echo "$NEW_NAMES" >> {log}
        paste -d ':' <(echo "$NEW_NAMES" | cut -f 2 | sed -e "s/\(.*\)/'\\1'/") \
                     <(cut -f 2 -d ':' {input}) \
                     > {output}
        """

rule check_output_tax:
    """
   The taxonomy of the output from the classifiers will be determined by means of the provided taxdump.
   However, due to name changes, this can differ from the taxonomy of the ground truth.
   This rule checks if there are any older names in the output.
   It renames the taxonomy of the given output file and creates an updated version, if necessary
   """
    input:
        ROOT / 'output' / '{class_type}' / '{class_level}' / 'output_{classifier}.tsv'
    output:
        temp(ROOT / "logs" / "update_output_{class_type}_{classifier}_{class_level}.log")
    params:
        TAXONOMY = TAXONOMY_DB
    shell:
        """
        list_uc="no_hit|ambigious|unclassified|missing_genus|not_distributed|rounding_error"
        taxa=$(grep -vE $list_uc {input})
        uc=$(grep -E $list_uc {input})
        ml taxonkit/0.15.0
        export TAXONKIT_DB={params.TAXONOMY}
        # Create table with old_name, found_taxid and new_name
        if [[ {wildcards.class_level} == "genus" ]]
        then
            NEW_NAMES=$(echo "$taxa" | cut -f 1 | taxonkit name2taxid -s | taxonkit reformat --format "{{g}}" -I 2)
        elif [[ {wildcards.class_level} == "species" ]]
        then
            NEW_NAMES=$(echo "$taxa" | cut -f 1 | taxonkit name2taxid -s | taxonkit reformat --format "{{s}}" -I 2)
        fi
        # Print if a name will be updated
        echo "level\tclassifier\ttaxid\told_name\tnew_name" > {output}
        echo "$NEW_NAMES" | awk -F "\t" '{{if ($1!=$3) print "{wildcards.class_level}\t{wildcards.classifier}\t"$1"\t"$2"\t"$3}}' >> {output}
        #paste <(echo "$NEW_NAMES" | cut -f 3) \
        #      <(echo "$taxa" | cut -f 2) \
        #      > {output}
        #echo "$uc" >> {output}
        """

rule merge_check_output_tax:
    input:
        expand(rules.check_output_tax.output, class_type=CLASS_TYPE, class_level=CLASS_LEVEL, classifier=use_classifiers)
    output:
        ROOT / 'logs' / 'merge_check_output_tax.log'
    shell:
        """
        echo "level\tclassifier\ttaxid\told_name\tnew_name" > {output}
        for file in {input}
        do
        tail -n+2 $file >> {output}
        done
        """

