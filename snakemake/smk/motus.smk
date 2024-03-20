MOTUS_ENV= CLASSIFIERS_LOCATIONS['motus']['path_env']
MOTUS_DB = CLASSIFIERS_LOCATIONS['motus']['db_path']
BWA = CLASSIFIERS_LOCATIONS['bwa']['path']
SAMTOOLS = CLASSIFIERS_LOCATIONS['samtools']['path']


rule motus_split_long_reads:
    """
    Splits the long reads into shorted reads
    """
    input:
        lambda wildcards: FILENAME[wildcards.class_type]
    output:
        converted_long_reads = ROOT / 'motus' / 'converted_long_reads.fasta.gz'
    params:
        db = MOTUS_DB,
    shell:
        """
        {MOTUS} prep_long -i {input} -o {output.converted_long_reads} -db {params.db}
        """

rule motus_classification:
    """
    Classifies reads/assemblies using motus.
    """
    input:
        converted_long_reads = rules.motus_split_long_reads.output.converted_long_reads
    output:
        tax_profile = ROOT / 'motus' / 'taxonomic_profile.txt'
    log:
        classification_log = ROOT / 'motus' / 'log.txt'
    params:
        db = MOTUS_DB,
    threads: 32
    shell:
        """
        export PATH=""$(dirname "{BWA}")":$PATH"
        export PATH=""$(dirname "{BOWTIE2}")":$PATH"
		. {MOTUS_ENV}/bin/activate
        motus profile -s {input.converted_long_reads} \
                      -t {threads} \
                      -p  \
                      -u \
                      -e \
                      -q \
                      -db {params.db} \
                         1> {output.tax_profile} \
                         2> {log}
        """

rule motus_clean:
    """
    Cleans the output of motus.
    Remove ref-mOTU if: Taxonomic lineage is 'NA'
    """
    input:
        tax_profile = rules.motus_classification.output.tax_profile
    output:
        cleaned_count = ROOT / 'motus' / 'cleaned_count.tsv',
        cleaned_tax_profile = ROOT / 'motus' / 'cleaned_taxonomic_profile.txt'

    run:
        import pandas as pd

        found_otus, found_otus_debug = [], []
        with open(input.tax_profile, 'r') as f:
            for line in f.read().splitlines():
                if line[0] != '#':
                    otu, lineage, tax_lineages, count = line.split('\t')
                    if otu == "unassigned": # keep unassigned count
                        # do not add unknown, they will be added in the last rule
                        #found_otus.append(['unassigned', float(count)])
                        #found_otus_debug.append([otu, tax_lineages, lca, count])
                        continue
                    if float(count) > 0: # skip if count is 0
                        lca = ''
                        for tax_lineage in tax_lineages.split('|')[::-1]: # backtrack lca until one if not NA
                            if tax_lineage != 'NA':
                                lca = tax_lineage
                                break
                            else:
                                continue
                        if lca == '': # check if an lca was found
                            print(f'Dropped {otu}; unknown taxid')
                            lca = 'NA'
                            found_otus_debug.append([otu, tax_lineages, lca, count])
                        found_otus.append([lca, float(count)])
                        found_otus_debug.append([otu, tax_lineages, lca, count])

        df = pd.DataFrame(found_otus, columns=['lca_tax', 'count']).groupby('lca_tax').sum()
        df_debug = pd.DataFrame(found_otus_debug, columns=['otu', 'lca_tax', 'lca', 'count'])

        df.to_csv(output.cleaned_count, sep='\t', header=False)
        df_debug.to_csv(output.cleaned_tax_profile, sep='\t', index=False)


rule motus_taxonomy:
    input:
        cleaned_count = rules.motus_clean.output.cleaned_count
    output:
        cleaned_count_tax = ROOT / 'motus' / 'cleaned_count_tax.tsv'
    params:
        db = TAXONOMY_DB
    shell:
        """
        ml taxonkit;
        paste <(cut -f 1 {input.cleaned_count} \
               | taxonkit lineage --data-dir {params.db} \
               | taxonkit reformat --data-dir {params.db} --format "{{g}}\t{{s}}" --miss-rank-repl "unclassified" \
               | cut -f 3,4 \
               | awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($1 == "unclassified" && $2 != "unclassified") {{print "missing_genus", $2}} else {{print}}}}' \
               ) \
              <(cut -f 2 {input.cleaned_count}) \
        > {output.cleaned_count_tax};
        """

rule motus_output:
    input:
        cleaned_count_tax = rules.motus_taxonomy.output.cleaned_count_tax,
        tax_profile = rules.motus_classification.output.tax_profile
    output:
        TSV_motus = ROOT / 'output' / '{class_level}' / 'output_motus.tsv'
    params:
        level = lambda wildcards: wildcards.class_level,
        duplicate_names= lambda wildcards: ROOT / 'output' / wildcards.class_level / 'corrected_duplicates_motus.log'
    run:
        import pandas as pd

        with open(input.tax_profile, 'r') as f:
            for line in f:
                if line.split('\t')[0] == 'unassigned':
                    unassigned = float(line.split('\t')[-1])
        motus_input = pd.read_table(input.cleaned_count_tax,
                        names=['genus', 'species', 'count']
                    )

        organism_count = pd.Series(motus_input['count'].values, index=motus_input[params.level])
        organism_count['no_hit'] = unassigned

        # Sum for the same organisms. This is important with genera as there can be multiple species with the same genus.
        organism_count = organism_count.groupby(organism_count.index).sum()

        if any(organism_count.index.duplicated(keep=False)):
            organism_count[organism_count.index.duplicated()].to_csv(params.duplicate_names,sep='\t',index=True,header=False)
            organism_count = organism_count.groupby(organism_count.index,sort=False).sum()

        organism_count.to_csv(output.TSV_motus,sep='\t',index=True,header=False)

