import pandas as pd
import os
from datetime import datetime

chr_list = [str(chromosome) for chromosome in range(1, 23)] + ['X', 'Y']


def prepare_input(input_file, annovar_dir, spliceai_dir):
    """
    Prepare input file for ANNOVAR and SpliceAI.
    :param input_file: Original input filename of MAGPIE. Input file should be tab or comma separated files containing five headers: Chr, Start, End, Ref, Alt
    :param annovar_dir: Dictionary of ANNOVAR input file. (e.g. /data/output/annovar/)
    :param spliceai_dir: Dictionary of SpliceAI input file. (e.g. /data/output/spliceai/)
    :return: None
    """
    filename = os.path.splitext(os.path.basename(input_file))[0]
    data = pd.read_csv(input_file, low_memory=False)
    data.Chr = [str(i).replace('chr', '') if 'chr' in str(i) else str(i) for i in data.Chr.tolist()]
    data = data[['Chr', 'Start', 'End', 'Ref', 'Alt']]
    data = data[data.Chr.isin(chr_list)]
    data = data.fillna('-')
    data.insert(data.shape[-1], 'info', '.')
    data.to_csv(f'{os.path.join(annovar_dir, filename)}.avinput', sep='\t', index=False, header=False)
    data.to_csv(f'{os.path.join(spliceai_dir, filename)}.vcf', sep='\t', index=False, header=False)
    now = datetime.now()
    str_to_insert = f'''##fileformat=VCFv4.1
##fileDate={"{:02d}-{:02d}-{:02d}".format(now.year, now.month, now.day)}
##source=magpie
##reference=GRCh38
##contig=<ID=1,length=248956422>
##contig=<ID=2,length=242193529>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##ID=<Description="ClinVar Variation ID">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
'''
    temp_filename = f"{os.path.join(spliceai_dir, filename)}_temp.txt"
    spliceai_file = f'{os.path.join(spliceai_dir, filename)}.vcf'
    with open(temp_filename, 'w') as temp_file:
        temp_file.write(str_to_insert)
        with open(spliceai_file, 'r') as file:
            temp_file.write(file.read())
    os.remove(spliceai_file)
    os.rename(temp_filename, spliceai_file)


def prepare_training_data(data):
    data = data.rename(columns={'Gene.refGene': 'gene'})
    data = pd.get_dummies(data, columns=['func', 'omim'])
    for added_feature in ['func_frameshift', 'func_nonframeshift', 'func_nonsynonymous SNV',
                          'func_startloss', 'func_stopgain', 'func_stoploss',
                          'omim_Autosomal_dominant', 'omim_Autosomal_recessive',
                          'omim_X_linked_dominant', 'omim_X_linked_recessive', 'omim_other']:
        if added_feature not in data.columns:
            data.insert(data.shape[-1], added_feature, 0)

    return data
