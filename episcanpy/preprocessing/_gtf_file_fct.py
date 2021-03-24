import pandas as pd

def check_gtf_composition(gtf_file, annotation=None, feature_type='gene'):
    """
    annotation can be either HAVANA or ENSEMBL
    feature_type can be gene, transcript, exon 
    """
    # loading the gtf file
    mtx = []
    with open(gtf_file) as f:
        for line in f:
            if line[0] != '#':
                mtx.append(line.split('\t'))
            
    mtx_df = pd.DataFrame(mtx)

    if annotation != None:
        is_true =  mtx_df[1]==annotation
        mtx_df2 = mtx_df[is_true]
    else: 
        mtx_df2 = mtx_df.copy()

    is_true =  mtx_df2[2]==feature_type
    mtx_df3 = mtx_df2[is_true]

    list_genetype = [x.split(';')[1].split('"')[1] for x in mtx_df3[8].tolist()]

    count_list = {}
    for key in list(set(list_genetype)):
        count_list[key] = 0
    for n in list_genetype:
        count_list[n] += 1
    return(count_list)
    