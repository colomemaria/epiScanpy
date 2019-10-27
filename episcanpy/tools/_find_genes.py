import anndata as ad

def top_feature_genes(adata, gtf_file, extension=5000):
    
    # extract_top_markers
    #print(adata.uns['rank_genes_groups'].keys())
    windows = [list(w) for w in adata.uns['rank_genes_groups']['names'].tolist()]
    windows_all = []
    for w in windows:
        windows_all += w

    # load the gtf file
    with open(gtf_file) as f:
        gtf_raw = f.readlines()
    gtf = []
    for line in gtf_raw:
        if line[0] != '#':
            gtf.append(line[:-2].split('\t'))
    del gtf_raw
    
    gtf_dic = {}
    for line in gtf:
        if line[0] not in gtf_dic.keys():
            gtf_dic[line[0]] = [line]
        else:
            gtf_dic[line[0]].append(line)
    del gtf

    markers = []
    for w in windows_all:
        curr_m = []
        w2 = w.split('_')
        chrom = w2[0][3:]
        w2 = [int(x) for x in w2[1:]]
        for gene in gtf_dic[chrom]:
            start = int(gene[3])-extension
            end = int(gene[4])+extension
            if (w2[0] < start < w2[1]) or (w2[0] < end < w2[1]):
                gene_name = gene[-1]
                for n in gene[-1].split(';'):
                    if 'gene_name' in n:
                        gene_name = n
                curr_m.append([w, gene_name])
        if curr_m != []:
            markers += curr_m
    markers = [list(x) for x in set(tuple(x) for x in markers)]
    markers_dict={}
    for n in markers:
        markers_dict[n[1].split(' "')[1][:-1]] = n[0]
    return(markers_dict)


### function to optimize running time

def find_genes(adata, gtf_file_name, path='', extension=5000,
               key_added='gene_name', feature_coordinates=None, copy=True):
    """
    Given a gtf file, you can match peak coordinates (stored in adata.var_names or
    in a var annotation) to genes.
    The peak annotation has to be written as chr1:20000-20500 or chr1_20000_20500.
    the corresponding gene (if any) will be sotred in a var annotation 
    It extend the search to match a gene to an window of + and - extensions size(5kb
    for example).
    """
    #start = time.time()

    # load the gtf file
    gtf_file = []
    with open(gtf_file_name) as f:
        for line in f:
            if line[0] != '#':
                gtf_file.append(line)
    gtf_file = pd.DataFrame([l.split('\t') for l in gtf_file])
    gtf_file.columns = ['Chromosome', 'source', 'gene_type', 'Start', 'End',
                   'NA', 'Strand', 'NA2', 'extra_info'] 
    
    del gtf_file['NA'], gtf_file['NA2']
    
    # extract the variable names
    feature_annot = adata.var
    if feature_coordinates ==None:
        feature_names = adata.var_names.tolist()
    else:
        feature_names = adata.var[feature_coordinates]
    
    # format the feature name
    start_feature = []
    end_feature = []
    chrom_feature = []
    for feature in feature_names:
        if ':' in feature: # if the feature is a Granger coordinate.
            feature2 = feature.split(':')
            w = [int(x) for x in feature2[1].split('-')]
        else:
            feature2 = feature.split('_')
            w = [int(x) for x in feature2[1:]]
        chrom_feature.append(feature2[0][3:])
        start_feature.append(w[0])
        end_feature.append(w[1])
    
    adata.var['Index'] = range(0,len(chrom_feature))
    adata.var['Chromosome'] = chrom_feature
    adata.var['Start'] = start_feature
    adata.var['End'] = end_feature
    adata.var['name_feature'] = adata.var_names.tolist()
    #adata.var['chrom_feature'] = chrom_feature
    adata.var['start_ext'] = [x-extension for x in start_feature]
    adata.var['end_ext'] = [x+extension for x in end_feature]
    #adata.var['Strand'] = len(end_feature)*['+']
    
    
    # match the feature with 
    gtf = pr.PyRanges(gtf_file)
    del gtf_file
    adata_var = pr.PyRanges(chromosomes=adata.var.loc[:,'Chromosome'], #strands=adata.var.loc[:,'strand_feature'],
               starts=adata.var.loc[:,'start_ext'], ends=adata.var.loc[:,'end_ext'])
    
    merge = gtf.join(adata_var, suffix="_ext")
    merge = merge.dfs
    overlap3 = pd.concat([merge[key] for key in merge.keys()])
    overlap3['Index'] = overlap3.index
    overlap4 = overlap3.sort_values(['Chromosome', 'Start_ext', 'End_ext', 'Index'])
     
    #print(time.time()-start)
    
    adata.var = adata.var.sort_values(['Chromosome', 'start_ext', 'end_ext'])
    adata_var = pr.PyRanges(adata.var)
    tot_gene_annot = []
    for chrom in list(set(adata.var['Chromosome'])):
        index_gtf = 0
        #next_index = 0
        #curr_adata = adata_var[chrom].df
        overlap3 = pd.concat([merge[key] for key in [(chrom, '+'), (chrom, '-')]])
        overlap3['Index'] = overlap3.index
        overlap3 = overlap3.sort_values(['Chromosome', 'Start_ext', 'End_ext', 'Index'])
        overlap_chrom = overlap3['Start_ext'].tolist()
        #for line_adata in curr_adata[['start_ext']].iterrows():
        for line_adata in adata_var[chrom].df[['start_ext']].iterrows():
            gene_annot = []
            for start_gtf in overlap_chrom[index_gtf:]:
                if start_gtf == line_adata[1][0]:
                    gene_annot.append(overlap3.iloc[index_gtf])
                    index_gtf += 1
                    continue
                else:
                    #index_gtf = next_index
                    break
                
            if gene_annot == []:
                tot_gene_annot.append(('NA'))
            else:
                tot_gene_annot.append(tuple(gene_annot))
        #print(chrom, time.time()-start)
    
    
    adata.var[key_added] = tot_gene_annot
    adata.var.sort_values(['Index'])
    #print(time.time()-start)
    return(tot_gene_annot, overlap4)
