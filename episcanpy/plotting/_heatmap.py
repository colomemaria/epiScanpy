import pandas as pd
import seaborn as sns
import numpy as np


def prct_overlap(adata, key_1, key_2, color="Blues", norm=False, ax_norm="line", sort_index=False):
    """
    % of overlap between cell types. 
    """
    
    data_1 = adata.obs[key_1].tolist()
    data_2 = adata.obs[key_2].tolist()
    
    count = {k:[] for k in list(set(data_1))}
    #count = {k:[] for k in sorted(list(set(data_1)))}
    i = 0
    for index in data_1:
        count[index].append(data_2[i])
        i += 1
    
    total_matrix = []
    for key, value in count.items():
        value = sorted(value)
        curr_key_list = []
        for element in sorted(list(set(data_2))):
            curr_count = 0
            for v in value:
                if element == v: 
                    curr_count += 1
            curr_key_list.append(curr_count)
        curr_sum = sum(curr_key_list) 
        #total_matrix.append([x/curr_sum for x in curr_key_list])
        total_matrix.append(curr_key_list)
        
    if norm and ax_norm == "line":
        total_matrix = []
        for key, value in count.items():
            value = sorted(value)
            curr_key_list = []
            for element in sorted(list(set(data_2))):
                curr_count = 0
                for v in value:
                    if element == v: 
                        curr_count += 1
                curr_key_list.append(curr_count)
            curr_sum = sum(curr_key_list) 
            total_matrix.append([x/curr_sum for x in curr_key_list])
        
    elif norm:
        print("""error in the argument ax_norm or it is col and 
              I haven't figure out how to make it for mow.
              , here is the heatmap with no normalisation""")
            
    if sort_index:
        data_heatmap = pd.DataFrame(data=np.matrix(total_matrix),
                                    index=list(set(data_1)),
                                    columns=sorted(list(set(data_2)))).sort_index()
    else:
        data_heatmap = pd.DataFrame(data=np.matrix(total_matrix),
                                index=list(set(data_1)),
                                columns=sorted(list(set(data_2))))
    
    return(data_heatmap)


def overlap_heatmap(adata, key_1, key_2, color="Blues", norm=False, ax_norm="line", sort_index=False):
    """
    If you want to normalize the ocunt matrix you put norm = True and you put in  
    ax_norm either "line" or "col"
    I am not 100 percent convinced with the normalisation right now. I prefer the visual
    as simple cell count. Despite the bias torward large cell clusters.
    
    Also, you need pandas and seaborn 
    
    Options for colors: "YlGnBu", "BuPu", "Greens", "Blues" etc.
    """
   
    data_heatmap = prct_overlap(adata, key_1, key_2, color, norm, ax_norm, sort_index)
    
    return(sns.heatmap(data_heatmap, xticklabels=True, yticklabels=True, cmap=color))
