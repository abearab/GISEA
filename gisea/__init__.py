import pandas as pd
import blitzgsea as blitz


def run_gsea(data,gmt,processes=6):
    FullResults = {}

    for i,_ in enumerate(data.index):
        signature = data.iloc[i,:].reset_index()

        result = blitzgsea.gsea(signature=signature,library=gmt,processes=processes)

        FullResults[i] = result.to_dict()
        del result, signature
    return FullResults


def gsea_to_matrix(data, results, info, pivot=True):
    """
    'es', 'nes', 'pval', 'sidak', 'fdr', 'geneset_size', 'leading_edge', 'leading_edge_size'
    """
    obs_dict = dict([(str(i),gene) for i,gene in enumerate(data.index)])
    mat = [(obs,term,val) for obs in results for term,val in results[obs][info].items()]
    mat = pd.DataFrame(mat,columns=['Gene','Term',info])
    if pivot:
        mat = mat.pivot(index=['Gene'],columns=['Term'],values=info)
    mat.columns.name = None
    
    mat = mat.loc[obs_dict.keys(),]
    mat.index = obs_dict.values()
    
    return mat
