import h5py
import numpy as np
import pandas as pd
import json
import re
import functools
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist
import time

class h5Reader:
    def __init__(self, filename):
        with h5py.File(filename,"r") as h5_data:
            self.L =list(np.array(h5_data['celltype']['gene_expression_average']['axis0'].asstr()))
        self.filename = filename

    def read_file(self, df_type,genes=None):
        """
        read in a .h5 file from the source directory,return a dataframe based on the user's specified dataset type and genes of interest

        :param df_type: value taken from celltype/celltype_dataset/celltype_dataset_timepoint
        :type df_type: string
        :param genes: a list of gene names, defaults to None
        :type genes: string, optional
        :return: a table(dataframe) with rows being genes and columns being all celltypes,values in the table are the expression level of genes in different cell types.
        :rtype: dataframe
        """
        genes = genes.split(",")
        with h5py.File('./static/scData/condensed_lung_atlas_in_cpm.h5',"r") as h5_data:
            # List of genes
            indexs = []
            if isinstance(genes,list):
                for gene in genes:
                    indexs.append(self.L.index(gene.capitalize()))
            else:
                indexs.append(self.L.index(genes.capitalize()))
            # sort the indexs
            indexs.sort()
            new_genes = []
            for index in indexs:
                new_gene = self.L[index]
                new_genes.append(new_gene)
            # expression table (only numbers,no index nor columns)
            data = np.array(h5_data[df_type]['gene_expression_average']['block0_values'][:, indexs]).astype(np.float32)
            columns=np.array(h5_data[df_type]['gene_expression_average']['axis1'].asstr())
            df = pd.DataFrame(
                data = data.T,
                index = new_genes,
                columns=columns,
            )

        return df

    def get_marker_genes_list(self,celltype=None):
        """
        For each available cell types in the whole dataset, their marker genes are pre-computed and stored in a large dictionary.
        By specifing the cell type, this function return a list of marker genes which belong to it. 

        :param celltype: cell type of interest, defaults to None
        :type celltype: string, optional
        :return: list of gene names
        :rtype: list
        """
        with open('./static/scData/celltypeMarkerGeneList.txt') as f:
            data = json.load(f)
            
        return data[celltype]

    def marker_genes_expression(self,celltype):
        """_summary_

        :param celltype: _description_
        :type celltype: _type_
        :return: _description_
        :rtype: _type_
        """
        
        gene_list = self.get_marker_genes_list(celltype)
        gene_list = ','.join(gene_list)
        df = self.read_file("celltype",gene_list)
        df['current'] = df[celltype]
        for column in df.columns:
            df[column] = (df[column] / df['current']).round(3)
        df.drop(['current'], axis=1,inplace=True)
        return json.loads(df.to_json())


    def timepoint_reorder(self,tp1, tp2):
        """
        sort the order of two given timepoints, such that timepoints are in order:
        E18.5, P1 P3 P7 P14 P21, 3m 18m 24m

        :param tp1: the first timepoint
        :type tp1: string
        :param tp2: the second timepoint
        :type tp2: string
        :return: 0 if tp_1 == tp_2,
                -1 if tp_1 comes before tp_2,
                1 if tp_2 comes before tp_1,
        :rtype: int
        """
        if '_' in tp1:
            tp1 = tp1.split('_')[1]
        
        if '_' in tp2:
            tp2 = tp2.split('_')[1]

        if tp1 == tp2:
            return 0
        
        if 'P' in tp1 and 'P' in tp2:
            value_1 = int(tp1.split('P')[1])
            value_2 = int(tp2.split('P')[1])
            if value_1 < value_2:
                return 1
            else:
                return -1
        
        if 'm' in tp1 and 'm' in tp2:
            value_1 = int(tp1.split('m')[0])
            value_2 = int(tp2.split('m')[0])
            if value_1 < value_2:
                return 1
            else:
                return -1
        
        # E always before others
        if 'E' in tp1 and 'E' not in tp2:
            return 1
        
        if 'E' in tp2 and 'E' not in tp1:
            return -1
        
        # m always after others
        if 'm' in tp1 and 'm' not in tp2:
            return -1

        if 'm' in tp2 and 'm' not in tp1:
            return 1 


    #  function that get the cell type dataset timepoint as a dictionary
    # user input a gene name.
    #  for each unique dataset of this gene, we plot a heatmap of all timepoint vs celltypes 
    def dataset_by_dataset(self,genename,df_type):
        """
        Generate a dictionary that contains data for plotting the dataset Heatmap of a selected gene

        :param genename: gene of interest
        :type genename: string
        :param df_type: value taken from celltype/celltype_dataset/celltype_dataset_timepoint
        :type df_type: string
        :return: a dictionary
        :rtype: dictionary
        """
        genename = genename.capitalize()
        # select and pre-preprocessing data
        df_tp = self.read_file(df_type,genename)

        # select data for a given gene name
        df_filtered=df_tp.loc[[genename]]
        # for each dataset name, we construct a new dataframe for it (celltypes x timepoints)
        datasets = set([name.split("_")[1] for name in df_filtered.columns])
        dic_per_dataset = {}
        new_cell_type_order = {}
        # ACZ,Hur.2021,TMS
        for i in datasets:
            columns_with_this_dataset = []
            for col_name in df_filtered.columns:
                if i in col_name:
                    columns_with_this_dataset.append(col_name)
            dataset_name = df_filtered[columns_with_this_dataset] 
            timepoints = list(set([name.split("_")[2] for name in dataset_name.columns]))
            celltypes = set([name.split("_")[0] for name in dataset_name.columns])
            
            # re-arrange the dataframe based on timepoint order:
            # e.g: P21, P3
            # sorted function depends on key(0,1,-1), output of the sorting function
            timepoints = sorted(timepoints, key=functools.cmp_to_key(self.timepoint_reorder))

            # create a new empty dataframe for each dataset   
            gene_exp_df = pd.DataFrame(np.eye(len(celltypes),len(timepoints)), columns=timepoints, index=celltypes)
            for tp in timepoints:
                for ct in celltypes:
                    col_name = '_'.join([ct,i,tp])
                    # celltype Car4+ cappilaries can not be recognized, as it contains a + sign. which means one or more 4 in regex
                    # replace regex '+' with string '\+'
                    col_name = re.sub('\+', '\+', col_name)
                    try:
                        # try to see if the col_name is in df
                        gene_exp_df[tp][ct] = dataset_name[col_name].values[0]
                    except:
                        # col_name not in df
                        pass

            # celltype order for hierarchical clustered set
            distance = pdist(gene_exp_df.values)
            Z = linkage(distance,optimal_ordering=True)
            new_order = leaves_list(Z)
            clustered_gene_exp_df = gene_exp_df.iloc[new_order]
            new_cell_type_order[i] = [ct for ct in clustered_gene_exp_df.index]
            dic_per_dataset[i] = json.loads(gene_exp_df.to_json())
        return {"result":dic_per_dataset, "hierarchicalCelltypeOrder":new_cell_type_order}
        

    def dataset_unified(self,genename):

        """
        Generate a json object that contains data for the Unified heatmap

        :param genename: gene of interest
        :type df_type: string
        :return: a dictionary contains expression value of a gene in different timepoints and celltypes
        :rtype: json object
        """
        genename = genename.capitalize()
        df = self.read_file('celltype_dataset_timepoint',genename)
        filtered_df = df.filter(items=[genename],axis=0)
        all_celltypes = []
        dt_combinations = []  # store all the existing dataset and timepoint combinations

        for column_name in filtered_df.columns:
            celltype = column_name.split("_")[0]
            dataset_timepoint = column_name.split(celltype+"_")[1]
            if celltype not in all_celltypes:
                all_celltypes.append(celltype)
            if dataset_timepoint not in dt_combinations:
                dt_combinations.append(dataset_timepoint) 

        gene_exp_df = pd.DataFrame(np.eye(len(all_celltypes),len(dt_combinations)), columns=dt_combinations, index=all_celltypes)
        for dt in dt_combinations:
            for ct in all_celltypes:
                name = "_".join([ct,dt])
                if name not in filtered_df.columns:
                    exp_value = -1
                else:
                    exp_value = float(filtered_df[name].values[0])

                gene_exp_df[dt][ct] = exp_value
        
        cell_type_label = all_celltypes
        distance = pdist(gene_exp_df.values)
        Z = linkage(distance,optimal_ordering=True)
        new_order = leaves_list(Z)
        clustered_gene_exp_df = gene_exp_df.iloc[new_order]
        new_cell_type_order = [ct for ct in clustered_gene_exp_df.index]
        dt_sorted = sorted([key for key in gene_exp_df.columns], key=functools.cmp_to_key(self.timepoint_reorder))

        return {'expression': json.loads(gene_exp_df.to_json()), 'dataset_timepoint': dt_sorted, 'cell_type': cell_type_label, 'gene': genename, 'hierarchicalCelltypeOrder':new_cell_type_order}
        
