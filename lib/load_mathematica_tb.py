import numpy as np
import pandas as pd
import  matplotlib.pyplot as plt

class cavity:
    '''This class loads the cavity results computed with Mathematica script "cavity.wls".
        mathematica script creates curves with percolation p=0.95
    '''

    def __init__(self, path,gene_on=True):

        self.gene_on = gene_on
        if gene_on:
            self.gene_cavity_tb = self.load_theoy_from_mathematica(path+"/genes_on_cavity.txt")
            self.TF_cavity_tb = self.load_theoy_from_mathematica(path+"TFs_on_cavity.txt")

        else:
            self.gene_cavity_tb = self.load_theoy_from_mathematica(path+"/genes_off_cavity.txt")
            self.TF_cavity_tb = self.load_theoy_from_mathematica(path+"TFs_on_cavity.txt")

    def get_transition_line(self):
        '''This method gives the transition line of the discontinuous transition. Returns the parameters used by the mathematica file.
        Returns :
            column id
                (<cout>-1 for shifted Poisson)
            row id (<dout>-1 for shifted Poisson)
        '''

        cond = np.where(self.gene_cavity_tb > 0.5, True, False)
        col_id = np.array(list(map(lambda row: min(np.arange(len(row))[row]) if any(row) else None,
                              cond))) # index of the row where the transition happens
        col_id_2 = col_id[col_id != None].astype(int) #filter the valid index position

        cs_thr = self.gene_cavity_tb.columns[col_id_2]
        ds_thr = self.gene_cavity_tb.index[col_id!=None]
        return cs_thr , ds_thr
    def load_theoy_from_mathematica(self, file):
        '''Create a Pandas dataframe, which contains "d" as index of rows, "c" as column name, content is the number of solutions'''

        def spacchetta(row):
            ds, root = zip(*[[a for a in row_list.strip(" {").split(",")] for row_list in row.strip("}").split("},")])
            ds = np.array(ds, dtype=float)
            root = np.array(root, dtype=float)
            return ds, root

        tb = pd.read_table(file, header=None, names=["c", "ds"], index_col=0)
        A = []
        for row in tb["ds"]:
            ds_cavity, root = spacchetta(row)
            A += [root]
        A = np.transpose(A)
        num_cavity = pd.DataFrame(A, index=list(map(lambda x: round(x, 3), ds_cavity)),
                                  columns=list(map(lambda x: round(x, 3), tb.index)))
        num_cavity.index.name = "d"
        return num_cavity

    def heatmap(self,**kwargs):
        fig = plt.figure()
        cavity_tb = self.gene_cavity_tb
        '''
        if self.TF_on:
            title("Number of solution. TF on")
        else:
            title("Number of solutions. TF off")
        '''
        #cmap = plt.cm.get_cmap('viridis')
        # bounds = np.linspace(0.5, 3.5, 4)
        # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        plt.imshow(cavity_tb, origin="lower", extent=(
        min(list(cavity_tb.index)) + 1, max(list(cavity_tb.index)) + 1, min(list(cavity_tb.columns)) + 1,
        max(list(cavity_tb.columns)) + 1), aspect="auto", **kwargs)
        # colorbar()
        plt.xlabel("Mean gene in-degree")
        plt.ylabel("Mean TF in-degree")
        plt.tight_layout()
        return fig

