import pandas as pd
import numpy as np
from IPython.display import display, HTML



class OrBi():
    def __init__(self, mat):
        self.mat = mat
        self.df = pd.DataFrame(mat)  
        
        
    def gini_mu(self, rows,cols):
        _, r_mu = self.r_par(cols)
        _, c_mu = self.c_par(rows)
        return r_mu+c_mu-r_mu*c_mu
        
        
    def r_par(self,cols):    
#         display(self.df.loc[:,cols])
        print('Columns')
        display(self.df.style.applymap(self._bg_map,
                               subset=pd.IndexSlice[:,cols]))
        
        
        mat_pat,index, counts = np.unique(self.mat[:,cols],axis=0,return_index=1, return_counts=1)
#         print(mat_pat,index,counts)
        print(f'Partitions of columns on row')
        display(pd.DataFrame(mat_pat, index=[counts,index]))
        
        
        gini = 1 - np.sum(np.square(counts/np.sum(counts)))        
        gini_m = 1 - np.dot(np.square(counts/np.sum(counts)), 
                            np.sum(mat_pat,axis=1)/mat_pat.shape[1])
        print(f'For row partitions on cols of{cols}:\n\n \
                gini : {gini}\n   \
              gini_m : {gini_m}           ')
        return gini, gini_m
        
        
        
    def c_par(self,rows):    
#         display(self.df.loc[:,rows])
        display(self.df.style.applymap(self._bg_map,
                               subset=pd.IndexSlice[rows,:]))
        
        print(f'Partitions of rows on columns')

        mat_pat,index, counts = np.unique(self.mat[rows,:],axis=1,return_index=1, return_counts=1)
        display(pd.DataFrame(mat_pat, columns=[counts,index]))
        
        
        gini = 1 - np.sum(np.square(counts/np.sum(counts)))        
        gini_m = 1 - np.dot(np.square(counts/np.sum(counts)), 
                            np.sum(mat_pat,axis=0)/mat_pat.shape[0])
        print(f'For column partitions on cols of{rows}:\n\n \
                gini : {gini}\n   \
              gini_m : {gini_m}           ')
        return gini, gini_m
        
    def proj(self, row, cols):
        
        print(f'Projections of columns {cols} on row {row}')
        return self.df.style.applymap(self._bg_map,
                               subset=pd.IndexSlice[row,cols])
    def rstr(self, col, rows):
        
        print(f'Restriction of rows {rows} on column {col}')        
        return self.df.style.applymap(self._bg_map,
                               subset=pd.IndexSlice[rows,col])
    
        
    def _bg_map(self,val, color='yellow'):
        return f'background-color: {color}'