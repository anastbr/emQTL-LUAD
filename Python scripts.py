#########################################################################################################################
############################ Spectral Co-Clustering Heatmap (in Python) ##################################################################
##################################################################################################################
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.cluster import SpectralCoclustering

# import emQTL matrix
r = pd.read_csv(my_directory + "/r_final_discovery_forClustering_SexFiltered_47296x6338.csv")

m = r.set_index('Unnamed: 0')
m = m.astype("float32")
model = SpectralCoclustering(n_clusters=5, random_state=0)
model.fit(m)

fit_data = m.values[np.argsort(model.row_labels_)]
fit_data = fit_data[:, np.argsort(model.column_labels_)]

plt.figure(figsize=(10,10))
fig = sns.heatmap(fit_data, vmin = -1, vmax=1, center = 0, 
		xticklabels = False, yticklabels = False, cmap="RdBu")
plt.xlabel("6,338 Genes", fontsize=15)
plt.ylabel("47,296 CpG", fontsize=15)
plt.rcParams["ytick.labelsize"]=35
plt.savefig(my_directory + "/Spectral_CoClustering_5_Biclusters", dpi=600)
##################################################################################################################


#########################################################################################################################
############################ Elbow method based on Distortion Score and K-Means (in Python) ##################################################################
##################################################################################################################

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sns


#from sklearn.cluster import SpectralCoclustering
r = pd.read_csv(my_directory + "/r_final_discovery_forClustering_SexFiltered_47296x6338.csv")
m = r.set_index('Unnamed: 0')
from sklearn.cluster import KMeans
from yellowbrick.cluster import KElbowVisualizer

model = KMeans(random_state=0)
# k is range of number of clusters.
visualizer = KElbowVisualizer(model, k=(2,20), timings= False)
visualizer.fit(m)        # Fit data to visualizer

visualizer.show(outpath = my_directory + "/Elbow_method_distortion_score.png")
#########################################################################################


#########################################################################################################################
############################ Elbow method based on MSR Score (in Jupyter Notebook) ##################################################################
##################################################################################################################
import pandas as pd
import numpy as np
from sklearn.cluster import SpectralCoclustering
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib import pyplot as plt
import glob
import csv

r = pd.read_csv(my_directory + "/r_final_discovery_forClustering_SexFiltered_47296x6338.csv")
r = pd.read_csv(my_directory + "/r_final_discovery_forClustering_SexFiltered_47296x6338.csv")
m = r.set_index('Unnamed: 0')

def mean_squared_residue(bicluster):
    nrows, ncols = np.shape(bicluster)
    bic_avg = np.mean(bicluster.values)
    avg_across_cols = np.mean(bicluster, axis=0)[np.newaxis, :]
    avg_across_rows = np.mean(bicluster, axis=1)[:, np.newaxis]
       
    avg_bic = (avg_across_cols + avg_across_rows) - bic_avg
            # Sanity check.
    assert np.shape(avg_bic) == np.shape(bicluster)

    msr_values = (bicluster - avg_bic) ** 2
    return np.sum(msr_values.values) / (nrows * ncols)


average_per_clustering = {}
for n_cluster in range(2,21):
    mean = []
    model = SpectralCoclustering(n_cluster, random_state=0)
    model.fit(m)
    bicluster = [model.get_submatrix(i, m) for i in range(n_cluster)]
    #print("bicluster ", len(bicluster)) check point
    #bicluster = pd.DataFrame(model.get_submatrix(i, m))
    for j in range(0,n_cluster):
       #print(bicluster)
        biclust_mrs = mean_squared_residue(pd.DataFrame(bicluster[j]))
        mean.append(biclust_mrs)
        #print("n_cluster ",n_cluster," j ",j)
    #print(mean) check point
    average_clust = sum(mean) / len(mean)
    average_per_clustering[n_cluster] = average_clust
    print(average_per_clustering)
	
#SAVING THE DICTIONARY with headers
myDic = average_per_clustering
with open('n_clusters_all_and_MRS_LUAD_gendered_filtered.csv', 'w', newline='') as csvfile:
    fieldnames = ['n_clusters', 'msr']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for key in myDic:
        writer.writerow({'n_clusters': key, 'msr': myDic[key]})		
my_dict = pd.read_csv("n_clusters_all_and_MRS_LUAD_gendered_filtered.csv")

plt.figure(figsize=(10,5))
sns.lineplot(data=my_dict, x=my_dict.n_clusters,y=my_dict.msr)
plt.xticks(np.arange(min(my_dict.n_clusters),max(my_dict.n_clusters)+1,1.0))
plt.title("Elbow method LUAD")
plt.ylabel("Average MSR score")
plt.xlabel("Number of clusters")
plt.savefig(my_directory + "/Elbow_MSR_LUAD.png",dpi=600)
###################################################################################################