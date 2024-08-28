# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
import os
import logging
from jmct.nbinom import analyzeSampleForMutationClusters
from jmct import csvDataFrameImport, drawSNPMap, getGenomeSize,  time_elapsed
from jmct.settings import settings
import datetime

@time_elapsed
def main():
    DATA_PATH = settings['data_path']
    PVAL_LIST = settings['pval_list']
    MAX_CLUST_DIST = settings['max_clust_distance']
    MIN_CLUST_SIZE = settings['min_clust_mutations']

    time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    logging.basicConfig(
        handlers=[logging.FileHandler(
            filename=f'logs/run_{time_now}.log',
            encoding='utf-8',
            mode='w')],
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:.%(funcName)s: %(message)s')

    if not os.path.exists(f'figures'):
        os.makedirs(f'figures')
    if not os.path.exists(f'outputs'):
        os.makedirs(f'outputs')

    sample_list = csvDataFrameImport(data_path=DATA_PATH, sample_name_len=settings['sample_name_len'])
    for df_sample in sample_list:
        sample_pval_dict = {}
        for pval in PVAL_LIST:
            df_clust = analyzeSampleForMutationClusters(
                df_sample[0],
                genome_size=getGenomeSize(), 
                min_p_value=pval,
                max_distance_between_mutations=MAX_CLUST_DIST,
                min_cluster_mutations=MIN_CLUST_SIZE
                )
            
            sample_pval_dict[pval] = df_clust

            logging.info(f"Found {len(df_clust)} clusters with p-value threshold {pval} for sample {df_sample[1]}.")

            df_clust.to_csv(f"outputs/{df_sample[1]}_{pval}_clusters.csv", index=False)
        
        drawSNPMap(
            df_SNVs=df_sample[0],
            pval_dict=sample_pval_dict,
            sample_name=df_sample[1],
            saveMap=True,
            showMap=False)
        
    logging.info(f"Analysis complete. Results saved in outputs/ and figures/ directories. Check logs/ for details.")

if __name__ == '__main__':
    main()
