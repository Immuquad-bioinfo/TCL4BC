import os
import pandas as pd
import numpy as np
import argparse
import datetime

def argparse_analyse():
    parser = argparse.ArgumentParser(description='ImmuScan-BC predict program')
    parser.add_argument('-i', type=str, required=True, help='sample_path')
    parser.add_argument('-o', type=str, required=True, help='output_path')
    args = parser.parse_args()
    return args.i, args.o

def index_score_cal(sample_path,output_path):
    Rscript_cmd="Rscript ./predict/cal_divindx.R "+sample_path+" "+output_path
    os.system(Rscript_cmd)

    while True:
        if os.path.exists(os.path.join(output_path,'stats_divIndx.csv')):
            break

    index_df=pd.read_csv(os.path.join(output_path,'stats_divIndx.csv'))
    index_df['index_score'] = 6.70860728 + 2.60231325 * index_df[
            "Clonality"] + 0.25264915 * index_df[
                "Shannon.Index"] - 1.97338766 * index_df[
                    "Hvj.Index"] - 1.65621167 * index_df[
                        "Singleton"] + 12.20182554 * index_df["Simpson.Index"]
    index_df=index_df[['FileName','index_score']]
    index_df.columns=['sample_id','index_score']
    return index_df

def cal_vdj_score(input_path):
    vdj_df = pd.read_csv(
        './predict/HD_penalty.list',
        sep='\t')
    penalty_df=pd.DataFrame(columns=['sample_id', 'vdj_score'])
    for sample_name in os.listdir(input_path):
        if not sample_name.endswith('.txt'):
            continue
        input_file=os.path.join(input_path, sample_name)
        df=pd.read_csv(input_file,sep='\t')
        df=pd.merge(df, vdj_df, on=['v', 'd', 'j'], how='inner')
        if df.shape[0]==0:
            vdj_score=0
        else:
            vdj_score=np.log10(df['freq'].sum()/0.005)*100 
        penalty_df = pd.concat([
        penalty_df, 
        pd.DataFrame([{'sample_id': sample_name.split('.')[0], 'vdj_score': vdj_score}])
], ignore_index=True)        
        penalty_df['label'] = input_path.split('/')[-1]
    return penalty_df

def cal_key_item(sample_path):
    feature_list='./predict/feature.txt'
    standard_list = ["cdr3aa"]

    feature_df = pd.read_csv(feature_list, sep="\t")
    feature_df = feature_df[standard_list]

    sample_list = [
        sample_name for sample_name in os.listdir(sample_path)
        if sample_name.endswith(".txt")
    ]
    clonotypes_list = []
    associated_clonotype_list = []

    for sample_name in sample_list:
        sample_file_name = os.path.join(sample_path, sample_name)
        sample_df = pd.read_csv(sample_file_name, sep="\t")

        merged_df = pd.merge(feature_df,
                             sample_df,
                             on=standard_list,
                             how="inner")

        associated_clonotype = merged_df[standard_list].shape[0]

        clonotypes = sample_df[standard_list].shape[0]

        clonotypes_list.append(clonotypes)
        associated_clonotype_list.append(associated_clonotype)

    result_df = pd.DataFrame({
        "sample_id": [sample_name.split('.')[0] for sample_name in sample_list],
        "clonotypes": clonotypes_list,
        "associated_clonotype": associated_clonotype_list
    })
    return result_df


def cal_tcl_score(clonotypes, associated_clonotype):
    if clonotypes == 0 or associated_clonotype == 0:
        return 20, '阴性'

    ratio = associated_clonotype / clonotypes
    max_value_ratio = -2.5
    high_value_ratio = -3
    max_threshold_ratio = -3.3
    min_value_ratio = -4.0
    ratio_log = np.log10(ratio)

    if ratio_log<min_value_ratio:
        clonotype_score_ratio = 20
    elif ratio_log < max_threshold_ratio:
        clonotype_score_ratio = 20 + (ratio_log - min_value_ratio) * 40 / (
            max_threshold_ratio - min_value_ratio)
    elif ratio_log < high_value_ratio:
        clonotype_score_ratio = 60 + (ratio_log - max_threshold_ratio) * 20 / (
            high_value_ratio - max_threshold_ratio)

    elif ratio_log < max_value_ratio:
        clonotype_score_ratio = 80 + (ratio_log - high_value_ratio) * 20 / (
            max_value_ratio - high_value_ratio)
    else:
        clonotype_score_ratio = 100

    return  clonotype_score_ratio


    
def run_sample_predict(sample_path, output_path):
    ## calculate index score
    index_df=index_score_cal(sample_path,output_path)

    ## calculate vdj score
    vdj_df=cal_vdj_score(sample_path)

    ## calculate tcl score
    tcl_df=cal_key_item(sample_path)
    tcl_df['tcl_score']=tcl_df.apply(lambda x: cal_tcl_score(x['clonotypes'], x['associated_clonotype']), axis=1)
    tcl_df=tcl_df[['sample_id','tcl_score']]

    ## calculate final score
    result_df=pd.merge(index_df, vdj_df, on='sample_id', how='left')
    result_df=pd.merge(result_df, tcl_df, on='sample_id', how='left')

    result_df['final_score']=0.3*result_df['index_score']+0.6*result_df['tcl_score']+0.1*result_df['vdj_score']
    result_df['predict_label']=result_df['final_score'].apply(lambda x: '阳性' if x>=80 else '警戒' if x>=60 else '阴性')   
    result_df=result_df[result_df['final_score'].notna()]

    timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    result_df.to_csv(os.path.join(output_path, f'TBCP_{timestamp}_score.csv'), index=False)


if __name__ == '__main__':
    sample_path, output_path = argparse_analyse()
    os.makedirs(output_path, exist_ok=True)

    run_sample_predict(sample_path, output_path)
