import numpy as np
import pandas as pd

from scipy import stats
from scipy.stats import zscore

import statsmodels.api as sm
from statsmodels.formula.api import ols

from itertools import chain

from . models import DataAnalysis , Contaminant
from . utils import (expandNCleanColumns,seperatesamples , expandcols , get_sna_cna, get_matrix_limma , calculate_adjusted_pval ,
                        impute_with_lowest5, missforest_impute)

from . R_utils import limma_diff_API

def cleanAcc(x):
    xs = list(x)
    xs = [x for x in xs if (str(x) != 'nan')]
    if len(xs) > 0:
        return str(xs[0]).strip()
    else:
        return ""

def caluculate_irs_norm(df, sample_columns , control_columns):
    irs_each_sample = []
    irs_each_control = []
    irs_samps = 0
    irs_samps_list = []
    irs_factor_list = []
    irs_cna = []
    irs_sna = [] 

    for samples in sample_columns:
        df['_irs_sum_of_'+str(irs_samps)] = df[samples].sum(axis = 1 )
        irs_samps_list.append('_irs_sum_of_'+str(irs_samps))
        irs_samps += 1

    for controls in control_columns:
        df['_irs_sum_of_'+str(irs_samps)] = df[samples].sum(axis = 1 )
        irs_samps_list.append('_irs_sum_of_'+str(irs_samps))
        irs_samps += 1

    df["_irs_gmean_"] = stats.gmean(df[irs_samps_list],axis=1)



    for itmes in irs_samps_list:
        df['factor_'+itmes] = df['_irs_gmean_'].div(df[itmes])
        irs_factor_list.append('factor_'+itmes)

    irs_samps = 0
    for samples in sample_columns:
        for replicates in samples:
            df['NORM_'+replicates] = df[replicates] * df[irs_factor_list[irs_samps]]
            irs_each_sample.append('NORM_'+replicates)

        irs_sna.append(irs_each_sample)
        irs_each_sample = []
        irs_samps += 1
    
    for controls in control_columns:
        for replicates in controls:
            df['NORM_'+replicates] = df[replicates] * df[irs_factor_list[irs_samps]]
            irs_each_control.append('NORM_'+replicates)
        
        irs_cna.append(irs_each_control)
        irs_each_control = []
        irs_samps += 1
    

    df.drop(irs_samps_list, axis=1,  inplace = True)
    df.drop(irs_factor_list, axis=1,  inplace = True)
    df.drop('_irs_gmean_',axis=1, inplace = True)
    return df ,  irs_sna, irs_cna 

def caluculate_zscore_norm(df):
    df_zscore= zscore(df)
    df_zscore.rename(lambda x: "NORM_"+x, axis='columns', inplace = True)
    return df_zscore


def normaliz_data(job_id,sample_columns,control_columns,norm_method,missing_val_rep,tmmpr , impute_method):

    data = DataAnalysis.objects.get(id = job_id)
    datafile = data.file.path
    df = pd.read_csv(datafile)

    df , contaminants_df = deletemultizero(df,sample_columns,control_columns)
    missing_val = float(missing_val_rep)

    if impute_method== 'impute-lowest5':
        df = impute_with_lowest5(df)

    elif impute_method== 'miss-forest':
        df = missforest_impute(df)
    else:
        df.fillna(missing_val, inplace = True)

    exp_samp , exp_cont , exp_sna, exp_cna = expandcols(sample_columns,control_columns)

    mediun_list = {}

    if (norm_method == 'Median'):

        for controls in control_columns:
            for replicates in controls:
                mediun_list[replicates] = df[replicates].median()

        for samples in sample_columns:
            for samp_replicates in samples:
                mediun_list[samp_replicates] = df[samp_replicates].median()


    elif (norm_method == 'TMM'):
        prcount = int(tmmpr)/100
        for controls in control_columns:
            for replicates in controls:
                mediun_list[replicates] = stats.trim_mean(df[replicates], prcount)

        for samples in sample_columns:
            for samp_replicates in samples:
                mediun_list[samp_replicates] = stats.trim_mean(df[samp_replicates], prcount)


    elif (norm_method == 'Sum'):
        for controls in control_columns:
            for replicates in controls:
                mediun_list[replicates] = df[replicates].mean()
        for samples in sample_columns:
            for samp_replicates in samples:
                mediun_list[samp_replicates] = df[samp_replicates].mean()



    elif (norm_method == "Quntail"):

        df_for_qunt = df[ exp_samp + exp_cont ]
        df_PCA_before = df_for_qunt
        quant_df = quantile_normalize(df_for_qunt)
        df_PCA_after = quant_df
        df = df.join(quant_df)
        sna , cna = get_sna_cna(sample_columns,control_columns)
        return df,df_PCA_before, df_PCA_after , cna, sna , contaminants_df
    
    elif(norm_method == "zscore"):
        
        df_for_zscore = df[ exp_samp + exp_cont ]

        df_PCA_before = df_for_zscore

        zscore_df = caluculate_zscore_norm(df_for_zscore)
        
        df_PCA_after = zscore_df
        df = df.join(zscore_df)
        sna , cna = get_sna_cna(sample_columns,control_columns)
       
        return df,df_PCA_before, df_PCA_after , cna, sna , contaminants_df


    elif(norm_method == "irs"):

        df, irs_sna, irs_cna =  caluculate_irs_norm( df , sample_columns , control_columns )
        df_PCA_before_irs = df[ exp_samp + exp_cont ]
        df_PCA_After_irs = df[ exp_sna + exp_cna ]
        return df,df_PCA_before_irs, df_PCA_After_irs , irs_cna, irs_sna , contaminants_df



    minn = min(mediun_list.values())
    #deviding each value with multiplication factor
    multiplication_fact_list = {}
    for key,value in mediun_list.items():
        multiplication_fact_list[key] = (minn/value)

    cna = []
    for controls in control_columns:
        each_control = []
        for replicates in controls:
            df['NORM_'+replicates] = df[replicates] * multiplication_fact_list[replicates]
            each_control.append('NORM_'+replicates)
        cna.append(each_control)


    sna = []
    for samples in sample_columns:
        each_sample = []
        for samp_replicates in samples:
            df['NORM_'+samp_replicates] = df[samp_replicates] * multiplication_fact_list[samp_replicates]
            each_sample.append('NORM_'+samp_replicates)
        sna.append(each_sample)

    
    df_PCA_before = df[exp_samp + exp_cont]
    df_PCA_after = df[exp_sna + exp_cna]
    return df,df_PCA_before, df_PCA_after ,cna, sna , contaminants_df


def batch_correct(df, sam_col,con_col ):
    batch_list = []
    df_list_sample = []
    i = 1
    for batch in sam_col:
        for sample in batch:
            batch_list.append(i)
            df_list_sample.append(sample)
        i +=1

    i = 1
    df_list_control = []
    for batch in con_col:
        for control in batch:
            batch_list.append(i)
            df_list_control.append(control)
        i +=1

    df_columns_for_bc = df_list_sample + df_list_control

    df_before_bc = df[df_columns_for_bc]
  

    df_after_bc = combat.pycombat(df_before_bc,batch_list)

    df.drop(df_before_bc.columns, axis = 1, inplace = True )

    return df , df_after_bc 



def normaliz_data_bio(job_id,sample_columns,control_columns,norm_method,missing_val_rep,tmmpr , impute_method):
    
    data = DataAnalysis.objects.get(id = job_id)
    datafile = data.file.path
    df = pd.read_csv(datafile)

    missing_val = float(missing_val_rep)
    exp_samp , exp_cont , exp_sna, exp_cna = expandcols(sample_columns,control_columns)

    new_samp_array,new_ctl_aaray,new_sna, new_cna = seperatesamples(sample_columns,control_columns)
    

    df, contaminants_df = deletemultizero(df,new_samp_array,new_ctl_aaray)

    if impute_method== 'impute-lowest5':
        df = impute_with_lowest5(df)
    elif impute_method== 'miss-forest':
        df = missforest_impute(df)
    else:
        df.fillna(missing_val, inplace = True)

    mediun_list = {}

    if (norm_method == 'Median'):

        for controls in control_columns:
            for replicates in controls:
                mediun_list[replicates] = df[replicates].median()

        for samples in sample_columns:

            for samp_replicates in samples:
                mediun_list[samp_replicates] = df[samp_replicates].median()

    elif (norm_method == 'Sum'):
        for controls in control_columns:
            for replicates in controls:
                mediun_list[replicates] = df[replicates].mean()

        for samples in sample_columns:
            for samp_replicates in samples:
                mediun_list[samp_replicates] = df[samp_replicates].mean()

    elif (norm_method == 'TMM'):
        prcount = int(tmmpr)/100

        for controls in control_columns:
            for replicates in controls:
                mediun_list[replicates] = stats.trim_mean(df[replicates], prcount)
        for samples in sample_columns:
            for samp_replicates in samples:
                mediun_list[samp_replicates] = stats.trim_mean(df[samp_replicates], prcount)

    elif (norm_method == "Quntail"):

        df_for_qunt = df[ exp_samp + exp_cont ]
        
        df_PCA_before = df_for_qunt

        quant_df = quantile_normalize(df_for_qunt)

        df_PCA_after = quant_df
        
        df = df.join(quant_df)
        sna, cna = get_sna_cna(sample_columns,control_columns)

        return df, df_PCA_before, df_PCA_after , cna , sna , contaminants_df

    elif(norm_method == "zscore"):
        df_for_zscore = df[ exp_samp + exp_cont ]

        df_PCA_before = df_for_zscore

        zscore_df = caluculate_zscore_norm(df_for_zscore)
        
        df_PCA_after = zscore_df
        df = df.join(zscore_df)
        sna , cna = get_sna_cna(sample_columns,control_columns)
       
        return df,df_PCA_before, df_PCA_after , cna, sna , contaminants_df
    
    elif(norm_method == "irs"):
        
        df, irs_sna, irs_cna = caluculate_irs_norm( df , new_samp_array , new_ctl_aaray )

        sna, cna = get_sna_cna(sample_columns,control_columns)

        df_PCA_before_irs = df[ exp_samp + exp_cont ]
        df_PCA_After_irs = df[ exp_sna + exp_cna ]
        

        return df ,df_PCA_before_irs, df_PCA_After_irs ,cna, sna , contaminants_df


    minn = min(mediun_list.values())
    multiplication_fact_list = {}
    for key,value in mediun_list.items():
        multiplication_fact_list[key] = (minn/value)

    cna = []
    for controls in control_columns:
        each_control = []
        for replicates in controls:
            df['NORM_'+replicates] = df[replicates] * multiplication_fact_list[replicates]
            each_control.append('NORM_'+replicates)
        cna.append(each_control)

    sna = []
    for samples in sample_columns:
        each_sample = []
        for samp_replicates in samples:
            df['NORM_'+samp_replicates] = df[samp_replicates] * multiplication_fact_list[samp_replicates]
            each_sample.append('NORM_'+samp_replicates)
        sna.append(each_sample)

    snaa, cnaa = get_sna_cna(sample_columns,control_columns)

    df_PCA_before = df[ exp_samp + exp_cont ]
    df_PCA_after = df[ exp_sna + exp_cna ]
    
    return df, df_PCA_before, df_PCA_after , cnaa, snaa , contaminants_df

def quantile_normalize(df):
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0),
                             index=df.index,
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    df_qn.rename(lambda x: "NORM_"+x, axis='columns', inplace = True)
    return df_qn                

def two_way_annova(x, all_cols, reps, s_and_c_col):

    mydict = dict()
    mydict['abundance'] = [ x[i] for i in all_cols ]
    mydict['sample'] = s_and_c_col
    mydict['reps'] = reps
    mydf = pd.DataFrame(mydict)
    model = ols('abundance ~ C(sample) + C(reps) + C(sample):C(reps)', data=mydf).fit()
    two_annovaa = sm.stats.anova_lm(model, typ=2)
    pvalue = two_annovaa.loc['C(sample):C(reps)']['PR(>F)']
    return pvalue

def pvalAndRatio(df, cna,sna, pvalue, key,fc_left,fc_right, lg2cut,both,pv_cutoff,snames, cnames, accessionKey, selected_control_key, job_id , adj_pval_method):

    average_normalized_sample_array = []
    average_normalized_control_array = []
    controlcols = []

    for index , controls in enumerate(cna):
        if index == int(selected_control_key):
            for control in controls:
                controlcols.append(control)


    foranova = cna + sna

    avrg_norm_array = []

    minuslog10array = list()
    log2fcarray = list()
    pvalue_array = list()
    # do this in normalization page
    df = df.loc[df[accessionKey] != 'sp'] 

    if pvalue == 'weltch':
        i = 0
        for samples in sna:
            df_sample = df[[y for y in samples]]
    
            average_normalized_sample_array.append("AVG_NORM"+snames[str(i)])
            df["AVG_NORM"+snames[str(i)]] = df_sample.mean(axis = 1)
            avrg_norm_array.append("AVG_NORM"+snames[str(i)])

            oneAnnova = False

            _,df["P VALUE of"+snames[str(i)]]= stats.ttest_ind(df[controlcols],df_sample,axis=1, equal_var = False)

            pvalue_array.append("P VALUE of"+snames[str(i)])
            df["Minus_log10_pvalue"+snames[str(i)]] = abs(np.log10(df["P VALUE of"+snames[str(i)]]))
            minuslog10array.append("Minus_log10_pvalue"+snames[str(i)])
            i += 1

    elif pvalue == 'ttest':

        i = 0
        for samples in sna:
            df_sample = df[[y for y in samples]]

            average_normalized_sample_array.append("AVG_NORM"+snames[str(i)])
            df["AVG_NORM"+snames[str(i)]] = df_sample.mean(axis = 1)
            avrg_norm_array.append("AVG_NORM"+snames[str(i)])

            oneAnnova = False

            _,df["P VALUE of"+snames[str(i)]]= stats.ttest_ind(df[controlcols],df_sample,axis=1, equal_var = True)

            pvalue_array.append("P VALUE of"+snames[str(i)])
            df["Minus_log10_pvalue"+snames[str(i)]] = abs(np.log10(df["P VALUE of"+snames[str(i)]]))
            minuslog10array.append("Minus_log10_pvalue"+snames[str(i)])
            i += 1
    # ==================================================================LIMMA=========================================
   
    elif pvalue == 'limma':
        oneAnnova = False
        rename_list = []
        log10array_limma = []
        log2fcarray_limma = []       
        
        matrix, fordf, pval_array, fc_array, rename_list = get_matrix_limma(sna,controlcols,snames)
        
        limma_df = limma_diff_API(df[fordf], matrix, job_id, rename_list)
        
        df = pd.concat([df.reset_index(drop=True),limma_df.reset_index(drop=True)], axis = 1)

        for k,v in snames.items():

            df["Minus_log10_pvalue"+v] = abs(np.log10(df[pval_array[int(k)]] ))
            log10array_limma.append("Minus_log10_pvalue"+v)

            df['LOG2 foldchange of'+v] = df[fc_array[int(k)]]
            log2fcarray_limma.append('LOG2 foldchange of'+v)

        df, dif_df_final, expression_array  = expression_calc(df,fc_array,log2fcarray_limma,log10array_limma,fc_left,fc_right,
            lg2cut, oneAnnova, pv_cutoff,key, both, accessionKey)
        
        forvolcano = list()
        i = 0
        for fc in log2fcarray_limma:
            volcano = []
            volcano.append(fc)
            volcano.append(log10array_limma[i])
            i+=1
            volcano.append(key)
            forvolcano.append(volcano)

        
        return df,forvolcano , dif_df_final, expression_array 
    
    #  =================================================limma ends here=================================================================================

    elif pvalue == '2anova':
        i = 0
        for samples in sna:
            df_sample = df[[y for y in samples]]

            average_normalized_sample_array.append("AVG_NORM"+snames[str(i)])
            df["AVG_NORM"+snames[str(i)]] = df_sample.mean(axis = 1)
            avrg_norm_array.append("AVG_NORM"+snames[str(i)])
            minuslog10array.append("Minus_log10_pvalue")

            i += 1

        oneAnnova = True
        sna_for_anova = list(chain.from_iterable(sna))
        cna_for_anova = list(chain.from_iterable(cna))

        all_cols_anova =  sna_for_anova + cna_for_anova

        reps_list = list()
        reps = list()

        # for 2 way annova
        s_and_c_col = list()
        for s in sna_for_anova:
            s_and_c_col.append('sample')
        for s in cna_for_anova:
            s_and_c_col.append('control')


        for r in range(0,len(sna[0])):
            reps_list.append( "rep"+str(r+1))

        for n in range (0,len(cna)+len(sna)):
            reps.append(reps_list)

        reps = list(chain.from_iterable(reps))

        df['p_value using TWO-way-ANOVA'] = df.apply(lambda x: two_way_annova(x, all_cols_anova, reps, s_and_c_col), axis=1)
        df["Minus_log10_pvalue"] = abs(np.log10(df["p_value using TWO-way-ANOVA"]))

    else:

        i = 0
        for samples in sna:
            df_sample = df[[y for y in samples]]
            average_normalized_sample_array.append("AVG_NORM"+snames[str(i)])
            df["AVG_NORM"+snames[str(i)]] = df_sample.mean(axis = 1)
            avrg_norm_array.append("AVG_NORM"+snames[str(i)])
            minuslog10array.append("Minus_log10_pvalue")


            i += 1
        oneAnnova = True

        _,df["P VALUE using One-Way-ANOVA"]= stats.f_oneway(*exapndd(foranova,df), axis = 1)
        df["Minus_log10_pvalue"] = abs(np.log10(df["P VALUE using One-Way-ANOVA"]))

    df["AVG_NORM_CONTROL"] = df[controlcols].mean(axis =1)

    #calculating foldchange
    foldchange_array = []
    for avg_sample in avrg_norm_array:
        sample_name  = avg_sample.replace("AVG_NORM" ,'')
        foldchange_array.append('FOLDCHANGE_'+ sample_name)
        df['FOLDCHANGE_'+ sample_name ] = df[avg_sample].div(df["AVG_NORM_CONTROL"])

    #caluclating log2 of foldchange
    for foldchange in foldchange_array:
        name = foldchange.replace("FOLDCHANGE_",'')
        df['LOG2 foldchange of'+ name ] = np.log2(df[foldchange])
        log2fcarray.append('LOG2 foldchange of'+ name )

    df.drop(avrg_norm_array, axis = 1,errors = 'ignore', inplace = True)
    df.drop('AVG_NORM_CONTROL', axis = 1,errors = 'ignore', inplace = True)



    df , dif_df_final, expression_array  = expression_calc(df,foldchange_array,log2fcarray,minuslog10array,fc_left,fc_right,
        lg2cut, oneAnnova, pv_cutoff,key, both, accessionKey)

    forvolcano = list()
    i = 0
    for fc in log2fcarray:
        volcano = []
        volcano.append(fc)
        volcano.append(minuslog10array[i])
        i+=1
        volcano.append(key)
        forvolcano.append(volcano)

    df = calculate_adjusted_pval(df,adj_pval_method)
    return df,forvolcano , dif_df_final, expression_array

def expression_calc(df,foldchange_array,log2fcarray,minuslog10array,fc_left,fc_right, lg2cut, oneAnnova, pv_cutoff,key, both, accessionKey):


    expressions = ("Upregulated", "Downregulated", "Not-sig")
    expression_array = list()
    pv_cutoff = -(np.log10(pv_cutoff))

    dif_df_final = pd.DataFrame()

    dif_df = pd.DataFrame()

    if both:

        for fc in log2fcarray:
            sample_name = fc.replace('LOG2 foldchange of','')
            sample_name = sample_name.strip()
            expr_samp = 'LOG2FC-Expression_'+sample_name

            if oneAnnova:
                pvalcol = "Minus_log10_pvalue"
            else:
                pvalcol = "Minus_log10_pvalue"+sample_name

            df.loc[(df[fc] >= lg2cut ) & (df[pvalcol] > pv_cutoff),expr_samp ] = expressions[0]
            df.loc[(df[fc] <= -lg2cut ) & (df[pvalcol] > pv_cutoff), expr_samp ] = expressions[1]

            dif_df = df.loc[(df[expr_samp] == 'Upregulated') | (df[expr_samp] == 'Downregulated'),[accessionKey,expr_samp,fc]]
            dif_df_final = pd.concat([dif_df_final,dif_df], axis = 1)

            df[expr_samp].fillna(expressions[2], inplace=True)
            expression_array.append(expr_samp)

    else:

        for fc in foldchange_array:
            sample_name = fc.replace('FOLDCHANGE_','')
            sample_name = sample_name.strip()

            expr_samp = 'FC-Expression_'+sample_name

            if oneAnnova:
                pvalcol = "Minus_log10_pvalue"
            else:
                pvalcol = "Minus_log10_pvalue"+sample_name

            df.loc[(df[fc] >= fc_right ) & (df[pvalcol] > pv_cutoff), expr_samp] = expressions[0]
            df.loc[(df[fc] <= fc_left ) & (df[pvalcol] > pv_cutoff), expr_samp] = expressions[1]

            dif_df = df.loc[(df[expr_samp] == 'Upregulated') | (df[expr_samp] == 'Downregulated'),[accessionKey,expr_samp,fc]]

            dif_df_final = pd.concat([dif_df_final,dif_df], axis = 1)

            df[expr_samp].fillna(expressions[2], inplace=True)
            expression_array.append(expr_samp)


    dif_df_final.set_index(accessionKey,inplace = True)
    dif_df_final.reset_index(inplace = True )
    dif_df_final[accessionKey] = dif_df_final[accessionKey].apply(cleanAcc)

    # incase if we do not convert gene symbol
    if (key != accessionKey):
        keydf = df[[key,accessionKey]]
        dif_df_final = dif_df_final.merge(keydf,on = accessionKey)

    dif_df_final.set_index(accessionKey,inplace = True)

    return df,dif_df_final,expression_array


def deletemultizero(df,sample_columns,control_columns):
    
    sc,cc = expandNCleanColumns(sample_columns,control_columns)
    all_cols = sc+cc
    df.replace(0, np.nan, inplace=True)
    contaminants = pd.DataFrame()
    no_of_repli = len(sample_columns[0])
    sample = sample_columns[0]
    nan_values = df[df[all_cols].isna().all(axis=1)]

    df = df.dropna(how='all', subset=all_cols)

    if no_of_repli == 3:
        indices_to_drop = list()
        for index, row in df[sample].iterrows():
            if ( pd.isnull(row[sample[0]])  and (  pd.isnull(row[sample[1]])  or  pd.isnull(row[sample[2]])  )  or ( pd.isnull(row[sample[1]])  and  pd.isnull(row[sample[2]]) )  ):
                indices_to_drop.append(index)
        
        contaminants =  df[df.index.isin(indices_to_drop)]

        df.drop(labels=indices_to_drop, inplace=True)
        
        df = df.dropna(how='all', subset=all_cols)
        
        contaminants = pd.concat([contaminants,nan_values])


    return df , contaminants


def exapndd(foranova,df):
    dflist = []
    for sam in foranova:
        dflist.append(df[sam])
    return dflist

