from django.shortcuts import render 
from django.http import HttpResponse
from django.http import Http404
from django.core.files.base import ContentFile

from django.conf import settings
from django.contrib import messages
from django.http import JsonResponse

import io
import os
import zipfile

import pandas as pd
import numpy as np
import matplotlib

from combat.pycombat import pycombat

from .import normaliz
from .import ibaq
from . import R_utils

from .utils import (abundances, savefile_csv , clean_coulumn_heading,intensities,lablesforbox, clean_custom_names,
    sort_name , columnsforbox, seperatesamples,cobmine_samp_control, calc_avg_exp, gene_ontology_calc, getAccesion,get_protien_interaction, clean_plot_df, get_gene_symbol,
     stringtolist, decode_svg_for_zip, getAccesionCol, getGeneCol, convert_acc_to_gene , save_plot_file, load_explot_file,  getbatches , z_convert_df)

from .myplots import (plot_pca,plot_heatmap,plot_volcano,box_plot_sns,plot_bar_grap, getcircbar,
    getElbowPlot,getkClusterMap,get_upset_plot, plot_heatmap_new, plot_volcano_new,get_bubble_plot,get_density_plot, get_histogram, get_pca_plot, plot_kmeans_new,
    get_raincloud_plot, get_violin_plot, get_scurve_plot, get_venn_plot,get_ma_plot,get_box_plot) 

from .keggpath import draw_pathway
from .models import DataAnalysis , Example, Contaminant , Ploting , Forzip, PlotExample


def home(request):
    matplotlib.rc_file_defaults()
    return render(request,'proteome/index.html')

def input(request):
    return render(request,'proteome/home.html')

def inputf(request): 

    lablled = True
    lablefree = False

    if (request.method == 'POST'):

        try:
            example_analysis = request.POST.get("wokring-with-ex")
            
            if example_analysis == "yes":
                
                q = Example.objects.filter(usethis = True).first()
                files = q.file
            else:
                files = request.FILES['file']

            
            
            if files.name.endswith('.xlsx') or files.name.endswith('.csv') or files.name.endswith('.txt') or files.name.endswith('.xls'):

                user="user"
                job_id , columns = savefile_csv(files , user ,lablled , lablefree, False)
                
                request.session['job_id'] = job_id

                accession_col = getAccesionCol(columns)
            
                gene_col = getGeneCol(columns)

                if (request.POST.get('rep_method')) == "techrep":

                    number_of_samples = int(request.POST.get('no_of_sample'))
                    number_of_control = int(request.POST.get('no_of_control'))
                    

                    abd_columns = abundances(columns)
                    
                    context = {'abd_columns':abd_columns,'number_of_samples':number_of_samples,
                    'number_of_control':number_of_control, 'accession_col':accession_col, 'gene_col':gene_col}
                    
                    return render(request,'proteome/pre_analyze.html',context)

                else:
                    number_of_batches = int(request.POST.get('no_of_batches'))
                    samples_in_bio = request.POST.get('samples_in_bio')
                    request.session['samples_in_bio'] = samples_in_bio

                    abd_columns = abundances(columns)
                    print(abd_columns)

                    context = {'abd_columns':abd_columns,'number_of_batches':number_of_batches,
                    'number_of_samples':samples_in_bio, 'accession_col':accession_col, 'gene_col':gene_col, 'example_analysis':example_analysis}

                    return render(request,'proteome/pre_anlz_bio.html',context)

            else:
                messages.error(request, 'Please upload only Excel , CSV or TXT file')
                return render(request, 'proteome/home.html')

        except:
            messages.error(request, 'ERROR:Something went wrong! Please check the input file')
            return render(request,'proteome/home.html')
        
    return render(request, 'proteome/home.html')



def analaze_cols(request):

    if (request.method == 'POST'):

        sample_data_columns = request.POST.get('final_sample_data')
        final_control_data = request.POST.get('final_control_data')

        final_sample_name = request.POST.get('final_sample_name')
        final_control_name = request.POST.get('final_control_name')

        final_sample_name = clean_custom_names(final_sample_name)
        final_control_name = clean_custom_names(final_control_name)

        job_id = request.session.get('job_id')
        
        missing_val_rep = request.POST.get('missing_val')
        norm_method = request.POST.get('norm_method')
        tmmpr = request.POST.get('tmmpr')

        impute_method = request.POST.get('impute')        

        sample_columns = clean_coulumn_heading(sample_data_columns)
        control_columns = clean_coulumn_heading(final_control_data)

        accession_col = request.POST.get('accession-col')
        convert = request.POST.get('convert_acc')
        gene_col = request.POST.get('gene-col')
        
        if gene_col != None:
            final_key = gene_col

        df,df_before_norm, df_after_norm, cna, sna , contaminants_df= normaliz.normaliz_data(job_id,sample_columns,
            control_columns,norm_method,missing_val_rep,tmmpr, impute_method)
        
        contaminant_and_res = contaminants_df.shape[0]
        norm_data_count = df.shape[0]
    
        
        if (gene_col == None or convert == "convert_acc") and accession_col != None:
            con_df , gs_convert_success = convert_acc_to_gene(df[accession_col].tolist())

            if gs_convert_success:
                df[accession_col] = df[accession_col].apply(lambda x: x.split(';')[0] if ';' in x else x)
                df = df.merge(con_df,left_on = accession_col , right_on='Accesion_gf' , how = 'left' , suffixes = (None , '_y'))
                df.drop('Accesion_gf', axis=1, inplace=True)
                final_key = "_GENE_SYMBOL_"
            else:
                final_key = accession_col

        if gene_col != None and accession_col == None:
            accession_col = gene_col



        request.session['cna'] = cna
        request.session['sna'] = sna
        request.session['final_sample_name'] = final_sample_name
        request.session['final_control_name'] = final_control_name


        request.session['accession_col'] = accession_col
        request.session['final_key'] = final_key

        new_df = df.to_csv(index = False)
        updated_file = ContentFile(new_df)
        updated_file.name = "normalized_data.csv"

        resedue_df = contaminants_df.to_csv(index = False)
        resedue_file = ContentFile(resedue_df)
        resedue_file.name = "resedue.csv"

        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = updated_file)
        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = resedue_file)

        result_q = DataAnalysis.objects.get(id =job_id)
        result_q.resultData = updated_file
        result_q.dleteted_resd = resedue_file

        result_q.save()


        pca_before = plot_pca(df_before_norm,sample_columns,control_columns, title = "PCA plot [Before normalization]", before = True,
            sname = final_sample_name, cname = final_control_name)
        pca_after = plot_pca(df_after_norm,sample_columns,control_columns,title = "PCA plot [After normalization]", before = False,
            sname = final_sample_name, cname = final_control_name)

        before_norm_box_log = box_plot_sns(df = df_before_norm, samp_col = sample_columns ,control_col = control_columns,
         title = "Box plot [Before normalization]", isbio = False)
        after_norm_box_log = box_plot_sns(df = df_after_norm, samp_col = sample_columns ,control_col = control_columns,
         title = "Box plot [After normalization]", isbio = False)

        b4_norm_image = decode_svg_for_zip(pca_before , name =  str(job_id)+'before normalization.svg')
        after_norm_image = decode_svg_for_zip(pca_after , name = str(job_id)+'After normalization.svg')

        b4_box_image = decode_svg_for_zip(before_norm_box_log , name =  str(job_id)+'box plot before normalization.svg')
        after_box_image = decode_svg_for_zip(after_norm_box_log , name = str(job_id)+'box plot After normalization.svg')

        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = b4_norm_image)
        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = after_norm_image)
        
        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = b4_box_image)
        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = after_box_image)
        
        forzip_qs.save()

        
        no_of_control = len(final_control_name)
        control_to_select =dict()
        if no_of_control > 1:
            control_to_select = final_control_name



        context = { 'pca_before':pca_before,'pca_after':pca_after,'before_norm_box_log':before_norm_box_log,
            'after_norm_box_log':after_norm_box_log, 'control_to_select':control_to_select, 'norm_data_count':norm_data_count,'contaminant_and_res':contaminant_and_res}

        return render(request, 'proteome/normalized.html', context)
    return render(request, 'proteome/home.html')


def downloadfile(request):
    job_id = request.session.get('job_id')
    q = DataAnalysis.objects.get(id = job_id)
    file = q.resultData.path

    download_path = os.path.join(settings.MEDIA_ROOT, file)

    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404


def download_ddf(request):
    job_id = request.session.get('job_id')
    q = DataAnalysis.objects.get(id = job_id)
    file = q.expression_data.path

    download_path = os.path.join(settings.MEDIA_ROOT, file)

    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404



def analaze_cols_bio(request):
    
    pd.options.mode.chained_assignment = None

    if (request.method == 'POST'):
        sample_data_columns = request.POST.get('final_sample_data')
        final_control_data = request.POST.get('final_control_data')

        job_id = request.session.get('job_id')

        missing_val_rep = request.POST.get('missing_val')
        norm_method = request.POST.get('norm_method')
        tmmpr = request.POST.get('tmmpr')

        impute_method = request.POST.get('impute')
        
        final_sample_name = request.POST.get('final_sample_name')
        final_control_name = request.POST.get('final_control_name')

        final_sample_name = clean_custom_names(final_sample_name)
        final_control_name = clean_custom_names(final_control_name)

        samples_in_bio = request.session.get('samples_in_bio')

        sample_columns = clean_coulumn_heading(sample_data_columns)
        control_columns = clean_coulumn_heading(final_control_data)

        # must have
        accession_col = request.POST.get('accession-col')
        convert = request.POST.get('convert_acc')
        gene_col = request.POST.get('gene-col')
        
        if gene_col != None:
            final_key = gene_col
        
        df, df_PCA_before, df_PCA_after, cna, sna , contaminants_df = normaliz.normaliz_data_bio(job_id,
        sample_columns,control_columns,norm_method,missing_val_rep,tmmpr , impute_method)
        
        new_samp_array,new_ctl_aaray,new_sna, new_cna = seperatesamples(sample_columns,control_columns)

        contaminant_and_res = contaminants_df.shape[0]
        norm_data_count = df.shape[0]
        if (gene_col == None or convert == "convert_acc") and accession_col != None:
            con_df , gs_convert_success = convert_acc_to_gene(df[accession_col].tolist())

            if gs_convert_success:
                df[accession_col] = df[accession_col].apply(lambda x: x.split(';')[0] if ';' in x else x)
                df = df.merge(con_df,left_on = accession_col , right_on='Accesion_gf' , how = 'left' , suffixes = (None , '_y'))
                df.drop('Accesion_gf', axis=1, inplace=True)
                final_key = "_GENE_SYMBOL_"
            else:
                final_key = accession_col

        if gene_col != None and accession_col == None:
            accession_col = gene_col

        request.session['cna'] = new_cna
        request.session['sna'] = new_sna

        request.session['cna_for_batch'] = cna
        request.session['sna_for_batch'] = sna

        request.session['final_sample_name'] = final_sample_name
        request.session['final_control_name'] = final_control_name

        request.session['accession_col'] = accession_col
        request.session['final_key'] = final_key


        new_df = df.to_csv(index = False)
        updated_file = ContentFile(new_df)
        updated_file.name = "normalized_data.csv"

        resedue_df = contaminants_df.to_csv(index = False)
        resedue_file = ContentFile(resedue_df)
        resedue_file.name = "resedue.csv"

        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = updated_file)
        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = resedue_file)

        result_q = DataAnalysis.objects.get(id =job_id)
        result_q.resultData = updated_file
        result_q.dleteted_resd = resedue_file
        result_q.save()

        pca_before = plot_pca(df_PCA_before,new_samp_array,new_ctl_aaray, title = "PCA plot [Before Normalization]", before = True,
            sname = final_sample_name, cname = final_control_name)
        pca_after = plot_pca(df_PCA_after,new_samp_array,new_ctl_aaray,title = "PCA plot [After Normalization]", before = False,
            sname = final_sample_name, cname = final_control_name)


        before_batch_box = box_plot_sns(df = df_PCA_before, samp_col = sample_columns ,control_col = control_columns,
         title = "Box plot [Before normalization]", isbio = True)
        after_batch_box = box_plot_sns(df = df_PCA_after, samp_col = sna ,control_col = cna,
         title = "Box plot [After normalization]", isbio = True)
        
        b4_norm_image = decode_svg_for_zip(pca_before , name =  str(job_id)+'before normalization.svg')
        after_norm_image = decode_svg_for_zip(pca_after , name = str(job_id)+'After normalization.svg')

        
        b4_box_image = decode_svg_for_zip(before_batch_box , name =  str(job_id)+'Before Batch correction.svg')
        after_box_image = decode_svg_for_zip(after_batch_box , name = str(job_id)+'After Batch correction.svg')

        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = b4_norm_image)
        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = after_norm_image)
        
        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = b4_box_image)
        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = after_box_image)
        
        forzip_qs.save()

        bio_data = True



        no_of_control = len(final_control_name)
        control_to_select =dict()
        if no_of_control > 1:
            control_to_select = final_control_name


        context = {'pca_before':pca_before,
            'pca_after':pca_after,'before_batch_box':before_batch_box,
             'after_batch_box':after_batch_box , 'bio_data': bio_data, 'control_to_select':control_to_select,'norm_data_count':norm_data_count, 'contaminant_and_res':contaminant_and_res}

        return render(request, 'proteome/normalized.html', context)
    return render(request, 'proteome/home.html')


def pvalues(request):
    
    pd.options.mode.chained_assignment = None

    if (request.method == 'POST'):
        job_id = request.session.get('job_id')

        pvalue = request.POST.get('pvalue')
        fc_left = float(request.POST.get('fc_left'))
        fc_right = float(request.POST.get('fc_right'))
        pv_cutoff = float(request.POST.get('pv_cutoff'))
        lg2cut = float(request.POST.get('lg2cut'))
        ratio_or_lg = request.POST.get('ratiolg2')
        adj_pval_method = request.POST.get('adj-pval-method')

        if ratio_or_lg == None:
            both = request.session.get('both')
            if both:
                ratio_or_lg = 'lg2'
            else:
                ratio_or_lg = 'ratio'

        if pvalue == None:
            pvalue = request.session.get('p_value_type')

        selected_control_key = request.POST.get('selected-control')
        # if there is only one control
        if selected_control_key == None:
            selected_control_key = '0'

        both = False
        
        if ratio_or_lg == 'lg2':
            both = True

        
        data = DataAnalysis.objects.get(id = job_id)
        datafile = data.resultData.path
        df = pd.read_csv(datafile)

        isdif_ex = request.POST.get('difex-flow') 
        
        if isdif_ex == 'difex':
            accessionKey = request.POST.get('accession-col')
            final_key = accessionKey
            sample_data_columns = request.POST.get('final_sample_data')
            final_control_data = request.POST.get('final_control_data')

            sna = clean_coulumn_heading(sample_data_columns)
            cna = clean_coulumn_heading(final_control_data)

            
            final_sample_name = request.POST.get('final_sample_name')
            final_control_name = request.POST.get('final_control_name')

            snames = clean_custom_names(final_sample_name)
            cnames = clean_custom_names(final_control_name)

            snames = {str(k):v for k,v in snames.items()}
            cnames = {str(k):v for k,v in cnames.items()}

        else:

            accessionKey = request.session.get('accession_col')
            final_key = request.session.get('final_key')
            cna = request.session.get('cna')
            sna = request.session.get('sna')
            snames = request.session.get('final_sample_name')
            cnames = request.session.get('final_control_name')

        no_of_control = len(cnames)

        control_to_select =dict()

        if no_of_control > 1:
            control_to_select = cnames

        request.session['fc_left'] = fc_left
        request.session['fc_right'] = fc_right
        request.session['lg2cut'] = lg2cut
        request.session['both'] = both
        request.session['p_value_type'] = pvalue    

        expression_present = request.session.get('expression_array')

        if expression_present:
            df.drop(expression_present, axis=1, inplace=True, errors = 'ignore')

        df,forvolcano , dif_df_final, expression_array  = normaliz.pvalAndRatio(df,cna,sna, pvalue, final_key,
              fc_left,fc_right, lg2cut, both, pv_cutoff, snames, cnames, accessionKey, selected_control_key, job_id , adj_pval_method)


        # do something here
        # df.replace([np.inf, -np.inf], 0 , inplace = True)
        # dif_df_final.replace([np.inf, -np.inf], np.nan , inplace = True)
        # forheatmap.replace([np.inf, -np.inf], np.nan , inplace = True)


        df_size = df.shape
        total_protiens = df_size[0]

        request.session['expression_array'] = expression_array
  
        total_up = dif_df_final.shape[0]
        
        volcanoplotlist = dict()
        
        zip_vlc_name = []

        for volcanocols in forvolcano:
            volcanodf = df[volcanocols]
            volcanodf[final_key] = df[final_key]

            sample_name = volcanocols[0].replace('LOG2 foldchange of','').strip()
            
            get_title = sample_name+"_vs._"+cnames[selected_control_key]
            
            zip_vlc_name.append(get_title)
            
            get_volacno , difex_genes = plot_volcano( volcanodf,lfc = volcanocols[0], pv = volcanocols[1], genes = final_key,
                lg2cut = lg2cut ,pvalue_cut_off = pv_cutoff ,title = get_title , both = both)

            volcanoplotlist[get_title] = [ get_volacno , difex_genes ]

        for k,v in volcanoplotlist.items():
            print(v[1])



        df_for_bar = dif_df_final[expression_array]

        # uplist_bar , dowlist_bar , samples_bar = plot_bar_grap(df_for_bar)

        bar_plot = plot_bar_grap(df_for_bar)
        
        new_df = df.to_csv(index = False)
        updated_file = ContentFile(new_df)
        updated_file.name = "finalresult.csv"

        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = updated_file)

        diff_express_df = dif_df_final.to_csv()
        dif_df = ContentFile(diff_express_df)
        dif_df.name = "diffrential expressed.csv"

        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = dif_df)

        result_q = DataAnalysis.objects.get(id =job_id)
        result_q.resultData = updated_file
        
        result_q.expression_data = dif_df
        
        result_q.save()

        request.session['key'] = final_key
        
        i = 0
        for vplt in volcanoplotlist:
            vplt_image = decode_svg_for_zip(volcanoplotlist[vplt][0] , name = zip_vlc_name[i] +'.svg')
            forzip_qs = Forzip.objects.create(analysis_id = job_id, files = vplt_image)
            i += 1

        # result_image = decode_svg_for_zip( bar_graph , name =  str(job_id)+'result_stat.svg')
        # forzip_qs = Forzip.objects.create(analysis_id = job_id, files = result_image)


        forzip_qs.save()

        context = {'volcanoplotlist': volcanoplotlist, 'total_up': total_up,
                    # 'uplist_bar': uplist_bar, 'dowlist_bar':dowlist_bar,'samples_bar':samples_bar,
                    'total_protiens':total_protiens, 'bar_plot':bar_plot, 
        'fc_left':fc_left, 'fc_right':fc_right, 'pv_cutoff':pv_cutoff, 'lg2cut':lg2cut, 'ratio_or_lg':ratio_or_lg,
        'key':final_key, 'accessionKey':accessionKey, 'control_to_select':control_to_select }

        return render(request, 'proteome/pvalandratio.html', context)


def exampledown(request):
    q = Example.objects.get(usethis = True)
    file = q.file.path
    download_path = os.path.join(settings.MEDIA_ROOT, file)

    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404


def exampledownbio(request):
    q = Example.objects.filter(biological_rep = True).first()
    file = q.file.path
    download_path = os.path.join(settings.MEDIA_ROOT, file)

    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404



def download_geneontology(request):

    job_id = request.session.get('job_id')
    q = DataAnalysis.objects.get(id = job_id)
    file = q.geneOntology.path

    download_path = os.path.join(settings.MEDIA_ROOT, file)

    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404


def perform_btch_cr(request):
    job_id = request.session.get('job_id')
    sna = request.session.get('sna_for_batch')
    cna = request.session.get('cna_for_batch')

    btch_method = request.POST.get('batchMethod')

    q = DataAnalysis.objects.get(id = job_id)
    
    file = q.resultData.path
    df = pd.read_csv(file)

    batch_list, btc_columns = getbatches(sna, cna)
    
    bf_btch_correct = df[btc_columns]

    if btch_method == "combat":
        
        aftr_bth_df = pycombat(bf_btch_correct,batch_list)

    elif btch_method == "limma":
        aftr_bth_df = R_utils.batch_correct_limma(bf_btch_correct,batch_list, job_id)
  
    df.drop(btc_columns, axis=1, inplace=True)
    df = pd.concat([df,aftr_bth_df], axis = 1 )
    new_df = df.to_csv(index = False)
    updated_file = ContentFile(new_df)
    updated_file.name = "normalized_data.csv"

    forzip_qs = Forzip.objects.create(analysis_id = job_id, files = updated_file)
    result_q = DataAnalysis.objects.get(id =job_id)
    result_q.resultData = updated_file
    result_q.save()

    before_batch_box = box_plot_sns(df = bf_btch_correct, samp_col = sna ,control_col = cna,
        title = "Before Batch Correction", isbio = True)
    after_batch_box = box_plot_sns(df = aftr_bth_df, samp_col = sna ,control_col = cna,
        title = "After Batch Correction", isbio = True)
    
    b4_batch_image = decode_svg_for_zip(before_batch_box , name =  str(job_id)+'Before Batch Correction.svg')
    after_batch_image = decode_svg_for_zip(after_batch_box , name = str(job_id)+'After Batch Correction.svg')

    forzip_qs = Forzip.objects.create(analysis_id = job_id, files = b4_batch_image)
    forzip_qs = Forzip.objects.create(analysis_id = job_id, files = after_batch_image)

    forzip_qs.save()

    data = dict()
    data['before_batch_box'] = before_batch_box
    data['after_batch_box'] = after_batch_box

    return JsonResponse(data)


def inputf_dexp(request):

    if (request.method == 'POST'):
        isexample = request.POST.get("wokring-with-ex")

        if isexample == 'yes':
            q = Example.objects.filter(name = "differential example").first()
            files = q.file
            number_of_samples = q.number_of_sample
            number_of_control = 1

        else:        
            files = request.FILES['file']
            number_of_samples = int(request.POST.get('no_of_sample'))
            number_of_control = int(request.POST.get('no_of_control'))
    
        user="user"

        job_id , all_column = savefile_csv(files , user ,lablled = False, lablefree = False , diffex = True)

        request.session['job_id'] = job_id
        accession_col = getAccesionCol(all_column)
        gene_col = getGeneCol(all_column)
        print(all_column)

        context = {'all_column':all_column,'number_of_samples':number_of_samples,
        'number_of_control':number_of_control,'accession_col':accession_col, 'gene_col':gene_col, 'isexample':isexample}
        return render(request,'proteome/preanalyze_dexp.html',context)

    else:
        return render(request,'proteome/input_dexp.html')


def kmeans(request):
    data = {}
    if request.method == 'POST':
        job_id = request.session.get('job_id')
        key = request.session.get("key")

        q = DataAnalysis.objects.get(id = job_id)

        datafile = q.expression_data.path

        df = pd.read_csv(datafile)

        elbowmap = getElbowPlot(df,key)

        data = dict()
        data['kmean'] = elbowmap
        return JsonResponse(data)


def kClusterPlot(request):

    if request.method == 'POST':

        job_id = request.session.get('job_id')
        key = request.session.get("key")
        clusters = int(request.POST.get("noOfCluster"))

        q = DataAnalysis.objects.get(id = job_id)

        datafile = q.expression_data.path

        df = pd.read_csv(datafile)
        expression_array = request.session.get('expression_array')
        df.drop( expression_array, axis=1, inplace = True )
        
        fc_left = request.session.get('fc_left')
        fc_right = request.session.get('fc_right')
        lg2cut = request.session.get('lg2cut')
        both = request.session.get('both')

        request.session['cluster_key'] = key
        kmap, cluster_df = getkClusterMap(df,key,clusters,fc_left,fc_right,lg2cut,both)

        go_for_cluster = cluster_df.loc[cluster_df['clusters'] == 0 ]

        go_for_cluster.reset_index(inplace  = True)
        genes_for_go = go_for_cluster[key].to_list()

        cluster_df = cluster_df.to_csv()
        updated_file = ContentFile(cluster_df)
        updated_file.name = "clustered_data.csv"
        q.clusterdb = updated_file
        q.save()

        total_clusters = []

        for n in range(0,clusters):
            total_clusters.append(n+1)

        data = dict()
        data['total_clusters'] = total_clusters
        data['kmap'] = kmap
        # data['circbar'] = circbar
        data['cluster'] = 1

        return JsonResponse(data)


def heirarchial_heatmap(request):
    if request.method == 'POST':

        job_id = request.session.get('job_id')
        key = request.session.get("key")
        fc_left = request.session.get('fc_left')
        fc_right = request.session.get('fc_right')
        lg2cut = request.session.get('lg2cut')
        both = request.session.get('both')
        z_convert = request.POST.get('zcoreConvert')
        
        print(z_convert)
        print(type(z_convert))
        q = DataAnalysis.objects.get(id = job_id)
        datafile = q.expression_data.path
        df = pd.read_csv(datafile)
        
        expression_array = request.session.get('expression_array')
        df.drop(expression_array, axis=1, inplace = True )

        df.set_index(key, inplace = True)
        df = df.select_dtypes(exclude=['object'])

        if z_convert == "true":
            print("converting")
            df = z_convert_df(df)
            fc_left = -2
            fc_right = 2
            lg2cut = 2

        get_heatmap = plot_heatmap(df,both,fc_left,fc_right,lg2cut , z_convert)
        heatmap_image = decode_svg_for_zip( get_heatmap , name =  str(job_id)+'heatmap.svg')
        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = heatmap_image)
        forzip_qs.save()

        data = dict()
        data['get_heatmap'] = get_heatmap
        return JsonResponse(data)


def goForClusters(request):

    if request.method == 'POST':
        job_id = request.session.get('job_id')
        
        key = request.session.get("key")

        selectedCluster = int(request.POST.get("clusterSelected"))

        q = DataAnalysis.objects.get(id = job_id)
        cluster = q.clusterdb.path
        cluster_df = pd.read_csv(cluster)

        go_for_cluster = cluster_df.loc[cluster_df['clusters'] == selectedCluster-1 ]
        go_for_cluster.reset_index(inplace  = True)

        genes_for_go = go_for_cluster[key].to_list()
        gene_ont_df = gene_ontology_calc(genes_for_go, 'Homo sapiens')

 
        gene_ont_df.rename(columns = {"intersection_size": "value","source":"group"},inplace = True)

        circbar = getcircbar(gene_ont_df)
        data = dict()
        data['circbar'] = circbar
        return JsonResponse(data)


def string_ajax(request):
    if request.method == 'POST':
        job_id = request.session.get('job_id')
        score = int(request.POST.get("scoreString"))

        for_string = request.session.get('for_string')
        enrichment_df = get_protien_interaction(for_string)
        
        q = DataAnalysis.objects.get(id = job_id)

        new_df = enrichment_df.to_csv(index = False)
        updated_file = ContentFile(new_df)
        updated_file.name = "stringdb.csv"
        q.stringdb = updated_file
        q.save()
        data = dict()
        data['score'] = score
        data['for_string'] = for_string
        return JsonResponse(data)


def string_db_download(request):

    job_id = request.session.get("job_id")

    q = DataAnalysis.objects.get(id = job_id)
    file = q.stringdb.path

    download_path = os.path.join(settings.MEDIA_ROOT, file)

    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404



def download_cluster_df(request):
    
    job_id = request.session.get('job_id')
    q = DataAnalysis.objects.get(id = job_id)
    file = q.clusterdb.path

    download_path = os.path.join(settings.MEDIA_ROOT, file)

    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404


def contaminant_down(request):
    q = Contaminant.objects.get(name = "cont")
    file = q.file.path
    download_path = os.path.join(settings.MEDIA_ROOT, file)
    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404


def download_resedue(request):

    job_id = request.session.get('job_id')
    q = DataAnalysis.objects.get(id = job_id)
    file = q.dleteted_resd.path
    download_path = os.path.join(settings.MEDIA_ROOT, file)
    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404


def gene_ontology_ajax(request):

    data = {}

    if request.method == 'POST':

        job_id = request.session.get('job_id')
        
        key = request.session.get("key")
        species = request.POST.get("species")

        pvcutOff = float(request.POST.get("pvcutOff"))

        q = DataAnalysis.objects.get(id = job_id)

        datafile = q.expression_data.path

        df = pd.read_csv(datafile)

        accession = df[key].tolist()

        request.session['for_string'] = accession

        gene_ont_df = gene_ontology_calc(accession, species , pvcutOff)

        new_df = gene_ont_df.to_csv(index = False)
        updated_file = ContentFile(new_df)
        updated_file.name = "geneOntology.csv"
        q.geneOntology = updated_file
        q.save()

        pathway_list = gene_ont_df[['native','name']].loc[gene_ont_df['source'] == 'KEGG']
        
        pathway_list = pathway_list.loc[pathway_list['native'] != 'KEGG:00000' ]
        
        pathway_list.set_index('native', inplace = True )

        pathway_list = pathway_list.to_json() 

        gene_ont_df = gene_ont_df[['intersection_size','source','native','p_value']]

        gene_ont_df.rename(columns = {"intersection_size": "value","source":"group",'native':'name'},inplace = True)

        filt = gene_ont_df['p_value'] < pvcutOff

        gene_ont_df = gene_ont_df.loc[filt]

        circbar = getcircbar(gene_ont_df)

        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = updated_file)
        circbar_decoded = decode_svg_for_zip(circbar, name = 'gene ontology.svg')
        forzip_qs = Forzip.objects.create(analysis_id = job_id, files = circbar_decoded)
        forzip_qs.save()

        data = dict()
        data['circbar'] = circbar
        data['job_id'] = job_id
        data['pathway_list'] = pathway_list
        return JsonResponse(data)


def kegg_pathway(request):
    pathway_color = ('red', 'green')
    data = {}

    if request.method == 'POST':

        job_id = request.session.get('job_id')  
        key = request.session.get("key")
        req_pathway = request.POST.get("selectedPathway")

        expression_array = request.session.get("expression_array")

        expression_array.append(key)

        q = DataAnalysis.objects.get(id = job_id)
        expression_df = q.expression_data.path
        expression_df = pd.read_csv(expression_df)

        go_df = q.geneOntology.path
        go_df = pd.read_csv(go_df, usecols = ['source', 'native', 'name', 'intersections'])

        go_df = go_df.loc[go_df['native'] == req_pathway]

        expression_df = expression_df[expression_array]
        expression_df.set_index(key , inplace = True )
        expression_df['avg_expression'] = expression_df.apply(calc_avg_exp, axis = 1)
        expression_df = expression_df[['avg_expression']]

        expression_df.loc[(expression_df['avg_expression'] == "Upregulated") , 'color'] = pathway_color[0]
        expression_df.loc[(expression_df['avg_expression'] == "Downregulated") , 'color'] = pathway_color[1]

        go_df['intersections'] = go_df['intersections'].apply(stringtolist)
        go_df = go_df.explode('intersections')

        pathway_df = go_df.merge(expression_df, left_on='intersections' , right_on=key)

        color_dict = dict(zip(pathway_df["intersections"], pathway_df["color"]))
        
        gene_color_df = pathway_df[["intersections","color"]]
        gene_color_df.set_index('intersections',inplace = True)
        gene_color = gene_color_df.to_json()

        pathway_image = draw_pathway(req_pathway , color_dict , request)

        data = dict()
        data['pathway_image'] = pathway_image
        data['gene_color'] = gene_color

        return JsonResponse(data)


def input_label_free(request):
    return render(request,'proteome/home_labelfree.html')


def download_ibaq(request):
    job_id = request.session.get('job_id')
    q = DataAnalysis.objects.get(id = job_id)
    file = q.ibaq.path
    download_path = os.path.join(settings.MEDIA_ROOT, file)
    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404


def inputfile_lblfree(request):
    lablefree = True
    if (request.method == 'POST'):
        try:
            exData = request.POST.get("wokring-with-ex")
            if exData == 'yes':
                q = Example.objects.filter(lableFreeData = True).first()
                fastafile = q.fastafile
                files = q.file
            else:
                files = request.FILES['file']

                fastafile = request.FILES['fasta-file']

            prot_ident = request.POST.get("prot_ident")
            enzyme = request.POST.get("enzyme")
            missed_clevage = int(request.POST.get("missed_clevage"))
            pep_min_len = int(request.POST.get("pep_min_len"))
            pep_max_len = int(request.POST.get("pep_max_len"))
            fastsdb = request.POST.get("fastsdb")

            if files.name.endswith('.xlsx') or files.name.endswith('.csv') or files.name.endswith('.txt') or files.name.endswith('.xls'):
                user="user"

                data_als = DataAnalysis.objects.create(file = files, user = user, labledData = False,lableFreeData = lablefree , fastafile = fastafile)
                
                job_id = data_als.id
                request.session['job_id'] = job_id
                data_als.save()

                df = pd.DataFrame()
                if files.name.endswith('.xlsx'):
                    df = pd.read_excel(files,engine='openpyxl')
                elif files.name.endswith('.xls'):
                    df = pd.read_excel(files,engine='xlrd')
                elif files.name.endswith('.csv'):
                    df = pd.read_csv(data_als.file.path)
                else:
                    data = DataAnalysis.objects.get(id = job_id)
                    df = pd.read_csv(data_als.file.path, delimiter = '\t')
                
                df = ibaq.protien_identify(df,prot_ident, missed_clevage,pep_min_len,pep_max_len, enzyme , job_id, fastsdb)
                new_df = df.to_csv()
                updated_file = ContentFile(new_df)
                updated_file.name = "ibaqresult.csv"
                data_als.ibaq = updated_file
                data_als.save()
                ibaq_head = df.head(10).to_html()
            
            else:
                messages.error(request, 'ERROR:Please check the input file format')
                return render(request,'proteome/home_labelfree.html')

            context = { 'ibaq_head':ibaq_head,}
            return render(request,'proteome/ibaqdone.html', context)
            
        except:
            messages.error(request, 'ERROR:Please check the input file format, input columns must contain Annotated Sequence, Gene symbol, Intensity')
            return render(request,'proteome/home_labelfree.html')
        
    return render(request,'proteome/home_labelfree.html')


def download_zip(request):
    job_id = request.session.get('job_id')

    query = Forzip.objects.filter(analysis_id = job_id)
    filelist = []
    for q in query:
        download_path = os.path.join(settings.MEDIA_ROOT, q.files.path)
        filelist.append(download_path)

    byte_data = io.BytesIO()
    zip_file = zipfile.ZipFile(byte_data, "w")

    for file in filelist:
        filename = os.path.basename(os.path.normpath(file))
        zip_file.write(file, filename)
    zip_file.close()

    response = HttpResponse(byte_data.getvalue(), content_type='application/zip')
    response['Content-Disposition'] = 'attachment; filename=Proteoark_results.zip'
    return response



def getplot(request):
    plot_type = request.POST.get('plot-radio')
    user="user"

    if len(request.FILES) != 0:
        files = request.FILES['file']
        
        if files.name.endswith(('.xlsx','.xls','.csv','.txt','.tsv')):
            plot_id , columns_temp  , column_dtype = save_plot_file(files, plot_type , user)
            request.session['plot_id'] = plot_id
        
        else:
            raise Http404

    else:
        plot_id , columns_temp , column_dtype = load_explot_file(plot_type , user)
        request.session['plot_id'] = plot_id

    print(column_dtype)

    if plot_type == "circbar":

        context = {'columns': column_dtype}
        return render(request,'proteome/plot_circbar.html', context)

    elif plot_type == "upset":

        context = {'columns': column_dtype}
        return render(request,'proteome/plot_upset_plot.html', context)


        # df=getdf(plot_id)
        # upset_plot = get_upset_plot(df)
        # context = {'upset_plot': upset_plot}
        # return render(request,'proteome/plot_upset_plot.html', context)

    elif plot_type == "histogram":
        context = {'columns': column_dtype}
        return render(request,'proteome/plot_histogram.html', context)
    
    elif plot_type == "density":
        context = {'columns': column_dtype}
        return render(request,'proteome/plot_density_plot.html', context)

    elif plot_type == "box":
        context = {'columns': column_dtype}
        return render(request,'proteome/plot_box_plot.html', context)

    elif plot_type == "heatmap":
        context = {'columns': column_dtype}
        return render(request,'proteome/plot_heatmap_plot.html', context)

    elif plot_type == "volcano":
        context = {'columns': column_dtype}
        return render(request,'proteome/plot_volcano.html', context)

    
    elif plot_type == "bubble":
        context = {'columns': column_dtype}
        return render(request,'proteome/plot_bubble.html', context)
    
    
    elif plot_type == "scurve":
        context = {'columns': column_dtype}
        return render(request,'proteome/plot_Scurve_plot.html', context)

    
    elif plot_type == "pca":
        context = {'columns': column_dtype}
        return render(request,'proteome/plot_pca.html', context)        

    elif plot_type == "rain":

        context = {'columns': column_dtype}
        return render(request,'proteome/plot_rain_cloud.html', context)  

    elif plot_type == "maplot":

        context = {'columns': column_dtype}
        return render(request,'proteome/plot_ma_plot.html', context)  

    
    elif plot_type == "venn":

        context = {'columns': column_dtype}
        return render(request,'proteome/plot_venn_plot.html', context)  

    elif plot_type == "violin":

        context = {'columns': column_dtype}
        return render(request,'proteome/plot_violin.html', context) 
    
    
    elif plot_type == "kmean":

        context = {'columns': column_dtype}
        return render(request,'proteome/plot_kmean_plot.html', context) 


def getdf(id):

    q = Ploting.objects.get(id = id)
    data = q.file.path

    df = pd.read_csv(data)

    return df


def plot_circbar(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')
        
        value = request.POST.get("values")
        termName = request.POST.get("termName")
        group = request.POST.get("group")

        q = Ploting.objects.get(id = plot_id)
        data = q.file.path
        df = pd.read_csv(data)
        df.rename(columns = {value: "value",group:"group", termName:"name"},inplace = True)
        circbar = getcircbar(df)
        data = dict()
        data['circbar'] = circbar
        return JsonResponse(data)


def plot_heatmap_plot(request):

    if request.method == 'POST':
        plot_width = request.POST.get("plotWidth")
        plot_height = request.POST.get("plotHeight")
        plot_id = request.session.get('plot_id')
        cutoff = request.POST.get("cutOff")
        index = request.POST.get("index")
        left = None
        right = None
        if cutoff == 'true':
            left = float(request.POST.get("left"))
            right = float(request.POST.get("right"))
        
        

        df = getdf(plot_id)
        heatmap ,figsize= plot_heatmap_new(df,index, left, right, cutoff, plot_width , plot_height)

        data = dict()
        data['heatmap'] = heatmap
        data['figsize'] = figsize
        response_data = {
        "plotWidth": plot_width,
        "plotHeight": plot_height,
        "figsize": figsize,
        "heatmap": heatmap

    }
    return JsonResponse(response_data)


def plot_bubble_plot(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')

        xAxis = request.POST.get("xAxis").strip()
        yAxis = request.POST.get("yAxis").strip()
        sizeCol = request.POST.get("sizeCol").strip()
        category = request.POST.get("category").strip()

        df = getdf(plot_id)
        newcols = [x.strip() for x in list(df.columns)]
        df.columns = newcols
        bubble = get_bubble_plot(df,xAxis, yAxis,sizeCol, category)

        data = dict()
        data['bubble'] = bubble
        return JsonResponse(data)


def plot_histogram_plot(request):
    print(request.user)
    if request.method == 'POST':
        plot_id = request.session.get('plot_id')
        colors = request.POST.get('colorValues')
        print(colors)

        fieldsHistogram = request.POST.get("fields-histogram")
        
        fieldsHistogram = fieldsHistogram.split(',')
        
        df = getdf(plot_id)
        histogram = get_histogram(df,fieldsHistogram)

        data = dict()
        data['histogram'] = histogram

        return JsonResponse(data)


def plot_ma_plot(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')

        fcCutoff = request.POST.get("fcCutoff").strip()
        sampleColumn = request.POST.get("sampleColumn").strip()
        controlColumn = request.POST.get("controlColumn").strip()
        log2fcColumn = request.POST.get("log2fcColumn").strip()

        fcCutoff=float(fcCutoff)

        df = getdf(plot_id)
        newcols = [x.strip() for x in list(df.columns)]
        df.columns = newcols
        ma= get_ma_plot(df,fcCutoff,sampleColumn,controlColumn,log2fcColumn)

        data = dict()
        data['ma'] = ma
        return JsonResponse(data)


def plot_pca_plot(request):
    
    if request.method == 'GET':
        plot_id = request.session.get('plot_id')
        df = getdf(plot_id)
        pca = get_pca_plot(df)
        data = dict()
        data['pca'] = pca
        return JsonResponse(data)


def plot_elbow_plot(request):

    data = {}
    if request.method == 'POST':
        plot_id = request.session.get('plot_id')
        key = request.POST.get("key")

        df = getdf(plot_id)
        
        index = request.POST.get("index")

        elbowmap = getElbowPlot(df,index)

        data = dict()
        data['elbowmap'] = elbowmap
        return JsonResponse(data)


def plot_kmeans_plot(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')

        index = request.POST.get("index")
        left = float(request.POST.get("left"))
        right = float(request.POST.get("right"))
        noOfCluster = int(request.POST.get("noOfCluster"))

        df = getdf(plot_id)
        heatmap = plot_kmeans_new(df,index, left, right,noOfCluster)
        
        data = dict()
        data['heatmap'] = heatmap

        return JsonResponse(data)


def plot_raincloud_plot(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')
        colors = request.POST.get('colorValues')
        print(colors)
        # index = request.POST.get("index")
        # category = request.POST.get("category")
        # dataColumn = request.POST.get("dataColumn")

        fieldsRain = request.POST.get("fieldsRain")
        fieldsRain = fieldsRain.split(',')

        df = getdf(plot_id)
        raincloud = get_raincloud_plot(df,fieldsRain)


        data = dict()
        data['raincloud'] = raincloud

        return JsonResponse(data)


def plot_box_plot(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')
        colors = request.POST.get('colorValues')

        print(colors)
        

        fieldsBox = request.POST.get("fields-box")

        fieldsBox = fieldsBox.split(',')

        df = getdf(plot_id)
        box = get_box_plot(df,fieldsBox , colors)

        data = dict()
        data['box'] = box

        return JsonResponse(data)


def plot_violin_plot(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')
        colors = request.POST.get('colorValues')



        fieldsViolin = request.POST.get("fields-violin")
        fieldsViolin = fieldsViolin.split(',')

        df = getdf(plot_id)

        violin = get_violin_plot(df,fieldsViolin , colors)
        data = dict()
        data['violin'] = violin

        return JsonResponse(data)


def plot_upset_plot(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')

        fieldsUpset = request.POST.get("fieldsUpset")
        fieldsUpset = fieldsUpset.split(',')

        df = getdf(plot_id)
        df=df.filter(fieldsUpset)

        upset_plot = get_upset_plot(df)

        data = dict()
        data['upset'] = upset_plot
        return JsonResponse(data)


def plot_scurve_plot(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')

        axis = request.POST.get("axis")

        df = getdf(plot_id)

        scurve = get_scurve_plot(df,axis)

        data = dict()
        data['scurve'] = scurve

        return JsonResponse(data)


def plot_volcano_plot(request):

    if request.method == 'POST':
        

        plot_id = request.session.get('plot_id')
        index = request.POST.get("index")
        pvalCutoff = float(request.POST.get("pvalCutoff"))
        fcCutoff = float(request.POST.get("fcCutoff"))

        pvColumn = request.POST.get("pvColumn")
        fcColumn = request.POST.get("fcColumn")

        new_colours = request.POST.get("colorValues")    

        lableProtiens = request.POST.get("lableProtiens")
        
        topInput = 5
        if lableProtiens == "top":
            topInput = int(request.POST.get("topInput"))

        df = getdf(plot_id)
        volcano = plot_volcano_new(df,index, pvalCutoff,fcCutoff, pvColumn,fcColumn, lableProtiens, topInput , new_colours)
        data = dict()
        data['volcano'] = volcano
        return JsonResponse(data)



def plot_venn_plot(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')
        colors = request.POST.get('colorValues')
        print(colors)
        sets = request.POST.get("sets")
        sets = sets.split(',')        
        df = getdf(plot_id)
        venn = get_venn_plot(df,sets)
        data = dict()
        data['venn'] = venn
        return JsonResponse(data)


def plot_density_plot(request):

    if request.method == 'POST':
        plot_id = request.session.get('plot_id')
        colors = request.POST.get('colorValues')
        print(colors)
        fieldsDensity = request.POST.get("fields-density")
        fieldsDensity = fieldsDensity.split(',')
        df = getdf(plot_id)
        density = get_density_plot(df,fieldsDensity)

        data = dict()
        data['density'] = density
        return JsonResponse(data)


def load_example_lablled(request):
    q = Example.objects.get(usethis = True)
    file = q.file.path
    download_path = os.path.join(settings.MEDIA_ROOT, file)

    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return JsonResponse(response)
    raise Http404


def load_plot_example(request):
    if request.method == 'GET':
        plottype = request.GET.get("plotType")
        request.session["plottype"] = plottype
        df = R_utils.getplot_example_s3(plottype)
        df = df.head()
        df = df.to_html()
        data = dict()
        data['exData'] = df
        return JsonResponse(data)


def load_example_lbl_free(request):
    if request.method == 'GET':
        q = Example.objects.filter(lableFreeData = True).first()
        file = q.file.path
        df = pd.read_excel(file, engine='openpyxl', nrows = 5)
        df = df.head().to_html()
        data = dict()
        data['exData'] = df
        return JsonResponse(data)


def load_example_diffex(request):
    if request.method == 'GET':
        q = Example.objects.filter(name = "differential example").first()

        file = q.file.path
        samples = q.number_of_sample

        df = pd.read_csv(file, nrows = 5)
        df = df.head().to_html()
        data = dict()
        data['exData'] = df
        data['samples'] = samples
        return JsonResponse(data)



def load_example_technical(request):
    q = Example.objects.filter(usethis = True).first()
    file = q.file.path
    samples = q.number_of_sample
    if file.endswith(".csv"):
        df = pd.read_csv(file, nrows = 5)
    else:
        df = pd.read_excel(file, engine='openpyxl', nrows = 5)

    df = df.head().to_html()
    data = dict()
    data['exData'] = df
    data['samples'] = samples
    return JsonResponse(data)    


def load_example_biological(request):
    q = Example.objects.filter(biological_rep = True).first()
    file = q.file.path
    bioReps = q.number_of_batches
    samples = q.number_of_sample

    if file.endswith(".csv"):
        df = pd.read_csv(file, nrows = 5)
    else:
        df = pd.read_excel(file, engine='openpyxl', nrows = 5)

    df = df.head().to_html()
    data = dict()
    data['exData'] = df
    print(samples)
    data['bioReps'] = bioReps
    data['samples'] = samples

    return JsonResponse(data)    


def download_example_diffex(request):

    q = Example.objects.filter(name = "differential example").first()
    file = q.file.path

    download_path = os.path.join(settings.MEDIA_ROOT, file)

    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    raise Http404


def down_plot_ex_data(request):
    plottype = request.session.get("plottype")
    q = PlotExample.objects.get(plot_type = plottype)
    file = q.file.path
    
    download_path = os.path.join(settings.MEDIA_ROOT, file)
    
    
    if os.path.exists(download_path):
        with open(download_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="text/csv")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(download_path)
            return response
    
    raise Http404
