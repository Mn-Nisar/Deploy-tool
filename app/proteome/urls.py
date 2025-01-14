from django.urls import path
from django.contrib.auth.decorators import login_required

from django.views.generic import TemplateView
from . import views

app_name = 'proteome'

urlpatterns = [
    path('',views.home, name = 'home'),
    path('input/', views.input, name = 'input'),
    path('inputfile/', views.inputf, name='inputfile'),
    path('inputf_dexp/', views.inputf_dexp, name='inputf_dexp'),
    path('download/',views.downloadfile, name='download_file'),
    path('analaze_cols/',views.analaze_cols, name = 'analaze_cols'),
    path('analaze_cols_bio/',views.analaze_cols_bio, name = 'analaze_cols_bio'),
    path('pvalue/',views.pvalues, name = 'pvalues'),
    path('exampledown/',views.exampledown, name = 'exampledown'),
    path('exampledownbio/',views.exampledownbio, name = 'exampledownbio'),
    
    path('faqs/',TemplateView.as_view(template_name = 'proteome/faqs.html'), name = 'faqs'),
    path('contact/',TemplateView.as_view(template_name = 'proteome/contact.html'), name = 'contact'),
    path('team/',TemplateView.as_view(template_name = 'proteome/team.html'), name = 'team'),
    path('download_ddf',views.download_ddf, name='download_ddf'),
    path('download_example_diffex',views.download_example_diffex, name='download_example_diffex'),
    path('inputDexp/', TemplateView.as_view(template_name='proteome/input_dexp.html'), name = 'inputDexp'),
    path('download_geneontology',views.download_geneontology, name='download_geneontology'),
    path('kmeans/', views.kmeans, name='kmeans'),
    path('kClusterPlot/', views.kClusterPlot, name='kClusterPlot'),
    path('heirarchial_heatmap/', views.heirarchial_heatmap, name='heirarchial_heatmap'),
    
    path('load_example_technical',views.load_example_technical, name='load_example_technical'),
    path('load_example_biological',views.load_example_biological, name='load_example_biological'),

    path('string_db_download', views.string_db_download, name = 'string_db_download'),
    path('download_cluster_df', views.download_cluster_df, name = 'download_cluster_df'),
    path('contaminant_down', views.contaminant_down, name = 'contaminant_down'),
    path('gene_ontology_ajax', views.gene_ontology_ajax, name = 'gene_ontology_ajax'),
    path('string_ajax', views.string_ajax, name = 'string_ajax'),
    path('goForClusters', views.goForClusters, name = 'goForClusters'),
    path('kegg_pathway', views.kegg_pathway, name = 'kegg_pathway'),
    path('load_example_lablled', views.load_example_lablled, name = 'load_example_lablled'),
    path('load_example_diffex', views.load_example_diffex, name = 'load_example_diffex'),
    path('download_resedue', views.download_resedue, name = 'download_resedue'),
    # lable free
    path('input_label_free',views.input_label_free, name = 'input_label_free'),
    path('inputfile_lblfree/', views.inputfile_lblfree, name='inputfile_lblfree'),
    path('download_ibaq/',views.download_ibaq, name = 'download_ibaq'),
    path('download_zip/',views.download_zip, name = 'download_zip'),
    # allplts
    #we have add login required here
    path('plots/',login_required(TemplateView.as_view(template_name='proteome/allplots.html')), name = 'allplots'),

    path('getplot/',views.getplot, name = 'getplot'),
    path('load_plot_example',views.load_plot_example, name = 'load_plot_example'),
    path('plot_circbar/',views.plot_circbar, name = 'plot_circbar'),
    path('plot_heatmap_plot/',views.plot_heatmap_plot, name = 'plot_heatmap_plot'),
    path('plot_volcano_plot/',views.plot_volcano_plot, name = 'plot_volcano_plot'),
    path('plot_bubble_plot/',views.plot_bubble_plot, name = 'plot_bubble_plot'),
    path('plot_pca_plot/',views.plot_pca_plot, name = 'plot_pca_plot'),
    path('plot_elbow_plot/',views.plot_elbow_plot, name = 'plot_elbow_plot'),
    path('plot_kmeans_plot/',views.plot_kmeans_plot, name = 'plot_kmeans_plot'),
    path('plot_raincloud_plot/',views.plot_raincloud_plot, name = 'plot_raincloud_plot'),
    path('plot_violin_plot/',views.plot_violin_plot, name = 'plot_violin_plot'),
    path('plot_scurve_plot/',views.plot_scurve_plot, name = 'plot_scurve_plot'),
    path('plot_venn_plot/',views.plot_venn_plot, name = 'plot_venn_plot'),
    path('plot_ma_plot/',views.plot_ma_plot, name = 'plot_ma_plot'),
    path('plot_density_plot/',views.plot_density_plot, name = 'plot_density_plot'),
    path('plot_box_plot/',views.plot_box_plot, name = 'plot_box_plot'),
    path('plot_histogram_plot/',views.plot_histogram_plot, name = 'plot_histogram_plot'),
    path('plot_upset_plot/',views.plot_upset_plot, name = 'plot_upset_plot'),
    path('perform_btch_cr/',views.perform_btch_cr, name = 'perform_btch_cr'),
    path('down_plot_ex_data',views.down_plot_ex_data, name = 'down_plot_ex_data'),
    path('load_example_lbl_free',views.load_example_lbl_free, name = 'load_example_lbl_free'),
    path('dockerpage',TemplateView.as_view(template_name = 'proteome/dockerpage.html'), name = 'dockerpage'),
    path('cyto',TemplateView.as_view(template_name = 'proteome/cyto.html'), name = 'cyto'),

]


