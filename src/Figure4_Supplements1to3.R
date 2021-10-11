# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#move up a directory
setwd("..")

#install.packages('pdftools')

pdftools::pdf_convert(pdf='outputs/sem/PlotSEMpathCompProdInvasion_PartialMediation_Diversity.lavfull.pdf',  filenames='outputs/figures/Figure4_Supplement1.tiff', dpi=300, format='tiff')
pdftools::pdf_convert(pdf='outputs/sem/PlotSEMpathCompProdInvasion_CompleteMediation_Diversity.lavfull.pdf',  filenames='outputs/figures/Figure4_Supplement2.tiff', dpi=300, format='tiff')
pdftools::pdf_convert(pdf='outputs/sem/PlotSEMpathCompProdInvasion_NoMediation_Diversity.lavfull.pdf',  filenames='outputs/figures/Figure4_Supplement3.tiff', dpi=300, format='tiff')