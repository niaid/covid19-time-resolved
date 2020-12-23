library(circlize) #for colorRamp2
severity.color <- c("Critical-Alive"="#374E55FF","Critical-Deceased"="#DF8F44FF","Moderate-Alive"="#00A1D5FF" ,"Severe-Alive"="#B24745FF","HC" = "#79AF97FF")
PC1.color <- colorRamp2(c(-3, 3.6), c("white", "#e97171"))
Class.color <- c("HC" = "#F8766D", "COVID" = "#eebb4d")
PC1class.color <- c("HC"="#00BA38", "PC1_low"="#619CFF", "PC1_high"="#F8766D")

severity.shape <- c("Critical-Alive"=15,"Critical-Deceased"=7,"Moderate-Alive"=17,"Severe-Alive"=18,"HC" = 3)
