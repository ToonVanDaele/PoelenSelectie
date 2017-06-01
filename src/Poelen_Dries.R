# init
library(plyr)
library(dplyr)
library(tidyr)

# Inlezen data en enkele aanpassingen

poelbio <- read.csv("bio_poelen_id.txt")
colnames(poelbio)[1] <- "BioID"
poeldist <- read.csv("poelen_dist.txt")
colnames(poeldist)[1] <- "ID"
poeldist <- poeldist[,c(1,3,4,5,6,7)]
lgb <- read.csv("lgb_agg_buff.txt")
colnames(lgb)[1] <- "ID"
lndbw  <- read.csv("lndbw_agg_buffer_biopct.txt")
colnames(lndbw)[1] <- "ID"

# Oppervlaktes voor elk landgebruik per poel buffer = 500m
lgbspr <- lgb %>%
  dplyr::filter(BUFF == 500) %>%
  dplyr::select(BW_ID, LGB_ORCA, buff_pct) %>%
  tidyr::spread(LGB_ORCA, buff_pct, fill = 0)


# Een functie die voor elke biopoel, de 10 meest nabije (gelijkaardige)
# conventionele poelen teruggeeft
getSimilarLandUse <- function(mypoel = "BW_40557") {

  # Maak een vector met het langebruik van de biopoel
  lgb_mypoel <- dplyr::filter(lgbspr, BW_ID == mypoel)
  v_lgb_mypoel <- as.numeric(lgb_mypoel[1,-c(1)])
  
  # Selectie van de conventionele poelen nabij (< 3km) de biopoel
  conv_sel <- dplyr::filter(poeldist, BW_ID == mypoel)
  conIDs <- conv_sel[,"BW_ID_ConPoel"]
  
  # Maak een matrix met het landgebruik van de conventionele poelen
  lgb_conv <- dplyr::filter(lgbspr, BW_ID %in% conIDs)
  mat_lgb_conv <- as.matrix(lgb_conv[,-c(1)])
  
  # Bereken afstanden
  #MD <- mahalanobis(mat_lgb_conv, v_lgb_mypoel, cov(mat_lgb_conv)) # eventueel later
  MD <- mahalanobis(mat_lgb_conv, v_lgb_mypoel, diag(9))
  
  # afstanden samenvoegen met IDs conventionele poelen en id biopoel
  df_dist <- data.frame(bio_ID = mypoel, conv_ID = lgb_conv[,"BW_ID"], dist = MD)
  
  # Sorteren van laag naar hoog en de 10 eerst 10 records nemen
  df_dist_top10 <- df_dist %>%
    dplyr::arrange(dist) %>%
    dplyr::slice(1:10)
  
  return(df_dist_top10)
}  

# Vector met bio poelen
v_poelbio <- as.character(poelbio[,"BW_ID"])

# Lege  data frame
output <- data.frame()

# Loop door alle biopoelen, voer de functie getSimilarLanduse uit
# en voeg de resultaten onderaan toe in de data frame output.
for (mypoel in v_poelbio) {
  
  top10 <- getSimilarLandUse(mypoel)
  output <- rbind(output, top10)
}

head(output)

write.csv(x = output, file = "Biopoel_top10.csv")
