# init
library(plyr)
library(dplyr)
library(tidyr)

# Inlezen data en enkele aanpassingen

poelbio <- read.csv("bio_poelen_id.txt", stringsAsFactors = FALSE)
colnames(poelbio)[1] <- "BioID"
poeldist <- read.csv("poelen_dist.txt", stringsAsFactors = FALSE)
colnames(poeldist)[1] <- "ID"
poeldist <- poeldist[,c(1,3,4,5,6,7)]
lgb <- read.csv("lgb_agg_buff.txt", stringsAsFactors = FALSE)
colnames(lgb)[1] <- "ID"
lndbw  <- read.csv("lndbw_agg_buffer_biopct.txt")
colnames(lndbw)[1] <- "ID"

# Data frame met de procentuele oppervlaktes voor elk landgebruik per poel (buffer = 500m)
lgbspr <- lgb %>%
  dplyr::filter(BUFF == 500) %>%                    # selectie buffer
  dplyr::select(BW_ID, LGB_ORCA, buff_pct) %>%      # Selectie variabelen
  tidyr::spread(LGB_ORCA, buff_pct, fill = 0) %>%   # voor elk landgebruik een kolom met de oppervlakte (%)
  dplyr::select(-overig)                            # de variabele 'overig' laten vallen (komt weinig voor)

#lgbspr <- dplyr::select(lgbspr, BW_ID, akker, cultuurgrasland) # Als check, de berekening met slecht 2 variabelen

# Landgebruik als matrix
m_lgbspr <- lgbspr %>%
  dplyr::select(-BW_ID) %>%
  as.matrix()

# Covariantiematrix voor de landgebruiken van de poelen
cov_lgb <- lgbspr %>%
  dplyr::select(-BW_ID) %>%
  as.matrix() %>%
  cov()

# Vector met de ID's van de bio poelen
v_poelbio <- as.character(poelbio[,"BW_ID"])

# Een functie die voor elke biopoel, de 10 meest nabije (gelijkaardige)
# conventionele poelen teruggeeft
getSimilarLandUse <- function(mypoel = "BW_40557") {

  # Een vector met het langebruik van de biopoel
  v_lgb_mypoel <- lgbspr %>%
    dplyr::filter(BW_ID == mypoel) %>%  # Selecteer de rij met 'mypoel'
    dplyr::select(-BW_ID) %>%           # Selecteer alle kolommen behalve BW_ID
    as.numeric()                        # in vector formaat

  # Bereken afstanden
  d_maha <- mahalanobis(m_lgbspr, v_lgb_mypoel, cov_lgb)
  d_norm <- mahalanobis(m_lgbspr, v_lgb_mypoel, diag(diag(cov_lgb)))
  d_eucl <- mahalanobis(m_lgbspr, v_lgb_mypoel, diag(dim(cov_lgb)[1]))
  
  # Afstanden samenvoegen met IDs poelen en ID biopoel
  df <- data.frame(bio_ID = mypoel, BW_ID = lgbspr[,"BW_ID"], 
                   d_maha, d_norm, d_eucl, stringsAsFactors = FALSE)
  
  # Lijst met conv poelen < 3km van de biopoel 'mypoel'
  conv_sel <- poeldist %>%
    dplyr::filter(BW_ID == mypoel) %>%
    dplyr::select(BW_ID_ConPoel) %>%
    .[["BW_ID_ConPoel"]]
  
  # Sorteren van laag naar hoog
  df <- df %>%
    dplyr::filter(BW_ID %in% conv_sel) %>%  # Selecteer enkel conv poelen in omgeving
    dplyr::arrange(d_maha)              # Sorteren volgens afstandsmaat mahalanobis

  # Ranking berekenen voor elke afstandsmaat
  df$r_maha <- rank(df$d_maha)
  df$r_norm <- rank(df$d_norm)
  df$r_eucl <- rank(df$d_eucl)
  
  # Toevoegen van de geografische afstand tussen de biopoel en conventionele poelen
  df <- df %>%
    dplyr::left_join(poeldist, by = c("bio_ID" = "BW_ID", "BW_ID" = "BW_ID_ConPoel"),) %>%   # geografische afstand
    dplyr::left_join(lgbspr, by = "BW_ID")     # landgebruik conventionele poelen 

  return(df)
}  


output <- data.frame()  # Lege  data frame

# Loop door alle biopoelen, voer de functie getSimilarLanduse uit
# en voeg de resultaten aan de data frame 'output'.
for (mypoel in v_poelbio) {
  
  df <- getSimilarLandUse(mypoel)
  output <- rbind(output, df)
}

head(output)

write.csv(x = output, file = "Biopoel_all.csv")

output_top10 <- output %>%
  dplyr::group_by(bio_ID) %>%             # groepeer op bio_ID
  dplyr::top_n(-10, r_maha)                # Selecteer enkel de top10

write.csv(x = output_top10, file = "Biopoel_top10.csv")


# Enkele regels code om het verschil tussen de drie afstandsmaten
# te illlustreren.
#
# De plot geeft slechts twee assen weer (akker & cultuurgrasland)
#
# Omwille van de indormatie uit de andere assen kan de ranking soms
# verrassend overkomen.
#
# Om de berekening te illustreren met slecht twee variabelen 
# kan je regel 25 toevoegen:
# "lgbspr <- dplyr::select(lgbspr, BW_ID, akker, cultuurgrasland)"


library(ggplot2)

mypoel <- "BW_40557"    # eender welke ID van een biopoel

df_plot <- output %>%
  dplyr::filter(bio_ID == mypoel) %>%
  dplyr::top_n(-20, r_maha)
df_bio <- dplyr::filter(lgbspr, BW_ID == mypoel)

ggplot(df_plot, aes(x = akker, y = cultuurgrasland)) + geom_point() +
  geom_text(aes(label = r_maha), hjust = 2, vjust = 0) +
  geom_point(aes(x = df_bio$akker, y = df_bio$cultuurgrasland), color = "red") +
  coord_fixed()
  
ggplot(df_plot, aes(x = akker, y = cultuurgrasland)) + geom_point() +
  geom_text(aes(label = r_eucl), hjust = 2, vjust = 0) +
  geom_point(aes(x = df_bio$akker, y = df_bio$cultuurgrasland), color = "red") +
  coord_fixed()

ggplot(df_plot, aes(x = akker, y = cultuurgrasland)) + geom_point() +
  geom_text(aes(label = r_norm), hjust = 2, vjust = 0) +
  geom_point(aes(x = df_bio$akker, y = df_bio$cultuurgrasland), color = "red") +
  coord_fixed()


