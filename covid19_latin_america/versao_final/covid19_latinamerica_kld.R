###############
# PACKAGE SETUP
###############

pacman::p_load(benford.analysis, BenfordTests, corrplot, data.table, dplyr, formattable, ggplot2, ggpubr, gridExtra, LaplacesDemon, leaflet, kableExtra, philentropy, rgdal, scales, scater, tidyr, tidyverse, writexl)

###############
# DATA WRANGLING
###############

general_data<-fread("https://covid.ourworldindata.org/data/owid-covid-data.csv")

countries<-c("Argentina", "Brazil", "Bolivia", "Chile", "Colombia", "Costa Rica", "Cuba", "Dominican Republic", "Ecuador", "El Salvador", "Guatemala", "Haiti", "Honduras", "Mexico", "Nicaragua", "Panama", "Paraguay", "Peru", "Uruguay", "Venezuela")

latin_america<-general_data %>% filter(location %in% countries)

###############
# KLD ANALYSIS
###############

# probabilities comparison: Kullback–Leibler divergence

latin_america <- latin_america %>% filter(new_cases>=0)

mlatin <- latin_america %>% group_by(location) %>% mutate(row = row_number()) %>% select(location, new_cases, row)

# filters data set to make sure that all countries have the same length of data
result <- mlatin %>% group_by(location) %>% filter(row == max(row))
mlatin <- mlatin %>% filter(row<=min(result$row)) 

# pivots the data.frame to wider for use with KL algorithm
mlatinw <- mlatin %>% pivot_wider(names_from = row, values_from = new_cases) %>% remove_rownames %>% tibble::column_to_rownames(var="location") 

mlatinwaux <- mlatin %>% pivot_wider(names_from = location, values_from = new_cases) %>% select(-row)

#estimates KL matrix
mlatinKL <- KL(as.matrix(mlatinw), est.prob = "empirical")

# changes the variables for country names
latin_america_names <- c('Argentina', 'Brazil', 'Bolivia', 'Chile', 'Colombia', 'Costa Rica', 'Cuba', 'Dominican Republic', 'Ecuador', 'El Salvador', 'Guatemala', 'Haiti', 'Honduras', 'Mexico', 'Nicaragua', 'Panama', 'Paraguay', 'Peru', 'Uruguay', 'Venezuela')
rownames(mlatinKL) <- latin_america_names
colnames(mlatinKL) <- latin_america_names

# generates table
round(mlatinKL, 3) %>% kbl() %>% kable_minimal()

mlatinKLDF <- as.data.frame(mlatinKL)

###############
# FIGURE 2
###############
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(round(mlatinKL, 2), method="color", col=col(200),  
         type="lower", order="hclust", 
         number.cex= 15/ncol(mlatinKL),
         is.corr=FALSE,
         addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         diag=FALSE 
)

###############
# FIGURE 3
###############
# creates a heatmap with dendogram
heatmap(mlatinKL)

# calculates the mean divergence for each country and prints styled table
mlatinKLMean <- data.frame(KLDmean = apply(mlatinKL, 1, mean))
mlatinKLMedian <- data.frame(KLDmedian = apply(mlatinKL, 1, median))

ggplot(data = mlatinKLMean, aes(x=names(mlatinKLMean), y = KLDmean)) + 
  geom_boxplot(outlier.color = "red", outlier.size = 3) + 
  stat_summary(aes(label = paste(row.names(mlatinKLMean)[which.max(mlatinKLMean$KLDmean)],"(", 
                                 round(stat(y), 1),")")),
               geom = "text", 
               fun.y = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
               hjust = -1)+
  scale_x_discrete() +
  labs(title = "", y = "", x="")+
  theme_classic()+
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        axis.text = element_text(size = 9),
        legend.text=element_text(size=9),
        legend.position = "bottom")

round(mlatinKLMean, 3) %>% kbl()
mlatinKLMean2 <- round(mlatinKLMean, 3) %>% arrange(desc(KLDmean), .by_group = TRUE)
mlatinKLMean2$KLDmean <- ifelse(
  mlatinKLMean2$KLDmean > 2,
  cell_spec(mlatinKLMean2$KLDmean, color = "red", bold = T),
  cell_spec(mlatinKLMean2$KLDmean, color = "green", italic = T)
)
kbl(mlatinKLMean2, escape = F) %>% kable_minimal() 

# filters only countries with KLD mean >= 2
mlatinKLMeanOut <- mlatinKLMean %>% filter (KLDmean >= 2)
countriesOut <- row.names(mlatinKLMeanOut)
latin_americaOut <- latin_america%>%filter(location %in% countriesOut)
latin_americaOut$fatorZero <- ifelse(latin_americaOut$new_cases == 0, 1, 0)

###############
# FIGURE 4
###############
# plots time series of cases for all countries on KLD anomalies list
ggplot(data = latin_americaOut,
       mapping = aes(x = as.Date(date), y = new_cases, color = fatorZero)) + 
  geom_jitter(alpha = 0.3) + scale_color_gradient(low="blue", high="red") +
  scale_x_date(date_breaks = "3 month", date_labels = "%b-%Y")+
  theme_classic()+
  facet_wrap(facets = vars(location), scales = "free") + 
  xlab("") + ylab("New Cases")+
  theme(strip.text.x = element_text(face = "bold", size = 10, hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 90),
        legend.text=element_text(size=10),
        legend.position = "bottom") 

###############
# DATA WRANGLING MAP
###############

america_map <- readOGR("Americas.shp", GDAL1_integer64_policy = TRUE)

latin_america_map <- subset(america_map, america_map$COUNTRY %in% c("Argentina", "Brazil", "Bolivia", "Chile", "Colombia", "Costa Rica", "Cuba", "Dominican Republic", "Ecuador", "El Salvador", "Guatemala", "Haiti", "Honduras", "Mexico", "Nicaragua", "Panama", "Paraguay", "Peru", "Uruguay", "Venezuela"))

mlatinKLMeanMap <- data.frame(names = row.names(mlatinKLMean), mlatinKLMean)

latin_america_map <- cbind(latin_america_map, mlatinKLMeanMap)

pal <- colorNumeric( palette = "YlOrRd", domain = latin_america_map$KLDmean
)

###############
# MAP
###############

leaflet(latin_america_map) %>%  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5, opacity = 1.0, fillOpacity = 0.5, fillColor = ~colorQuantile("YlOrRd", KLDmean)(KLDmean), highlightOptions = highlightOptions(color = "white", weight = 2, bringToFront = TRUE)) %>% addLegend("bottomright", pal = pal, values = ~KLDmean, opacity = 1, title = "KLD (Mean)",)  %>% addProviderTiles("Esri.WorldImagery")
