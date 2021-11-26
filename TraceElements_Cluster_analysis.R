##--------------------------------------------------------------------------------------------------------
## SCRIPT : Regroupment of capture fisheries species from the Seychelles, according to their trace element profiles (i.e. mean trace element concentrations)
## Specific content : - Selection of trace element with detection frequencies > 70%
##                    - Substitution of values < LOQ with values drawn randomly from the range ]0;LOQ[
##                    - Scaling of mean trace element concentrations
##                    - Hierarchical cluster analysis (Ward's clustering method)
##                    - Testing significant difference (p < 0.05) in trace element concentrations among clusters (Kruskal-Wallis followed by post-hoc Dunn test)
##                    - Representation (plot) of dendrogram and heatmap plot
## As part of :
##    Sabino et al. "The role of tropical capture fisheries in trace element delivery for a Small Island Developing State community, the Seychelles"
##
## Author : Magali Sabino
## First publication on GitHub : 2021-11-XX
## Last update : -
##
## For more information on the clustering method, see:
##    - Factoextra R Package: Easy Multivariate Data Analyses and Elegant Visualization (http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization)
##    - K-means Cluster Analysis (https://uc-r.github.io/kmeans_clustering?fbclid=IwAR0XKVVCFGk83qEPkWcpbARJhcoEpuWQNZFwSe62QojO-I7J7YliUW4duo0)
##
####
##
## R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
## Copyright (C) 2021 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

### 0 // Packages ##########################################################################################

## Open libraries
lapply(c("tidyverse", "vegan", "ggdendro", "dendextend", "NbClust", "factoextra", "FSA", "ggpubr"),
       library, character.only = TRUE)



### I // Data ##############################################################################################

## 1 / Dataset importation

# The dataset used in the paper associated to this script is available from the corresponding author on
# reasonable request. Thus, here we use the type dataset given with this script (Type_Dataset_Test_Script.csv)
# containing fictitious trace element concentrations (µg.g-1 dw) and moisture contents (%) by species. In this dataset,
# each line corresponds to one individual while trace element concentrations and moisture content are in columns.

TraceElement_content <- read.table(file = "C:/Users/msabin02/Nextcloud/PUBLI - TE patterns in Seychelles seafood/Script/Type_Dataset_Test_Script.csv", # Replace by your own path or orking directory
                             header = T, dec = ".",sep = ";")


## 2 / Substitution of values < LOQ with values drawn randomly from the range ]0;LOQ[

# Creation of a function for counting decimals
count_decimals = function(x) {
  #length zero input
  if (length(x) == 0) return(numeric())

  #count decimals
  x_nchr = x %>% abs() %>% as.character() %>% nchar() %>% as.numeric()
  x_int = floor(x) %>% abs() %>% nchar()
  x_nchr = x_nchr - 1 - x_int
  x_nchr[x_nchr < 0] = 0

  x_nchr
}

# Creation of database with LOQ
LOQ_data <- TraceElement_content %>%
  filter(!is.na(TE1)) %>%
  mutate(nb_TE1 = count_decimals(TE1),
         nb_TE2 = count_decimals(TE2),
         nb_TE3 = count_decimals(TE3),
         nb_TE4 = count_decimals(TE4),
         nb_TE5 = count_decimals(TE5),
         nb_TE6 = count_decimals(TE6),
         nb_TE7 = count_decimals(TE7),
         nb_TE8 = count_decimals(TE8),
         nb_TE9 = count_decimals(TE9),
         nb_TE10 = count_decimals(TE10),
         nb_TE11 = count_decimals(TE11),
         nb_TE12 = count_decimals(TE12),
         nb_TE13 = count_decimals(TE13),
         nb_TE14 = count_decimals(TE14)) %>%
  mutate(TE1 = ifelse(nb_TE1 <= 3, TE1, NA),
         TE2 = ifelse(nb_TE2 <= 3, TE2, NA),
         TE3 = ifelse(nb_TE3 <= 3, TE3, NA),
         TE4 = ifelse(nb_TE4 <= 3, TE4, NA),
         TE5 = ifelse(nb_TE5 <= 3, TE5, NA),
         TE6 = ifelse(nb_TE6 <= 3, TE6, NA),
         TE7 = ifelse(nb_TE7 <= 3, TE7, NA),
         TE8 = ifelse(nb_TE8 <= 3, TE8, NA),
         TE9 = ifelse(nb_TE9 <= 3, TE9, NA),
         TE10 = ifelse(nb_TE10 <= 3, TE10, NA),
         TE11 = ifelse(nb_TE11 <= 3, TE11, NA),
         TE12 = ifelse(nb_TE12 <= 3, TE12, NA),
         TE13 = ifelse(nb_TE13 <= 3, TE13, NA),
         TE14 = ifelse(nb_TE14 <= 3, TE14, NA)) %>%
  select(sample_identifier,TE1,TE2,TE3,TE4,TE5,TE6,TE7,TE8,TE9,TE10,TE11,TE12,TE13,TE14) %>%
  gather(TraceElement,LOQ,-sample_identifier)

# Replace LOQ by NA in main database with trace element concentrations
TraceElement_content_LOQ <- TraceElement_content %>%
  mutate(nb_TE1 = count_decimals(TE1),
         nb_TE2 = count_decimals(TE2),
         nb_TE3 = count_decimals(TE3),
         nb_TE4 = count_decimals(TE4),
         nb_TE5 = count_decimals(TE5),
         nb_TE6 = count_decimals(TE6),
         nb_TE7 = count_decimals(TE7),
         nb_TE8 = count_decimals(TE8),
         nb_TE9 = count_decimals(TE9),
         nb_TE10 = count_decimals(TE10),
         nb_TE11 = count_decimals(TE11),
         nb_TE12 = count_decimals(TE12),
         nb_TE13 = count_decimals(TE13),
         nb_TE14 = count_decimals(TE14)) %>%
  mutate(TE1 = ifelse(nb_TE1 <= 3, NA, TE1),
         TE2 = ifelse(nb_TE2 <= 3, NA, TE2),
         TE3 = ifelse(nb_TE3 <= 3, NA, TE3),
         TE4 = ifelse(nb_TE4 <= 3, NA, TE4),
         TE5 = ifelse(nb_TE5 <= 3, NA, TE5),
         TE6 = ifelse(nb_TE6 <= 3, NA, TE6),
         TE7 = ifelse(nb_TE7 <= 3, NA, TE7),
         TE8 = ifelse(nb_TE8 <= 3, NA, TE8),
         TE9 = ifelse(nb_TE9 <= 3, NA, TE9),
         TE10 = ifelse(nb_TE10 <= 3, NA, TE10),
         TE11 = ifelse(nb_TE11 <= 3, NA, TE11),
         TE12 = ifelse(nb_TE12 <= 3, NA, TE12),
         TE13 = ifelse(nb_TE13 <= 3, NA, TE13),
         TE14 = ifelse(nb_TE14 <= 3, NA, TE14)) %>%
  select(sample_identifier,TE1,TE2,TE3,TE4,TE5,TE6,TE7,TE8,TE9,TE10,TE11,TE12,TE13,TE14) %>%
  gather(TraceElement,value,-sample_identifier) %>%
  left_join(LOQ_data, by = c("sample_identifier","TraceElement"))

# Random generation of values from the range ]0;LOQ[
for (i in 1:nrow(TraceElement_content_LOQ)){
  if (is.na(TraceElement_content_LOQ[i,]$value)){
    TraceElement_content_LOQ[i,]$value <- runif(n = 1, min = 0, max = TraceElement_content_LOQ[i,]$LOQ)
  }
}

# Creation of dataframe with all information related to samples
sample_information <- TraceElement_content %>% 
  select(sample_identifier,Species,moist)

# Cleansing database and adding information for each sample
TraceElement_content_LOQ <- TraceElement_content_LOQ %>%
  select(-LOQ) %>%
  spread(TraceElement, value) %>% 
  left_join(sample_information, by = "sample_identifier")

rm(LOQ_data,i,sample_information)


## 3 / Selection of trace element with detection frequency > 70%

# Calculate the detection frequency (DF) for each trace element
n_row = 100 # Replace by the total number of sample in your database
TraceElement_content_DF <- TraceElement_content %>%
  mutate(nb_TE1 = count_decimals(TE1),
         nb_TE2 = count_decimals(TE2),
         nb_TE3 = count_decimals(TE3),
         nb_TE4 = count_decimals(TE4),
         nb_TE5 = count_decimals(TE5),
         nb_TE6 = count_decimals(TE6),
         nb_TE7 = count_decimals(TE7),
         nb_TE8 = count_decimals(TE8),
         nb_TE9 = count_decimals(TE9),
         nb_TE10 = count_decimals(TE10),
         nb_TE11 = count_decimals(TE11),
         nb_TE12 = count_decimals(TE12),
         nb_TE13 = count_decimals(TE13),
         nb_TE14 = count_decimals(TE14)) %>%
  mutate(TE1 = ifelse(nb_TE1 <= 3, NA, TE1),
         TE2 = ifelse(nb_TE2 <= 3, NA, TE2),
         TE3 = ifelse(nb_TE3 <= 3, NA, TE3),
         TE4 = ifelse(nb_TE4 <= 3, NA, TE4),
         TE5 = ifelse(nb_TE5 <= 3, NA, TE5),
         TE6 = ifelse(nb_TE6 <= 3, NA, TE6),
         TE7 = ifelse(nb_TE7 <= 3, NA, TE7),
         TE8 = ifelse(nb_TE8 <= 3, NA, TE8),
         TE9 = ifelse(nb_TE9 <= 3, NA, TE9),
         TE10 = ifelse(nb_TE10 <= 3, NA, TE10),
         TE11 = ifelse(nb_TE11 <= 3, NA, TE11),
         TE12 = ifelse(nb_TE12 <= 3, NA, TE12),
         TE13 = ifelse(nb_TE13 <= 3, NA, TE13),
         TE14 = ifelse(nb_TE14 <= 3, NA, TE14)) %>% 
  select(TE1,TE2,TE3,TE4,TE5,TE6,TE7,TE8,TE9,TE10,TE11,TE12,TE13,TE14) %>%
  gather(TraceElement,value) %>%
  filter(!is.na(value)) %>%
  group_by(TraceElement) %>%
  summarise(n = n()) %>%
  mutate(n_tot = n_row,
         percent_above_LOQ = n*100/n_tot) %>%
  select(TraceElement,percent_above_LOQ) %>%
  filter(percent_above_LOQ >= 70)

TraceElement_content_LOQ <- TraceElement_content_LOQ %>% 
  select(sample_identifier,Species,moist,TraceElement_content_DF$TraceElement)

rm(TraceElement_content_DF, n_row)



### III // Cluster analysis ##############################################################################################

## 1 / Calculation of mean trace element concentrations for each species and each trace element

TraceElement_profiles <- TraceElement_content_LOQ %>% 
  select(Species,TE2,TE3,TE6,TE7,TE8,TE9,TE10,TE11,TE14) %>% 
  gather(TraceElement, value, -Species) %>% 
  group_by(Species, TraceElement) %>% 
  summarise(mean = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  spread(TraceElement, mean)


## 2 / Scaling data to reduce effect of value gaps between trace elements
names(TraceElement_profiles)
clustering <- as.matrix(TraceElement_profiles[,2:10])
rownames(clustering) <- TraceElement_profiles$Species
clustering <- scale(clustering)


## 3 / Data translation to avoid negative data

# Creation of a function to achieve linear translation
Linear_positiv_translation <- function(matrix){
  result <- matrix + abs(min(matrix, na.rm = TRUE)) + 0.1
}

# Data translation
clustering <- Linear_positiv_translation(clustering)
clustering_mtx <- clustering
clustering <- as.data.frame(clustering)


## 4 / Cluster analysis
dist_matrix <- vegdist(clustering, method = "euclidean") # Calculation of distance matrix with Euclidean distance
hc <- hclust(dist_matrix, method = "ward.D2") # Clustering using Ward's method


## 5 / determination of the optimal number of cluster

res.nbclust <- NbClust(clustering_mtx, distance = "euclidean",
          min.nc = 2, max.nc = 10,
          method = "ward.D2", index = "all")

fviz_nbclust(res.nbclust, ggtheme = theme_minimal()) # Here, optimal number of cluster is 3

# Compute silhouette plot
res.hc <- clustering_mtx %>%
  eclust("hclust", k = 3, graph = FALSE)
fviz_silhouette(res.hc)


## 6 / Data formatting for further plotting

graph.cluster <- as.dendrogram(hc) # Build dendogram object from hclust results
graph.cluster <- dendro_data(graph.cluster, type = "rectangle") # Extract the data for rectangle lines
clusters <- TraceElement_profiles %>% 
  select(Species) %>% 
  rename(label = Species) # Extract infos on species
graph.cluster[["labels"]] <- merge(graph.cluster[["labels"]], clusters, by = "label") # Merge infos with labels of dendogram

# Cut data into the optimal number of cluster (here, n = 3)
clust_TE <- as.data.frame(cutree(hc, k = 3)) # k is the final number of cluster
clust_TE$label <- row.names(clust_TE)
colnames(clust_TE)[1] <- c("cluster")
clust_TE$cluster <- as.character(clust_TE$cluster)

# The numbering of clusters does not always correspond to the reading direction of the dendrogram :
# sometimes need to renumber the clusters.
# Here, need to invert clusters 2 and 3.

clust_TE <- clust_TE %>%
  mutate(new_cluster = cluster,
         new_cluster = ifelse(cluster == "2", "3", new_cluster),
         new_cluster = ifelse(cluster == "3", "2", new_cluster)) %>% 
  select(-cluster)

# "new_cluster" becomes the column with the correct numbering of clusters (in the order of display on the dendrogram),
# it is thus the column that will be used for further plotting.

# With the following line, we add cluster numbers in the graph.cluster object, that will be used for plotting the dendrogram
graph.cluster[["labels"]] <- merge(graph.cluster[["labels"]], clust_TE, by = "label")


## 7 / Computing plots

dendro.plot <- ggplot()+
  geom_segment(data = segment(graph.cluster), aes(x = y, y = x, xend = yend, yend = xend))+
  geom_text(data = label(graph.cluster), aes(x = y, y = x, label = label, hjust = 0,
                                             color = new_cluster,
                                             angle = 0), size = 3)+
  scale_color_manual(values = c("#ED7340","#88BB9A","#D2062E"))+ # All html codes for colors used in the paper: "#ED7340","#FBAB19","#88BB9A","#D2062E","#048C7F"
  #coord_flip()+
  scale_x_reverse()+
  scale_y_reverse()+
  labs(x = "Weight", color = "Cluster")+
  theme(#legend.position = "left",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "italic"),
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank())

rm(clustering, hc, res.hc, dist_matrix, clusters, clustering_mtx, res.nbclust)



### IV // Testing significant difference in trace element concentrations among clusters ##############################################################################################

## 1 / Kruskal-Wallis followed by post-hoc Dunn test when significant difference (p < 0.05)

# Creating dataset with trace element concentrations to test for significant difference and the different groups (here, the clusters)
data_for_testing <- TraceElement_content_LOQ %>% 
  rename(label = Species) %>% 
  left_join(clust_TE, by = "label") %>% 
  select(-sample_identifier, -label, -moist) %>% # Removing columns with no interest; here "sample_identifier", "label" and "moist"
  mutate(new_cluster = factor(new_cluster, levels = c("1","2","3"))) # We change "new_cluster" into a factor and we precise the levels we want

# Creating datasets for the loop
TraceElement <- colnames(data_for_testing)[1:(ncol(data_for_testing)-1)]
output_all_N <- as.data.frame(unique(data_for_testing$new_cluster))
names(output_all_N) <- "new_cluster"
output_all_test <- NULL
output_all_posthoc <- NULL

# The following loop allows an automated testing of significant difference among more than 2
# groups using either ANOVA or Kruskal-Wallis, followd by  their associated post-hoc tests
# (Tukey's HSD or Dunn test, respectively). For this, Kruskal-Wallis is used when the number
# of sample is not high enough (n <= 5) or when one condition among homoscedasticity and
# normality is not respected. Homoscedasticity is tested using Fligner test and normality
# of residus is tested using Shapiro test.

for (ii in 1:length(TraceElement)){
  #ii = 2
  # Création des variables
  TE_name <- TraceElement[ii]
  data_test <- data_for_testing %>% 
    select(new_cluster,TraceElement[ii])
  data_test <- data_test[!is.na(data_test[2]),]
  data_test <- as.data.frame(data_test)
  
  output_test <- data.frame(TE = TE_name,
                            Test = NA,
                            P_val = NA,
                            stat_val = NA)
  
  # Calcul du nombre de données par espèce
  output_N <- data_test %>% 
    group_by(new_cluster) %>% 
    summarise(N = n())
  colnames(output_N)[2] <- TE_name
  output_all_N <- output_all_N %>% 
    left_join(output_N, by = "new_cluster")
  
  # Tests stats en fonction du nbre de samples
  if (min(output_N[,2]) <= 5) {
    test_stat <- kruskal.test(data_test[,2] ~ data_test[,1])
    
    output_test$Test <- "Kruskal"
    output_test$P_val <- test_stat$p.value
    output_test$stat_val <- test_stat$statistic
    
    output_all_test <- rbind(output_all_test, output_test)
    
    # Post hoc
    if (output_test$P_val < 0.05) {
      test <- dunnTest(data_test[,2] ~ data_test[,1], method = "bh")
      
      test <- as.data.frame(test[[2]])
      #test <- test[!(test[,4]) >= 0.05,]
      test[,4] <- round(test[,4], 3)
      
      output_all_posthoc[[ii]] <- test
      names(output_all_posthoc)[ii] <- paste(TE_name,"Dunn", sep = ".")
      
      rm(test)
    } else {
      output_all_posthoc[[ii]] <- NULL
      names(output_all_posthoc)[ii] <- paste(TE_name,"Dunn", sep = ".")
    }
    
    rm(test_stat)
  } else {
    
    # Test homoscedasticité & normalité des résidus
    test_a <- fligner.test(data_test[,2] ~ data_test[,1]) # Homoscedasticity
    a1 <- aov(data_test[,2] ~ data_test[,1])
    test_b <- shapiro.test(resid(a1)) # Normality
    
    if (test_a$p.value < 0.05 | test_b$p.value < 0.05){
      # Kruskal
      test_stat <- kruskal.test(data_test[,2] ~ data_test[,1])
      
      output_test$Test <- "Kruskal"
      output_test$P_val <- test_stat$p.value
      output_test$stat_val <- test_stat$statistic
      
      output_all_test <- rbind(output_all_test, output_test)
      
      # Post hoc
      if (output_test$P_val < 0.05) {
        test <- dunnTest(data_test[,2] ~ data_test[,1], method = "bh")
        
        test <- as.data.frame(test[[2]])
        #test <- test[!(test[,4]) >= 0.05,]
        test[,4] <- round(test[,4], 3)
        
        output_all_posthoc[[ii]] <- test
        names(output_all_posthoc)[ii] <- paste(TE_name,"Dunn", sep = ".")
        
        rm(test)
      } else {
        output_all_posthoc[[ii]] <- NULL
        names(output_all_posthoc)[ii] <- paste(TE_name,"Dunn", sep = ".")
      }
      
      rm(test_stat)
    } else {
      # ANOVA
      lmFA <- lm(data_test[,2] ~ data_test[,1])
      test_stat <- anova(lmFA)
      rm(lmFA)
      
      output_test$Test <- "ANOVA"
      output_test$P_val <- test_stat[1,5]
      output_test$stat_val <- test_stat[1,4]
      
      output_all_test <- rbind(output_all_test, output_test)
      
      # Post-hoc
      if (output_test$P_val < 0.05) {
        test <- TukeyHSD(a1, 'data_test$new_cluster', conf.level = 0.95)
        output_all_posthoc[[ii]] <- as.data.frame(test[[2]])
        names(output_all_posthoc)[ii] <- paste(TE_name,"Tuke", sep = ".")
        
        rm(test)
      } else {
        output_all_posthoc[[ii]] <- NULL
        names(output_all_posthoc)[ii] <- paste(TE_name,"Tuke", sep = ".")
      }
      
      rm(test_stat)
    }
    
    rm(a1,test_a,test_b)
  }
  
  rm(TE_name,data_test,output_test,output_N)
}

output_all_test$P_val <- round(output_all_test$P_val, 3)

# After the loop, the object "output_all_N" contains the number of values used for the tests,
# for each cluster and each trace element tested.
# The object "output_all_test" contains all results of statistical tests (either ANOVA or
# Kruskal-Wallis).
# The object "output_all_posthoc" contains all results of post-hoc tests (either Tukey's HSD
# or Kruskal-Wallis.

rm(ii, data_for_testing, TraceElement)

## 2 / Plotting the mean (+/- SD) concentration for each trace element in each cluster

# This step allows to determine which cluster has the highest mean value and which has the
# lowest mean value, in combination with post-hoc test results.

TraceElement_content_LOQ %>% 
  rename(label = Species) %>% 
  left_join(clust_TE, by = "label") %>% 
  select(-sample_identifier, -label, -moist) %>% 
  gather(TraceElement, value, -new_cluster) %>% 
  group_by(TraceElement,new_cluster) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(TraceElement = factor(TraceElement, levels = c("TE2","TE3","TE6","TE7","TE8","TE9","TE10","TE11","TE14"))) %>% 
  ggplot()+
  geom_errorbar(aes(x = new_cluster, y = mean, ymin = mean-(sd/5), ymax = mean+sd), color = "#404041")+
  geom_bar(aes(x = new_cluster, y = mean, fill = new_cluster), stat = "identity", color = "#404041")+
  #geom_text(aes(x = cluster, y = mean+se+(1.3*se), label = letter), color = "#404041")+
  labs(y = "Concentrations (µg.g-1 ww)")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")+
  facet_wrap(TraceElement ~ ., ncol = 4, nrow = 3, scales = "free")+
  theme(strip.text.x = element_text(face = "bold"))

mean_values <- TraceElement_content_LOQ %>% 
  rename(label = Species) %>% 
  left_join(clust_TE, by = "label") %>% 
  select(-sample_identifier, -label, -moist) %>% 
  gather(TraceElement, value, -new_cluster) %>% 
  group_by(TraceElement,new_cluster) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE)) %>% 
  ungroup()

# To display post-hoc test results, from the object "output_all_posthoc"
names(output_all_posthoc)
output_all_posthoc[[3]] # Change the number to display the results you want

# According to the plot and to post-hoc test results:
# - for TE2, Cluster 1 has the highest mean value while Clusters 2 and 3 have the lowest
# - for TE3, Cluster 1 has the highest mean value, Cluster 2 has an intermediate value and Cluster 3 has the lowest
# - for TE6, Cluster 1 has the highest mean value while Clusters 2 and 3 have the lowest
# - for TE7, Cluster 2 has the highest mean value while Clusters 1 and 3 have the lowest
# - for TE8, Cluster 3 has the highest mean value, Cluster 2 has an intermediate value and Cluster 1 has the lowest
# - for TE9, Cluster 3 has the highest mean value while Clusters 1 and 2 have the lowest
# - for TE10, Cluster 1 and 2 have the highest mean value while Clusters 3 has the lowest
# - for TE11, Cluster 3 has the highest mean value, Cluster 2 has an intermediate value and Cluster 1 has the lowest
# - for TE14, Cluster 1 has the highest mean value while Clusters 2 and 3 have the lowest


## 3 / Plotting the heatmap for rapid visualisation of differences among clusters

heatmap.plot <- clust_TE %>% 
  # We give a specific value for each trace element in each cluster:
  # - if highest mean concentration for a given trace element, give 1
  # - if lowest mean concentration, give -1
  # - if intermediate concentration, give 0
  mutate(TE2 = ifelse(new_cluster == "1", "1", "-1"), # The ifelse function is written this way: ifelse(condition, what to do if "yes", what to do if "no")
         TE3 = ifelse(new_cluster == "1", "1", "-1"),
         TE3 = ifelse(new_cluster == "2", "0", TE3),
         TE6 = ifelse(new_cluster == "1", "1", "-1"),
         TE7 = ifelse(new_cluster == "2", "1", "-1"),
         TE8 = ifelse(new_cluster == "3", "1", "-1"),
         TE8 = ifelse(new_cluster == "2", "0", TE8),
         TE9 = ifelse(new_cluster == "3", "1", "-1"),
         TE10 = ifelse(new_cluster == "3", "-1", "1"),
         TE11 = ifelse(new_cluster == "3", "1", "-1"),
         TE11 = ifelse(new_cluster == "2", "0", TE11),
         TE14 = ifelse(new_cluster == "1", "1", "-1")) %>% 
  select(-new_cluster) %>% 
  gather(TraceElement, value, -label) %>% 
  # Define the levels (in factors) for the labels (i.e. species names) and the values defined above
  # This step is necessary to have the same orders than for the dendrogram.
  # /!\ To have the same ordre for the species' names, give them in the reverse order (from bottom to top)
  mutate(label = factor(label, levels = c("SP12","SP10","SP11","SP13","SP6","SP18","SP17","SP19","SP15",
                                          "SP9","SP5","SP8","SP7","SP20","SP14","SP16","SP3","SP1","SP2","SP4")),
         value = factor(value, levels = c("1","0","-1")),
         TraceElement = factor(TraceElement, levels = c("TE2","TE3","TE6","TE7","TE8","TE9","TE10","TE11","TE14"))) %>% 
  ggplot()+
  geom_tile(aes(x = TraceElement, y = label, fill = value))+
  scale_fill_manual(values = c("#3B6078","#F1F1F2","#99BFC8"))+ # Colors given to highest, intermediate and lowest mean concentrations
  scale_x_discrete(position = "top")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")


## 3 / Plotting the dendrogram and the heatmap on the same plot (side by side)

ggarrange(dendro.plot, heatmap.plot,
          align = "hv", widths = c(0.7, 1))

