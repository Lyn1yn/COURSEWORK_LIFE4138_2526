rm(list = ls())

#read in data
pumpkins <- read.csv("pumpkins_11.csv")

#heaviest pumpkin:2947.032. Variety:Pastorelli, Bruna. From: Radda in Chianti,Tuscany,Italy. when: 2021.
library(tidyverse) #Loading tidyverse
pumpkins %>% #Use pipes to pass the pumpkins data to the next step.
  arrange(desc(weight_lbs))%>% #Sort by weight_lbs in descending order
  head() #Only the first 6 rows are displayed

#change the weight in pounds (lbs) to kilograms(kg)
weight_kg <- pumpkins %>%
  summarise(weight_kg = weight_lbs * 0.453592) # do the calculate
pumpkins$weight_kg = weight_kg #add into a new column

#weight class
weight_class <- ifelse(pumpkins$weight_lbs < 500,"light", 
                       ifelse(pumpkins$weight_lbs<1000, "medium", "heavy")) #Conditional judgment, classifying by weight
pumpkins$weight_class = weight_class #add into a new column

#plot
ggplot(aes(x = weight_lbs, y = est_weight, color = weight_class), data = pumpkins)+ #aesthetic mapping
  geom_point() #Draw a scatter plot

#select 3 countries
pumpkins %>% #Use pipes to pass the pumpkins data to the next step.
  count(country) #Count the number of rows for each country.
pumpkins_filtered <- pumpkins %>% #Save the filtered results into a new data frame.
  filter(country == "Australia" |country == "Austria" | country == "Belgium") #selected 3 countries
write.csv(pumpkins_filtered, "pumpkins_filtered.csv", row.names = FALSE) #save as csv file.


#summaries:
#a: mean weight: Australia-179, Austria-608, Belgium-1077. Belgium is the highest.
pumpkins_filtered %>% #Pass the pumpkins_filtered data frame as input to the next step ggplot().
  group_by(country) %>% #group according to country
  summarise(mean_weight = mean(weight_lbs), na.rm = TRUE) #Calculate the average weight (lbs) for each country.
#b: lowest weight: Australia Autumn Gold
mean_value <- pumpkins_filtered %>% 
  group_by(country, variety) %>% #group according to both country and variety.
  summarise(mean_weight = mean(weight_lbs), na.rm = TRUE) %>% #Calculate the average weight for each country and variety.
  arrange(mean_weight) #Sort by mean_weight in ascending order
head(mean_value) #Only the first 6 rows are displayed, which are the lightest 6 country and variety combinations.


#boxplot
pumpkins_filtered %>% #Pass the pumpkins_filtered data frame as input to the next step ggplot().
  ggplot(aes(x = country, y = weight_lbs))+ #Initialize the ggplot layer and set the default aesthetic mapping.
  geom_boxplot() #Add a box plot layer


#facet plot
pumpkins_filtered %>% #Pass the pumpkins_filtered data frame as input to the next step ggplot().
  ggplot(aes(x = variety, y = weight_lbs))+ #Initialize the ggplot layer and set the default aesthetic mapping.
  geom_boxplot()+ #Add a box plot layer
  facet_wrap(~country, ncol=1)+ #facet according to country
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #Adjust the display method of X-axis text